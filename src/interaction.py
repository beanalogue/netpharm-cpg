"""
단백질-리간드 상호작용 분석 + 2D 다이어그램
Vina 도킹 PDBQT → 접촉 잔기 검출 → Wheel 다이어그램 + rdkit 리간드 2D

접촉 분류 기준 (거리 규칙, 문헌 표준):
  H-Bond      : N/O ↔ N/O  ≤ 3.5 Å
  Salt Bridge : 하전 잔기 ↔ 반대극성 원자 ≤ 4.0 Å
  Hydrophobic : C↔C (비극성 잔기) ≤ 5.0 Å
  Pi-Stacking : 방향족 잔기 ↔ 리간드 방향족 C ≤ 5.5 Å
  Van der Waals: 기타 ≤ 4.0 Å

ProLIF + MDAnalysis 설치 시 자동 전환:
  pip install prolif mdanalysis
"""

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

RESULTS_DIR = Path(__file__).parent.parent / "results"
RESULTS_DIR.mkdir(exist_ok=True)

# ── ProLIF 선택적 로드 ────────────────────────────────────────────────────────
try:
    import prolif
    import MDAnalysis as mda
    PROLIF_OK = True
    log.info("[Interaction] ProLIF + MDAnalysis 사용 가능")
except ImportError:
    PROLIF_OK = False

# ── rdkit 선택적 로드 ─────────────────────────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdDepictor, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False

# ── 잔기 분류 테이블 ──────────────────────────────────────────────────────────
HYDROPHOBIC_RES  = {"ALA","VAL","ILE","LEU","MET","PHE","PRO","TRP","TYR","CYS"}
AROMATIC_RES     = {"PHE","TYR","TRP","HIS"}
CHARGED_POS_RES  = {"LYS","ARG","HIS"}
CHARGED_NEG_RES  = {"ASP","GLU"}
HBOND_RES        = {"SER","THR","ASN","GLN","TYR","TRP","HIS",
                     "ASP","GLU","LYS","ARG","CYS"}

INTERACTION_COLORS = {
    "H-Bond":       "#2196F3",
    "Salt Bridge":  "#F44336",
    "Hydrophobic":  "#FF9800",
    "Pi-Stacking":  "#9C27B0",
    "Van der Waals":"#78909C",
}


# ── 데이터 구조 ───────────────────────────────────────────────────────────────
@dataclass
class Atom3D:
    serial:   int
    atom_name:str
    res_name: str
    chain_id: str
    res_seq:  int
    x: float; y: float; z: float
    element:  str = ""

    @property
    def res_id(self) -> str:
        return f"{self.res_name}{self.res_seq}{self.chain_id}"


@dataclass
class ContactResidue:
    res_id:     str
    res_name:   str
    res_seq:    int
    chain_id:   str
    interaction:str         # "H-Bond" / "Hydrophobic" / ...
    min_dist:   float       # closest contact distance (Å)
    lig_atoms:  list = field(default_factory=list)
    prot_atoms: list = field(default_factory=list)


# ─────────────────────────────────────────────────────────────────────────────
# 1. PDB / PDBQT 파서
# ─────────────────────────────────────────────────────────────────────────────
def _parse_atom_line(line: str) -> Optional[Atom3D]:
    """ATOM / HETATM 레코드 → Atom3D (실패 시 None)"""
    if len(line) < 54:
        return None
    try:
        return Atom3D(
            serial    = int(line[6:11]),
            atom_name = line[12:16].strip(),
            res_name  = line[17:20].strip(),
            chain_id  = line[21].strip() or "A",
            res_seq   = int(line[22:26]),
            x         = float(line[30:38]),
            y         = float(line[38:46]),
            z         = float(line[46:54]),
            element   = line[76:78].strip() if len(line) > 76 else line[12:14].strip().lstrip("0123456789"),
        )
    except (ValueError, IndexError):
        return None


def parse_pdb_protein(pdb_path: Path) -> list[Atom3D]:
    """수용체 PDB → ATOM 레코드만 파싱"""
    atoms = []
    for line in pdb_path.read_text(errors="replace").splitlines():
        if line.startswith("ATOM"):
            a = _parse_atom_line(line)
            if a:
                atoms.append(a)
    return atoms


def extract_first_pose_pdbqt(pdbqt_path: Path, out_dir: Path) -> Optional[Path]:
    """Vina 출력 PDBQT → 첫 번째 포즈(MODEL 1) PDB로 저장"""
    lines = pdbqt_path.read_text(errors="replace").splitlines()
    pose_lines, in_model = [], False
    for line in lines:
        if line.startswith("MODEL"):
            in_model = True
            continue
        if line.startswith("ENDMDL"):
            break
        if in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
            # PDBQT의 리간드 잔기명을 UNL로 통일
            if len(line) >= 20:
                line = line[:17] + "UNL" + line[20:]
            pose_lines.append(line)

    if not pose_lines:
        log.warning(f"[Interaction] {pdbqt_path.name}: 첫 번째 포즈 없음")
        return None

    out = out_dir / (pdbqt_path.stem + "_pose1.pdb")
    out.write_text("\n".join(pose_lines))
    return out


def parse_ligand_pdb(lig_pdb_path: Path) -> list[Atom3D]:
    """리간드 PDB (HETATM 또는 ATOM) 파싱"""
    atoms = []
    for line in lig_pdb_path.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            a = _parse_atom_line(line)
            if a:
                atoms.append(a)
    return atoms


# ─────────────────────────────────────────────────────────────────────────────
# 2. 거리 기반 접촉 검출
# ─────────────────────────────────────────────────────────────────────────────
def _dist(a: Atom3D, b: Atom3D) -> float:
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)


def _classify(prot: Atom3D, lig: Atom3D, d: float) -> Optional[str]:
    """거리 + 원소 + 잔기 유형 기반 상호작용 분류"""
    pe = prot.element.upper()[:1]
    le = lig.element.upper()[:1]
    rn = prot.res_name.upper()

    # H-Bond: N/O ↔ N/O ≤ 3.5Å
    if d <= 3.5 and pe in ("N","O") and le in ("N","O"):
        return "H-Bond"

    # Salt Bridge: 하전 잔기 ↔ 반대극성 N/O ≤ 4.0Å
    if d <= 4.0:
        if rn in CHARGED_POS_RES and le == "O":
            return "Salt Bridge"
        if rn in CHARGED_NEG_RES and le == "N":
            return "Salt Bridge"

    # Pi-Stacking: 방향족 잔기 ↔ 리간드 C ≤ 5.5Å
    if d <= 5.5 and rn in AROMATIC_RES and pe == "C" and le == "C":
        return "Pi-Stacking"

    # Hydrophobic: 비극성 잔기 C↔C ≤ 5.0Å
    if d <= 5.0 and pe == "C" and le == "C" and rn in HYDROPHOBIC_RES:
        return "Hydrophobic"

    # Van der Waals: 기타 ≤ 4.0Å
    if d <= 4.0:
        return "Van der Waals"

    return None


def detect_contacts(
    prot_atoms: list[Atom3D],
    lig_atoms:  list[Atom3D],
    cutoff:     float = 5.5,
) -> list[ContactResidue]:
    """단백질 ↔ 리간드 접촉 잔기 목록 반환"""
    residues: dict[str, ContactResidue] = {}

    for la in lig_atoms:
        for pa in prot_atoms:
            d = _dist(la, pa)
            if d > cutoff:
                continue
            itype = _classify(pa, la, d)
            if itype is None:
                continue

            rid = pa.res_id
            if rid not in residues:
                residues[rid] = ContactResidue(
                    res_id    = rid,
                    res_name  = pa.res_name,
                    res_seq   = pa.res_seq,
                    chain_id  = pa.chain_id,
                    interaction = itype,
                    min_dist  = d,
                )
            else:
                cr = residues[rid]
                # 더 강한 상호작용으로 업데이트
                priority = list(INTERACTION_COLORS.keys())
                if priority.index(itype) < priority.index(cr.interaction):
                    cr.interaction = itype
                cr.min_dist = min(cr.min_dist, d)

            residues[rid].lig_atoms.append(la.atom_name)
            residues[rid].prot_atoms.append(pa.atom_name)

    return sorted(residues.values(), key=lambda x: (
        list(INTERACTION_COLORS.keys()).index(x.interaction), x.min_dist
    ))


def contacts_to_dataframe(contacts: list[ContactResidue]) -> pd.DataFrame:
    rows = []
    for c in contacts:
        rows.append({
            "residue_id":   c.res_id,
            "res_name":     c.res_name,
            "res_seq":      c.res_seq,
            "chain":        c.chain_id,
            "interaction":  c.interaction,
            "min_dist_A":   round(c.min_dist, 2),
            "lig_atoms":    ", ".join(sorted(set(c.lig_atoms))[:4]),
            "prot_atoms":   ", ".join(sorted(set(c.prot_atoms))[:4]),
        })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# 3. ProLIF 강화 분석 (선택적)
# ─────────────────────────────────────────────────────────────────────────────
def run_prolif(receptor_pdb: Path, ligand_pdb: Path) -> Optional[pd.DataFrame]:
    """ProLIF 기반 정밀 분석 (MDAnalysis + ProLIF 설치 필요)"""
    if not PROLIF_OK:
        return None
    try:
        import MDAnalysis as mda
        import prolif

        u = mda.Universe(str(receptor_pdb))
        lig_u = mda.Universe(str(ligand_pdb))
        prot = u.select_atoms("protein")
        lig  = lig_u.select_atoms("all")

        fp = prolif.Fingerprint()
        fp.run_from_atoms(lig, prot)
        df = fp.to_dataframe()
        log.info(f"[ProLIF] 분석 완료: {len(df.columns)}개 상호작용")
        return df
    except Exception as e:
        log.warning(f"[ProLIF] 오류 (거리 분석으로 폴백): {e}")
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 4. 시각화
# ─────────────────────────────────────────────────────────────────────────────
def draw_interaction_wheel(
    contacts:    list[ContactResidue],
    ligand_name: str,
    prefix:      str = "analysis",
    filename:    str = None,
    dpi:         int = 200,
) -> Optional[Path]:
    """
    리간드 중심 Wheel 다이어그램
    - 중앙: 리간드 이름
    - 주변: 접촉 잔기 (색상 = 상호작용 유형)
    - 연결선: 상호작용 유형별 색상 + 거리 표기
    """
    if not contacts:
        log.warning("[Interaction] 접촉 잔기 없음 — 다이어그램 생략")
        return None

    n = len(contacts)
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect("equal")
    ax.axis("off")

    # ── 중앙 리간드 박스 ─────────────────────────────────────────
    lig_box = mpatches.FancyBboxPatch(
        (-0.18, -0.08), 0.36, 0.16,
        boxstyle="round,pad=0.02",
        facecolor="#E3F2FD", edgecolor="#1565C0", linewidth=2,
    )
    ax.add_patch(lig_box)
    ax.text(0, 0, ligand_name, ha="center", va="center",
            fontsize=10, fontweight="bold", color="#1565C0",
            wrap=True)

    # ── 잔기 배치 ────────────────────────────────────────────────
    radius = 0.55
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)

    for i, (contact, angle) in enumerate(zip(contacts, angles)):
        rx = radius * math.cos(angle)
        ry = radius * math.sin(angle)
        color = INTERACTION_COLORS.get(contact.interaction, "#90A4AE")

        # 연결선
        # 리간드 박스 테두리 방향점 (근사)
        lx = 0.18 * math.cos(angle) if abs(math.cos(angle)) > abs(math.sin(angle)/2.25) \
             else 0.0
        ly = 0.08 * math.sin(angle) / (abs(math.sin(angle)) or 1)
        lx = np.clip(0.18 * math.cos(angle), -0.18, 0.18)
        ly = np.clip(0.08 * math.sin(angle), -0.08, 0.08)

        ax.annotate(
            "", xy=(rx, ry),
            xytext=(lx, ly),
            arrowprops=dict(
                arrowstyle="-",
                color=color, lw=1.8,
                connectionstyle="arc3,rad=0.05",
            ),
        )
        # 거리 표기 (선 중간)
        mx, my = (lx+rx)/2, (ly+ry)/2
        ax.text(mx, my, f"{contact.min_dist:.1f}Å",
                fontsize=6, color=color, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                           edgecolor="none", alpha=0.8))

        # 잔기 원
        circle = plt.Circle((rx, ry), 0.09, color=color, alpha=0.2, zorder=3)
        ax.add_patch(circle)
        circle2 = plt.Circle((rx, ry), 0.09, fill=False, edgecolor=color, lw=1.5, zorder=4)
        ax.add_patch(circle2)

        # 잔기 라벨
        res_label = f"{contact.res_name}\n{contact.res_seq}{contact.chain_id}"
        ax.text(rx, ry, res_label, ha="center", va="center",
                fontsize=7.5, fontweight="bold", zorder=5)

    # ── 범례 ─────────────────────────────────────────────────────
    legend_patches = [
        mpatches.Patch(color=c, label=itype, alpha=0.8)
        for itype, c in INTERACTION_COLORS.items()
        if any(cont.interaction == itype for cont in contacts)
    ]
    ax.legend(handles=legend_patches, loc="lower center",
              bbox_to_anchor=(0.5, -0.05), ncol=len(legend_patches),
              fontsize=8, framealpha=0.9)

    ax.set_xlim(-0.8, 0.8)
    ax.set_ylim(-0.8, 0.8)
    ax.set_title(f"Protein–Ligand Interaction: {ligand_name}",
                 fontsize=12, fontweight="bold", pad=14)

    fname = filename or f"{prefix}_interaction_wheel.png"
    path  = RESULTS_DIR / fname
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"[Interaction] Wheel 다이어그램 → {path}")
    return path


def draw_ligand_2d(
    smiles:      str,
    contacts:    list[ContactResidue],
    ligand_name: str,
    prefix:      str = "analysis",
    filename:    str = None,
    dpi:         int = 200,
) -> Optional[Path]:
    """rdkit으로 2D 리간드 구조 + 접촉 원자 하이라이팅"""
    if not RDKIT_OK or not smiles:
        return None
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        rdDepictor.Compute2DCoords(mol)

        # 간단한 2D 이미지 (접촉 원자 하이라이팅은 atom_name 매핑 없이 전체 표시)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()

        fname = filename or f"{prefix}_ligand_2d.svg"
        path  = RESULTS_DIR / fname
        path.write_text(svg)
        log.info(f"[Interaction] 리간드 2D SVG → {path}")
        return path
    except Exception as e:
        log.warning(f"[Interaction] 2D 리간드 생성 실패: {e}")
        return None


def draw_interaction_barplot(
    contacts: list[ContactResidue],
    prefix:   str = "analysis",
    filename: str = None,
    dpi:      int = 200,
) -> Optional[Path]:
    """상호작용 유형별 잔기 수 막대 + 세부 잔기 목록"""
    if not contacts:
        return None

    df = contacts_to_dataframe(contacts)
    counts = df["interaction"].value_counts()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, max(4, len(df)*0.35+2)))

    # ── 왼쪽: 유형별 막대 ────────────────────────────────────────
    colors = [INTERACTION_COLORS.get(t, "#90A4AE") for t in counts.index]
    ax1.barh(counts.index[::-1], counts.values[::-1], color=colors[::-1], alpha=0.85)
    ax1.set_xlabel("잔기 수", fontsize=10)
    ax1.set_title("Interaction Type Summary", fontsize=11, fontweight="bold")
    for i, v in enumerate(counts.values[::-1]):
        ax1.text(v + 0.05, i, str(v), va="center", fontsize=9)

    # ── 오른쪽: 잔기 목록 ────────────────────────────────────────
    ax2.axis("off")
    table_data = df[["residue_id","interaction","min_dist_A"]].head(20).values.tolist()
    col_labels  = ["Residue","Interaction","Dist (Å)"]
    tbl = ax2.table(
        cellText   = table_data,
        colLabels  = col_labels,
        loc        = "center",
        cellLoc    = "center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8)
    tbl.scale(1, 1.3)
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor("#1565C0")
            cell.set_text_props(color="white", fontweight="bold")
        else:
            itype = table_data[r-1][1] if r-1 < len(table_data) else ""
            cell.set_facecolor(INTERACTION_COLORS.get(itype, "#ECEFF1") + "44")
    ax2.set_title("Interacting Residues (Top 20)", fontsize=11, fontweight="bold", pad=12)

    plt.tight_layout()
    fname = filename or f"{prefix}_interaction_bar.png"
    path  = RESULTS_DIR / fname
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log.info(f"[Interaction] Bar plot → {path}")
    return path


# ─────────────────────────────────────────────────────────────────────────────
# 5. 통합 분석 실행
# ─────────────────────────────────────────────────────────────────────────────
def analyze_docking_result(
    receptor_pdb:    Path,
    vina_out_pdbqt:  Path,
    smiles:          str   = "",
    ligand_name:     str   = "Ligand",
    prefix:          str   = "analysis",
    work_dir:        Path  = None,
) -> dict:
    """
    단일 도킹 결과 → 상호작용 분석 + 시각화 일괄 처리
    반환: {contacts_df, wheel_png, bar_png, ligand_svg, prolif_df}
    """
    if work_dir is None:
        work_dir = vina_out_pdbqt.parent

    # 리간드 첫 번째 포즈 추출
    lig_pdb = extract_first_pose_pdbqt(vina_out_pdbqt, work_dir)
    if lig_pdb is None:
        return {}

    # 원자 파싱
    prot_atoms = parse_pdb_protein(receptor_pdb)
    lig_atoms  = parse_ligand_pdb(lig_pdb)

    if not prot_atoms or not lig_atoms:
        log.warning(f"[Interaction] 원자 파싱 실패: prot={len(prot_atoms)}, lig={len(lig_atoms)}")
        return {}

    log.info(f"[Interaction] {ligand_name} vs {receptor_pdb.stem}: "
             f"prot={len(prot_atoms)} atoms, lig={len(lig_atoms)} atoms")

    # ProLIF (선택적)
    prolif_df = run_prolif(receptor_pdb, lig_pdb)

    # 거리 기반 접촉
    contacts = detect_contacts(prot_atoms, lig_atoms)
    log.info(f"[Interaction] 접촉 잔기: {len(contacts)}개")

    contacts_df = contacts_to_dataframe(contacts)

    lbl = f"{prefix}_{ligand_name[:20]}"
    wheel_png  = draw_interaction_wheel(contacts, ligand_name, prefix=lbl)
    bar_png    = draw_interaction_barplot(contacts, prefix=lbl)
    ligand_svg = draw_ligand_2d(smiles, contacts, ligand_name, prefix=lbl) if smiles else None

    return {
        "contacts_df": contacts_df,
        "wheel_png":   wheel_png,
        "bar_png":     bar_png,
        "ligand_svg":  ligand_svg,
        "prolif_df":   prolif_df,
        "ligand_name": ligand_name,
        "n_contacts":  len(contacts),
    }


def batch_analyze(
    docking_df:   pd.DataFrame,
    receptor_pdbs: dict,      # {gene_symbol: receptor_pdb_path}
    ligand_smiles: dict,      # {compound_name: smiles}
    docking_dir:   Path,
    prefix:        str = "analysis",
) -> list[dict]:
    """
    BatchDocking 결과 DataFrame 전체 일괄 분석
    docking_df: compound, target, affinity, out_pdbqt
    """
    results = []
    for _, row in docking_df[docking_df["success"] == True].iterrows():
        compound = str(row.get("compound",""))
        target   = str(row.get("target",""))
        pdbqt    = row.get("out_pdbqt","")

        if not pdbqt or not Path(pdbqt).exists():
            continue
        rec_pdb = receptor_pdbs.get(target)
        if rec_pdb is None:
            # 수용체 PDB 경로 추정
            rec_pdb = docking_dir / "receptors" / f"{target}.pdb"
            if not rec_pdb.exists():
                rec_pdb = docking_dir / "receptors" / f"AF_{target}.pdb"
            if not rec_pdb.exists():
                log.warning(f"[Interaction] 수용체 PDB 없음: {target}")
                continue

        smiles = ligand_smiles.get(compound, "")
        r = analyze_docking_result(
            receptor_pdb   = rec_pdb,
            vina_out_pdbqt = Path(pdbqt),
            smiles         = smiles,
            ligand_name    = compound,
            prefix         = f"{prefix}_{compound[:15]}_vs_{target}",
            work_dir       = Path(pdbqt).parent,
        )
        if r:
            r["compound"] = compound
            r["target"]   = target
            r["affinity"] = row.get("affinity")
            results.append(r)

    log.info(f"[Interaction] 배치 분석 완료: {len(results)}개 쌍")
    return results
