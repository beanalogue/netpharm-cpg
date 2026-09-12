"""
자동 분자 도킹 모듈 — AutoDock Vina 래퍼
핵심 성분(Ligand) × 허브 단백질(Receptor) 자동 배치 실행

설치:
  wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
  chmod +x vina_1.2.5_linux_x86_64 && mv vina_1.2.5_linux_x86_64 ~/bin/vina

의존 도구:
  - vina       : 도킹 실행
  - obabel     : 포맷 변환 (SMILES → PDBQT)
  - curl/wget  : PDB 수동 다운로드

PDB 수신 전략:
  AlphaFold DB 또는 RCSB PDB API에서 단백질 구조 자동 다운로드
"""

import logging
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import requests
import pandas as pd

BASE_DIR    = Path(__file__).parent.parent
DOCKING_DIR = BASE_DIR / "docking"
DOCKING_DIR.mkdir(exist_ok=True)

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# 도구 확인
# ─────────────────────────────────────────────────────────────────────────────
def check_tools() -> dict:
    status = {}
    for tool in ["vina", "obabel"]:
        status[tool] = shutil.which(tool) is not None
    return status


# ─────────────────────────────────────────────────────────────────────────────
# 단백질 구조 다운로더
# ─────────────────────────────────────────────────────────────────────────────
class StructureDownloader:

    from collect import make_retry_session as _mks
    session = _mks()

    @classmethod
    def get_pdb_id(cls, uniprot_id: str) -> str | None:
        """UniProt ID → 대표 PDB ID 조회 (UniProt cross-reference API)"""
        if not uniprot_id:
            return None
        try:
            r = cls.session.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot_id}",
                params={"fields": "xref_pdb", "format": "json"},
                timeout=15,
            )
            dbs = r.json().get("uniProtKBCrossReferences", [])
            pdb_ids = [x["id"] for x in dbs if x.get("database") == "PDB"]
            return pdb_ids[0] if pdb_ids else None
        except Exception as e:
            log.debug(f"[RCSB] {uniprot_id}: {e}")
            return None

    @classmethod
    def download_pdb(cls, pdb_id: str, out_dir: Path) -> Path | None:
        """RCSB에서 PDB 파일 다운로드"""
        path = out_dir / f"{pdb_id.upper()}.pdb"
        if path.exists():
            return path
        try:
            r = cls.session.get(
                f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb",
                timeout=30,
            )
            r.raise_for_status()
            path.write_bytes(r.content)
            log.info(f"[PDB] {pdb_id} 다운로드 완료")
            return path
        except Exception as e:
            log.error(f"[PDB] {pdb_id} 다운로드 실패: {e}")
            return None

    @classmethod
    def download_alphafold(cls, uniprot_id: str, out_dir: Path) -> Path | None:
        """AlphaFold DB에서 구조 다운로드 (실험 구조 없을 때 대안) — 최신 버전 자동 조회"""
        path = out_dir / f"AF_{uniprot_id}.pdb"
        if path.exists():
            return path
        try:
            # 최신 버전 URL을 AlphaFold API에서 조회
            meta = cls.session.get(
                f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}",
                timeout=15,
            )
            meta.raise_for_status()
            pdb_url = meta.json()[0].get("pdbUrl", "")
            if not pdb_url:
                return None
            r = cls.session.get(pdb_url, timeout=60)
            r.raise_for_status()
            path.write_bytes(r.content)
            log.info(f"[AlphaFold] {uniprot_id} 구조 다운로드 완료")
            return path
        except Exception as e:
            log.error(f"[AlphaFold] {uniprot_id}: {e}")
            return None

    @classmethod
    def gene_to_uniprot(cls, gene_symbol: str) -> str | None:
        """Gene symbol → UniProt canonical accession (human, reviewed)"""
        if not gene_symbol:
            return None
        try:
            r = cls.session.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
                    "fields": "accession",
                    "format": "json",
                    "size": "1",
                },
                timeout=10,
            )
            results = r.json().get("results", [])
            if results:
                return results[0]["primaryAccession"]
        except Exception as e:
            log.debug(f"[UniProt] {gene_symbol}: {e}")
        return None

    @classmethod
    def get_structure(cls, gene_symbol: str, uniprot_id: str, out_dir: Path) -> Path | None:
        """PDB → AlphaFold 순으로 단백질 구조 확보 (uniprot_id 없으면 gene symbol로 조회)"""
        if not uniprot_id and gene_symbol:
            uniprot_id = cls.gene_to_uniprot(gene_symbol) or ""
            if uniprot_id:
                log.info(f"[UniProt] {gene_symbol} → {uniprot_id}")
        pdb_id = cls.get_pdb_id(uniprot_id)
        if pdb_id:
            path = cls.download_pdb(pdb_id, out_dir)
            if path:
                return path
        if uniprot_id:
            return cls.download_alphafold(uniprot_id, out_dir)
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 리간드 준비 (SMILES → PDBQT)
# ─────────────────────────────────────────────────────────────────────────────
class LigandPrep:

    @staticmethod
    def smiles_to_pdbqt(smiles: str, name: str, out_dir: Path) -> Path | None:
        """obabel로 SMILES → PDBQT 변환 (3D 생성 포함)"""
        if not shutil.which("obabel"):
            log.error("[도킹] obabel 미설치. 'sudo apt install openbabel' 필요")
            return None

        path = out_dir / f"{name}.pdbqt"
        try:
            result = subprocess.run(
                ["obabel", f"-:{smiles}", "-O", str(path),
                 "--gen3d", "-h", "--partialcharge", "gasteiger"],
                capture_output=True, text=True, timeout=60,
            )
            if path.exists() and path.stat().st_size > 0:
                log.info(f"[Ligand] {name} PDBQT 생성 완료")
                return path
            log.error(f"[Ligand] {name} 변환 실패: {result.stderr[:200]}")
            return None
        except Exception as e:
            log.error(f"[Ligand] {name}: {e}")
            return None


# ─────────────────────────────────────────────────────────────────────────────
# 수용체 준비 (PDB → PDBQT)
# ─────────────────────────────────────────────────────────────────────────────
class ReceptorPrep:

    @staticmethod
    def clean_pdb(pdb_path: Path, out_dir: Path, chain: str = None) -> Path:
        """
        물 분자(HOH/WAT) 및 소분자 HETATM 제거
        ATOM 레코드(단백질 주쇄)만 남김
        chain: 특정 체인만 추출 (예: 'A') — 대형 복합체 경량화
        """
        out_path = out_dir / (pdb_path.stem + "_clean.pdb")
        kept = []
        for line in pdb_path.read_text(errors="replace").splitlines():
            record = line[:6].strip()
            if record == "HETATM":
                continue  # 물·리간드·이온 전부 제거
            if record == "ATOM":
                if chain and line[21:22].strip() != chain:
                    continue  # 지정 체인 외 제거
            kept.append(line)
        out_path.write_text("\n".join(kept))
        size_kb = out_path.stat().st_size // 1024
        log.info(f"[Receptor] 정제 완료 → {out_path.name} ({size_kb} KB)")
        return out_path

    @staticmethod
    def add_hydrogens_pdb(pdb_path: Path, out_dir: Path) -> Path:
        """obabel로 극성 수소 추가"""
        if not shutil.which("obabel"):
            return pdb_path
        out_path = out_dir / (pdb_path.stem + "_h.pdb")
        size_kb = pdb_path.stat().st_size // 1024
        timeout = max(120, size_kb // 5)  # 파일 크기 비례 타임아웃 (최소 120초)
        try:
            subprocess.run(
                ["obabel", str(pdb_path), "-O", str(out_path), "-h"],
                capture_output=True, text=True, timeout=timeout,
            )
            return out_path if out_path.exists() and out_path.stat().st_size > 0 else pdb_path
        except Exception as e:
            log.warning(f"[Receptor] 수소 추가 실패 ({e}), 원본 사용")
            return pdb_path

    @staticmethod
    def pdb_to_pdbqt(pdb_path: Path, out_dir: Path, chain: str = None) -> Path | None:
        """
        PDB → PDBQT 3단계 파이프라인
          1) 물·HETATM 제거 + 체인 필터  (clean_pdb)
          2) 극성 수소 추가               (add_hydrogens_pdb)
          3) PDBQT 변환                   (obabel --partialcharge gasteiger)
        chain: 특정 체인만 사용 (대형 복합체는 'A' 권장)
        """
        if not shutil.which("obabel"):
            return None

        clean  = ReceptorPrep.clean_pdb(pdb_path, out_dir, chain=chain)
        with_h = ReceptorPrep.add_hydrogens_pdb(clean, out_dir)

        out_path = out_dir / (pdb_path.stem + "_receptor.pdbqt")
        size_kb  = with_h.stat().st_size // 1024
        timeout  = max(180, size_kb // 3)  # 크기 비례 타임아웃
        try:
            r = subprocess.run(
                ["obabel", str(with_h), "-O", str(out_path),
                 "--partialcharge", "gasteiger",
                 "-xr"],   # -xr: rigid receptor (no torsion tree) — Vina 필수
                capture_output=True, text=True, timeout=timeout,
            )
            if out_path.exists() and out_path.stat().st_size > 0:
                log.info(f"[Receptor] PDBQT 생성 → {out_path.name} ({out_path.stat().st_size//1024} KB)")
                return out_path
            log.error(f"[Receptor] PDBQT 변환 실패: {r.stderr[:300]}")
            return None
        except Exception as e:
            log.error(f"[Receptor] {pdb_path.name}: {e}")
            return None

    # 물·이온 등 non-ligand HETATM 제외 목록
    _SKIP_HETATM = {
        "HOH","WAT","DOD","H2O",           # 물
        "SO4","PO4","GOL","EDO","PEG",      # 결정화 시약
        "MG","ZN","CA","NA","FE","MN","CU","CO","NI","CD","HG","PT",  # 금속
        "CL","BR","IOD","F",               # 할로겐
        "ACT","ACE","IMD","DMS","MPD",     # 용매/additive
    }

    @classmethod
    def get_box_center(cls, pdb_path: Path) -> tuple[float, float, float]:
        """
        Binding box 중심 결정 (우선순위):
          1) PDB 내 공결정 리간드 HETATM 중심  (실험적 결합 부위)
          2) Cα 무게중심 폴백                  (구조 전체)
        """
        # 1) 공결정 리간드 탐색
        lig_x, lig_y, lig_z = [], [], []
        try:
            for line in pdb_path.read_text(errors="replace").splitlines():
                if not line.startswith("HETATM"):
                    continue
                resname = line[17:20].strip().upper()
                if resname in cls._SKIP_HETATM:
                    continue
                try:
                    lig_x.append(float(line[30:38]))
                    lig_y.append(float(line[38:46]))
                    lig_z.append(float(line[46:54]))
                except ValueError:
                    pass
        except Exception:
            pass

        if lig_x:
            cx = round(sum(lig_x) / len(lig_x), 2)
            cy = round(sum(lig_y) / len(lig_y), 2)
            cz = round(sum(lig_z) / len(lig_z), 2)
            log.info(f"[Box] 공결정 리간드 중심 사용: ({cx}, {cy}, {cz})")
            return (cx, cy, cz)

        # 2) Cα 무게중심 폴백
        xs, ys, zs = [], [], []
        try:
            for line in pdb_path.read_text(errors="replace").splitlines():
                if not line.startswith("ATOM"):
                    continue
                if line[12:16].strip() != "CA":
                    continue
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
        except Exception:
            pass
        if not xs:
            try:
                for line in pdb_path.read_text(errors="replace").splitlines():
                    if line.startswith("ATOM"):
                        xs.append(float(line[30:38]))
                        ys.append(float(line[38:46]))
                        zs.append(float(line[46:54]))
            except Exception:
                pass
        if not xs:
            return (0.0, 0.0, 0.0)
        log.info(f"[Box] Cα 무게중심 사용 (리간드 없음)")
        return (
            round(sum(xs) / len(xs), 2),
            round(sum(ys) / len(ys), 2),
            round(sum(zs) / len(zs), 2),
        )


# ─────────────────────────────────────────────────────────────────────────────
# Vina 도킹 실행
# ─────────────────────────────────────────────────────────────────────────────
class VinaDocking:

    def __init__(
        self,
        vina_bin:      str  = "vina",
        exhaustiveness: int = 8,
        num_modes:     int  = 5,
        energy_range:  int  = 3,
        box_size:      tuple = (20, 20, 20),
    ):
        self.vina_bin       = vina_bin
        self.exhaustiveness = exhaustiveness
        self.num_modes      = num_modes
        self.energy_range   = energy_range
        self.box_size       = box_size

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt:   Path,
        center:         tuple,
        out_dir:        Path,
        label:          str = "dock",
    ) -> dict:
        """
        단일 도킹 실행
        반환: {label, affinity (kcal/mol), out_pdbqt, success}
        """
        if not shutil.which(self.vina_bin):
            return {"label": label, "affinity": None, "success": False,
                    "error": "vina 미설치"}

        out_pdbqt = out_dir / f"{label}_out.pdbqt"

        # Vina 1.2.5 이상: --log 옵션 없음 → stdout 직접 캡처
        cmd = [
            self.vina_bin,
            "--receptor", str(receptor_pdbqt),
            "--ligand",   str(ligand_pdbqt),
            "--center_x", str(round(center[0], 3)),
            "--center_y", str(round(center[1], 3)),
            "--center_z", str(round(center[2], 3)),
            "--size_x",   str(self.box_size[0]),
            "--size_y",   str(self.box_size[1]),
            "--size_z",   str(self.box_size[2]),
            "--out",      str(out_pdbqt),
            "--exhaustiveness", str(self.exhaustiveness),
            "--num_modes",      str(self.num_modes),
            "--energy_range",   str(self.energy_range),
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600
            )
            stdout = proc.stdout + proc.stderr
            # 로그를 파일로 저장 (디버그용)
            (out_dir / f"{label}.log").write_text(stdout)

            affinity = self._parse_affinity_stdout(stdout)
            if proc.returncode != 0 and affinity is None:
                log.error(f"[Vina] {label} 오류: {stdout[-500:]}")
            else:
                log.info(f"[Vina] {label}: {affinity} kcal/mol")
            return {
                "label":    label,
                "affinity": affinity,
                "out_pdbqt": str(out_pdbqt) if out_pdbqt.exists() else None,
                "success":  affinity is not None,
                "stdout":   stdout,
            }
        except subprocess.TimeoutExpired:
            return {"label": label, "affinity": None, "success": False,
                    "error": "timeout"}
        except Exception as e:
            return {"label": label, "affinity": None, "success": False,
                    "error": str(e)}

    @staticmethod
    def _parse_affinity_stdout(stdout: str) -> float | None:
        """Vina 1.2.5 stdout에서 최저 결합 에너지 파싱
        출력 형식: '   1       -7.5      0.000      0.000'
        """
        for line in stdout.splitlines():
            stripped = line.strip()
            if stripped.startswith("1 ") or stripped.startswith("1\t"):
                parts = stripped.split()
                try:
                    return float(parts[1])
                except (IndexError, ValueError):
                    pass
        return None


# ─────────────────────────────────────────────────────────────────────────────
# 배치 도킹 파이프라인
# ─────────────────────────────────────────────────────────────────────────────
class BatchDocking:

    def __init__(self, max_ligands: int = 20, max_receptors: int = 5):
        self.max_ligands   = max_ligands
        self.max_receptors = max_receptors
        self.downloader    = StructureDownloader()
        self.vina          = VinaDocking()

    def run(
        self,
        compounds:   list,          # [{chembl_id, pref_name, smiles}, ...]
        hub_targets: list,          # [{label(gene), uniprot_id}, ...]
        run_dir:     Path = DOCKING_DIR,
        box_size:    tuple = (20, 20, 20),   # (x, y, z) Å — 수동 지정 가능
        box_center:  tuple | None = None,    # None → 수용체 Cα 무게중심 자동 계산
    ) -> pd.DataFrame:
        """
        화합물 × 허브 단백질 배치 도킹
        반환: DataFrame (compound, target, affinity)
        """
        tools = check_tools()
        if not tools["vina"]:
            log.error("[도킹] AutoDock Vina 미설치. 아래 명령어로 설치하세요:\n"
                      "  wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/"
                      "download/v1.2.5/vina_1.2.5_linux_x86_64\n"
                      "  chmod +x vina_1.2.5_linux_x86_64\n"
                      "  sudo mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina")
            return pd.DataFrame()

        if not tools["obabel"]:
            log.error("[도킹] OpenBabel 미설치: sudo apt install openbabel")
            return pd.DataFrame()

        run_dir.mkdir(exist_ok=True)
        lig_dir = run_dir / "ligands"
        rec_dir = run_dir / "receptors"
        out_dir = run_dir / "results"
        for d in [lig_dir, rec_dir, out_dir]:
            d.mkdir(exist_ok=True)

        # 리간드 준비
        ligands = []
        for comp in compounds[:self.max_ligands]:
            smiles    = comp.get("smiles", "")
            chembl_id = comp.get("chembl_id", "")
            pref_name = comp.get("pref_name", "")
            # 파일명용 식별자 (영문 안전)
            file_id   = chembl_id or pref_name or "lig"
            # 논문 표시용 이름: 성분명 우선, 없으면 CHEMBL ID
            display_name = pref_name or chembl_id or "lig"
            if not smiles:
                continue
            pdbqt = LigandPrep.smiles_to_pdbqt(smiles, file_id, lig_dir)
            if pdbqt:
                ligands.append({
                    "name":         display_name,
                    "chembl_id":    chembl_id,
                    "pdbqt":        pdbqt,
                    "smiles":       smiles,
                })

        # 수용체 준비
        receptors = []
        for tgt in hub_targets[:self.max_receptors]:
            gene     = tgt.get("label") or tgt.get("gene_symbol", "")
            uniprot  = tgt.get("uniprot_id") or tgt.get("uniprot", "")
            pdb_path = self.downloader.get_structure(gene, uniprot, rec_dir)
            if not pdb_path:
                continue
            rec_pdbqt = ReceptorPrep.pdb_to_pdbqt(pdb_path, rec_dir)
            # box_center 수동 지정 없으면 수용체별 Cα 중심 자동 계산
            center = box_center if box_center else ReceptorPrep.get_box_center(pdb_path)
            log.info(f"[Receptor] {gene} grid center: {center}, size: {box_size}")
            if rec_pdbqt:
                receptors.append({"gene": gene, "pdbqt": rec_pdbqt,
                                   "center": center, "box_size": box_size})

        if not ligands or not receptors:
            log.warning(f"[도킹] 리간드 {len(ligands)}개 / 수용체 {len(receptors)}개 — 중단")
            return pd.DataFrame()

        log.info(f"[도킹] {len(ligands)}개 리간드 × {len(receptors)}개 수용체 = "
                 f"{len(ligands)*len(receptors)}회 도킹 예정")

        # VinaDocking에 box_size 반영
        self.vina.box_size = box_size

        rows = []
        for rec in receptors:
            for lig in ligands:
                label  = f"{lig['name']}_vs_{rec['gene']}"
                result = self.vina.dock(
                    rec["pdbqt"], lig["pdbqt"], rec["center"], out_dir, label
                )
                rows.append({
                    "compound":   lig["name"],
                    "chembl_id":  lig.get("chembl_id", ""),
                    "target":     rec["gene"],
                    "affinity":   result.get("affinity"),
                    "success":    result.get("success", False),
                    "out_pdbqt":  result.get("out_pdbqt", ""),
                })
                time.sleep(0.1)

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values("affinity")
            csv_path = run_dir / "docking_results.csv"
            df.to_csv(csv_path, index=False)
            log.info(f"[도킹] 결과 저장: {csv_path}")

        return df
