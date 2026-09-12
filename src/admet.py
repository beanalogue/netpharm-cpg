"""
ADMET 약물성 스크리닝 필터
Lipinski's Rule of 5 + OB/DL 기준으로 성분 필터링

기준 (네트워크약리학 논문 표준):
  OB  >= 30%        (oral bioavailability, TCMSP 기준)
  DL  >= 0.18       (drug-likeness, TCMSP 기준)
  MW  <= 500 Da     (Lipinski)
  LogP <= 5         (Lipinski)
  HBD <= 5          (H-bond donors, Lipinski)
  HBA <= 10         (H-bond acceptors, Lipinski)
  TPSA <= 140 Å²    (topological polar surface area, oral absorption)
  RotBonds <= 10    (rotatable bonds)
"""

import logging
from dataclasses import dataclass, field, asdict
from typing import Optional

log = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
    RDKIT_OK = True
except ImportError:
    RDKIT_OK = False
    log.warning("[ADMET] rdkit 미설치 — SMILES 기반 계산 불가")


# ─────────────────────────────────────────────────────────────────────────────
# 기본 기준값
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class ADMETCriteria:
    ob_min:       float = 30.0    # Oral Bioavailability ≥ 30% (TCMSP)
    dl_min:       float = 0.18    # Drug-likeness ≥ 0.18 (TCMSP)
    mw_max:       float = 500.0   # Molecular Weight ≤ 500 Da
    logp_max:     float = 5.0     # LogP ≤ 5
    hbd_max:      int   = 5       # H-bond donors ≤ 5
    hba_max:      int   = 10      # H-bond acceptors ≤ 10
    tpsa_max:     float = 140.0   # TPSA ≤ 140 Å²
    rotbonds_max: int   = 10      # Rotatable bonds ≤ 10


@dataclass
class ADMETResult:
    compound_id:  str
    compound_name: str
    smiles:       str = ""
    # TCMSP 기반
    ob:           Optional[float] = None
    dl:           Optional[float] = None
    # 계산값 (rdkit)
    mw:           Optional[float] = None
    logp:         Optional[float] = None
    hbd:          Optional[int]   = None
    hba:          Optional[int]   = None
    tpsa:         Optional[float] = None
    rotbonds:     Optional[int]   = None
    # 판정
    pass_filter:  bool  = False
    fail_reasons: list  = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


# ─────────────────────────────────────────────────────────────────────────────
# 계산 함수
# ─────────────────────────────────────────────────────────────────────────────
def calc_properties(smiles: str) -> dict:
    """SMILES → rdkit 물성 계산. rdkit 없으면 빈 dict 반환."""
    if not RDKIT_OK or not smiles:
        return {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw":       round(Descriptors.MolWt(mol), 2),
        "logp":     round(Crippen.MolLogP(mol), 3),
        "hbd":      rdMolDescriptors.CalcNumHBD(mol),
        "hba":      rdMolDescriptors.CalcNumHBA(mol),
        "tpsa":     round(Descriptors.TPSA(mol), 2),
        "rotbonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


# ─────────────────────────────────────────────────────────────────────────────
# 스크리너
# ─────────────────────────────────────────────────────────────────────────────
class ADMETScreener:

    def __init__(self, criteria: ADMETCriteria = None):
        self.crit = criteria or ADMETCriteria()

    def screen_compound(self, compound_data: dict) -> ADMETResult:
        """
        단일 화합물 ADMET 평가
        compound_data: ChEMBL/TCMSP 수집 결과의 molecule dict
        """
        mol = compound_data.get("molecule", compound_data)
        cid   = mol.get("chembl_id") or compound_data.get("query", "")
        cname = mol.get("pref_name")  or compound_data.get("query", cid)
        smiles = mol.get("smiles", "")

        result = ADMETResult(compound_id=cid, compound_name=cname, smiles=smiles)

        # OB / DL (TCMSP 기반, 있으면 사용)
        result.ob = mol.get("ob") or mol.get("OB") or mol.get("OB (%)")
        result.dl = mol.get("dl") or mol.get("DL")

        # rdkit 계산 (SMILES 있을 때)
        props = calc_properties(smiles)
        result.mw       = props.get("mw")       or mol.get("mol_weight") or mol.get("full_mwt")
        result.logp     = props.get("logp")
        result.hbd      = props.get("hbd")
        result.hba      = props.get("hba")
        result.tpsa     = props.get("tpsa")
        result.rotbonds = props.get("rotbonds")

        # ChEMBL에서 mol_weight가 문자열로 올 때 변환
        if isinstance(result.mw, str):
            try:
                result.mw = float(result.mw)
            except ValueError:
                result.mw = None

        # 판정
        reasons = []
        c = self.crit

        if result.ob is not None and result.ob < c.ob_min:
            reasons.append(f"OB {result.ob:.1f}% < {c.ob_min}%")
        if result.dl is not None and result.dl < c.dl_min:
            reasons.append(f"DL {result.dl:.2f} < {c.dl_min}")
        if result.mw is not None and result.mw > c.mw_max:
            reasons.append(f"MW {result.mw:.0f} > {c.mw_max}")
        if result.logp is not None and result.logp > c.logp_max:
            reasons.append(f"LogP {result.logp:.2f} > {c.logp_max}")
        if result.hbd is not None and result.hbd > c.hbd_max:
            reasons.append(f"HBD {result.hbd} > {c.hbd_max}")
        if result.hba is not None and result.hba > c.hba_max:
            reasons.append(f"HBA {result.hba} > {c.hba_max}")
        if result.tpsa is not None and result.tpsa > c.tpsa_max:
            reasons.append(f"TPSA {result.tpsa:.1f} > {c.tpsa_max}")
        if result.rotbonds is not None and result.rotbonds > c.rotbonds_max:
            reasons.append(f"RotBonds {result.rotbonds} > {c.rotbonds_max}")

        result.fail_reasons = reasons
        result.pass_filter  = len(reasons) == 0
        return result

    def screen_all(self, herb_results: dict) -> dict:
        """
        collect_herb_compounds_chembl() 결과 전체 스크리닝
        반환: {
            "passed":  [ADMETResult],  # 통과
            "failed":  [ADMETResult],  # 탈락
            "summary": DataFrame
        }
        """
        import pandas as pd

        passed, failed = [], []
        for herb_name, herb_data in herb_results.items():
            for comp in herb_data.get("compounds", []):
                r = self.screen_compound(comp)
                r.compound_name = f"{r.compound_name} [{herb_name}]"
                if r.pass_filter:
                    passed.append(r)
                else:
                    failed.append(r)

        log.info(f"[ADMET] 통과: {len(passed)}개 / 탈락: {len(failed)}개")

        rows = [r.to_dict() for r in passed + failed]
        df = pd.DataFrame(rows) if rows else pd.DataFrame()

        return {"passed": passed, "failed": failed, "summary": df}

    def filter_herb_results(self, herb_results: dict) -> dict:
        """
        herb_results에서 ADMET 통과 화합물만 남긴 새 dict 반환
        pipeline과 연계: collect → admet_filter → network
        """
        screen = self.screen_all(herb_results)
        passed_ids = {r.compound_id for r in screen["passed"]}

        filtered = {}
        for herb_name, herb_data in herb_results.items():
            kept = [
                c for c in herb_data.get("compounds", [])
                if (c.get("molecule", {}).get("chembl_id") or c.get("query", "")) in passed_ids
            ]
            filtered[herb_name] = {**herb_data, "compounds": kept}

        total_before = sum(len(h.get("compounds", [])) for h in herb_results.values())
        total_after  = sum(len(h.get("compounds", [])) for h in filtered.values())
        log.info(f"[ADMET] 필터 결과: {total_before}개 → {total_after}개 성분")
        return filtered
