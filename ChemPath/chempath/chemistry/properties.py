"""
Lightweight molecular property estimation from SMILES strings.
No RDKit required — uses string pattern analysis.

These are approximations useful for risk scoring, not exact calculations.
For precise values, use RDKit when available.
"""

import math
import re
from dataclasses import dataclass


@dataclass
class MolecularProperties:
    smiles: str
    estimated_mw: float
    estimated_hbd: int  # hydrogen bond donors
    estimated_hba: int  # hydrogen bond acceptors
    estimated_rotatable: int
    estimated_logp: float
    lipinski_violations: int
    risk_score: float  # 0.0 (safe) to 1.0 (risky)
    risk_flags: list[str]


# Average atomic weights for common SMILES atoms
_ATOM_WEIGHTS = {
    "C": 12.0, "N": 14.0, "O": 16.0, "S": 32.1, "F": 19.0,
    "Cl": 35.5, "Br": 79.9, "I": 126.9, "P": 31.0, "B": 10.8,
    "Si": 28.1, "Se": 79.0,
}


def estimate_molecular_weight(smiles: str) -> float:
    """Estimate molecular weight from SMILES atom counts."""
    weight = 0.0
    i = 0
    atom_count = 0
    while i < len(smiles):
        ch = smiles[i]
        # Two-letter atoms
        if ch in ("C", "B", "S") and i + 1 < len(smiles) and smiles[i + 1] == "l":
            weight += _ATOM_WEIGHTS.get("Cl", 35.5)
            atom_count += 1
            i += 2
            continue
        if ch == "B" and i + 1 < len(smiles) and smiles[i + 1] == "r":
            weight += _ATOM_WEIGHTS.get("Br", 79.9)
            atom_count += 1
            i += 2
            continue
        if ch == "S" and i + 1 < len(smiles) and smiles[i + 1] == "e":
            weight += _ATOM_WEIGHTS.get("Se", 79.0)
            atom_count += 1
            i += 2
            continue
        if ch == "S" and i + 1 < len(smiles) and smiles[i + 1] == "i":
            weight += _ATOM_WEIGHTS.get("Si", 28.1)
            atom_count += 1
            i += 2
            continue
        # Single-letter atoms (uppercase = atom)
        if ch.isupper() and ch in _ATOM_WEIGHTS:
            weight += _ATOM_WEIGHTS[ch]
            atom_count += 1
        i += 1

    # Add implicit hydrogens (rough: ~1H per C, 1-2H per N, 1H per O)
    c_count = smiles.count("C") - smiles.count("Cl")
    # Very rough hydrogen estimate
    h_estimate = max(0, atom_count * 0.7)
    weight += h_estimate * 1.008

    return round(weight, 1)


def estimate_hbd(smiles: str) -> int:
    """Estimate hydrogen bond donors (OH, NH, NH2 groups)."""
    # Count N and O not in rings or double bonds (very rough)
    count = 0
    count += len(re.findall(r"[NH]", smiles))  # NH in brackets
    count += len(re.findall(r"(?<![a-z])O(?![a-z])", smiles))  # O not in aromatic
    return min(count, 15)


def estimate_hba(smiles: str) -> int:
    """Estimate hydrogen bond acceptors (N, O atoms)."""
    n_count = smiles.upper().count("N")
    o_count = smiles.upper().count("O")
    return n_count + o_count


def estimate_rotatable_bonds(smiles: str) -> int:
    """Estimate rotatable bonds from single bonds between non-ring heavy atoms."""
    # Count single bonds (not in rings) — rough: non-aromatic non-ring bonds
    single_bonds = smiles.count("-") + smiles.count("C") // 3
    return max(0, min(single_bonds, 30))


def estimate_logp(smiles: str) -> float:
    """Very rough LogP estimate from atom composition."""
    c_count = smiles.count("C") + smiles.count("c")
    n_count = smiles.count("N") + smiles.count("n")
    o_count = smiles.count("O") + smiles.count("o")
    f_count = smiles.count("F")
    cl_count = smiles.count("Cl")

    # Wildman-Crippen-like fragment approach (very simplified)
    logp = c_count * 0.12 - n_count * 0.7 - o_count * 0.5 + f_count * 0.4 + cl_count * 0.6
    return round(logp, 2)


def compute_properties(smiles: str) -> MolecularProperties:
    """Compute all estimated molecular properties and risk score."""
    mw = estimate_molecular_weight(smiles)
    hbd = estimate_hbd(smiles)
    hba = estimate_hba(smiles)
    rotatable = estimate_rotatable_bonds(smiles)
    logp = estimate_logp(smiles)

    # Lipinski Rule of 5 violations
    violations = 0
    risk_flags = []

    if mw > 500:
        violations += 1
        risk_flags.append(f"MW={mw:.0f} > 500")
    if hbd > 5:
        violations += 1
        risk_flags.append(f"HBD={hbd} > 5")
    if hba > 10:
        violations += 1
        risk_flags.append(f"HBA={hba} > 10")
    if logp > 5:
        violations += 1
        risk_flags.append(f"LogP={logp:.1f} > 5")

    # Additional risk factors
    if mw > 800:
        risk_flags.append("Very high MW (>800) — poor oral bioavailability")
    if rotatable > 10:
        risk_flags.append(f"High flexibility ({rotatable} rotatable bonds)")
    if len(smiles) > 200:
        risk_flags.append("Very large molecule")

    # Risk score: 0 to 1
    risk_score = min(1.0, violations * 0.2 + (mw > 800) * 0.15 + (rotatable > 10) * 0.1)

    return MolecularProperties(
        smiles=smiles,
        estimated_mw=mw,
        estimated_hbd=hbd,
        estimated_hba=hba,
        estimated_rotatable=rotatable,
        estimated_logp=logp,
        lipinski_violations=violations,
        risk_score=round(risk_score, 2),
        risk_flags=risk_flags,
    )


def estimate_tpsa(smiles: str) -> float:
    """Estimate topological polar surface area from N and O atom counts."""
    n_count = smiles.upper().count("N")
    o_count = smiles.upper().count("O")
    s_count = smiles.upper().count("S") - smiles.count("Si") - smiles.count("Se")
    # Rough: each N contributes ~26A^2, each O ~20A^2, each S ~25A^2
    return round(n_count * 26.0 + o_count * 20.0 + max(0, s_count) * 25.0, 1)


def count_aromatic_rings(smiles: str) -> int:
    """Count aromatic rings from lowercase atoms in SMILES."""
    aromatic_atoms = sum(1 for ch in smiles if ch.islower() and ch in "cnos")
    return max(0, aromatic_atoms // 5)  # ~5 aromatic atoms per ring


def estimate_qed(smiles: str) -> float:
    """
    Estimate Quantitative Estimate of Drug-likeness (QED).
    QED ranges 0-1, higher = more drug-like.
    Based on desirability functions for MW, LogP, HBA, HBD, TPSA, RotBonds, AromaticRings.
    """
    mw = estimate_molecular_weight(smiles)
    logp = estimate_logp(smiles)
    hba = estimate_hba(smiles)
    hbd = estimate_hbd(smiles)
    tpsa = estimate_tpsa(smiles)
    rotatable = estimate_rotatable_bonds(smiles)
    arom_rings = count_aromatic_rings(smiles)

    # Desirability functions (sigmoid-like, penalize extremes)
    def _desirability(value, optimal, width):
        """Gaussian desirability centered at optimal."""
        return math.exp(-0.5 * ((value - optimal) / width) ** 2)

    d_mw = _desirability(mw, 300, 150)
    d_logp = _desirability(logp, 2.5, 2.0)
    d_hba = _desirability(hba, 4, 3)
    d_hbd = _desirability(hbd, 2, 2)
    d_tpsa = _desirability(tpsa, 75, 50)
    d_rot = _desirability(rotatable, 4, 4)
    d_arom = _desirability(arom_rings, 2, 2)

    # Geometric mean of desirability scores
    scores = [d_mw, d_logp, d_hba, d_hbd, d_tpsa, d_rot, d_arom]
    product = 1.0
    for s in scores:
        product *= max(s, 0.001)
    qed = product ** (1.0 / len(scores))
    return round(min(1.0, max(0.0, qed)), 3)


def estimate_synthetic_accessibility(smiles: str) -> float:
    """
    Estimate synthetic accessibility score (1=easy, 10=hard).
    Based on molecular complexity heuristics: ring count, stereocenters, heteroatom diversity.
    """
    ring_count = smiles.count("1") + smiles.count("2") + smiles.count("3") + smiles.count("4")
    stereo = smiles.count("@") + smiles.count("/") + smiles.count("\\")
    length = len(smiles)
    hetero = sum(1 for ch in smiles.upper() if ch in "NOSPF")

    # Complexity score
    complexity = (
        ring_count * 0.5
        + stereo * 0.8
        + length * 0.02
        + hetero * 0.15
    )
    sa_score = min(10.0, max(1.0, 1.0 + complexity * 0.8))
    return round(sa_score, 1)


@dataclass
class EnrichedProperties(MolecularProperties):
    """Extended properties for multi-dimensional screening."""
    estimated_tpsa: float = 0.0
    aromatic_rings: int = 0
    qed_score: float = 0.0
    sa_score: float = 1.0  # synthetic accessibility


def compute_enriched_properties(smiles: str) -> EnrichedProperties:
    """Compute all molecular properties including QED, TPSA, SA for multi-objective optimization."""
    base = compute_properties(smiles)
    tpsa = estimate_tpsa(smiles)
    arom = count_aromatic_rings(smiles)
    qed = estimate_qed(smiles)
    sa = estimate_synthetic_accessibility(smiles)

    return EnrichedProperties(
        smiles=base.smiles,
        estimated_mw=base.estimated_mw,
        estimated_hbd=base.estimated_hbd,
        estimated_hba=base.estimated_hba,
        estimated_rotatable=base.estimated_rotatable,
        estimated_logp=base.estimated_logp,
        lipinski_violations=base.lipinski_violations,
        risk_score=base.risk_score,
        risk_flags=base.risk_flags,
        estimated_tpsa=tpsa,
        aromatic_rings=arom,
        qed_score=qed,
        sa_score=sa,
    )


def compute_risk_scores(compounds: list[dict], smiles_key: str = "smiles") -> dict[str, float]:
    """Compute risk scores for a list of compounds. Returns {chembl_id: risk_score}."""
    scores = {}
    for c in compounds:
        smi = c.get(smiles_key, "")
        if not smi:
            scores[c["chembl_id"]] = 0.5  # unknown
            continue
        props = compute_properties(smi)
        scores[c["chembl_id"]] = props.risk_score
    return scores
