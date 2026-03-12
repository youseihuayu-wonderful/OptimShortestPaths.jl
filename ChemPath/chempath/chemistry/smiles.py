"""
SMILES validation — prevents pipeline crashes from invalid molecular strings.

Three levels:
1. Basic syntax check (no dependencies)
2. RDKit validation (if rdkit installed)
3. Batch validation with summary report
"""

import re
from dataclasses import dataclass


@dataclass
class ValidationResult:
    smiles: str
    is_valid: bool
    level: str  # "basic" or "rdkit"
    error: str | None = None
    canonical_smiles: str | None = None
    mol_weight: float | None = None
    num_atoms: int | None = None


_SMILES_PATTERN = re.compile(
    r'^[A-Za-z0-9@+\-\[\]\(\)\\/=#$:.%~&!]+$'
)


def validate_smiles_basic(smiles: str) -> ValidationResult:
    """Basic SMILES validation without external dependencies."""
    if not smiles or not isinstance(smiles, str):
        return ValidationResult(smiles=str(smiles), is_valid=False, level="basic",
                                error="Empty or non-string SMILES")

    smiles = smiles.strip()
    if len(smiles) == 0:
        return ValidationResult(smiles=smiles, is_valid=False, level="basic",
                                error="Empty SMILES after stripping whitespace")

    if not _SMILES_PATTERN.match(smiles):
        invalid_chars = set(smiles) - set(
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789@+-[]()\\/#=$:.%~&!"
        )
        return ValidationResult(smiles=smiles, is_valid=False, level="basic",
                                error=f"Invalid characters: {invalid_chars}")

    # Check matched brackets
    for open_ch, close_ch, name in [('[', ']', 'bracket'), ('(', ')', 'parenthesis')]:
        count = 0
        for ch in smiles:
            if ch == open_ch:
                count += 1
            elif ch == close_ch:
                count -= 1
            if count < 0:
                return ValidationResult(smiles=smiles, is_valid=False, level="basic",
                                        error=f"Unmatched closing {name} '{close_ch}'")
        if count != 0:
            return ValidationResult(smiles=smiles, is_valid=False, level="basic",
                                    error=f"Unmatched opening {name} '{open_ch}'")

    if not any(c.isalpha() for c in smiles):
        return ValidationResult(smiles=smiles, is_valid=False, level="basic",
                                error="No atoms found in SMILES")

    return ValidationResult(smiles=smiles, is_valid=True, level="basic")


def validate_smiles_rdkit(smiles: str) -> ValidationResult:
    """Full SMILES validation using RDKit. Falls back to basic if unavailable."""
    basic_result = validate_smiles_basic(smiles)
    if not basic_result.is_valid:
        return basic_result

    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except ImportError:
        basic_result.level = "basic (rdkit unavailable)"
        return basic_result

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ValidationResult(smiles=smiles, is_valid=False, level="rdkit",
                                error="RDKit failed to parse SMILES")

    canonical = Chem.MolToSmiles(mol)
    return ValidationResult(
        smiles=smiles,
        is_valid=True,
        level="rdkit",
        canonical_smiles=canonical,
        mol_weight=round(Descriptors.MolWt(mol), 2),
        num_atoms=mol.GetNumAtoms(),
    )


def validate_smiles(smiles: str) -> ValidationResult:
    """Auto-select best available validation."""
    try:
        from rdkit import Chem  # noqa: F401
        return validate_smiles_rdkit(smiles)
    except ImportError:
        return validate_smiles_basic(smiles)


def validate_batch(smiles_list: list[str]) -> dict:
    """Validate a batch of SMILES and return summary."""
    valid = []
    invalid = []

    for smi in smiles_list:
        result = validate_smiles(smi)
        if result.is_valid:
            valid.append(result)
        else:
            invalid.append(result)

    total = len(smiles_list)
    return {
        "valid": valid,
        "invalid": invalid,
        "summary": {
            "total": total,
            "valid_count": len(valid),
            "invalid_count": len(invalid),
            "valid_rate": len(valid) / total if total > 0 else 0,
        }
    }


def filter_valid_compounds(
    compounds: list[dict], smiles_key: str = "smiles"
) -> tuple[list[dict], list[dict]]:
    """
    Filter compound dicts, keeping only those with valid SMILES.
    Returns (valid_compounds, rejected_compounds).
    """
    valid = []
    rejected = []

    for compound in compounds:
        smi = compound.get(smiles_key, "")
        result = validate_smiles(smi)
        if result.is_valid:
            compound["_smiles_validation"] = "passed"
            if result.canonical_smiles:
                compound["_canonical_smiles"] = result.canonical_smiles
            valid.append(compound)
        else:
            compound["_smiles_validation"] = f"failed: {result.error}"
            rejected.append(compound)

    return valid, rejected
