from pathlib import Path

from rnapolis.adapter import ExternalTool, parse_external_output
from rnapolis.tertiary import Structure3D

from eltetrado.analysis import eltetrado


def test_5zev_tracts():
    """
    Analyse 5ZEV structure with DSSR-generated JSON and verify the detected
    G-quadruplex tracts.
    """
    cif_path = Path("tests/files/5zev-assembly1.cif.gz")
    json_path = Path("tests/files/5zev-assembly1.json")

    # Load the 3D structure (API differences handled for older/newer rnapolis)
    if hasattr(Structure3D, "from_file"):
        structure3d = Structure3D.from_file(cif_path, model=1)
    elif hasattr(Structure3D, "from_cif"):
        structure3d = Structure3D.from_cif(cif_path, model=1)  # type: ignore[attr-defined]
    else:
        raise RuntimeError("Unable to create Structure3D instance from mmCIF")

    # Parse external interactions produced by DSSR
    base_interactions = parse_external_output(
        [str(json_path)], ExternalTool.DSSR, structure3d
    )

    # Run ElTetrado analysis
    analysis = eltetrado(
        base_interactions,
        structure3d,
        strict=False,
        no_reorder=False,
        stacking_mismatch=0,
    )

    # Expected tracts (in 5'→3' order within each tract)
    expected_tracts = {
        ("A.DG9", "A.DG10", "A.DG11"),
        ("A.DG21", "A.DG22", "A.DG23"),
        ("A.DG13", "A.DG15", "A.DG16"),
        ("A.DG2", "A.DG1", "A.DG20"),
    }

    # Collect tracts detected by ElTetrado
    detected_tracts = {
        tuple(nt.full_name for nt in tract.nucleotides)
        for helix in analysis.helices
        for quadruplex in helix.quadruplexes
        for tract in quadruplex.tracts
    }

    for tract in expected_tracts:
        assert tract in detected_tracts, (
            f"Expected tract {tract} not found; detected tracts: {detected_tracts}"
        )
