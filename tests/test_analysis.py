from pathlib import Path

import rnapolis.annotator
import rnapolis.parser
from rnapolis.adapter import ExternalTool, parse_external_output

from eltetrado.analysis import eltetrado
from eltetrado.cli import handle_input_file


def test_5zev_tracts():
    """
    Analyse 5ZEV structure with DSSR-generated JSON and verify the detected
    G-quadruplex tracts.
    """
    cif_path = Path("tests/files/5zev-assembly1.cif.gz")
    json_path = Path("tests/files/5zev-assembly1.json")

    # Load the 3D structure using the same helper routine as the CLI
    cif_or_pdb = handle_input_file(str(cif_path))
    structure3d = rnapolis.parser.read_3d_structure(
        cif_or_pdb, 1, nucleic_acid_only=False
    )

    # Parse external interactions produced by DSSR
    base_interactions = parse_external_output(
        [str(json_path)], ExternalTool.DSSR, structure3d
    )

    # Run ElTetrado analysis
    analysis = eltetrado(
        base_interactions,
        structure3d,
        no_reorder=False,
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


def test_6fc9_path():
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    assert [
        tuple(nt.full_name for nt in tract.nucleotides)
        for tract in analysis.helices[0].quadruplexes[0].tracts
    ] == [
        ("A.DG1", "A.DG2"),
        ("A.DG27", "A.DG26"),
        ("A.DG22", "A.DG23"),
        ("A.DG6", "A.DG5"),
    ]
    assert analysis.helices[0].quadruplexes[0].path == [
        "A1",
        "B1",
        "B4",
        "A4",
        "A3",
        "B3",
        "B2",
        "A2",
    ]


def test_2ms9_is_left_handed():
    cif = handle_input_file("tests/files/2ms9.cif")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.handedness is not None
    assert quadruplex.handedness.value == "left"
    assert [
        p.value if p is not None else None for p in quadruplex.tetrad_polarities
    ] == [
        "clockwise",
        "clockwise",
        "anticlockwise",
        "anticlockwise",
    ]
