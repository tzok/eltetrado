import math
from pathlib import Path

import numpy
import rnapolis.annotator
import rnapolis.parser
from rnapolis.adapter import ExternalTool, parse_external_output

import eltetrado.analysis as analysis_module
from eltetrado.analysis import calculate_best_fit_rotation_around_axis, eltetrado
from eltetrado.model import StrandPolarity
from eltetrado.cli import handle_input_file, read_secondary_structure_from_dssr
from eltetrado.g4composer import (
    canonical_dot_bracket,
    export_residues,
    generate_g4composer_entry,
)


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


def test_path_tetrad_letters_follow_5p_order_in_5v3f():
    cif = handle_input_file("tests/files/5v3f-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/5v3f-assembly-1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False)

    assert analysis.helices[0].quadruplexes[0].path == [
        "A1",
        "B1",
        "C1",
        "A4",
        "B4",
        "C4",
        "A3",
        "B3",
        "C3",
        "A2",
        "B2",
        "C2",
    ]


def test_path_starts_with_a1_in_5zev():
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/5zev-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    path = analysis.helices[0].quadruplexes[0].path

    assert path[0] == "A1"


def test_path_starts_with_a1_in_5dea():
    cif = handle_input_file("tests/files/5dea-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    path = analysis.helices[0].quadruplexes[0].path

    assert path[0] == "A1"


def test_g4composer_export_in_5dea():
    cif = handle_input_file("tests/files/5dea-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(analysis, quadruplex, "5dea-assembly1")

    assert entry.path == "A1;A4;B4;C4;B1;C1;A2;B2;C2;B3;C3;A3"
    assert entry.orient == "A-;B+;C+"
    assert entry.rise == "-6.9;3.5"
    assert entry.twist == "11.1;26.8"


def test_loop_signs_follow_columns_in_1i34():
    cif = handle_input_file("tests/files/1i34-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert [loop.loop_type.value for loop in quadruplex.loops] == [
        "diagonal",
        "propeller+",
        "diagonal",
    ]


def test_tetrad_polarity_inferred_for_non_canonical_tetrads_in_2km3():
    cif = handle_input_file("tests/files/2km3-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/2km3-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert len(quadruplex.tetrads) == 3
    assert all(p is not None for p in quadruplex.tetrad_polarities)


def test_g4composer_exports_multiple_quadruplexes_in_2rsk(tmp_path):
    cif = handle_input_file("tests/files/2rsk-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/2rsk-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    from eltetrado.g4composer import write_g4composer

    output_base = str(tmp_path / "2rsk-assembly1.inp")
    write_g4composer(analysis, "tests/files/2rsk-assembly1.cif.gz", output_base)

    outputs = sorted(tmp_path.glob("2rsk-assembly1*.inp"))
    assert len(outputs) == 2
    assert any("-A" in p.name for p in outputs)
    assert any("-B" in p.name for p in outputs)


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


def test_best_fit_rotation_around_axis_reports_signed_rotation():
    axis = numpy.array([0.0, 0.0, 1.0])
    points1 = numpy.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ]
    )
    angle = math.radians(30.0)
    rotation = numpy.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    points2 = points1 @ rotation.T + numpy.array([0.0, 0.0, 2.5])

    twist = calculate_best_fit_rotation_around_axis(points1, points2, axis)

    assert math.isclose(twist, 30.0, abs_tol=1.0e-6)


def test_best_fit_rotation_around_axis_handles_distorted_rectangle():
    axis = numpy.array([0.0, 0.0, 1.0])
    points1 = numpy.array(
        [
            [2.0, 0.1, 0.0],
            [0.3, 1.2, 0.0],
            [-1.8, -0.2, 0.0],
            [-0.4, -1.0, 0.0],
        ]
    )
    angle = math.radians(-45.0)
    rotation = numpy.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    distortion = numpy.array(
        [
            [0.05, -0.03, 0.0],
            [-0.02, 0.04, 0.0],
            [0.01, 0.02, 0.0],
            [-0.04, -0.01, 0.0],
        ]
    )
    points2 = points1 @ rotation.T + distortion + numpy.array([0.0, 0.0, -3.0])

    twist = calculate_best_fit_rotation_around_axis(points1, points2, axis)

    assert math.isclose(twist, -45.0, abs_tol=2.0)


def test_tetrad_pair_twist_is_nan_when_c1prime_is_missing(monkeypatch):
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)

    original = analysis_module.residue_c1prime_coordinates

    def patched_residue_c1prime_coordinates(residue):
        if residue.full_name == "A.DG1":
            return None
        return original(residue)

    monkeypatch.setattr(
        analysis_module,
        "residue_c1prime_coordinates",
        patched_residue_c1prime_coordinates,
    )

    analysis = eltetrado(base_interactions, structure3d, False)

    assert math.isnan(analysis.helices[0].quadruplexes[0].tetrad_pairs[0].twist)


def test_incomplete_tetrad_residue_does_not_crash_geometry(monkeypatch):
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)

    original = analysis_module.residue_topology_point

    def patched_residue_topology_point(residue):
        if residue.full_name == "A.DG1":
            return None
        return original(residue)

    monkeypatch.setattr(
        analysis_module,
        "residue_topology_point",
        patched_residue_topology_point,
    )

    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.path
    assert len(quadruplex.tetrad_polarities) == len(quadruplex.tetrads)


def test_g4composer_non_linear_intervals_follow_build_order():
    """Document a nonlinear build-order case used by the g4composer exporter.

    In 5ZEV the adjacent stack order is B-A-C, while the g4composer build order
    inferred from the path is A-B-C. The first exported interval therefore
    traverses the stack graph in reverse, which flips the sign of rise/twist.
    """
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/5zev-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(analysis, quadruplex, "5zev-assembly1")

    assert entry.orient == "A-;B+;C-"
    assert entry.rise == "-2.9;6.3"
    assert entry.twist == "-31.8;63.5"


def test_g4composer_structure_keeps_canonical_brackets_for_flanks():
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    residues = export_residues(analysis, quadruplex)

    assert canonical_dot_bracket(analysis, residues) == "......((((((...))))))......"

    entry = generate_g4composer_entry(analysis, quadruplex, "6fc9-assembly-1")

    assert entry.structure == "^^..^^((((((...))))))^^..^^"


def test_6fc9_strand_polarities_are_consistent_within_tracts():
    """
    In a canonical antiparallel quadruplex (6fc9) every tract should have
    a uniform strand polarity (no inverted polarity inside a single tract).
    """
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    assert len(quadruplex.strand_polarities) == len(quadruplex.tracts)

    for tract_polarities in quadruplex.strand_polarities:
        non_none = [p for p in tract_polarities if p is not None]
        if non_none:
            assert all(p == non_none[0] for p in non_none), (
                "Expected uniform strand polarity within a tract, got mixed signs"
            )


def test_5de5_g20_has_minus_strand_polarity():
    """
    In 5DE5 (published inverted strand polarity example) every tract should
    show inverted polarity and A.G20 in particular must be MINUS.
    """
    cif = handle_input_file("tests/files/5de5-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    assert len(quadruplex.strand_polarities) == len(quadruplex.tracts)

    # Every tract should have mixed signs (inverted strand polarity)
    for tract_polarities in quadruplex.strand_polarities:
        non_none = [p for p in tract_polarities if p is not None]
        if non_none:
            has_plus = any(p.value == "plus" for p in non_none)
            has_minus = any(p.value == "minus" for p in non_none)
            assert has_plus and has_minus, (
                "Expected inverted strand polarity in 5DE5 tract"
            )

    # Find G20 specifically and assert it is MINUS
    g20_polarity = None
    for tract, tract_polarities in zip(quadruplex.tracts, quadruplex.strand_polarities):
        for nt, polarity in zip(tract.nucleotides, tract_polarities):
            if nt.full_name == "A.G20":
                g20_polarity = polarity
                break
        if g20_polarity is not None:
            break

    assert g20_polarity is not None, "A.G20 not found in strand polarities"
    assert g20_polarity.value == "minus", (
        f"Expected A.G20 to be minus, got {g20_polarity.value}"
    )


def test_5v3f_most_5prime_guanine_is_plus():
    """
    In 5V3F the most 5' guanine (G8) should be the reference and drawn
    upright (PLUS). Guanines with opposite strand direction (G16, G21, G26)
    should be MINUS (inverted), matching the published convention.
    """
    cif = handle_input_file("tests/files/5v3f-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    # First quadruplex (chain A)
    quadruplex = analysis.helices[0].quadruplexes[0]

    # Build a mapping from full_name to polarity
    polarity_map = {}
    for tract, tract_polarities in zip(quadruplex.tracts, quadruplex.strand_polarities):
        for nt, polarity in zip(tract.nucleotides, tract_polarities):
            polarity_map[nt.full_name] = polarity

    # Most 5' guanine should be PLUS (normal/upright)
    assert polarity_map.get("A.G8") == StrandPolarity.PLUS
    assert polarity_map.get("A.G9") == StrandPolarity.PLUS
    assert polarity_map.get("A.G10") == StrandPolarity.PLUS

    # These should be MINUS (inverted relative to the reference)
    assert polarity_map.get("A.G16") == StrandPolarity.MINUS
    assert polarity_map.get("A.G21") == StrandPolarity.MINUS
    assert polarity_map.get("A.G26") == StrandPolarity.MINUS

    # These should be PLUS (same direction as reference)
    assert polarity_map.get("A.G14") == StrandPolarity.PLUS
    assert polarity_map.get("A.G19") == StrandPolarity.PLUS
    assert polarity_map.get("A.G24") == StrandPolarity.PLUS
