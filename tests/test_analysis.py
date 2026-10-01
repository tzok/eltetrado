import math
from pathlib import Path

import numpy
import rnapolis.annotator
import rnapolis.parser
from rnapolis.adapter import ExternalTool, parse_external_output

import eltetrado.analysis as analysis_module
from eltetrado.analysis import calculate_best_fit_rotation_around_axis, eltetrado
from eltetrado.cli import handle_input_file, read_secondary_structure_from_dssr
from eltetrado.g4composer import (
    G4ComposerMode,
    canonical_dot_bracket,
    export_residues,
    generate_g4composer_entry,
)
from eltetrado.model import SnapbackOrientation, StrandPolarity


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
    assert entry.rise == "6.9;-3.5"
    assert entry.twist == "-11.1;-26.8"


def test_g4composer_consistent_mode_in_5dea():
    """Consistent mode shifts the classic numbering cyclically per tetrad.

    In 5DEA the chain visits tetrad A in tract order 1, 4, 2, 3 (so A keeps
    its classic numbering, anchored at A1) and tetrad B starts at tract 4, so
    classic ``B4`` becomes consistent ``B1`` and the remaining B numbers keep
    their cyclic order (B1->2, B2->3, B3->4). Orient never changes. Rise is
    positive-first (the shared axis is oriented accordingly). The consistent
    twist is the torsion signed by the interval's rise: the first interval
    (rise +6.9) keeps the torsion 77.1, the second one (rise -3.5) flips it.
    """
    cif = handle_input_file("tests/files/5dea-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(
        analysis, quadruplex, "5dea-assembly1", G4ComposerMode.CONSISTENT
    )

    assert entry.path == "A1;A4;B1;C1;B2;C2;A2;B3;C3;B4;C4;A3"
    assert entry.orient == "A-;B+;C+"
    assert entry.rise == "6.9;-3.5"
    assert entry.twist == "77.1;-26.9"


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

    output_path = str(tmp_path / "2rsk-assembly1.inp")
    write_g4composer(analysis, "tests/files/2rsk-assembly1.cif.gz", output_path)

    outputs = sorted(tmp_path.glob("2rsk-assembly1*.inp"))
    assert len(outputs) == 1
    assert outputs[0] == tmp_path / "2rsk-assembly1.inp"

    content = outputs[0].read_text()
    blocks = content.split("\n\n")
    assert len(blocks) == 2
    assert all(b.startswith("name        2rsk-assembly1\n") for b in blocks)


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


def test_tags_two_block_and_interlocked():
    """two_block applies to any molecularity: a column visited by two or more
    chain segments that each stack at least two guanines, with a clean row
    partition (2ms9 unimolecular; the bimolecular 6jje/6jjf family and 6sx6
    own whole columns or split cleanly by rows). 1+2 / 2+1 tetrad caps
    (5dea, 5de5, 5v3f Mango) are NOT two-blocks. Staggered multi-chain
    arrangements (1jjp, locked tetramolecular 1nyd/1xce) have no row
    partition and are tagged interlocked instead."""

    def tags(path: str) -> list[list[str]]:
        cif = handle_input_file(path)
        structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
        base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
        analysis = eltetrado(base_interactions, structure3d, False)
        return [
            quadruplex.tags
            for helix in analysis.helices
            for quadruplex in helix.quadruplexes
        ]

    assert all("two_block" in quad_tags for quad_tags in tags("tests/files/2ms9.cif"))

    for path in ["6jjf-assembly1.cif.gz", "6sx6-assembly1.cif.gz"]:
        multi_chain = tags(f"tests/files/{path}")
        assert all("two_block" in quad_tags for quad_tags in multi_chain), path
        assert all("interlocked" not in quad_tags for quad_tags in multi_chain), path

    # 1T+2T / 2T+1T aptamers are not two-blocks
    for path in [
        "tests/files/5dea-assembly1.cif.gz",
        "tests/files/5de5-assembly1.cif.gz",
        "tests/files/5v3f-assembly1.cif.gz",
    ]:
        assert all("two_block" not in quad_tags for quad_tags in tags(path)), path

    for path in ["1jjp", "1nyd", "1xce"]:
        interlocked = tags(f"tests/files/{path}.cif.gz")
        assert all("interlocked" in quad_tags for quad_tags in interlocked), path
        assert all("two_block" not in quad_tags for quad_tags in interlocked), path

    # chain-pure multi-chain quadruplexes are neither two-block nor interlocked
    for path in [
        "tests/files/2awe-assembly-1.cif.gz",
        "tests/files/6a85-assembly1.cif.gz",
    ]:
        for quad_tags in tags(path):
            assert "two_block" not in quad_tags
            assert "interlocked" not in quad_tags


def test_v_loop_and_non_g4_tag_rules():
    """V-loops require an isolated tetrad nucleotide at a propeller-type
    connector; same-tetrad V-BRAKs, missing-residue breaks and non-G4
    quadruple interactions are excluded."""
    cif = handle_input_file("tests/files/4p1d-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    interactions = parse_external_output(
        ["tests/files/4p1d-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(interactions, structure3d, False)
    assert all("v_loop" not in q.tags for h in analysis.helices for q in h.quadruplexes)

    cif = handle_input_file("tests/files/2kpr-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    analysis = eltetrado(
        rnapolis.annotator.extract_base_interactions(structure3d),
        structure3d,
        False,
    )
    quadruplex = analysis.helices[0].quadruplexes[0]
    assert "v_loop" in quadruplex.tags
    assert any(
        s.orientation == SnapbackOrientation.FIVE_PRIME for s in quadruplex.snapbacks
    )

    cif = handle_input_file("tests/files/7sxp-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    analysis = eltetrado(
        rnapolis.annotator.extract_base_interactions(structure3d),
        structure3d,
        False,
    )
    assert all("v_loop" not in q.tags for h in analysis.helices for q in h.quadruplexes)

    cif = handle_input_file("tests/files/8gxc-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    interactions = parse_external_output(
        ["tests/files/8gxc-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(interactions, structure3d, False)
    assert all(q.tags == ["non_g4"] for h in analysis.helices for q in h.quadruplexes)


def test_reentry_snapbacks_in_asc_g4_examples():
    expected = {
        "2mgn": SnapbackOrientation.THREE_PRIME,
        "6jcd": SnapbackOrientation.THREE_PRIME,
        "6xt7": SnapbackOrientation.FIVE_PRIME,
    }
    for pdb, orientation in expected.items():
        cif = handle_input_file(f"tests/files/{pdb}-assembly1.cif.gz")
        structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
        analysis = eltetrado(
            rnapolis.annotator.extract_base_interactions(structure3d),
            structure3d,
            False,
        )
        quadruplexes = [
            quadruplex
            for helix in analysis.helices
            for quadruplex in helix.quadruplexes
        ]
        assert quadruplexes
        assert all(
            any(snapback.orientation == orientation for snapback in q.snapbacks)
            for q in quadruplexes
        ), pdb


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
    traverses the stack graph in reverse; the shared axis is oriented so the
    first rise is positive.
    """
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/5zev-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(analysis, quadruplex, "5zev-assembly1")

    assert entry.orient == "A+;B-;C+"
    assert entry.rise == "2.9;-6.3"
    assert entry.twist == "31.8;-63.5"


def test_g4composer_consistent_mode_in_5zev():
    """Consistent mode numbers each tetrad geometrically: the first guanine
    visited by the chain becomes 1 and the remaining guanines follow in the
    direction that is clockwise along the g4 axis (validated against the
    expert's canonical 5ZEV entry). The orient follows the H-bond successor
    cycle direction around the same axis; the torsion is signed by the
    interval's rise."""
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/5zev-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(
        analysis, quadruplex, "5zev-assembly1", G4ComposerMode.CONSISTENT
    )

    assert entry.path == "A1;B1;B3;A3;C1;B4;A4;C2;C3;B2;A2;C4"
    assert entry.orient == "A+;B-;C+"
    assert entry.rise == "2.9;-6.3"
    assert entry.twist == "31;116.2"


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


def test_g4composer_consistent_mode_in_6fc9():
    """In 6FC9 the chain visits each tetrad's tracts in ascending order, so
    the consistent path is identical to the classic one; the torsion twist of
    this right-handed quadruplex stays positive and close to the classic one."""
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]
    entry = generate_g4composer_entry(
        analysis, quadruplex, "6fc9-assembly-1", G4ComposerMode.CONSISTENT
    )

    assert entry.path == "A1;B1;B2;A2;A3;B3;B4;A4"
    assert entry.twist == "17.9"


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
    for tract, tract_polarities in zip(
        quadruplex.tracts, quadruplex.strand_polarities, strict=True
    ):
        for nt, polarity in zip(tract.nucleotides, tract_polarities, strict=True):
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
    for tract, tract_polarities in zip(
        quadruplex.tracts, quadruplex.strand_polarities, strict=True
    ):
        for nt, polarity in zip(tract.nucleotides, tract_polarities, strict=True):
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


def test_8psi_has_four_loops_and_reentry_snapback():
    """
    8psi forms a fold in which the chain visits one of the tracts twice,
    giving four cross-tract loops and a 5' re-entry snapback. The loop
    classification requires exactly three loops, so it stays unset.
    """
    cif = handle_input_file("tests/files/8psi-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert [loop.loop_type.value for loop in quadruplex.loops] == [
        "lateral+",
        "lateral-",
        "propeller-",
        "lateral-",
    ]
    assert all(loop.nucleotides for loop in quadruplex.loops)
    assert len(quadruplex.snapbacks) == 1
    assert quadruplex.snapbacks[0].orientation == SnapbackOrientation.FIVE_PRIME
    assert len(quadruplex.snapbacks[0].nucleotides) == 8
    assert quadruplex.loop_class is None


def test_7wgw_reports_three_prime_snapback():
    """
    7wgw ends with a chain re-entry into an already visited tract: after the
    third propeller loop the chain leaves the column and re-enters it three
    nucleotides later (8-oxo-G, A, A), injecting the last guanosine into the
    column - a snapback 3'.
    """
    cif = handle_input_file("tests/files/7wgw-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert [loop.loop_type.value for loop in quadruplex.loops] == [
        "propeller-",
        "propeller-",
        "propeller-",
    ]
    assert len(quadruplex.snapbacks) == 1
    snapback = quadruplex.snapbacks[0]
    assert snapback.orientation == SnapbackOrientation.THREE_PRIME
    assert [nt.full_name for nt in snapback.nucleotides] == [
        "X.8OG18",
        "X.DA19",
        "X.DA20",
    ]


def test_5zev_emits_zero_length_loop():
    """
    5zev contains a 0-nt cross-tract connector: two tetrad guanines directly
    adjacent in the chain but stacked in different tracts. Such a loop
    carries no nucleotides and brings the loop count to four, so the loop
    classification stays unset.
    """
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert len(quadruplex.loops) == 4
    zero_length = [loop for loop in quadruplex.loops if not loop.nucleotides]
    assert len(zero_length) == 1
    assert zero_length[0].loop_type.value == "propeller-"
    assert len(quadruplex.snapbacks) == 1
    assert quadruplex.snapbacks[0].orientation == SnapbackOrientation.FIVE_PRIME
    assert len(quadruplex.snapbacks[0].nucleotides) == 17
    assert quadruplex.loop_class is None


def test_5de5_reports_zero_length_snapback():
    """
    5de5 re-enters the same tract with no nucleotides in between (the chain
    leaves the column and immediately returns at a non-adjacent layer),
    producing a 0-nt snapback - contrast with 5zev's 0-nt loop, which joins
    two different tracts.
    """
    cif = handle_input_file("tests/files/5de5-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert len(quadruplex.loops) == 4
    assert len(quadruplex.snapbacks) == 2
    assert quadruplex.snapbacks[0].orientation == SnapbackOrientation.THREE_PRIME
    assert quadruplex.snapbacks[0].nucleotides == []
    assert quadruplex.snapbacks[1].orientation == SnapbackOrientation.FIVE_PRIME
    assert len(quadruplex.snapbacks[1].nucleotides) == 8


def test_6fc9_tags_chair_and_quadruplex_duplex():
    """
    6fc9's long lateral loop forms a Watson-Crick duplex stem, so the
    quadruplex is tagged with the common topology name and
    quadruplex_duplex.
    """
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.tags == ["chair", "quadruplex_duplex"]


def test_2awe_tags_bulge_and_inverted_strand_polarity():
    """
    2awe carries one bulged U per tract and tracts whose nucleotides point
    in both directions along the helix axis.
    """
    cif = handle_input_file("tests/files/2awe-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.tags == ["has_bulge", "has_inverted_strand_polarity"]


def test_7wgw_tags_snapback_and_parallel():
    cif = handle_input_file("tests/files/7wgw-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.tags == [
        "has_snapback",
        "has_inverted_strand_polarity",
        "parallel",
    ]


def test_6rs3_tags_v_loop():
    """
    6rs3 contains a seemingly missing loop: two tetrad guanines directly
    adjacent in the chain but stacked in different tracts (a V-loop).
    """
    cif = handle_input_file("tests/files/6rs3-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert "v_loop" in quadruplex.tags
    assert "quadruplex_duplex" not in quadruplex.tags


def test_2ms9_tags_two_block():
    """
    2ms9 is unimolecular but built from two separate two-tetrad blocks:
    three of its columns are each visited by two chain segments stacking
    two guanines apiece.
    """
    cif = handle_input_file("tests/files/2ms9.cif")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert "two_block" in quadruplex.tags


def test_5zev_tags_zero_length_loop_without_duplex():
    """
    5zev's V-loop (its 0-nt cross-tract connector) is tagged, while its
    loops form no canonical duplex.
    """
    cif = handle_input_file("tests/files/5zev-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert "v_loop" in quadruplex.tags
    assert "quadruplex_duplex" not in quadruplex.tags


def test_8psi_has_snapback_and_v_loop_tags():
    """
    8psi has a re-entry snapback and a V-loop at the broken first column.
    """
    cif = handle_input_file("tests/files/8psi-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    assert quadruplex.tags == ["has_snapback", "v_loop"]


class _StubPair:
    def __init__(self, rise, twist):
        self.rise = rise
        self.twist = twist


class _StubBulgeGroup:
    def __init__(self, nucleotides):
        self.nucleotides = [None] * nucleotides


class _StubQuadruplex:
    def __init__(self, tetrads=2, bulges=0, pairs=()):
        self.tetrads = [None] * tetrads
        self.bulges = [_StubBulgeGroup(bulges)] if bulges else []
        self.tetrad_pairs = [_StubPair(*pair) for pair in pairs]


def test_spurious_reasons_accepts_realistic_quadruplex():
    from eltetrado.analysis import quadruplex_spurious_reasons

    quadruplex = _StubQuadruplex(tetrads=3, bulges=15, pairs=[(3.4, 30.0)])
    assert quadruplex_spurious_reasons(quadruplex) == []


def test_spurious_reasons_rejects_single_tetrad():
    from eltetrado.analysis import quadruplex_spurious_reasons

    reasons = quadruplex_spurious_reasons(_StubQuadruplex(tetrads=1))
    assert len(reasons) == 1
    assert "tetrad" in reasons[0]


def test_spurious_reasons_rejects_huge_bulge():
    from eltetrado.analysis import quadruplex_spurious_reasons

    reasons = quadruplex_spurious_reasons(_StubQuadruplex(bulges=105))
    assert len(reasons) == 1
    assert "bulge" in reasons[0]


def test_spurious_reasons_rejects_unrealistic_rise():
    from eltetrado.analysis import quadruplex_spurious_reasons

    reasons = quadruplex_spurious_reasons(
        _StubQuadruplex(pairs=[(0.5, 30.0), (15.0, 30.0)])
    )
    assert len(reasons) == 2
    assert "rise" in reasons[0]
    assert "rise" in reasons[1]


def test_spurious_reasons_allows_extreme_twist():
    from eltetrado.analysis import quadruplex_spurious_reasons

    # extreme twist occurs in real structures (7d5d, 2la5) and must not
    # be rejected
    quadruplex = _StubQuadruplex(pairs=[(3.3, 141.19), (3.23, -170.74)])
    assert quadruplex_spurious_reasons(quadruplex) == []


def test_spurious_reasons_ignores_nan_geometry():
    from eltetrado.analysis import quadruplex_spurious_reasons

    quadruplex = _StubQuadruplex(pairs=[(math.nan, math.nan)])
    assert quadruplex_spurious_reasons(quadruplex) == []


def test_g4composer_8tns_matches_expert_canonical():
    """8TNS is the expert-validated canonical entry: the orient is judged
    from the H-bond successor cycle direction around the g4 axis (the axis
    oriented so the first build-order rise is positive), and the consistent
    numbering starts at the first guanine visited by the chain going
    clockwise along the same axis. Both modes must keep this output."""
    cif = handle_input_file("tests/files/8tns-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1, nucleic_acid_only=False)
    base_interactions = parse_external_output(
        ["tests/files/8tns-assembly1.json"], ExternalTool.DSSR, structure3d
    )
    analysis = eltetrado(base_interactions, structure3d, False)

    quadruplex = analysis.helices[0].quadruplexes[0]

    classic = generate_g4composer_entry(analysis, quadruplex, "8tns-assembly1")
    assert classic.orient == "A-;B+;C-;D+"
    assert classic.rise == "3.3;6.7;-3.1"
    assert classic.path == "A1;B1;C1;D1;A2;B2;C2;D2;A3;B3;C3;D3;A4;B4;C4;D4"

    consistent = generate_g4composer_entry(
        analysis, quadruplex, "8tns-assembly1", G4ComposerMode.CONSISTENT
    )
    assert consistent.orient == "A-;B+;C-;D+"
    assert consistent.twist == "-13.6;-18.1;-38.7"
    assert consistent.path == classic.path
