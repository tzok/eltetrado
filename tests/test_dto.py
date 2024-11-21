import rnapolis.annotator
import rnapolis.parser

from eltetrado.analysis import eltetrado
from eltetrado.cli import handle_input_file, read_secondary_structure_from_dssr
from eltetrado.dto import convert_nucleotides, convert_tetrads, generate_dto


def test_convert_nucleotides():
    """
    Load structure in nucleic_acid_only=False mode, but still do not serialize non-nucleotides in convert_nucleotides()
    """
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False, False, 2)
    nucleotides = convert_nucleotides(analysis)
    assert len(nucleotides) == 27


def test_ions():
    """
    Load structure in nucleic_acid_only=False mode and get the ions analysis
    """
    cif = handle_input_file("tests/files/7dju-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    base_interactions = rnapolis.annotator.extract_base_interactions(structure3d)
    analysis = eltetrado(base_interactions, structure3d, False, False, 2)
    tetrads = convert_tetrads(analysis.helices[0].quadruplexes[0])
    assert len(tetrads) > 0
    assert len(tetrads[0].ionsChannel) > 0
    assert tetrads[0].ionsChannel[0] == "Pt"


def test_7zko():
    """
    In 7zko there are two helices, but the second one has only one tetrad
    so it should be omitted from the output
    """
    cif = handle_input_file("tests/files/7zko-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/7zko-assembly-1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    assert len(analysis.helices) == 2
    assert len(analysis.helices[1].quadruplexes) == 1
    assert len(analysis.helices[1].quadruplexes[0].tetrads) == 1
    dto = generate_dto(analysis)
    assert len(dto.helices) == 1


def test_1v3p():
    """
    In 1v3p there is one helix with two quadruplexes, but each has only
    one tetrad so altogether we do not want to serialize it
    """
    cif = handle_input_file("tests/files/1v3p-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/1v3p-assembly-1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    assert len(analysis.helices) == 1
    assert len(analysis.helices[0].quadruplexes) == 2
    assert len(analysis.helices[0].quadruplexes[0].tetrads) == 1
    assert len(analysis.helices[0].quadruplexes[1].tetrads) == 1

    dto = generate_dto(analysis)
    assert len(dto.helices) == 0


def test_2awe():
    """
    In 2awe there is a residue G.U4 which only has the phosphorus group
    and is therefore not recognized as a nucleotide, but it was present in the
    loop section of the output JSON
    """
    cif = handle_input_file("tests/files/2awe-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/2awe-assembly-1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    dto = generate_dto(analysis)
    assert "G.U4" not in [nt.fullName for nt in dto.nucleotides]
    assert "G.U4" not in [
        nt
        for h in dto.helices
        for q in h.quadruplexes
        for l in q.loops
        for nt in l.nucleotides
    ]


def test_5v3f():
    """
    In 5V3F there are O+ and O- tetrads and the two-line dot-bracket has to take that into account
    """
    cif = handle_input_file("tests/files/5v3f-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/5v3f-assembly-1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    dto = generate_dto(analysis)
    assert (
        dto.dotBracket.sequence == "GUGCGAAGGGACGGUGCGGAGAGGAGAGCA-CGGGACGGUGCGGAGAGGAG"
    )
    assert dto.dotBracket.line1 == ".......([{..)].(.[{.).]}.}....-.([{..)].(.[{.).]}.}"
    assert dto.dotBracket.line2 == ".......([(..[{.).]}.{.)].}....-.([(..[{.).]}.{.)].}"


def test_6a85():
    """
    In 6A85 the two-line dot-bracket notation requires more than 10 different bracket levels
    """
    cif = handle_input_file("tests/files/6a85-assembly1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, 1)
    structure2d = read_secondary_structure_from_dssr(
        structure3d, 1, "tests/files/6a85-assembly1.json"
    )
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    dto = generate_dto(analysis)
    assert (
        dto.dotBracket.sequence
        == "AGAGAGATGGGTGCGT-TAGAGAGATGGGTGCGT-TAGAGAGATGGGTGCGT-TAGAGAGATGGGTGCGTT"
    )
    assert (
        dto.dotBracket.line1
        == ".([[{{.<ABCDE.F.-..)(][}.>abcde.f.-..{)(]<.ABCDEF.G.-..}])}>.abcdef.g.."
    )
    assert (
        dto.dotBracket.line2
        == ".(([[{.<ABCDE.F.-..{)(][.<ABCDE.F.-..}{)G].>abcde.f.-..)}]g}.>abcde.f.."
    )
