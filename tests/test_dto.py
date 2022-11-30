import pytest
from eltetrado.analysis import eltetrado

from eltetrado.cli import handle_input_file
import rnapolis.annotator
import rnapolis.parser

from eltetrado.dto import convert_nucleotides, convert_quadruplexes, convert_tetrads


def test_convert_nucleotides():
    """
    Load structure in nucleic_acid_only=False mode, but still do not serialize non-nucleotides in convert_nucleotides()
    """
    cif = handle_input_file("tests/files/6fc9-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    structure2d = rnapolis.annotator.extract_secondary_structure(structure3d)
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    nucleotides = convert_nucleotides(analysis)
    assert len(nucleotides) == 27


def test_ions():
    """
    Load structure in nucleic_acid_only=False mode and get the ions analysis
    """
    cif = handle_input_file("tests/files/7dju-assembly-1.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    structure2d = rnapolis.annotator.extract_secondary_structure(structure3d)
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    tetrads = convert_tetrads(analysis.helices[0].quadruplexes[0])
    assert len(tetrads) > 0
    assert len(tetrads[0].ionsChannel) > 0
    assert tetrads[0].ionsChannel[0] == "PT"
