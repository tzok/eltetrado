import pytest
from eltetrado.analysis import eltetrado

from eltetrado.cli import handle_input_file
import rnapolis.annotator
import rnapolis.parser

from eltetrado.dto import convert_nucleotides


def test_convert_nucleotides():
    cif = handle_input_file("tests/files/6fc9.cif.gz")
    structure3d = rnapolis.parser.read_3d_structure(cif, nucleic_acid_only=False)
    structure2d = rnapolis.annotator.extract_secondary_structure(structure3d)
    analysis = eltetrado(structure2d, structure3d, False, False, 2)
    nucleotides = convert_nucleotides(analysis)
    assert len(nucleotides) == 27
