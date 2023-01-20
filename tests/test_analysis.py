from eltetrado.analysis import eltetrado

from eltetrado.cli import handle_input_file
import rnapolis.annotator
import rnapolis.parser

from eltetrado.cli import read_secondary_structure_from_dssr
from eltetrado.model import ONZM
from eltetrado.dto import generate_dto


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
    dto = generate_dto(analysis)
    assert len(dto.helices) == 1
