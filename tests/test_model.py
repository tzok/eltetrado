import pytest

from eltetrado.model import ONZM


def test_onzm_invalid_value():
    with pytest.raises(RuntimeError):
        ONZM.from_value("testing")
