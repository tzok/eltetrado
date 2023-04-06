from eltetrado.cli import eltetrado_cli, has_tetrad_cli


def test_has_tetrad_with_dssr(capfd):
    has_tetrad_cli(
        [
            "--input",
            "tests/files/1v3p-assembly-1.cif.gz",
            "--dssr-json",
            "tests/files/1v3p-assembly-1.json",
        ]
    )
    out, err = capfd.readouterr()
    assert out == "True\n"


def test_has_tetrad_without_dssr(capfd):
    has_tetrad_cli(
        [
            "--input",
            "tests/files/1v3p-assembly-1.cif.gz",
        ]
    )
    out, err = capfd.readouterr()
    assert out == "True\n"


def test_eltetrado_with_dssr(capfd):
    eltetrado_cli(
        [
            "--input",
            "tests/files/2awe-assembly-1.cif.gz",
            "--dssr-json",
            "tests/files/2awe-assembly-1.json",
            "--no-image",
        ]
    )
    out, err = capfd.readouterr()
    with open("tests/files/2awe-assembly-1.out.json") as f:
        assert out == f.read()


def test_eltetrado_without_dssr(capfd):
    eltetrado_cli(
        [
            "--input",
            "tests/files/6fc9-assembly-1.cif.gz",
            "--no-image",
        ]
    )
    out, err = capfd.readouterr()
    with open("tests/files/6fc9-assembly-1.out.json") as f:
        assert out == f.read()
