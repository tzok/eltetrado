from eltetrado.cli import eltetrado_cli, has_tetrad_cli


def test_has_tetrad_with_dssr(capfd):
    has_tetrad_cli(
        [
            "--input",
            "tests/files/1v3p-assembly-1.cif.gz",
            "--external-files",
            "tests/files/1v3p-assembly-1.json",
            "--tool",
            "dssr",
        ]
    )
    out, err = capfd.readouterr()
    assert out == "True\n"


def test_has_tetrad_with_fr3d(capfd):
    has_tetrad_cli(
        [
            "--input",
            "tests/files/1v3p-assembly-1.cif.gz",
            "--external-files",
            "tests/files/1v3p-assembly-1.txt",
            "--tool",
            "fr3d",
        ]
    )
    out, err = capfd.readouterr()
    assert out == "True\n"


def test_has_tetrad(capfd):
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
            "--external-files",
            "tests/files/2awe-assembly-1.json",
            "--tool",
            "dssr",
            "--no-image",
        ]
    )
    out, err = capfd.readouterr()
    with open("tests/files/2awe-assembly-1.out") as f:
        assert out == f.read()


def test_eltetrado_with_fr3d(capfd):
    eltetrado_cli(
        [
            "--input",
            "tests/files/2awe-assembly-1.cif.gz",
            "--external-files",
            "tests/files/2awe-assembly-1.txt",
            "--tool",
            "fr3d",
            "--no-image",
        ]
    )
    out, err = capfd.readouterr()
    with open("tests/files/2awe-assembly-1-fr3d.out") as f:
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
    with open("tests/files/6fc9-assembly-1.out") as f:
        assert out == f.read()
