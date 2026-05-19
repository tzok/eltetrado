import pytest

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
        ]
    )
    out, err = capfd.readouterr()
    with open("tests/files/6fc9-assembly-1.out") as f:
        assert out == f.read()


def test_g4composer_export_for_6fc9(tmp_path, capfd):
    output_path = tmp_path / "6fc9.g4c"

    eltetrado_cli(
        [
            "--input",
            "tests/files/6fc9-assembly-1.cif.gz",
            "--g4composer-output",
            str(output_path),
        ]
    )

    out, err = capfd.readouterr()
    with open("tests/files/6fc9-assembly-1.out") as f:
        assert out == f.read()

    assert (
        output_path.read_text()
        == "name        6fc9-assembly-1\n"
        "sequence    ggttggcgcgaagcattcgcgggttgg\n"
        "structure   ^^..^^...............^^..^^\n"
        "chi         S...S................S...S.\n"
        "sugar       SSNNSSSSNSSSSSSNNSSSSNSNNSS\n"
        "orient      A+;B+\n"
        "rise        3.4\n"
        "twist       28.6\n"
        "path        A1;B1;B4;A4;A3;B3;B2;A2\n"
    )


def test_g4composer_export_rejects_multiple_quadruplexes(tmp_path):
    output_path = tmp_path / "2awe.g4c"

    with pytest.raises(SystemExit):
        eltetrado_cli(
            [
                "--input",
                "tests/files/2awe-assembly-1.cif.gz",
                "--external-files",
                "tests/files/2awe-assembly-1.json",
                "--tool",
                "dssr",
                "--g4composer-output",
                str(output_path),
            ]
        )

    assert not output_path.exists()
