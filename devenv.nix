{ pkgs, ... }:
{
  languages.python = {
    enable = true;
    poetry.enable = true;
  };
  packages = [
    pkgs.graphviz
    pkgs.highs
    pkgs.ruff
    pkgs.zlib
  ];
  enterShell = ''
    export PYTHONPATH=src/
  '';
}
