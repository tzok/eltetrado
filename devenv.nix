{ pkgs, ... }:
{
  languages.python = {
    enable = true;
    poetry.enable = true;
  };
  packages = with pkgs; [
    graphviz
    highs
    pandoc-include
    ruff
    zlib
  ];
  enterShell = ''
    export PYTHONPATH=src/
  '';
}
