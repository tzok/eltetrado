{ pkgs, ... }:
{
  languages.python.uv = {
    enable = true;
    sync.enable = true;
  };
  packages = with pkgs; [
    graphviz
    highs
    pandoc-include
    ruff
    zlib
  ];
}
