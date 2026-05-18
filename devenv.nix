{ pkgs, ... }:
{
  languages.python.uv = {
    enable = true;
    sync = {
      enable = true;
      arguments = [ "--locked" ];
    };
  };
  packages = with pkgs; [
    graphviz
    highs
    pandoc-include
    ruff
    zlib
  ];
}
