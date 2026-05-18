{ inputs, pkgs, ... }:
{
  overlays = [
    (final: prev: { ruff = (import inputs.nixpkgs-unstable { system = prev.stdenv.system; }).ruff; })
  ];
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
