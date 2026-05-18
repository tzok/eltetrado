{
  inputs,
  lib,
  pkgs,
  ...
}:
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
  env.LD_LIBRARY_PATH = "${lib.makeLibraryPath [
    pkgs.stdenv.cc.cc.lib
    pkgs.zlib
  ]}";
}
