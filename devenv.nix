{
  inputs,
  lib,
  pkgs,
  ...
}:
{
  languages.python.uv = {
    enable = true;
    sync.enable = true;
  };
  packages = with pkgs; [
    graphviz
    highs
    pandoc-include
    zlib
  ];
  env.LD_LIBRARY_PATH = "${lib.makeLibraryPath [
    pkgs.stdenv.cc.cc.lib
    pkgs.zlib
  ]}";
  tasks = {
    "venv:patchelf-ruff" = {
      description = "Patch ruff binary in virtualenv to work on NixOS";
      exec = "${pkgs.auto-patchelf}/bin/auto-patchelf --libs ${pkgs.gcc} --paths $VIRTUAL_ENV/bin/ruff";
      after = [ "devenv:python:uv" ];
      before = [ "devenv:enterShell" ];
    };
  };
}
