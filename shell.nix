let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  mach-nix = import (builtins.fetchGit {
    url = "https://github.com/DavHau/mach-nix/";
    ref = "refs/tags/3.2.0";
  }) {
    pkgs = pkgs;
    python = "python38";
  };
  customPython = mach-nix.mkPython rec {
    requirements = ''
      numpy
      ase
      PyYAML
      pytest
      pytest-datadir
      sh
    '';
  };
  mkShellNewEnv = pkgs.mkShell.override { stdenv = pkgs.gcc10Stdenv; };
  eigenClang337 = pkgs.eigen.overrideAttrs(old: rec {
    stdenv = pkgs.clangStdenv;
  });
  eigen339 = pkgs.eigen.overrideAttrs(old: rec {
    version = "3.3.9";
  });
in mkShellNewEnv {
  nativeBuildInputs = [ pkgs.cmake ];
  buildInputs = with pkgs; [
    customPython
    gtest
    bashInteractive
    which
    gcc10Stdenv
    gfortran
#   valgrind
    gdb
    eigen339

    zstd
    zlib
    lzma
    bzip2

    fmt
    abseil-cpp
    boost175
  ];
  shellHook = ''
  export PYTHONPATH=$(pwd):$PYTHONPATH
  export PATH=$(pwd)/bin:$PATH
'';
}
