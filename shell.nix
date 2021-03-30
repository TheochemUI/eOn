let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  mach-nix = import (builtins.fetchGit {
    url = "https://github.com/DavHau/mach-nix/";
    ref = "refs/tags/3.1.1";
  }) {
    pkgs = pkgs;

    # optionally specify the python version
    # python = "python27";

    # optionally update pypi data revision from https://github.com/DavHau/pypi-deps-db
    # pypiDataRev = "some_revision";
    # pypiDataSha256 = "some_sha256";
  };
  tsase = mach-nix.buildPythonPackage {
    src = "http://theory.cm.utexas.edu/code/tsase.tgz";
    extras = "ase";
    packagesExtra = "perl";
  };
  customPython = mach-nix.mkPython rec {
    requirements = ''
      ase
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
    # tsase
    gtest
    bashInteractive
    which
    gcc10Stdenv
    gfortran
    valgrind
    gdb
    eigen339
  ];
}
