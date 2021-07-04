let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
  mach-nix = import (builtins.fetchGit {
    url = "https://github.com/DavHau/mach-nix.git";
    ref = "refs/tags/3.3.0";
  }) {
    pkgs = pkgs;
    python = "python38";
  };
  customPython = mach-nix.mkPython {
    requirements = builtins.readFile ./requirements.txt;
    packagesExtra = [
      # "https://github.com/psf/requests/tarball/2a7832b5b06d"   # from tarball url
      ./.                                     # from local path
      # mach-nix.buildPythonPackage { ... };                     # from package
    ];  };
 #   mkShellNewEnv = pkgs.mkShell.override { stdenv = pkgs.gcc10Stdenv; };
    eigenClang337 = pkgs.eigen.overrideAttrs(old: rec {
      stdenv = pkgs.clangStdenv;
    });
    eigen339 = pkgs.eigen.overrideAttrs(old: rec {
      version = "3.3.9";
    });
  myCmop = pkgs.wrapCC (pkgs.gcc10.cc.override {
    langFortran = true;
    langCC = true;
    langC = true;
    enableShared = true;
    enableMultilib = false;
    staticCompiler = false;
    profiledCompiler = false;
  });
  mycompiler = myCmop.overrideAttrs(old: rec {
  hardeningEnable = ["pic"];
  });
in pkgs.mkShell {
  nativeBuildInputs = with pkgs; [ cmake blas mycompiler mycompiler.cc.lib lapack ];
  buildInputs = with pkgs; [
    gtest
    bashInteractive
    which
    customPython
  #  gcc10Stdenv
  #  gfortran
    #   valgrind
    lldb
    graphviz
   # gfortran
   # gfortran.cc
   # mycompiler.cc

    zstd
    zlib
    lzma
    bzip2
    openblas

    fmt
    abseil-cpp
    boost175
    eigen339
  ];
  shellHook = ''
    # For eonclient
    export PATH=$(pwd)/client/build:$PATH
  '';
}
