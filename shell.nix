# Generates a development environment for eonclient
#
# Usage examples:
#
# To create a dev environment:
#
#   nix-shell
#
# With clang:
#
#   nix-shell --argstr compilerUsed clang
#
{ withEonclient ? false, compilerUsed ? "gcc" }:
let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs {};
  mach-nix = import (
    builtins.fetchGit {
      url = "https://github.com/DavHau/mach-nix.git";
      ref = "refs/tags/3.3.0";
    }
  ) {
    pkgs = pkgs;
    python = "python38";
  };
  myGems = pkgs.bundlerEnv {
    name = "eonc-gems";
    # Kinda fragile
    gemdir = builtins.path { path = ./.; name = "EONgit"; };
  };
  customPython = mach-nix.mkPython {
    requirements = builtins.readFile ./requirements.txt;
    packagesExtra = [
      # "https://github.com/psf/requests/tarball/2a7832b5b06d"   # from tarball url
      ./. # from local path
      # mach-nix.buildPythonPackage { ... };                     # from package
    ];
  };
  compilerEnv = (
    if compilerUsed == "gcc" then pkgs.gcc10Stdenv
    else if compilerUsed == "clang" then pkgs.clang12Stdenv
    else pkgs.stdenv
  );
  mkShellNewEnv = pkgs.mkShell.override { stdenv = compilerEnv; };
  eigen339 = (pkgs.eigen.override{stdenv = compilerEnv; }).overrideAttrs (
    old: rec {
      version = "3.3.9";
      src = pkgs.fetchFromGitLab {
        owner = "libeigen";
        repo = "eigen";
        rev = "${version}";
        sha256 = "0m4h9fd5s1pzpncy17r3w0b5a6ywqjajmnr720ndb7fc4bn0dhi4";
      };
      # From https://github.com/foolnotion/aoc2020/blob/master/eigen_include_dir.patch
      patches = [ ./eigen_include_dir.patch ];
    }
  );
  fmt881 = (pkgs.fmt.override{stdenv = compilerEnv; }).overrideAttrs (
    old: rec {
      version = "8.1.1";
      src = pkgs.fetchFromGitHub {
        owner = "fmtlib";
        repo = "fmt";
        rev = "${version}";
        sha256 = "sha256-leb2800CwdZMJRWF5b1Y9ocK0jXpOX/nwo95icDf308=";
      };
      doCheck = false;
    }
  );
  macHook = ''
    # eonclient
    export PATH=$(pwd)/client/builddir:$PATH
  '';
  linuxHook = ''
    # eonclient
    export PATH=$(pwd)/client/builddir:$PATH
    # Locale
    export LOCALE_ARCHIVE=${pkgs.glibcLocales}/lib/locale/locale-archive
  '';
  myCmop = (
    if compilerUsed == "gcc" then (
      pkgs.wrapCC (
        pkgs.gcc10.cc.override {
          langFortran = true;
          langCC = true;
          langC = true;
          enableShared = true;
          enableMultilib = false;
          staticCompiler = false;
          profiledCompiler = false;
        }
      )
    )
    else if compilerUsed == "clang" then (
      compilerEnv
    ) else pkgs.stdEnv
  );
  mycompiler = (
    if compilerUsed == "gcc" then (
      myCmop.overrideAttrs (
        old: rec {
          hardeningEnable = [ "pic" ];
        }
      )
    ) else myCmop
  );
  # Can't use nix lldb on macos, must use native
  notDarwin = pkgs.lib.optionals (! pkgs.stdenv.isDarwin ) [ pkgs.lldb pkgs.mold ];
  isDarwin = pkgs.lib.optionals (pkgs.stdenv.isDarwin) (with pkgs.darwin; [ cctools ]);
 # CoreFoundation CoreServices Foundation
  compilerpkg = (if compilerUsed == "clang" then pkgs.clang_12 else pkgs.gcc10);
  eonclient = pkgs.callPackage ./default.nix {};
in
mkShellNewEnv {
  # nativeBuildInputs = with pkgs; [ cmake blas mycompiler (if compiler == "gcc" then mycompiler.cc.lib else null) openblas ninja ];
  buildInputs = with pkgs; [
    notDarwin
    isDarwin
    compilerpkg
    gtest
    bashInteractive
    which
    customPython
    meson ninja pkg-config cmake
    # mycompiler
    (if withEonclient then (if pkgs.stdenv.isDarwin then null else eonclient) else null)
    graphviz

    zstd
    zlib
    lzma
    bzip2
    openblas blas

    fmt881
    abseil-cpp
    boost175
    eigen339

    # For eonc.rb
    myGems
    (lowPrio myGems.wrappedRuby)
  ];
  shellHook = if pkgs.stdenv.isDarwin then macHook else linuxHook;
}
