{ sources ? import ./nix/sources.nix }:
let
  pkgs = import sources.nixpkgs { };
    eigen339 = pkgs.eigen.overrideAttrs (old: rec {
    version = "3.3.9";
    stdenv = pkgs.gcc10Stdenv;
    src = pkgs.fetchFromGitLab {
      owner = "libeigen";
      repo = "eigen";
      rev    = "${version}";
      sha256 = "0m4h9fd5s1pzpncy17r3w0b5a6ywqjajmnr720ndb7fc4bn0dhi4";
    };
    # From https://github.com/foolnotion/aoc2020/blob/master/eigen_include_dir.patch
    patches = [ ./eigen_include_dir.patch ];
  });
in  pkgs.stdenv.mkDerivation rec {
  name = "eonclient";
  src = ./client;

  stdenv = pkgs.gcc10Stdenv;
  nativeBuildInputs = with pkgs; [ cmake gfortran ninja ];
  buildInputs = with pkgs; [ eigen339 gtest ];

  cmakeFlags = [
    "-DCMAKE_BUILD_TYPE=Release"
    "-DPACKAGE_TESTS=OFF"
    "-DNO_WARN=TRUE"
    "-DFIND_EIGEN=TRUE"
    # "-DWITH_GPRD=TRUE"
    # "-DWITH_AMS=TRUE"
    "-DUSE_SYSTEM_GTEST=ON"
    "-GNinja"
    ];

  # cmakeDir = "${src}/client";
  dontUseCmakeBuildDir="";
  buildPhase = ''
    ninjaBuildPhase
  '';

  preConfigure = ''
    ./version.sh > version.h
    mkdir build
  '';

    meta = with pkgs.lib; {
    description = "EON C++ client";
    homepage = "https://github.com/HaoZeke/eongit";
    maintainers = [ maintainers.HaoZeke ];
    license = licenses.bsd3;
    platforms = platforms.linux++ platforms.darwin;
  };

}
