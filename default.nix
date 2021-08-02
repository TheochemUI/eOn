{ sources ? import ./nix/sources.nix }:
let
  pkgs = import sources.nixpkgs { };
    eigen339 = pkgs.eigen.overrideAttrs(old: rec {
      version = "3.3.9";
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
    "-DWITH_GPRD=TRUE"
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
