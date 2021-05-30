let
  sources = import ./nix/sources.nix;
  pkgs = import sources.nixpkgs { };
    eigen339 = pkgs.eigen.overrideAttrs(old: rec {
      version = "3.3.9";
    });
in  pkgs.stdenv.mkDerivation rec {
  name = "eonclient";
  src = ./client;
  stdenv = pkgs.gcc10Stdenv;

  cmakeFlags = [
    "-DCMAKE_BUILD_TYPE=Release"
    "-DPACKAGE_TESTS=OFF"
    "-DNO_WARN=TRUE"
    "-DFIND_EIGEN=TRUE"
    ];

  preConfigure = ''
    pwd
    ls -R
    ./version.sh > version.h
  '';

  buildInputs = with pkgs; [ cmake gfortran eigen339 ];

    meta = with pkgs.lib; {
    description = "EON C++ client";
    homepage = "https://github.com/HaoZeke/eongit";
    maintainers = [ maintainers.HaoZeke ];
    license = licenses.bsd3;
    platforms = platforms.linux++ platforms.darwin;
  };

}
