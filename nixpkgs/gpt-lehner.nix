{ stdenv
, config
, lib
, fetchFromGitHub
, python
, buildPythonPackage
, mpi
, grid-lehner
, numpy
, pkg-config
, git
, which
, flock
, rsync
, cudaSupport ? config.cudaSupport
, cudaPackages ? {}
}:

let
  orig-stdenv = stdenv;
in

buildPythonPackage rec {

  pname = "gpt-lehenr";
  version = "28efde3d7ef7ef6417abbf4bc2740120f7064aee";

  pyproject = false;

  src = fetchFromGitHub {
    owner = "lehner";
    repo = "gpt";
    rev = version;
    hash = "sha256-3w1D7C3VSrtJNVPfZHHT1mOT5U1EYuajvaWc4LhCDYw=";
  };

  enableParallelBuilding = true;

  stdenv = if cudaSupport then cudaPackages.backendStdenv else orig-stdenv;

  build-system = [
    pkg-config
    numpy
  ];

  nativeBuildInputs = [
    mpi
    git
    grid-lehner
    which
    flock
    rsync
  ]
  ++ lib.optionals cudaSupport (with cudaPackages; [ cuda_nvcc ])
  ;

  propagatedBuildInputs = [
    grid-lehner
  ];

  dependencies = [
    mpi
    numpy
  ];

  preConfigure = ''
	export OMPI_CXX=c++
	export OMPI_CC=cc
    #
    export
    #
    cd "$NIX_BUILD_TOP/source/lib/cgpt"
    pwd
    #
    echo "BASIS_SIZE(4)" > lib/basis_size.h
    echo "BASIS_SIZE(10)" >> lib/basis_size.h
    echo "BASIS_SIZE(25)" >> lib/basis_size.h
    #
    echo "SPIN(4)" > lib/spin_color.h
    echo "COLOR(3)" >> lib/spin_color.h
    echo "COLOR(2)" >> lib/spin_color.h
    echo "COLOR(1)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,3)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,2)" >> lib/spin_color.h
    echo "SPIN_COLOR(4,1)" >> lib/spin_color.h
    #
    cd "$NIX_BUILD_TOP/source"
    pwd
  '';

  preBuild = ''
    cd "$NIX_BUILD_TOP/source/lib/cgpt"
    pwd
    ls -l
    #
    patchShebangs --build ./clean ./update
    #
    echo '-- clean source...'
    bash ./clean
    echo '-- update source...'
    sed -i 's/git add /echo /' ./update
    bash ./update
    #
    echo '-- update source...'
    which grid-config
    bash ./make %grid-config "$NIX_BUILD_CORES"
    #
    cat build/logs/basis.err
    #
    cd "$NIX_BUILD_TOP/source"
    pwd
  '';

  preInstall = ''
    cd "$NIX_BUILD_TOP/source"
    pwd
    mkdir -pv "$out"/${python.sitePackages}
    mkdir -pv "$out"/src
    rsync -a --delete lib/gpt "$out"/${python.sitePackages}/
    rsync -a --delete lib/cgpt/build/cgpt.so "$out"/${python.sitePackages}/
    rsync -a --delete tests "$out"/src/
    rsync -a --delete applications "$out"/src/
    rsync -a --delete benchmarks "$out"/src/
    rsync -a --delete documentation "$out"/src/
    rsync -a --delete docker "$out"/src/
  '';

}
