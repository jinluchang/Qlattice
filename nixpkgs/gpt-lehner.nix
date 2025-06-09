{ stdenv
, config
, lib
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
, use-gitee ? null
}:

let
  orig-stdenv = stdenv;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
in

buildPythonPackage rec {

  pname = "gpt-lehenr";
  version = "861e8ee38beb6c97638cf8fdc1a024b04618961f";

  pyproject = false;

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/gpt" else "https://github.com/lehner/gpt";
    # url = "https://github.com/jinluchang/gpt";
    ref = "master";
    rev = version;
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
    rsync -a --delete lib/gpt "$out"/${python.sitePackages}/
    rsync -a --delete lib/cgpt/build/cgpt.so "$out"/${python.sitePackages}/
    mkdir -pv "$out"/share/${pname}
    rsync -a --delete tests "$out"/share/${pname}/
    rsync -a --delete applications "$out"/share/${pname}/
    rsync -a --delete benchmarks "$out"/share/${pname}/
    rsync -a --delete documentation "$out"/share/${pname}/
    rsync -a --delete docker "$out"/share/${pname}/
  '';

  postFixup = ''
    mkdir -pv "$out"/share/version
    echo ${version} >"$out"/share/version/${pname}
  '';

}
