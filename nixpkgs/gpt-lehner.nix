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
  use-gitee-wd = if use-gitee == null then false else use-gitee;
  pname = "gpt-lehenr";
  # version = "9bb16fc9fe134623b94116b2a47f7b9288f360df"; # 2025/06/15
  version = "d69e4d0389015b338a2956de806a53bb0354efaf"; # 2025/12/12
in

buildPythonPackage {

  pname = pname;
  version = version;

  pyproject = false;

  src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/gpt" else "https://github.com/lehner/gpt";
    # url = "https://github.com/jinluchang/gpt";
    ref = "master";
    rev = version;
  };

  enableParallelBuilding = true;

  stdenv = stdenv;

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
    if [ "$(uname)" == "Darwin" ] ; then
      export CGPT_EXTRA_LDFLAGS="-undefined dynamic_lookup"
    fi
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
    which python3-config
    python3-config --prefix
    python3-config --includes
    python3-config --ldflags
    # echo
    # cat "$(python3-config --prefix)/bin/python3-config"
    echo
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
    head -n 100 build/logs/*.err
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
