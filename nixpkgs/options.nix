{
  # Path or fetchTarball result to nixpkgs. Auto-detected from `version` if null.
  # E.g. (fetchTarball "https://channels.nixos.org/nixos-26.05/nixexprs.tar.xz")
  nixpkgs ? null,
  # NixOS release version string. Used to fetch matching nixpkgs when `nixpkgs` is null.
  # E.g. "26.05"
  version ? null,
  # Number of NVIDIA GPUs as a string. Auto-detected from /dev/nvidia* if null.
  # E.g. "2"
  ngpu ? null,
  # CUDA compute capability string (major.minor). Auto-detected via nvidia-smi if null.
  # E.g. "8.6"
  cudaCapability ? null,
  # Enable CUDA forward compatibility (bool). Defaults to false if null.
  cudaForwardCompat ? null,
  # Use gitee mirrors instead of github (bool). Defaults to false if null.
  use-gitee ? null,
}:

let

  nixpkgs-ini = nixpkgs;
  version-ini = version;
  ngpu-ini = ngpu;
  cudaCapability-ini = cudaCapability;
  cudaForwardCompat-ini = cudaForwardCompat;
  use-gitee-ini = use-gitee;

in let

  version = if version-ini == null then "" else version-ini;
  use-gitee = if use-gitee-ini == null then false else use-gitee-ini;
  nixpkgs-default = import ./nixpkgs.nix {
    version = version;
    use-gitee = use-gitee;
  };
  nixpkgs = if nixpkgs-ini == null then nixpkgs-default else nixpkgs-ini;
  import-nixpkgs = import nixpkgs;

  force = x: builtins.deepSeq x x;

  pkgs = import-nixpkgs {};

  lib = pkgs.lib;

  nixpkgs-release = lib.trivial.release;

  nixgl-src = builtins.fetchGit {
    url = if use-gitee then "https://gitee.com/jinluchang/nixGL" else "https://github.com/jinluchang/nixGL";
    ref = "main";
    rev = "28366884d82a3f471b471eea70baa2e6668d6b4e";
  };

  import-nixgl = import nixgl-src;

  options-default = let
    #
    pkgs = import-nixpkgs {
      config.allowUnfree = true;
    };
    #
    runCommandLocal = pkgs.runCommandLocal;
    #
    cpuinfo-sys = builtins.readFile (runCommandLocal
    "impure-cpuinfo-file"
    {
      time = builtins.currentTime;
    }
    ''
      echo "cpuinfo="
      echo "$(grep '^flags' /proc/cpuinfo 2>/dev/null | head -n 1)" >$out
      cat $out
    ''
    );
    #
    nixgl = (import-nixgl { pkgs = pkgs; }).auto.nixGLDefault;
    #
    ngpu-sys = if ngpu-ini != null
    then ngpu-ini
    else builtins.head (builtins.match
    "(.*)\n"
    (builtins.readFile (runCommandLocal
    "impure-ngpu-file"
    {
      time = builtins.currentTime;
    }
    ''
      mkdir tmp
      cd tmp
      ls /dev/nvidia{?,??} 2>/dev/null | wc -l >$out 2>/dev/null || echo "0" >$out
      echo "ngpu=$(cat $out)"
    ''
    )));
    #
    nvidia_x11_bin = if ngpu-sys == "0"
    then null
    else pkgs.linuxPackages.nvidia_x11.bin;
    #
    cudaCapability-sys = if ngpu-sys == "0"
    then null
    else if cudaCapability-ini != null
    then null
    else builtins.head (builtins.match
    "(.*)\n"
    (builtins.readFile (runCommandLocal
    "impure-cuda-capability-file"
    {
      time = builtins.currentTime;
    }
    ''
      ${nixgl}/bin/nixGL ${nvidia_x11_bin}/bin/nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null \
        | head -n 1 >$out 2>/dev/null \
        || echo >$out
      echo "cudaCapability=$(cat $out)"
    ''
    )));
    #
    cudaCapabilities-sys = if cudaCapability-sys == null
    then []
    else [ cudaCapability-sys ];
    #
    cudaForwardCompat-sys = if cudaForwardCompat-ini != null
    then cudaForwardCompat-ini
    else false;
    #
    get-nvcc-arch-from-cudaCapability = cudaCapability:
    "sm_" + builtins.replaceStrings [ "." ] [ "" ] cudaCapability;
    #
    nvcc-arch-sys = if cudaCapability-sys == null
    then null
    else get-nvcc-arch-from-cudaCapability cudaCapability-sys;
    #
    opts = {
      qlat-name = "";
      cpuinfo-sys = cpuinfo-sys;
      ngpu = ngpu-sys;
      nvcc-arch = nvcc-arch-sys;
      cudaCapabilities = cudaCapabilities-sys;
      cudaForwardCompat = false;
      use-cuda-software = false;
      use-grid-gpt = true;
      use-cps = true;
      use-cuda = false;
      use-cudasupport = false;
      use-cubaquad = true;
      use-clang = false;
      use-ucx = true;
      use-pypi = null;
    };
    #
  in force opts;

  mk-options = options:
  # Order the precedence:
  # options from arguments of this function.
  # options from arguments of this entire `q-pkgs.nix` file.
  # options-default from system parameters.
  let
    opts-0 = options-default // options;
    corrections-1 = opts: opts // {
      ${if opts.use-grid-gpt then "use-cubaquad" else null} = true;
    };
    corrections-2 = opts: opts // {
      ${if opts.use-cudasupport then "use-cuda" else null} = true;
    };
    corrections-3 = opts: opts // {
      ${if opts.use-cuda then "use-cuda-software" else null} = true;
      ${if opts.use-cuda then "use-clang" else null} = false;
    };
    corrections-4 = opts: opts // {
      ${if ! opts.use-cuda-software then "ngpu" else null} = 0;
      ${if ! opts.use-cuda-software then "nvcc-arch" else null} = null;
      ${if ! opts.use-cuda-software then "cudaCapabilities" else null} = [];
    };
    corrections-5 = opts: opts // {
      ${if opts.ngpu == "0" then "use-cuda-software" else null} = false;
      ${if opts.ngpu == "0" then "use-cuda" else null} = false;
      ${if opts.ngpu == "0" then "use-cudasupport" else null} = false;
    };
    opts = corrections-5 (corrections-4 (corrections-3 (corrections-2 (corrections-1 opts-0))));
    opts-1 =
      assert (opts.use-cudasupport -> opts.use-cuda);
      assert (opts.use-cuda -> opts.use-cuda-software);
      assert (opts.use-cuda-software -> opts.nvcc-arch != null);
      assert (opts.nvcc-arch != null -> opts.cudaCapabilities != []);
      assert (opts.use-grid-gpt -> opts.use-cubaquad);
      opts;
  in opts-1;

  mk-qlat-name = options:
  let
    opts = options;
  in (opts.qlat-name
  + lib.optionalString (! opts.use-grid-gpt && ! opts.use-cps) "-std"
  + lib.optionalString (opts.use-grid-gpt && ! opts.use-cps) "-cpsless"
  + lib.optionalString (! opts.use-grid-gpt && opts.use-cps) "-gridless"
  + lib.optionalString opts.use-cuda-software "-cu"
  + lib.optionalString opts.use-cuda "da"
  + lib.optionalString opts.use-cudasupport "support"
  + lib.optionalString (! opts.use-cubaquad) "-cubaquadless"
  + lib.optionalString opts.use-clang "-clang"
  + lib.optionalString (! opts.use-ucx) "-ucxless"
  + lib.optionalString (opts.use-pypi != null) "-pypi"
  );

  version-pypi = "1.7";

  options-list = [
    {}
    { use-cps = false; use-grid-gpt = false; }
    { use-cps = false; }
    { use-grid-gpt = false; }
    { use-cuda-software = true; }
    { use-cuda = true; }
    { use-cudasupport = true; }
    { use-ucx = false; }
    { use-clang = true; }
    { use-pypi = version-pypi; }
    #
    { use-clang = true; use-ucx = false; }
    { use-cuda = true; use-ucx = false; }
    { use-grid-gpt = false; use-cubaquad = false; }
    { use-grid-gpt = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-ucx = false; }
    { use-cps = false; use-grid-gpt = false; use-clang = true; use-ucx = false; }
    #
    { use-cps = false; use-ucx = false; }
    { use-cps = false; use-clang = true; use-ucx = false; }
    { use-cps = false; use-clang = true; }
    #
    { use-cps = false; use-grid-gpt = false; use-cuda-software = true; }
    { use-cps = false; use-grid-gpt = false; use-cuda = true; }
    { use-cps = false; use-grid-gpt = false; use-cudasupport = true; }
  ];

  qlat-name-list = let r = lib.lists.unique (builtins.map mk-qlat-name (builtins.map mk-options options-list)); in builtins.deepSeq r r;

  qlat-name-list-file-from-str = builtins.toFile "qlat-name-list"
  (builtins.foldl' (s: v: s + "${v}\n") "" qlat-name-list);

  qlat-name-list-file = pkgs.runCommand
  "qlat-name-list"
  {}
  ''
    cp -v ${qlat-name-list-file-from-str} $out
  '';

in {
  inherit nixpkgs version ngpu cudaCapability cudaForwardCompat use-gitee;
  inherit options-default mk-options mk-qlat-name;
  inherit import-nixpkgs import-nixgl;
  inherit nixpkgs-release;
  inherit version-pypi;
  inherit options-list qlat-name-list qlat-name-list-file-from-str qlat-name-list-file;
}
