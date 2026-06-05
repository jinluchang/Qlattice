{
  nixpkgs ? null, # nixpkgs. E.g. (fetchTarball "https://channels.nixos.org/nixos-25.11/nixexprs.tar.xz")
  version ? null, # version of the nixpkgs. E.g. "25.11"
  ngpu ? null, # adjust with desired number of GPUs. E.g. "2"
  cudaCapability ? null, # adjust with desired cudaCapability. E.g. "8.6"
  cudaForwardCompat ? null, # adjust with desired cudaForwardCompat. E.g. false
  use-gitee ? null, # true or false (default false)
}:

let

  version-wd = if version == null then "" else version;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
  nixpkgs-default = import ./nixpkgs.nix {
    version = version-wd;
    use-gitee = use-gitee-wd;
  };
  nixpkgs-wd = if nixpkgs == null then nixpkgs-default else nixpkgs;
  import-nixpkgs-wd = import nixpkgs-wd;

  # TODO: copy to q-pkgs.nix (and use it there when needed)
  force = x: builtins.deepSeq x x;

  # TODO: rename to pkgs (and do not return it)
  o-pkgs = import-nixpkgs-wd {};

  lib = o-pkgs.lib;

  nixgl-src = builtins.fetchGit {
    url = if use-gitee-wd then "https://gitee.com/jinluchang/nixGL" else "https://github.com/jinluchang/nixGL";
    ref = "main";
    rev = "28366884d82a3f471b471eea70baa2e6668d6b4e";
  };

  # TODO: move to q-pkgs.nix
  is-linux = (lib.lists.elem builtins.currentSystem lib.platforms.linux);

  options-default = let
    #
    pkgs = import-nixpkgs-wd {
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
    nixgl = (import nixgl-src { pkgs = pkgs; }).auto.nixGLDefault;
    #
    ngpu-sys = if ngpu != null
    then ngpu
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
    else if cudaCapability != null
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
    cudaForwardCompat-sys = if cudaForwardCompat != null
    then cudaForwardCompat
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

in {
  # TODO: For all args, return processed value with its original name (i.e., without -wd suffix)
  inherit version-wd use-gitee-wd nixpkgs-wd import-nixpkgs-wd;
  inherit nixgl-src is-linux lib;
  inherit o-pkgs;
  inherit options-default;
  inherit mk-options mk-qlat-name;
}
