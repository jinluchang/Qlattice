{
  # List of nixpkgs versions to build for, e.g. ["" "26.05" "25.11"].
  # "" means the default (unversioned) nixpkgs. null uses a built-in default.
  version-list ? null,
  # Optional filter: only build qlat variants whose names are in this list.
  # null means build all variants from q-pkgs.
  qlat-name-list ? null,

  # Which test packages to build, per selected qlat name:
  #   null         build the env and the test packages (no filter)
  #   "none"       build no test package (env packages only)
  #   "none-cuda"  build the test packages of every name but the '-cuda*' ones
  #   "only"       build only the test packages (no env packages)
  #   "only-cuda"  build only the test packages of the '-cuda*' names (no env)
  qlat-tests ? null,

  # Number of GPUs to target (passed to q-pkgs). null uses the q-pkgs default.
  ngpu ? null,
  # CUDA compute capability, e.g. "8.6" (passed to q-pkgs).
  cudaCapability ? null,
  # Enable CUDA forward compatibility (passed to q-pkgs).
  cudaForwardCompat ? null,
  # Use gitee mirrors instead of GitHub (passed to q-pkgs).
  use-gitee ? null,
}:

let

  version-list-ini = version-list;

in

  assert builtins.elem qlat-tests [ null "none" "none-cuda" "only" "only-cuda" ];

  let

  is-cuda-name = name: builtins.match ".*-cuda.*" name != null;
  include-tests = if qlat-tests == null then (name: true)
    else if qlat-tests == "none" then (name: false)
    else if qlat-tests == "none-cuda" then (name: ! is-cuda-name name)
    else if qlat-tests == "only" then (name: true)
    else if qlat-tests == "only-cuda" then (name: is-cuda-name name)
    else builtins.throw "qlat-tests must be null, \"none\", \"none-cuda\", \"only\" or \"only-cuda\", got: ${builtins.toString qlat-tests}";
  include-env = qlat-tests != "only" && qlat-tests != "only-cuda";

  version-list = if version-list-ini != null
  then version-list-ini
  else [
    ""
    "26.05"
    "25.11"
  ];

  mk-version-entries = version: let
    q-pkgs = import ./q-pkgs.nix {
      inherit version ngpu cudaCapability cudaForwardCompat use-gitee;
    };
    name-list = builtins.map (n: "q-pkgs${n}") (
      if qlat-name-list != null
      then builtins.filter (n: builtins.elem n qlat-name-list) q-pkgs.qlat-name-list
      else q-pkgs.qlat-name-list
    );
    ver-suffix = if version == "" then "" else "-${builtins.replaceStrings ["."] ["-"] version}";
    mk-name-entries = name: (if include-env then {
      "${name}${ver-suffix}-qlat-env" = q-pkgs.${name}.qlat-env;
    } else {}) // (if include-tests name then {
      "${name}${ver-suffix}-qlat-tests" = q-pkgs.${name}.qlat-tests;
    } else {});
  in builtins.foldl' (s: v: s // v) {} (builtins.map mk-name-entries name-list);

in builtins.foldl' (s: v: s // v) {} (builtins.map mk-version-entries version-list)
