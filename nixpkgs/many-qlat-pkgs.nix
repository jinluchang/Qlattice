{
  # List of nixpkgs versions to build for, e.g. ["" "26.05" "25.11"].
  # "" means the default (unversioned) nixpkgs. null uses a built-in default.
  version-list ? null,
  # Optional filter: only build qlat variants whose names are in this list.
  # null means build all variants from q-pkgs.
  qlat-name-list ? null,

  # Controls which CUDA-related test variants are built.
  # null (default): build all. "none": exclude CUDA. "only": only CUDA (no env).
  qlat-cuda-tests ? null,

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

  assert builtins.elem qlat-cuda-tests [ null "none" "only" ];

  let

  is-cuda-name = name: builtins.match ".*-cuda.*" name != null;
  include-tests = if qlat-cuda-tests == null then (name: true)
    else if qlat-cuda-tests == "none" then (name: ! is-cuda-name name)
    else if qlat-cuda-tests == "only" then (name: is-cuda-name name)
    else builtins.throw "qlat-cuda-tests must be null, \"none\", or \"only\", got: ${builtins.toString qlat-cuda-tests}";
  include-env = qlat-cuda-tests != "only";

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
