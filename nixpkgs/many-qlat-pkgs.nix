{
  version-list ? null,
  ngpu ? null,
  cudaCapability ? null,
  cudaForwardCompat ? null,
  use-gitee ? null,
}:

let

  version-list-ini = version-list;

in let

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
    name-list = builtins.map (n: "q-pkgs${n}") q-pkgs.qlat-name-list;
    ver-suffix = if version == "" then "" else "-${version}";
    mk-name-entries = name: {
      "${name}${ver-suffix}-qlat-tests" = q-pkgs.${name}.qlat-tests;
      "${name}${ver-suffix}-qlat-env" = q-pkgs.${name}.qlat-env;
    };
  in builtins.foldl' (s: v: s // v) {} (builtins.map mk-name-entries name-list);

in builtins.foldl' (s: v: s // v) {} (builtins.map mk-version-entries version-list)
