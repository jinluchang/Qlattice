{
  version ? null,
  use-gitee ? null,
}:

let
  version-wd = if version == null then "25.05" else version;
  use-gitee-wd = if use-gitee == null then false else use-gitee;
  url = if use-gitee-wd then
  "https://mirrors.ustc.edu.cn/nix-channels/nixos-${version-wd}/nixexprs.tar.xz"
  # "https://mirror.nju.edu.cn/nix-channels/nixos-${version-wd}/nixexprs.tar.xz"
  # "https://mirrors.tuna.tsinghua.edu.cn/nix-channels/nixos-${version-wd}/nixexprs.tar.xz"
  else
  "https://channels.nixos.org/nixos-${version-wd}/nixexprs.tar.xz";
in
  fetchTarball url
