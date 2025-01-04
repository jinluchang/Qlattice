Common command:

To build without cuda

$ time nix-build -j8 --cores 24

or

$ time nix-build -j8 --cores 24 default.nix

To build with cuda

$ time nix-build -j8 --cores 24 default.nix --argstr name "many-qlat-pkgs-core-w-cuda"

To start a shell without cuda

$ time nix-shell

or

$ time nix-shell shell.nix

To start a shell with cuda

$ time nix-shell shell.nix --argstr name "qlat-jhub-sh-cuda"

To start a shell with only qlat

$ time nix-shell shell.nix --argstr name "qlat-sh"

To start a shell with Filesystem Hierarchy Standard (FHS)

$ time nix-shell shell-fhs.nix

To start a shell with Filesystem Hierarchy Standard (FHS) and cuda

$ time nix-shell shell-fhs.nix --argstr name "qlat-jhub-fhs-cuda"
