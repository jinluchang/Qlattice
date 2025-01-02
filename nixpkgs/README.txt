Common command:

To build without cuda

$ time nix-build -j8 --cores 24

or

$ time nix-build -j8 --cores 24 default.nix

To build with cuda

$ time nix-build -j8 --cores 24 default-cuda.nix

To start a shell without cuda

$ time nix-shell

or

$ time nix-shell shell.nix

To start a shell with cuda

$ time nix-shell shell-cuda.nix

To start a shell with many other packages

$ time nix-shell shell-jhub.nix
