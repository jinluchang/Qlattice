{
  nixpkgs ? import ./nixpkgs.nix { version = "24.05"; },
  ngpu ? null, # adjust with desired number of GPUs. E.g. "2"
  cudaCapability ? null, # adjust with desired cudaCapability. E.g. "8.6"
  cudaForwardCompat ? null, # adjust with desired cudaForwardCompat. E.g. false
}@args:

import ./q-pkgs.nix { inherit nixpkgs ngpu cudaCapability cudaForwardCompat; }
