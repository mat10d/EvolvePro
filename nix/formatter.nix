{ inputs, ... }:
{
  imports = [ inputs.treefmt-nix.flakeModule ];
  perSystem.treefmt = {
    projectRootFile = "flake.nix";
    programs = {
      nixfmt.enable = true;
      deadnix.enable = true;
      statix.enable = true;
      keep-sorted.enable = true;
      ruff-check.enable = true;
      ruff-format.enable = true;
    };
  };
}
