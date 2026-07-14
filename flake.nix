{
  description = "Stitch microscope tiles and generate a standalone HTML viewer";

  inputs = {
    nixpkgs.url = "nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };

        python = pkgs.python312.withPackages (ps: with ps; [
          numpy
          tifffile
          pillow
          pip
        ]);
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [
            python
            pkgs.fiji
          ];

          shellHook = ''
            echo ""
            echo "  microscope2html development shell (Python 3.12)"
            echo ""
            echo "  Run the CLI directly:"
            echo "    python -m microscope2html.cli --help"
            echo ""
            echo "  Or create a virtual environment (recommended for pip install):"
            echo "    python -m venv .venv && source .venv/bin/activate && pip install -e ."
            echo ""
          '';
        };
      }
    );
}