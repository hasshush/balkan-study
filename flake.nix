{
  description = "poreklo - genetic structure of Western Balkan populations";

  inputs.nixpkgs.url = "nixpkgs";

  outputs =
    { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };
      admixture = pkgs.callPackage ./nix/admixture.nix { };
      python = pkgs.python3.withPackages (
        ps: with ps; [
          numpy
          scipy
          pandas
          scikit-learn
          matplotlib
          seaborn
          statsmodels
          jupyterlab
          ipykernel
          openpyxl
        ]
      );
    in
    {
      packages.${system} = {
        inherit admixture;
        default = admixture;
      };

      devShells.${system}.default = pkgs.mkShell {
        packages = [
          python
          admixture
          pkgs.plink2
          pkgs.bcftools
          pkgs.htslib
          pkgs.curl
          pkgs.wget
          pkgs.unzip
          pkgs.jq
        ];

        shellHook = ''
          export POREKLO_ROOT="$PWD"
          echo "poreklo dev shell"
          echo "  plink2:   $(plink2 --version 2>&1 | head -1)"
          echo "  admixture: $(admixture --version 2>&1 | grep -o 'Version [0-9.]*' | head -1)"
          echo "  python:   $(python --version 2>&1)"
        '';
      };
    };
}
