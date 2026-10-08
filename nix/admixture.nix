{ lib, stdenv, fetchurl }:

stdenv.mkDerivation (finalAttrs: {
  pname = "admixture";
  version = "1.4.0";

  src = fetchurl {
    url = "https://dalexander.github.io/admixture/binaries/admixture_linux-${finalAttrs.version}.tar.gz";
    hash = "sha256-VSDE8YUl9JqyxxXXEbQPSggY7xSFNhqU7IlU0TZObxk=";
  };

  sourceRoot = "admixture_linux-${finalAttrs.version}";
  dontConfigure = true;
  dontBuild = true;

  installPhase = ''
    runHook preInstall
    install -Dm755 admixture $out/bin/admixture
    runHook postInstall
  '';

  meta = {
    description = "Maximum likelihood estimation of individual ancestries from SNP genotype datasets";
    homepage = "https://dalexander.github.io/admixture/";
    license = lib.licenses.gpl2Only;
    platforms = [ "x86_64-linux" ];
    mainProgram = "admixture";
  };
})
