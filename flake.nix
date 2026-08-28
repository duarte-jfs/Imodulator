{
  description = "imodulator — reproducible Python environment (uv2nix, built from uv.lock)";

  # CI populates this cache, so `nix develop` fetches the prebuilt env instead of
  # rebuilding gmsh/solcore locally. Only applied for nix trusted-users; everyone
  # else silently falls back to building from source.
  nixConfig = {
    extra-substituters = [ "https://imodulator.cachix.org" ];
    extra-trusted-public-keys = [
      "imodulator.cachix.org-1:nNgtAAmkkOFHmwd9uOEFlraIqvtkIsRvJnTVwSms5wY="
    ];
  };

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

    pyproject-nix = {
      url = "github:pyproject-nix/pyproject.nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    uv2nix = {
      url = "github:pyproject-nix/uv2nix";
      inputs.pyproject-nix.follows = "pyproject-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    pyproject-build-systems = {
      url = "github:pyproject-nix/build-system-pkgs";
      inputs.pyproject-nix.follows = "pyproject-nix";
      inputs.uv2nix.follows = "uv2nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    {
      self,
      nixpkgs,
      uv2nix,
      pyproject-nix,
      pyproject-build-systems,
      ...
    }:
    let
      inherit (nixpkgs) lib;

      # x86_64-linux is the tested target; aarch64-darwin (Apple Silicon) reuses the
      # macOS wheels, which are self-contained (no autoPatchelf/libGLU needed there).
      systems = [
        "x86_64-linux"
        "aarch64-darwin"
      ];
      forAllSystems = lib.genAttrs systems;

      # Load the uv workspace (pyproject.toml + uv.lock) from this repo.
      workspace = uv2nix.lib.workspace.loadWorkspace { workspaceRoot = ./.; };

      # Prefer prebuilt wheels: reuse the manylinux wheels that ship the hard bits already
      # compiled (solcore's Fortran PDD ddModel.so, gmsh's libgmsh, scipy, shapely), so nix
      # never has to build gfortran/f2py/etc. Native .so's that dlopen system libs are fixed
      # up per-package below with autoPatchelfHook.
      overlay = workspace.mkPyprojectOverlay {
        sourcePreference = "wheel";
      };

      # Everything system-specific lives in this one builder, so the file stays flat.
      mkEnv =
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};

          # pyproject.toml pins requires-python = ">=3.11,<3.12".
          python = pkgs.python311;

          # Native libs the prebuilt wheels dlopen at runtime but do NOT bundle. gmsh's
          # libgmsh.so lists all of these as NEEDED (see: ldd on libgmsh.so). In buildInputs
          # they get baked into the wheel's RPATH by autoPatchelfHook, so the env works with
          # no nix-ld / LD_LIBRARY_PATH hacks.
          gmshRuntimeLibs = with pkgs; [
            libGLU
            libGL
            zlib
            fontconfig
            stdenv.cc.cc.lib # libstdc++.so.6, libgomp.so.1
            libx11
            libxcursor
            libxext
            libxfixes
            libxft
            libxinerama
            libxrender
          ];

          # Per-package fixups on top of the uv2nix-generated overlay.
          pyprojectOverrides = final: prev: {
            # femwell declares `meshwell` as a dependency, but neither femwell nor
            # imodulator imports it at runtime (verified by grep). meshwell drags in a
            # heavy, unused closure — cadquery -> {cadquery-ocp (VTK), numba (TBB),
            # trame -> wslink} — whose wheels need extensive native patching. Drop the
            # edge so mkVirtualEnv never pulls that subtree. (pyvista/vtk only enter via
            # nextnanopy, which this env already excludes.)
            #
            # femwell is pinned to a git commit (see the femwell extra in pyproject.toml),
            # so uv2nix builds it from the sdist rather than fetching a PyPI wheel. uv.lock
            # has no build-system metadata for git sources, so the hatchling build backend
            # femwell declares in its own pyproject.toml has to be supplied here explicitly.
            femwell = prev.femwell.overrideAttrs (old: {
              nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ final.resolveBuildSystem {
                hatchling = [ ];
              };
              passthru = (old.passthru or { }) // {
                dependencies = builtins.removeAttrs (old.passthru.dependencies or { }) [ "meshwell" ];
              };
            });

            # solcore seeds ~/.solcore on first import via shutil.copy2, which
            # preserves the source's mode bits. Copied out of the read-only nix store
            # that yields an unwritable ~/.solcore/solcore_config.txt, so the next
            # config write dies with EACCES on any machine that has no ~/.solcore yet.
            # Chmod the copy so the user's own config stays writable.
            solcore = prev.solcore.overrideAttrs (
              old:
              {
                postInstall = (old.postInstall or "") + ''
                  substituteInPlace $out/${python.sitePackages}/solcore/config_tools.py \
                    --replace-fail \
                      'shutil.copy2(self.default_config, self.user_config)' \
                      'shutil.copy2(self.default_config, self.user_config); os.chmod(self.user_config, 0o644)'
                '';
              }
              # solcore wheel: ddModel is a compiled Fortran extension; make libgfortran/
              # libquadmath and the C++ runtime resolvable for its .so. Linux-only --
              # on macOS the wheel resolves its own dylibs via @rpath.
              // lib.optionalAttrs pkgs.stdenv.hostPlatform.isLinux {
                nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [ pkgs.autoPatchelfHook ];
                buildInputs = (old.buildInputs or [ ]) ++ [
                  pkgs.gfortran.cc.lib # libgfortran.so, libquadmath.so
                  pkgs.stdenv.cc.cc.lib
                ];
              }
            );
          }
          # Linux-only native-lib fixups. On macOS the wheels resolve their own dylibs
          # via @rpath (and autoPatchelfHook doesn't exist there), so these are skipped.
          // lib.optionalAttrs pkgs.stdenv.hostPlatform.isLinux {
            # gmsh wheel: patch libGLU/X11/etc. into libgmsh.so's RPATH.
            gmsh = prev.gmsh.overrideAttrs (old: {
              nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [ pkgs.autoPatchelfHook ];
              buildInputs = (old.buildInputs or [ ]) ++ gmshRuntimeLibs;
            });
          };

          # Base set from pyproject.nix's build infra, then layer: build-system backends,
          # the uv2nix overlay from uv.lock, and our native-lib fixups.
          pythonSet = (pkgs.callPackage pyproject-nix.build.packages { inherit python; }).overrideScope (
            lib.composeManyExtensions [
              pyproject-build-systems.overlays.default
              overlay
              pyprojectOverrides
            ]
          );

          # imodulator with the FOSS solver extras (femwell + solcore). nextnanopy is
          # Windows/commercial and deliberately excluded.
          venv = pythonSet.mkVirtualEnv "imodulator-env" {
            imodulator = [
              "femwell"
              "solcore"
            ];
          };

          # devShell-only venv: same solver extras, plus "dev" (ruff, ipykernel — the
          # latter registers this env as the "imodulator_venv" Jupyter kernel the
          # docs/Tutorials notebooks expect). Kept separate from `venv` so `nix build`'s
          # output stays the minimal, production-shaped env.
          devVenv = pythonSet.mkVirtualEnv "imodulator-dev-env" {
            imodulator = [
              "femwell"
              "solcore"
              "dev"
            ];
          };
        in
        {
          inherit pkgs venv devVenv;
        };
    in
    {
      # `nix build`   -> a venv with imodulator + the FOSS solvers (femwell + solcore).
      packages = forAllSystems (system: { default = (mkEnv system).venv; });

      # `nix develop` -> devVenv (solver extras + dev tooling) on PATH, plus uv for lock edits.
      devShells = forAllSystems (
        system:
        let
          inherit (mkEnv system) pkgs devVenv;
        in
        {
          default = pkgs.mkShell {
            packages = [
              devVenv
              pkgs.uv
            ];
            env = {
              # Reproducible env from uv.lock — don't let uv sync/download a 2nd interpreter.
              UV_NO_SYNC = "1";
              UV_PYTHON = "${devVenv}/bin/python";
              UV_PYTHON_DOWNLOADS = "never";
            };
            shellHook = ''
              echo "imodulator uv2nix env — python: $(python --version)"
            '';
          };
        }
      );
    };
}
