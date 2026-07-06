# cpp/src/dependencies — HiGHS 1.15.0 (built from source)

The C++ MAF solver links the [HiGHS](https://highs.dev) LP/MIP library. It is **built from
source** here rather than depending on a machine-local Julia `HiGHS_jll` binary. This makes the
build self-contained, redistributable, and free of the `libblastrampoline` / BLAS dependency the
jll binary carried (HIPO/OpenBLAS are disabled, and our code uses only the simplex/MIP path).

**HiGHS 1.15.0 is the sole supported version.**

## Layout

```
dependencies/
  README.md                 this file
  .gitignore                ignores generated HiGHS-1.15.0/build|install
  HiGHS-1.15.0/             the vendored HiGHS source tree (contains all license files)
```

The HiGHS 1.15.0 source tree is **vendored in the repo** (committed under `HiGHS-1.15.0/`), so the
build needs no network.

## Verifying the vendored source

To confirm `HiGHS-1.15.0/` matches the official upstream release, run from `cpp/`:

```bash
bash cleanup.sh
```

It removes build artifacts, fetches `~/HiGHS-1.15.0.zip` (from the GitHub release tag `v1.15.0`)
if not already cached, extracts it to `/tmp`, and diffs it file-by-file against `HiGHS-1.15.0/`.

## Building

`cpp/CMakeLists.txt` builds HiGHS from `HiGHS-1.15.0/` as a subproject (`add_subdirectory`), into
the CMake build dir, and links the solver against it. `docker_setup.sh` instead builds HiGHS once
into `.highs/` (install prefix `HIGHS_INSTALL_DIR`) so the per-track builds link it without recompiling.

## Parity note

The original Julia reference used HiGHS **1.14.0**. Building **1.15.0** from source deliberately
relaxes bit-identical LP/MIP equivalence with that reference — a validated trade-off (same optima,
comparable speed) for a clean, dependency-light, redistributable build.

## Licensing

HiGHS is **MIT-licensed**; the full text ships inside the extracted tree at
`HiGHS-1.15.0/LICENSE.txt`. HiGHS bundles several third-party components under
`HiGHS-1.15.0/extern/`, each retaining its own license file. See
[`../../THIRD_PARTY_NOTICES.md`](../../THIRD_PARTY_NOTICES.md) for the consolidated acknowledgment.
