# Third-Party Notices

The MAF solver links the **HiGHS** optimization library, which it builds from
source (see `src/dependencies/README.md`). HiGHS and its bundled dependencies are redistributed
under their respective permissive licenses. Full, unmodified license texts ship inside the
extracted HiGHS source tree; the paths below are relative to
`src/dependencies/HiGHS-1.15.0/`.

## HiGHS

- **License:** MIT — `LICENSE.txt`
- **Upstream:** https://github.com/ERGO-Code/HiGHS
- **Version:** 1.15.0, vendored in the repo.
- **Build configuration:** `FAST_BUILD=ON`, `HIPO=OFF`, `BUILD_OPENBLAS=OFF`,
  `BUILD_TESTING=OFF`, `BUILD_EXAMPLES=OFF`. Interior-point/PDLP paths that need BLAS are
  disabled, so the resulting `libhighs` has **no `libblastrampoline` / BLAS dependency**.

## Third-party components bundled with and linked into HiGHS

These live under `extern/` in the HiGHS source tree and are compiled into `libhighs`:

| Component | License | License file (under `extern/`) |
|-----------|---------|--------------------------------|
| pdqsort   | zlib    | `pdqsort/license.txt` |
| AMD (SuiteSparse ordering) | BSD-3-Clause | `amd/License.txt` |
| METIS     | Apache-2.0 | `metis/LICENSE.txt` |
| RCM       | MIT     | `rcm/LICENSE` |
| zstr      | MIT     | `zstr/LICENSE` |

Also present in the source tree but **not** compiled into our build (test/CLI/optional-BLAS
only, disabled by the flags above): Catch2 (`extern/catch`, Boost Software License 1.0),
CLI11 (`extern/cli11`, BSD-3-Clause), and the reference BLAS shim (`extern/blas`).

## No BLAS runtime dependency

Building HiGHS from source with the flags above avoids a BLAS runtime dependency, so
`libblastrampoline` is not shipped or linked.
