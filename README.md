# MAF (C++) — Maximum Agreement Forest solver

A branch-and-price solver for the Maximum Agreement Forest problem (PACE 2026). One codebase builds
any of the three competition tracks (`exact`, `heuristic`, `lowerbound`) by selecting the `MAF_TRACK`
CMake option. It links the solver [HiGHS](https://highs.dev), built from vendored source.

## Build

**PACE 2026 organizers:** run `bash docker_setup.sh` (or extract the commands from it) inside a debian
13.5 docker that mounts the repository. This will automatically install the toolchain, build HiGHS once,
and build **all three** track binaries under `build/<track>/main_<track>`.

**Manual build** of a single track (`<track>` = `exact` | `heuristic` | `lowerbound`):

```bash
cmake -S . -B build/<track> -DCMAKE_BUILD_TYPE=Release -DMAF_TRACK=<track>
cmake --build build/<track> -j --target maf
```

HiGHS 1.15.0 is vendored in `src/dependencies/`, so the compile needs no network. CMake builds it
from that source; pass `-DHIGHS_INSTALL_DIR="$PWD/.highs"` to link a prebuilt copy (as
`docker_setup.sh` does) instead of recompiling it.

## Run

Each binary reads one instance from **stdin** (PACE Newick; `#`-lines are metadata) and writes the
forest to **stdout**. 

The solver can be called as follows:
```bash
./build/<track>/main_<track> < instance.nw
```


Using `MAF_TIMELIMIT=<seconds>` one can override the competition's default per-track budget:

```bash
MAF_TIMELIMIT=1800 ./build/exact/main_exact < instance.nw
```

| track | default budget | output | on interrupt / timeout |
|---|---|---|---|
| exact | 1800 s | Newick forest, proven optimal | emit **nothing** |
| heuristic | 300 s | Newick forest, best incumbent | on SIGTERM/SIGINT, flush the best valid incumbent |
| lowerbound | 600 s | Newick forest of size provenly not greater than ⌊a·k*⌋+b (reads the `#a a b` line) | emit **nothing** |

## Cleaning up

```bash
bash cleanup.sh                             # remove build artifacts, then verify the vendored HiGHS source
bash cleanup.sh --ignore-highs-verification # remove build artifacts only
```

## Licensing

The MAF code is MIT-licensed — see [`LICENSE.md`](LICENSE.md). Bundled HiGHS and its vendored
components are redistributed under their own permissive licenses — see
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) and
[`src/dependencies/README.md`](src/dependencies/README.md).
