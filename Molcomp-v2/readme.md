# Molcomp-v2 (Pure C11 / OpenMP Implementation)

> [!NOTE]
> A fast, standalone, dependency-free **C11 / OpenMP** implementation of **Molcomp** (pairwise 3P-KHG conformer comparison) and **Moldedupl** (high-throughput parallel conformer deduplication).  
> For exact numerical reproduction of the benchmarks published in the paper, please use the original prototype in the root directory.

---

## Key Improvements over v1

- **Zero External Dependencies:** Replaces GSL (`libgsl`, `libgslcblas`) with built-in stack-allocated Jacobi $3 \times 3$ and $4 \times 4$ eigensolvers.
- **Native C / OpenMP `moldedupl` Binary:** Eliminates the Python `ctypes` wrapper, temporary file splitting (`tmp/`), and process-forking overhead.
- **Precomputed Molecular Topology (`McPrecomp`):** Connectivity graphs, radial BFS equivalence matrices (`IDN`), and torsion lists are computed once per ensemble rather than per pair.
- **Memory-Efficient Pair Caching:** Uses an open-addressing 64-bit hash set instead of a dense $N \times N$ bit matrix.
- **Deterministic Exact Post-Pass & Reliability Fixes:**
  - Adds canonical sorting of duplicate pairs and a deterministic exact post-pass over remaining pairs within $\Delta E \le 2.0\text{ kcal/mol}$ (or all pairs when energies are absent).
  - Completes the final outer-layer comparison in radial BFS (`MC_RADIAL`) and remaps neighbor/torsion indices after Hungarian atom reordering.
  - Initializes atom structures cleanly and validates neighbor bounds and element sets.

---

## Benchmark Summary (3P-KHG `angle` mode)

| System | Input Conformers | v1 (Python + GSL) | `Molcomp-v2` (Pure C11 / OpenMP) | Speedup | Notes |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Vitamin A** | 18,245 | 21.79 s (`138` unique) | **2.52 s** (`138` unique) | **8.6×** | Identical unique count; `sys` time reduced from 4.24 s to 0.07 s |
| **Vitamin B2** | 12,148 | 1m 01.13s (`2,120` unique) | **16.71 s** (`2,119` unique) | **3.7×** | Exact chunk match (`2,120`) + 1 extra duplicate found in exact pass |
| **Chloroquine** | 32,048 | 3m 07.00s (`4,260` unique) | **1m 04.00s** (`4,249` unique) | **2.9×** | 11 extra duplicates found in deterministic exact pass |
| **Vitamin B9** | 51,857 | 7m 23.90s (`9,253` unique) | **2m 44.78s** (`9,177` unique) | **2.7×** | `sys` time reduced from 4m 02s to 0.31s; 76 extra duplicates found |

---

## Compilation

Requires only `gcc` (or `clang`) with OpenMP support:

```bash
chmod +x compile.sh
./compile.sh
```

This builds `molcomp` and `moldedupl` with native CPU optimizations (`-O3 -march=native -fopenmp`).

---

## Usage

### 1. Pairwise Comparison (`molcomp`)

```bash
./molcomp struct1.xyz struct2.xyz --compare=angle
```

Outputs `RMSD (Å)  MaxDist (Å)  MaxTors (°)`. All original CLI options (`--alg=std|adv|kabsch`, `--keepenantio=true|false`, `--connect=true|false`, `--readtopo=true`, `--mol2=true`, `--print=true`) are supported.

### 2. Ensemble Deduplication (`moldedupl`)

**Single multi-frame `.xyz` file:**
```bash
./moldedupl conformers.xyz angle [readtopo]
```
Generates `conformers-unique-angle.xyz`, `conformers-dup-angle.txt`, and `conformers-DUPXYZ-angle.txt`.

**Directory of `.xyz` files:**
```bash
./moldedupl angle [readtopo]
```
Generates `unique-angle/`, `dup-angle.txt`, and `DUPXYZ-angle.txt`.

**Available comparison modes:**
- `angle` — default 3P-KHG criterion (`RMSD < 0.5 Å`, `MaxDist < 1.5 Å`, `MaxTors < 30°`)
- `rmsd-graph` — RMSD with graph-based atom matching (`RMSD < 0.125 Å`)
- `rmsd-hungarian` — RMSD with Hungarian atom matching (`RMSD < 0.125 Å`)
- `rmsd-kabsch` — direct Kabsch RMSD without atom reordering (`RMSD < 0.125 Å`)

---

## Citation & Contact

When using this code, please cite:

> *"Molcomp: a software tool for the structural characterization and deduplication of conformers"*, A.M. Genaev, A.D. Moshchenkov, A.S. Ryzhako, A.A. Otlyotov, and Y. Minenkov, *Submitted*, 2026.

For inquiries, questions, or suggestions: `Yury.Minenkov@gmail.com`

---

## Development Note

The C11/OpenMP refactoring, memory-safety fixes, and precomputation architecture in `Molcomp-v2` were developed using AI-assisted coding tools and validated against the reference v1 implementation.
