# DMY Correctness Bug — Report for Decision

**Date:** 2026-06-09
**Found while:** building the Section 7 (Experimental Evaluation) benchmark — adding a
Graphs.jl baseline and a correctness check before timing.
**Severity:** High — contradicts the paper's central correctness claim and blocks Section 7.

---

## TL;DR

`dmy_sssp!` returns **incorrect shortest-path distances on general random directed
graphs** (errors are always over-estimates; some reachable vertices are dropped to ∞).
The bug is **reproducible with the repository's own `compare_with_dijkstra`** and is
**two layered problems**:

1. The DMY recursive **core is itself wrong** — it does not independently compute correct SSSP.
2. `dmy_sssp!` masks this with a post-hoc **Bellman-Ford "safety net" capped at 3 passes**,
   which only rescues the core when error chains are ≤ 3 hops.

The test suite passes only because **every random test graph is a forward tree/DAG**, where
the core's residual is short enough for the 3-pass cap to repair.

---

## Reproduction (≈10 s, uses only package functions)

```bash
julia --project=. docs/dmy_bug_repro.jl
```

Output (n = 1000, m = 2000, source = 1):

| backbone | `compare_with_dijkstra` | shipped `dmy_sssp!` errors | pure-core errors |
|----------|-------------------------|----------------------------|------------------|
| tree (forward spanning tree) | match=**true**, disc=0 | 0 | 23 |
| path (random Hamiltonian path) | match=**false**, disc=10 | **10** | 296 |

Ground truth is the package's own `simple_dijkstra`; independently confirmed against
Graphs.jl `dijkstra_shortest_paths` (the two agree exactly). The only difference between
the two rows is the spanning structure of the graph.

---

## Root-cause analysis

### Layer 1 — the recursive core is incorrect
Bypassing the safety net and calling `recursive_layer!` directly, the **pure core output is
wrong on both graph types** (23 errors on tree, 296 on path at n=1000). The frontier
sparsification / block partitioning / recursion in `src/dmy_algorithm.jl::recursive_layer!`
(lines ~28–133) drops vertices and leaves over-estimated distances. The core never
independently produces correct shortest paths.

### Layer 2 — the safety net hides it, but only sometimes
`dmy_sssp!` runs a Bellman-Ford relaxation after the core, capped at
`max_passes = min(n, 3)` (`src/dmy_algorithm.jl:180`; same in `dmy_sssp_with_parents!:238`
and `dmy_sssp_bounded!:302`). Bellman-Ford repairs the core's residual only if the error
propagates in ≤ 3 passes:

- **tree**: core residual needs **2** passes → 3-cap is enough → shipped result correct.
- **path**: core residual needs **4** passes → 3-cap insufficient → 10 errors survive.

### Why the 2,044-assertion suite never caught it
`test/test_correctness.jl` only exercises:
- hand-built graphs of n = 3–6,
- random graphs of **n = rand(5:20)** generated with the forward-tree generator
  (`parent = rand(1:i-1)`), and
- one structured n = 50 chain.

All random cases are tree/DAG-shaped, where the 3-pass safety net always succeeds. No test
uses a general random directed graph at the scale where the core's error chain exceeds 3 hops.

---

## Implications for the paper

- The headline **"DMY produces identical results to Dijkstra (47,000/47,000)"** likely holds
  only because Hetionet's structure keeps residual chains short; it is **not** true for
  general directed graphs.
- **"Fast *and* correct" does not hold as implemented.** Current behavior is *fast but wrong*
  (core + 3 BF passes). Making it correct requires Bellman-Ford to convergence
  (~graph-diameter passes): cheap on random graphs (~4) but **O(√n)–O(n) on grid / road
  networks** — exactly the graph families Section 7 needs — which degenerates toward
  O(nm), **slower than Dijkstra**. Simply raising the cap removes the performance story.
- The current speedup figures (38×–322×) were measured on the *fast-but-wrong* path and
  against a non-tuned Dijkstra; they need re-measurement against a tuned binary-heap /
  Graphs.jl baseline regardless.

---

## Decisions needed (Daniel)

1. **Was this known?** How was the 47,000/47,000 match verified — which Dijkstra, which graph,
   exact equality or tolerance?
2. **Do you have a known-good DMY core** (so we replace, not re-debug)?
3. **Direction:** fix the core so it is correct without the Bellman-Ford crutch, *or* narrow
   the paper's scope/claims to the graph classes where it holds?
4. **Your ready-to-run harness:** does it show DMY as *correct*? If so, how does it build/invoke
   the graph differently from the repro above?

---

## Pointers

- Algorithm: `src/dmy_algorithm.jl` (`recursive_layer!`, `dmy_sssp!`)
- Subroutines: `src/bmssp.jl`, `src/pivot_selection.jl`
- Repo's own comparator: `src/utilities.jl::compare_with_dijkstra` (vs `simple_dijkstra`)
- Repro: `docs/dmy_bug_repro.jl`
- Optional next step (not yet done): pin which of FindPivots / BMSSP / block-recursion drops
  the vertices, via read-only instrumentation of the core.
