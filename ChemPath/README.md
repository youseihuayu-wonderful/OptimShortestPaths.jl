# ChemPath

Drug discovery demonstration built on [OptimShortestPaths.jl](../README.md). ChemPath casts drug repurposing as a shortest-path problem on biological knowledge graphs, using the DMY algorithm as the computational backend.

## How It Works

**Core idea**: Drugs, targets, pathways, and diseases are vertices in a directed graph. Biological relationships (binding affinity, pathway participation, disease association) become weighted edges. The shortest path from a drug to a disease represents the most likely mechanism of action.

```
Drug ──[binding]──> Target ──[pathway]──> Pathway ──[association]──> Disease
  w = -log(P_bind)    w = -log(P_path)      w = -log(P_assoc)
```

Weight transformation: `w = -log(P)` converts probabilities to additive costs. Shortest path = maximum-likelihood path.

## Architecture

```
chempath/
├── data/                  # Data layer
│   ├── curated_data.py    # 18 validated drug-target pairs
│   ├── chembl_client.py   # ChEMBL REST API client (no API key needed)
│   ├── biological_network.py  # Hetionet knowledge graph loader
│   └── data_validation.py # 8 automated quality checks
├── graph/                 # Graph layer
│   ├── builder.py         # Data → NetworkX graph conversion
│   ├── optimizer.py       # Ranking, Pareto analysis, selectivity
│   ├── julia_bridge.py    # Python → Julia IPC (JSON + subprocess)
│   └── dmy_bridge.jl      # Julia-side DMY execution script
├── agent/                 # LLM layer
│   ├── engine.py          # Claude multi-turn conversation
│   └── tools.py           # 6 tool definitions for Claude Tool Use
└── ui/                    # Frontend layer
    ├── app.py             # Chainlit chatbot interface
    └── dashboard.py       # Streamlit visualization dashboard
```

## Datasets

| Dataset | Source | Scale | Usage |
|---|---|---|---|
| Curated | Manual curation | 18 drugs, 51 nodes | Development & unit tests |
| ChEMBL | REST API (real-time) | 1,553 drugs | Scale-up validation |
| Hetionet | CC0 knowledge graph | 47K nodes, 2M edges | Large-scale benchmark |

## Results

**Prediction quality** (shortest-path distance predicts drug-disease association):

| Dataset | AUROC | MRR | vs Random |
|---|---|---|---|
| Curated (18 drugs) | 0.806 | 0.691 | p=0.0001 |
| ChEMBL (1,553 drugs) | 0.909 | 0.733 | -- |
| Hetionet (47K nodes) | 0.794 | 0.210 | 5.2x MRR |

**Julia bridge verification**: DMY and Dijkstra produce identical distances (max discrepancy = 0.00) across all tested graphs.

## Setup

```bash
cd ChemPath
uv sync                  # Install dependencies
uv run pytest tests/     # Run test suite (89 tests)
```

Requires Python 3.12 and Julia 1.9+ with OptimShortestPaths.jl installed.

## Scripts

| Script | Purpose |
|---|---|
| `scripts/demo.py` | Basic pipeline demo on curated data |
| `scripts/dmy_vs_dijkstra.py` | DMY vs Dijkstra wall-clock benchmark |
| `scripts/scaleup_test.py` | ChEMBL 1,553-drug validation |
| `scripts/hetionet_benchmark.py` | Hetionet 47K-node benchmark |
| `scripts/generate_final_figures.py` | Publication figures |

## Key Design Decisions

1. **No RDKit dependency**: SMILES validation uses a 3-tier fallback (regex, basic checks, optional RDKit) so the package installs without conda.
2. **Julia bridge via JSON IPC**: Python exports the graph as JSON, calls Julia subprocess, reads results. Avoids fragile PyJulia/pyjulia bindings.
3. **Probability-based weights**: All edge weights derived from `w = -log(P)` where P is a biological probability (binding, efficacy, evidence). This ensures non-negative weights (required by DMY) and makes path costs interpretable as negative log-likelihoods.
4. **Hold-out validation**: AUROC computed by holding out known drug-disease edges and checking if shortest-path distance ranks them above random pairs.
