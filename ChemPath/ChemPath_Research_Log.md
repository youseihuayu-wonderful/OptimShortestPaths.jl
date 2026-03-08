# ChemPath 项目研究日志

> 记录项目构思、调研、决策的全过程

---

## 1. 项目背景：OptimShortestPaths.jl 仓库概览

**问题：** 请全面总结这个repo做了什么？

**回答摘要：**

OptimShortestPaths.jl 是一个 Julia 包，实现了 **DMY（Duan-Mao-Yin）算法**（STOC 2025），用于有向图单源最短路径，时间复杂度 O(m log^(2/3) n)，优于 Dijkstra 的 O(m log n)。

核心组件：
- **DMY 算法**：三阶段——FindPivots（前沿稀疏化）、BMSSP（有界多源最短路径）、递归分解
- **多目标扩展**：Pareto前沿计算（加权和、ε约束、字典序）
- **问题转换范式**：实体→顶点、关系→边、目标→权重、解→路径
- **领域应用**：药物-靶点网络、代谢通路、治疗方案、供应链优化

---

## 2. 项目方向选择

**问题：** 融合图优化与前沿AI技术，构想项目方向

**三个方向（从易到难）：**

| 方向 | 内容 | 难度 |
|------|------|------|
| 方向一 | LLM + 图优化的药物筛选助手（chatbot） | 最容易，2-3周 |
| 方向二 | 多智能体 + 图优化的研发决策系统 | 中等，方向一的扩展 |
| 方向三 | RAG + 知识图谱 + 路径优化（药物重定位） | 最有深度，时间最长 |

**决策：** 先做方向一，再扩展到方向二。

**技术栈：** Python + LangChain/CrewAI + NetworkX + Streamlit + ChEMBL API + Chemprop

---

## 3. 数据库选择与项目价值

**问题：** 我用的 database 是什么？做完后可以解决哪些问题？

### 数据库体系

**核心数据源（必选）：**
- **ChEMBL 36**（2025秋发布）：200万+化合物，2000万+活性数据，药物-靶点实验结合力（IC50, Ki）
- **DrugBank**：17,000+已批准药物的详细信息

**知识图谱（构建更大的图）：**
- **DRKG**（Amazon/AWS）：97,238节点，587万+边，107种关系
- **PrimeKG**（哈佛）：17,080种疾病，405万+关系，含适应症/禁忌/超说明书用药边
- **Hetionet**：47,031节点，225万+关系（经典但2017年后未更新）

**AI预测层：**
- **Chemprop v2.2.2**（MIT，2026年1月发布）：D-MPNN预测未知化合物活性，速度提升2x，内存降低3x

### 能解决的问题

1. **药物虚拟筛选**：从百万化合物中几秒筛出Top候选，缩小实验范围
2. **药物重定位（Drug Repurposing）**：发现已有药物治疗新疾病的隐藏路径
3. **多目标优化选药**：Pareto前沿展示效果/毒性/成本的完整trade-off
4. **降低使用门槛**：LLM自然语言接口，非编程人员也能用
5. **预测未知活性**：Chemprop扩大图覆盖范围

### ChEMBL vs Chemprop 区别
- **ChEMBL = 数据库（仓库）**：存真实实验数据，只能查已有的
- **Chemprop = 预测模型（工具）**：从数据学习，预测未做过实验的化合物性质
- 互补关系：ChEMBL提供训练数据 → Chemprop学习 → 预测新组合 → 扩大图

---

## 4. 竞品分析与差异化定位

**问题：** 网上现有的软件和我正要做的有什么不同？如何做一个不一样的？

### 现有竞品

| 工具 | 来源 | 做什么 | 局限 |
|------|------|--------|------|
| **Llamole** | MIT, ICLR 2025 | 自然语言→生成新分子+合成路径 | 只做分子生成，不做筛选和优化排序 |
| **DrugGPT** | Oxford, Nature BME 2025 | 临床用药推荐chatbot | 只针对已批准药物的临床决策 |
| **DrugAssist** | Briefings in Bioinfo 2024 | 对话式分子优化 | 不涉及靶点网络和图优化 |
| **DrugChat** | 2023 | 输入分子图→问答化合物性质 | 单分子分析，不做网络级优化 |
| **DRKG** | Amazon/AWS 2020 | 知识图谱+图嵌入drug repurposing | 没有LLM，没有自然语言接口 |
| **K-Paths** | KDD 2025 | 知识图谱上K条最短路径 | 纯算法论文，没有LLM和用户界面 |
| **ChatInvent** | AstraZeneca 2025 | 多智能体分子设计+合成 | 闭源，不公开 |
| **MADD** | 2026 | 多智能体de novo化合物生成 | 侧重生成新分子，不是已有化合物筛选 |

### 市场空白

```
                    有LLM自然语言    无LLM
                    ─────────────   ──────────
有图优化/最短路径  │  空白！        │ K-Paths, DRKG
                    │               │
无图优化            │ Llamole,       │ 传统筛选工具
(嵌入/GNN/生成)     │ DrugGPT,      │
                    │ DrugAssist    │
```

**结论：没有任何现有工具同时具备 LLM自然语言接口 + 图最短路径优化 + 多目标Pareto + 已有化合物筛选。**

### 五个差异化点

1. **最短路径可解释性**（vs 图嵌入的黑箱分数）：输出完整路径 Drug A → Target X → Pathway Y → Disease Z
2. **多目标Pareto优化**（独有）：展示所有非支配方案的trade-off
3. **自然语言→图查询自动翻译**（vs K-Paths的手动查询）
4. **Chemprop预测扩展覆盖**（vs 只用已有数据）
5. **聚焦"筛选"而非"生成"**（vs Llamole/MADD）：已有化合物安全性数据更完善

### 项目命名与定位

> **ChemPath**: An LLM-powered drug screening assistant that uses shortest-path graph optimization and multi-objective Pareto analysis on ChEMBL bioactivity data, with Chemprop-predicted activity expansion, to deliver explainable, natural-language-driven compound recommendations.

---

## 5. 项目架构（已搭建框架）

```
ChEMBL 36 (真实数据) → 训练集
                          ↓
                     Chemprop (训练模型)
                          ↓
                     预测未知组合的IC50
                          ↓
              真实数据 + 预测数据 → 构建更大的图
                          ↓
                     最短路径优化 → 推荐结果
                          ↓
                     LLM 自然语言界面
                          ↓
                     Streamlit 前端
```

已编写的模块（模拟数据版）：
- `data_loader.py` — 数据加载
- `graph_builder.py` — 数据→图结构（NetworkX）
- `optimizer.py` — 多目标优化、Pareto前沿
- `llm_interface.py` — 自然语言解析
- `chemprop_integration.py` — Chemprop集成
- `demo.py` — 完整pipeline演示

---

## 6. 时间规划

| 周次 | 任务 |
|------|------|
| 第1周 | ChEMBL数据获取+清洗 → 图构建 |
| 第2周 | 最短路径算法 + LLM接入（Claude API） |
| 第3周 | Streamlit前端 + Chemprop集成 |

---

## 7. 三个优先修复项（已实现）

**问题：** 建议优先修复三个问题：SMILES验证、置信度标签、sensitivity analysis

### 修复1：SMILES验证 → `src/smiles_validator.py`

**问题：** 无效SMILES会导致pipeline崩溃
**方案：** 三层验证
- Basic：正则匹配、括号配对、原子检查（零依赖）
- RDKit：完整化学验证 + canonical SMILES + 分子量（如果安装了rdkit）
- Batch：批量验证 + 通过率报告

**效果：** mock数据10个化合物中，8个通过，2个被拒（BadCompound: 非法字符, EmptySmiles: 空字符串），pipeline正常运行不崩溃。

### 修复2：置信度标签 → `src/optimizer.py`

**问题：** 面试官和用户会问"你多有信心？"
**方案：** 4维评分体系
- 数据来源（实验 vs AI预测）：权重最高
- 临床阶段（Phase 0-4）
- IC50强度（<10nM高效 / <100nM中效 / >100nM低效）
- 毒性水平（<0.2低 / <0.35中 / >0.35高）

**标签：**
- `HIGH [+++]`：总分≥8（实验数据 + 上市药物 + 高效 + 低毒）
- `MEDIUM [++ ]`：总分5-7
- `LOW [+  ]`：总分<5（AI预测 或 高毒性 或 无临床数据）

**效果：** Osimertinib→EGFR 获得 HIGH（实验数据+Phase 4+IC50=0.5nM+Tox=0.18）

### 修复3：Sensitivity Analysis → `src/sensitivity_analysis.py`

**问题：** 需要展示对模型局限性有清醒认识
**方案：** 三项敏感性分析

1. **毒性惩罚权重敏感性**：penalty从0.0到5.0
   - 结果：Osimertinib和Erlotinib在所有权重下都排前列（ROBUST）
   - Lapatinib和Gefitinib在高惩罚下排名变动（SENSITIVE）

2. **IC50测量噪声敏感性**：±0%到±100%扰动
   - 结果：Top 3排名在所有噪声水平下稳定（实验CV通常20-50%）

3. **数据缺失敏感性**：随机移除0-50%的边
   - 结果：移除30%+数据后排名开始变动——数据覆盖率很重要

**局限性声明：**
- 排名依赖毒性权重选择
- IC50有固有实验噪声（20-50% CV）
- 缺失数据可以改变推荐
- AI预测值增加额外不确定性
- 本系统是筛选工具，不能替代实验验证

### 当前代码文件（Phase 1重构后）

**新包结构** (`chempath/`)：

| 文件 | 功能 | 状态 |
|------|------|------|
| `chempath/__init__.py` | 包根，版本号 | 完成 |
| `chempath/data/mock_data.py` | Mock ChEMBL 36数据（8真实药物+2无效测试化合物） | 完成 |
| `chempath/data/loader.py` | 数据访问函数 | 完成 |
| `chempath/chemistry/smiles.py` | SMILES三层验证（basic→RDKit→batch） | 完成 |
| `chempath/graph/builder.py` | 数据→图（NetworkX），含SMILES验证+毒性惩罚 | 完成 |
| `chempath/graph/optimizer.py` | 排名+Pareto+置信度标签+knee point | 完成 |
| `chempath/graph/analysis.py` | 三项敏感性分析 | 完成 |
| `chempath/ml/__init__.py` | Chemprop集成占位（Phase 4） | 占位 |
| `chempath/agent/__init__.py` | Claude Tool Use/MCP占位（Phase 3） | 占位 |
| `chempath/ui/__init__.py` | Chainlit前端占位（Phase 3） | 占位 |
| `scripts/demo.py` | 完整pipeline演示（绝对import） | 完成，已验证 |
| `pyproject.toml` | uv项目配置（hatchling构建） | 完成 |

**旧文件** (`src/` — 已废弃，待清理)：

| 文件 | 状态 |
|------|------|
| `src/data_loader.py` | 已迁移到 `chempath/data/` |
| `src/smiles_validator.py` | 已迁移到 `chempath/chemistry/smiles.py` |
| `src/graph_builder.py` | 已迁移到 `chempath/graph/builder.py` |
| `src/optimizer.py` | 已迁移到 `chempath/graph/optimizer.py` |
| `src/sensitivity_analysis.py` | 已迁移到 `chempath/graph/analysis.py` |
| `src/demo.py` | 已迁移到 `scripts/demo.py` |

---

## 8. 深度技术审计与修订计划（2026-03-06）

**问题：** 深度复查所有技术选型，确保使用cutting-edge、stable、future-proof、modular的技术

### 审计发现：6个问题 + 修订方案

---

### 问题1：Python 3.8 太老了

**现状：** 系统Python是3.8.9（2021年），不支持`str | None`语法、match语句等。
**风险：** Python 3.8已于2024年10月EOL，不再接收安全更新。Chemprop v2要求Python ≥3.11。

**修订：** 使用 **`uv`**（Rust写的Python包管理器）管理项目。
- uv是2026年Python社区共识的最佳工具（10-100x快于pip）
- 自动管理Python版本（无需手动装3.11+）
- lockfile确保可复现性
- 替代pip + venv + pyenv的全部功能

```
# 一条命令搞定
uv init chempath --python 3.12
uv add networkx rdkit-pypi chemprop anthropic chainlit
```

---

### 问题2：项目结构不够模块化

**现状：** 所有.py文件平铺在`src/`里，互相用相对import。
**风险：** 无法`pip install`，无法被其他项目引用，测试不方便。

**修订：** 采用标准Python包结构 + `pyproject.toml`

```
ChemPath/
├── pyproject.toml              # uv/pip 项目配置（替代 requirements.txt）
├── README.md
├── chempath/                   # 正式Python包（可pip install -e .）
│   ├── __init__.py
│   ├── data/                   # 数据层
│   │   ├── __init__.py
│   │   ├── chembl_client.py    # ChEMBL 36 API
│   │   ├── loader.py           # 统一数据加载接口
│   │   └── mock_data.py        # 测试用mock数据
│   ├── chemistry/              # 化学计算层
│   │   ├── __init__.py
│   │   ├── smiles.py           # SMILES验证（RDKit 2025.09）
│   │   └── descriptors.py      # 分子描述符计算
│   ├── graph/                  # 图构建与优化层
│   │   ├── __init__.py
│   │   ├── builder.py          # 数据→图
│   │   ├── optimizer.py        # 最短路径+Pareto
│   │   └── analysis.py         # 敏感性分析
│   ├── ml/                     # 机器学习预测层
│   │   ├── __init__.py
│   │   └── predictor.py        # Chemprop v2集成
│   ├── agent/                  # LLM/Agent层
│   │   ├── __init__.py
│   │   ├── query_parser.py     # 自然语言→结构化查询
│   │   ├── tools.py            # Claude tool definitions (MCP-ready)
│   │   └── engine.py           # 查询执行引擎
│   └── ui/                     # 前端层
│       ├── __init__.py
│       └── app.py              # Chainlit聊天界面
├── tests/                      # 测试
│   ├── test_smiles.py
│   ├── test_graph.py
│   ├── test_optimizer.py
│   └── test_integration.py
├── notebooks/                  # 探索性分析
│   └── data_exploration.ipynb
└── scripts/                    # 独立脚本
    └── demo.py
```

---

### 问题3：前端选型应改为 Chainlit

**原计划：** Streamlit
**问题：** Streamlit是dashboard框架，每次交互重跑整个脚本。对chatbot场景不理想。

**修订：** 改用 **Chainlit** — 专为LLM chatbot设计的框架。
- 原生支持对话历史、流式输出、多轮交互
- 内置LangChain/LlamaIndex集成
- 支持文件上传（用户可上传自己的化合物列表）
- Streamlit保留用于数据可视化dashboard（Pareto图等）

---

### 问题4：LLM接入应使用 Claude Tool Use（而非关键词匹配）

**原计划：** 关键词匹配解析自然语言
**问题：** 脆弱，不能处理复杂查询。

**修订：** 使用 **Claude API Tool Use**（function calling）。
- Anthropic已将ChEMBL集成到Claude的healthcare能力中
- 定义structured tools让Claude自动调用：

```python
tools = [
    {
        "name": "screen_compounds",
        "description": "Screen compounds for a target protein",
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {"type": "string", "description": "Target protein name (e.g., EGFR)"},
                "max_ic50": {"type": "number", "description": "Maximum IC50 in nM"},
                "max_toxicity": {"type": "number", "description": "Maximum toxicity score (0-1)"},
                "strategy": {"enum": ["efficacy", "safety", "balanced"]},
            },
            "required": ["target"]
        }
    },
    {
        "name": "compute_pareto_front",
        "description": "Find Pareto-optimal compounds balancing efficacy and toxicity",
        ...
    },
    {
        "name": "run_sensitivity_analysis",
        "description": "Test robustness of recommendations",
        ...
    }
]
```

- **MCP-ready**：Tool定义符合Model Context Protocol标准，未来可被任何LLM调用

---

### 问题5：图库选型确认

**原计划：** NetworkX
**审计结论：** **保持NetworkX**，但做好迁移准备。

理由：
- NetworkX是Python图分析下载量最大的库（50M+下载/月），生态最好
- 我们的数据规模（ChEMBL几万化合物、几十万边）在NetworkX能力范围内
- igraph（C底层，快10x）和graph-tool（C++底层）更快，但API不如NetworkX友好

**迁移策略：** 将图操作封装在`graph/builder.py`中，不在其他模块直接调用nx。
如果未来数据量增长需要性能，只改一个文件。

**知识图谱扩展路径：**
- 当前：NetworkX（内存图，<100K节点）
- 未来：Neo4j（持久化图数据库，百万+节点，支持Cypher查询）
- DRKG/PrimeKG都有Neo4j部署方案

---

### 问题6：Chemprop vs 替代方案

**审计结论：** **保持Chemprop v2**，它是最佳选择。

对比：
| 模型 | 优势 | 劣势 |
|------|------|------|
| **Chemprop v2** | D-MPNN，端到端，多GPU扩展，2026.1最新发布 | 需要GPU训练 |
| MolPROP (ChemBERTa+GNN) | 临床毒性预测更好(95% vs 90% AUC) | 更复杂，非标准 |
| Uni-Mol+ | 3D构象精度高 | 专注量子化学，不适合IC50预测 |

Chemprop v2是MIT开发的标准工具，peer-reviewed (J. Chem. Inf. Model, 2026)，社区最大。

---

### 修订后的技术栈

| 层级 | 原计划 | 修订 | 理由 |
|------|--------|------|------|
| **Python版本** | 3.8 (系统) | 3.12 (via uv) | 3.8 EOL，Chemprop需3.11+ |
| **包管理** | pip + requirements.txt | **uv + pyproject.toml** | 10-100x快，lockfile可复现 |
| **项目结构** | 平铺src/ | **标准Python包** | 可安装、可测试、模块化 |
| **化学验证** | 自写basic validator | **RDKit 2025.09** + basic fallback | 行业标准，canonical SMILES |
| **图库** | NetworkX | **NetworkX**（封装隔离） | 生态最好，数据规模够用 |
| **LLM接口** | 关键词匹配 | **Claude Tool Use** (MCP-ready) | 结构化、可靠、可扩展 |
| **前端** | Streamlit | **Chainlit**（chatbot）+ Streamlit（dashboard） | 各用其长 |
| **Agent框架** | LangChain/CrewAI | **Claude API直接调用**（不过度依赖框架） | 减少依赖，MCP兼容 |
| **预测模型** | Chemprop v2 | **Chemprop v2.2.2** | peer-reviewed，最快最省内存 |
| **数据API** | chembl_webresource_client | **直接REST API** + 本地缓存 | 官方client维护停滞 |

---

### 修订后的执行计划

#### Phase 1：基础重构（当前 → 1周）
1. 安装uv，初始化项目为标准Python包
2. 重构现有6个文件到新目录结构
3. 添加pyproject.toml + uv.lock
4. 安装RDKit，升级SMILES验证器
5. 写基础测试（pytest）

#### Phase 2：真实数据（1周）
1. 写ChEMBL REST API client（不依赖官方Python client）
2. 拉取10个抗癌靶点的IC50数据
3. 数据清洗pipeline（单位标准化、去重、SMILES规范化）
4. 替换mock data，验证pipeline在真实数据上工作

#### Phase 3：LLM + Tool Use（1周）
1. 定义Claude Tool Use schema（screen, pareto, sensitivity）
2. 实现query_parser：自然语言 → tool call → 结果 → 自然语言回复
3. 多轮对话支持（追问、修改条件）
4. Chainlit前端集成

#### Phase 4：预测扩展 + 知识图谱（1-2周）
1. Chemprop v2训练（在ChEMBL数据上）
2. 预测未知compound-target组合的IC50
3. 预测边加入图（带uncertainty penalty和LOW confidence标签）
4. 集成DRKG/PrimeKG做multi-hop drug repurposing路径

#### Phase 5：打磨 + 展示（1周）
1. Streamlit dashboard（Pareto图、网络可视化、敏感性热力图）
2. README + 截图 + GIF demo
3. 写一个reproducible案例研究（如EGFR inhibitor selectivity分析）

---

### 关键参考文献（技术选型依据）

| 技术 | 来源 |
|------|------|
| uv包管理 | 2026年Python社区共识最佳工具 |
| Claude Tool Use | Anthropic官方文档 + Healthcare集成 |
| Chainlit | 专为LLM chatbot设计，原生支持Tool Use |
| MCP (Model Context Protocol) | Anthropic 2024年发布，Linux Foundation管理，OpenAI/Google已采用 |
| RDKit 2025.09 | C++20现代化，canonical SMILES标准化 |
| Chemprop v2 | J. Chem. Inf. Model, 2026年1月，peer-reviewed |
| K-Paths | KDD 2025，Yen's K最短路径做drug repurposing |
| NetworkX | 50M+月下载，Python图分析事实标准 |

---

## 9. Phase 1 重构完成（2026-03-06）

### 完成内容

Phase 1基础重构已全部完成，项目从平铺`src/`目录迁移到标准Python包结构。

**关键成果：**
- **uv v0.9.28** 初始化项目，自动下载Python 3.12.12
- **pyproject.toml** 配置：hatchling构建后端，NetworkX 3.6.1为核心依赖
- **6层模块结构**：`chempath/{data,chemistry,graph,ml,agent,ui}`
- **全部使用绝对import**：`from chempath.graph.builder import ...`
- **Python 3.12语法**：`str | None`、`list[dict]` 等现代类型注解
- **Demo验证通过**：`uv run python scripts/demo.py` 输出与旧pipeline完全一致

### 解决的问题

1. **rdkit-pypi不兼容Python 3.12**：移除pip依赖，改为conda安装（SMILES验证器有basic fallback）
2. **Python版本解析错误**：固定`requires-python = ">=3.12,<3.13"`
3. **模块找不到**：`uv sync --reinstall-package chempath`重建wheel
4. **敏感性分析重复警告**：`build_drug_target_graph()`增加`verbose`参数

### 待完成

- [x] 清理旧`src/`目录 ✓
- [x] 编写pytest测试（60个测试全部通过） ✓
- [x] Phase 2：ChEMBL REST API client + 真实数据 ✓

---

## 10. Phase 2 真实数据完成（2026-03-06）

### 完成内容

Phase 2已完成，项目从mock数据升级到真实ChEMBL数据。

**新文件：**
- `chempath/data/chembl_client.py` — ChEMBL REST API client
- `scripts/fetch_chembl_data.py` — 数据拉取脚本
- `scripts/demo_real.py` — 真实数据pipeline演示
- `data/chembl_real.json` — 缓存的真实数据（1.3MB）
- `tests/test_chembl_client.py` — 13个单元测试

**ChEMBL Client特性：**
- 直接HTTP请求（urllib，无外部依赖）
- 分页获取（自动翻页，可配置最大记录数）
- 本地JSON文件缓存（避免重复API调用）
- 速率限制（0.5秒/请求）
- 数据清洗：单位标准化（nM/uM/pM/M→nM）、去重、无效数据过滤

**数据规模：**
- 10个抗癌靶点：EGFR, HER2, ALK, BRAF, ABL1, VEGFR2, MET, FGFR1, PI3Kα, mTOR
- 3,595个unique化合物（100% SMILES验证通过）
- 4,032条bioactivity记录
- 56个化合物有clinical phase数据（包括已批准药物如Sirolimus）

**IC50→Weight公式改进：**
- 旧公式：`max(0, 10 - pIC50)` — 所有IC50<0.1nM的化合物weight=0
- 新公式：`clamp(14 - pIC50, 0.01, 10)` — 亚皮摩尔级化合物仍有distinct权重

**关键发现：**
- Sirolimus (Rapamycin) 在VEGFR2排名中获得 **HIGH** 置信度（Phase 4 + experimental + high potency）
- 无毒性数据时Pareto前沿退化为单点——需要Phase 4的毒性预测模块或外部毒性数据库

### 待完成

- [x] Phase 3：Claude Tool Use + Chainlit前端 ✓
- [ ] Phase 4：Chemprop预测 + 知识图谱
- [ ] 集成外部毒性数据（如ToxCast/Tox21）使Pareto优化有意义

---

## 11. Phase 3 LLM Agent + Chat UI 完成（2026-03-06）

### 完成内容

Phase 3实现了Claude Tool Use集成和两种前端界面（CLI + Chainlit）。

**新文件：**
- `chempath/agent/tools.py` — 6个Claude Tool Use定义 + `ChemPathToolExecutor`执行器
- `chempath/agent/engine.py` — `ChemPathAgent` 多轮对话引擎（同步 + 流式）
- `chempath/ui/app.py` — Chainlit聊天界面（异步Tool Use + Step可视化）
- `scripts/chat.py` — CLI终端聊天界面
- `tests/test_agent_tools.py` — 17个工具执行器测试

**6个Tool定义（Anthropic format，MCP兼容）：**

| Tool | 功能 | 参数 |
|------|------|------|
| `list_targets` | 列出所有可用靶点 | 无 |
| `screen_compounds` | 筛选排名化合物 | target, strategy, top_n, max_ic50 |
| `compute_pareto` | Pareto前沿 + knee point | target |
| `run_sensitivity` | 敏感性分析 | target, analysis_type |
| `get_compound_info` | 化合物详情 | compound |
| `graph_summary` | 图结构摘要 | 无 |

**Agent Engine特性：**
- 多轮对话：自动管理消息历史
- Tool Use循环：Claude调用工具 → 执行 → 返回结果 → Claude总结
- 流式输出：`chat_stream()` 方法支持逐字输出
- 错误处理：未知工具、未找到靶点等情况有友好提示

**启动方式：**
```bash
# CLI聊天
export ANTHROPIC_API_KEY=sk-ant-...
uv run python scripts/chat.py

# Chainlit Web UI
uv run chainlit run chempath/ui/app.py
```

**测试：77个测试全部通过**

---

## 12. Phase 4 完成：高级分析功能

**实现的新功能：**

### 4a. 选择性分析（Selectivity Analysis）
- `compute_selectivity(G, compound_id, primary_target)` — 计算化合物对主靶点vs脱靶的选择性比
- 选择性比 = IC50(脱靶) / IC50(主靶点)，比值越高越好
- Agent工具 `compare_selectivity` 输出带风险标签：HIGH RISK (<10x), MODERATE (10-100x), LOW RISK (>100x)

### 4b. 分子性质估算 & 毒性代理
- `chempath/chemistry/properties.py` — 纯SMILES字符串分析，无需RDKit
- 估算属性：分子量、HBD、HBA、可旋转键、LogP
- Lipinski五规则违规检测 → 风险评分 (0.0-1.0)
- 当ChEMBL缺少毒性数据时，风险评分作为毒性代理注入图权重
- 3,595个化合物中37%（1,351个）获得非零风险分

### 4c. 多靶点对比（Head-to-Head）
- `compare_compounds_across_targets()` — IC50矩阵构建
- Agent工具 `head_to_head` 展示top N化合物在所有靶点的活性谱
- 显示亚纳摩尔标记、分子风险分、置信度

**测试：89个测试全部通过**

---

## 13. Phase 5 完成：Streamlit可视化仪表板

### 5a. 交互式Dashboard
- `chempath/ui/dashboard.py` — Streamlit + Plotly
- **5个标签页**：
  1. **Compound Rankings** — 排名表格 + IC50分布直方图 + 置信度饼图
  2. **Pareto Front** — 散点图（所有化合物灰色背景 + Pareto前沿蓝色菱形 + 膝点红色星标）
  3. **Network View** — 化合物-靶点网络图（节点颜色=置信度，节点大小=排名，距离=权重）
  4. **Selectivity** — 选择性条形图（10x阈值线）+ 脱靶表格
  5. **Molecular Properties** — Lipinski雷达图 + 风险标签

- **侧边栏控件**：靶点选择、策略选择、Top N滑块、毒性惩罚权重
- **启动命令**：`uv run streamlit run chempath/ui/dashboard.py`

### 5b. Agent格式优化
- `graph_summary` 工具增加：每靶点化合物数、IC50范围、临床阶段统计
- 数据统计：3,595化合物、10靶点、4,032相互作用、IC50范围0.02-100,000 nM

### 项目文件结构（最终）

```
ChemPath/
├── chempath/
│   ├── __init__.py
│   ├── agent/
│   │   ├── __init__.py
│   │   ├── engine.py          # Claude Agent多轮对话引擎
│   │   └── tools.py           # 8个Claude Tool Use定义 + 执行器
│   ├── chemistry/
│   │   ├── __init__.py
│   │   ├── properties.py      # SMILES分子性质估算
│   │   └── smiles.py          # SMILES验证
│   ├── data/
│   │   ├── __init__.py
│   │   ├── chembl_client.py   # ChEMBL REST API客户端
│   │   ├── loader.py          # 数据查询工具
│   │   └── mock_data.py       # 测试用模拟数据
│   ├── graph/
│   │   ├── __init__.py
│   │   ├── analysis.py        # 灵敏度分析
│   │   ├── builder.py         # 图构建（IC50→权重 + 毒性惩罚）
│   │   └── optimizer.py       # 排名、Pareto前沿、选择性分析
│   └── ui/
│       ├── __init__.py
│       ├── app.py             # Chainlit聊天界面
│       └── dashboard.py       # Streamlit可视化仪表板
├── data/
│   └── chembl_real.json       # 缓存数据（3,595化合物）
├── scripts/
│   ├── chat.py                # CLI聊天脚本
│   ├── demo_real.py           # 真实数据演示
│   └── fetch_chembl_data.py   # 数据获取脚本
├── tests/                     # 89个测试
└── pyproject.toml
```

### 待完成

- [ ] Phase 4（原计划）：Chemprop v2 IC50预测（需GPU）
- [ ] Phase 4（原计划）：DRKG/PrimeKG知识图谱集成
- [ ] 案例研究文档 + 使用教程

---

*最后更新：2026-03-06*
