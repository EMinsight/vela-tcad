# BVmethods NMOS 接触平衡态与 TDR 节点映射审计（2026-08-02）

## 结论

固定 `edge 3457 (914 -> 915)` 后，净掺杂、OldSlotboom、`ni_eff`、ohmic
平衡载流子和接触内建势已逐项复算。最早的不一致不是 ohmic 公式，而是
Sentaurus TDR 区域节点数据的导入顺序：TDR 的区域 nodal values 按全局节点
ID 升序存储，旧导入器却按节点在三角形流中的首次出现顺序绑定。

该错误同时打乱了状态场和 Vela 使用的逐节点掺杂。修正导入器、重新导入
网格掺杂并重算后：

- 1、2、5、10 mV 的 Vela/Sentaurus 漏极电流比稳定在 `0.842--0.848`；
- 绝对误差为 `0.071--0.075 dex`，原先的 `4.58--4.66 dex` 差异消失；
- `edge 3457` 上原先约 `0.65 V/0.19 nm` 的假跳变消失；
- 剩余约 `35.6 mV` 接触势差来自 Sentaurus `Fermi` 统计与 Vela Boltzmann
  统计的明确模型差异。

## 1. 节点映射证据

以 10 mV Sentaurus 状态的硅区为例，分别用两种节点顺序重建三角形势梯度：

| 映射 | 势梯度中位数 (V/m) | P99 (V/m) | 最大值 (V/m) |
|---|---:|---:|---:|
| 旧：三角形首次出现 | 2.2298e7 | 1.4300e9 | 1.9941e10 |
| 新：全局节点 ID 升序 | 2.6996e6 | 4.1744e7 | 4.4574e7 |
| Sentaurus ElectricField | 2.6986e6 | -- | 4.4567e7 |

把势和电场都按升序映射后，`log10(|grad psi|)` 与 Sentaurus 电场的相关系数
为 `0.99709`，幅值比中位数为 `0.98824`。这给出了独立于 Vela 求解器的
映射判据。

## 2. edge 3457 的掺杂更正

单位为 `cm^-3`：

| node | 旧净掺杂 | 修正净掺杂 | 影响 |
|---:|---:|---:|---|
| 914 | +2.60298e18 | +3.25692e20 | 漏极接触掺杂低估约 125 倍 |
| 915 | -7.29443e17 | +3.25684e20 | 极性被错误翻成 p 型 |
| 1224 | +5.81220e18 | +2.59000e20 | 漏极接触掺杂低估约 44.6 倍 |

旧映射把 node 915 错误标成 p 型，正是先前逐边审计中出现巨大势差和 12 dex
局部通量差的直接原因。

## 3. OldSlotboom、ni_eff 与 ohmic 接触复算

node 914，300 K：

| 量 | 数值 |
|---|---:|
| donor | 3.266916861e20 cm^-3 |
| acceptor | 1.000000000e18 cm^-3 |
| net doping | 3.256916861e20 cm^-3 |
| total impurity | 3.276916861e20 cm^-3 |
| OldSlotboom `Delta Eg` | 0.145981282 eV |
| Vela `ni_eff` | 1.683405723e11 cm^-3 |
| Boltzmann `n_eq` | 3.256916861e20 cm^-3 |
| Boltzmann `p_eq=ni_eff^2/n_eq` | 87.01035 cm^-3 |
| `Vt*ln(n_eq/ni_eff)` | 0.552799203 V |
| Vela 0 V 接触 `psi` | 0.552799203 V |

计算值与 Vela 边界状态逐位一致，说明修正掺杂后 OldSlotboom、`ni_eff`、
电中性载流子和 ohmic 内建势公式内部自洽。

同一节点的 Sentaurus 0 V 状态为：

- `n = 3.256916861e20 cm^-3`，与净掺杂和 Vela 主载流子一致；
- `p = 190.4397 cm^-3`；
- `psi = 0.588423157 V`。

Sentaurus 官方 deck 启用了 `Fermi` 统计。高达 `3.26e20 cm^-3` 的退化电子
气体不能满足 Boltzmann 的 `n=ni_eff*exp((psi-phin)/Vt)`，因此其接触势比
Vela Boltzmann 值高 `35.624 mV`。这属于已识别的物理模型差异，不是符号或
接触掺杂错误。

## 4. 严格容差低偏压电流

Vela 使用 `reltol=1e-12`、`abstol=1e-10`，避免 10 mV 点按相对残差过早
验收：

| Vd (V) | Vela (A/um) | Sentaurus (A/um) | Vela/S | 差值 (dex) |
|---:|---:|---:|---:|---:|
| 0.001 | 8.72846e-11 | 1.029e-10 | 0.84825 | -0.07148 |
| 0.002 | 1.70877e-10 | 2.015e-10 | 0.84803 | -0.07159 |
| 0.005 | 4.00654e-10 | 4.737e-10 | 0.84580 | -0.07273 |
| 0.010 | 7.20994e-10 | 8.564e-10 | 0.84189 | -0.07475 |

低场电导已经进入 `<0.1 dex` 的可比范围。0 V 两边电流均处于数值零附近，
不使用其比值评价物理误差。

## 5. 实施内容与产物

- `src/io/SentaurusTdrReader.cpp`：区域节点顺序改为全局节点 ID 升序；
- `tests/test_sentaurus_tdr_reader.cpp`：新增“部分区域、置乱单元顺序”回归；
- `scripts/run_bvmethods_nmos_equilibrium_stages.py`：支持指定重新导入的掺杂；
- `scripts/run_bvmethods_nmos_bias_validation.py`：支持严格 `reltol/abstol`；
- 修正导入：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/mesh_import_sorted_node_order`；
- 严格低偏压结果：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sorted_node_order_fix_20260802/low_bias_strict`；
- 修正后的逐边结果：同目录的 `low_bias_strict_edge_compare`。

## 6. 后续边界

当前不应再调迁移率来解释原来的五数量级差异，因为该差异已经由导入错误
消除。下一层模型差异是 Fermi/Boltzmann 统计；若目标是严格复现官方 deck，
需要为载流子状态、ohmic 接触和广义 SG 输运整体加入一致的 Fermi-Dirac
统计，不能只给接触势增加经验偏移。
