# BVmethods 工作进展日报素材（2026-08-07 至 2026-08-08）

## 一句话进展

完成 BVmethods NMOS 外接电阻和 Voltage-to-Current 两种击穿边界方法的实现与
验收；两种方法相对 Sentaurus 的击穿电压误差均低于 `0.3%`。同时完成同电流
工作点的电势、电场、载流子、准费米势和 I-V 对比，并确认当前主要性能问题是
高场 Newton 更新过多、外层边界根重复调用完整 DD 求解，以及单次 Newton 更新
成本偏高。

## 2026-08-07（昨晚）

### 1. 修正 Sentaurus IIC 路径参考口径

- 确认此前用 Sentaurus Visual `ElectricField-V` 生成的流线不能等价替代
  `sdevice ComputeIonizationIntegrals` 的内部 `Eparallel` IIC 路径。
- 改用 WriteAll 最大场约束、`eCurrentDensity-V`/`hCurrentDensity-V` 方向和
  e/h/Mean 三项离化积分平台联合筛选候选路径。
- 在 `10.4482667308 V` 状态恢复出三条满足积分平台约束的候选：
  - rank-1：空穴电流方向，最大积分误差约 `0.053%`；
  - rank-2：电子电流方向，最大积分误差约 `0.0001%`；
  - rank-3：电子电流方向，最大积分误差约 `0.197%`。
- 完成 Sentaurus 与 Vela 路径几何比较：
  - rank-2 主体通道最接近，Vela 路径 `72.4%` 位于 Sentaurus 候选路径
    `10 nm` 范围内；
  - rank-3 种子仅相差 `1.25 nm`，但 Vela 低场尾部明显过长，问题集中在
    路径停止/通道保留语义；
  - rank-1 种子相差 `87.3 nm`，说明峰分组或物理通道排序尚未闭合。

### 2. 推进外接电阻与电流边界高场求解

- 固定已有 IIC 物理基线，没有调整迁移率或雪崩系数。
- 完成边界求解检查点、跨进程括区恢复、正负端点恢复和割线状态预测。
- 将局部 carrier-row 比率保留为诊断项，以全局电子/空穴连续性闭合作为硬验收；
  没有放松耦合 Newton 残差或自洽雪崩方程。
- 定位 `6.08709 -> 6.09959 V` 高场停滞：需要连续 Newton 轨迹和更完整的
  迭代预算，重启同一偏压会丢失单调收敛轨迹；该问题不是负载线方程符号错误。
- 完成 `0 -> 5.9 V` 的自洽雪崩预偏置链，为两种边界方法提供共同初始状态。

## 2026-08-08（今天）

### 1. 完成两种 BV 边界方法及最终验收

- 新增通用 `BoundaryControl` 括区标量闭合，DCSweep 支持互斥的：
  - `external_circuit.mode = series_resistor`；
  - `voltage_to_current`。
- 外接电阻负载线按 Sentaurus 二维单位闭合：

  `OuterVoltage = InnerVoltage + Iterminal[A/um] * R[ohm*um]`

- 增加 InnerVoltage、OuterVoltage、负载线残差、目标电流、边界残差、评估次数、
  检查点和恢复信息等输出。
- 最终结果：

| 方法 | Vela BV | Sentaurus BV | 相对误差 | 3% 验收 |
|---|---:|---:|---:|---|
| 外接电阻 `1e7 ohm*um` | `6.395887866 V` | `6.379791636 V` | `0.252300%` | PASS |
| Voltage-to-Current | `6.395904175 V` | `6.383184201 V` | `0.199273%` | PASS |

- 独立电流边界终点为 `6.395904200 V`、
  `1.0000000617e-4 A/um`，电流残差为 `6.17e-12 A/um`；电子和空穴全局
  连续性比率分别约 `1.22e-9` 和 `1.31e-9`。

### 2. 完成新 Sentaurus 同电流状态及关键物理量比较

- 在 Sentaurus O-2018.06-SP2 上重新执行 Voltage-to-Current，得到
  `1e-4 A/um` 时 `6.384111662 V`；Vela 相对该新结果误差为 `0.184717%`。
- 1909 个半导体节点与 Vela 网格按 node ID 精确对齐。
- 电势绝对误差中位数 `0.764 mV`、P95 `15.81 mV`，空间相关系数
  `0.9999957`。
- 载流子占据区比较：

| 物理量 | 中位误差 | P95 误差 | 相关系数 |
|---|---:|---:|---:|
| 电子浓度 | `0.003353 dex` | `0.078628 dex` | `0.999810` |
| 空穴浓度 | `0.023540 dex` | `0.325500 dex` | `0.996756` |
| 电子准费米势 | `0.002900 V` | `0.016477 V` | `0.999996` |
| 空穴准费米势 | `0.005291 V` | `0.026275 V` | `0.999990` |

- 高场空间形状保持一致。按匹配网格边投影，峰值场 Vela 比 Sentaurus 高
  `12.07%`；Sentaurus 峰值 10% 以上区域的中位相对误差 `2.76%`、P95
  `17.29%`，相关系数 `0.992835`。
- 生成电势、电场、I-V、运行时间和 Newton 累积图，以及可交互 HTML 报告。

### 3. 完成第一轮性能诊断

| 运行 | 墙钟时间 | Newton 更新 | 有效时间/更新 |
|---|---:|---:|---:|
| Sentaurus 完整预偏置和电流边界 | `103.24 s` | 370 | `0.279 s` |
| Vela Voltage-to-Current 恢复段 | `1032.37 s` | 279 | `3.700 s` |
| Vela 外接电阻恢复段 | `3885.10 s` | 1356 | `2.865 s` |

- 外接电阻最终括区需要 7 次完整 DD 评估；Voltage-to-Current 需要 2 次。
- 检查点恢复后约 `99.9%` 的时间位于新高场 Newton 求解中，I/O 不是主瓶颈。
- Sentaurus 已报告核心时间中，残差/RHS/装配占 `58.0%`，线性求解占
  `33.0%`，Jacobian 占 `7.2%`。
- Vela 当前尚无分阶段计时，无法可靠区分雪崩装配、Jacobian、SparseLU 数值
  分解和线搜索的内部占比；这是下一阶段首要任务。

### 4. 完成参考算例与模板归档

- 将 Sentaurus BVmethods NMOS 的 SDE、公共参数及 6 种 SDevice 方法脚本归档到
  `reference_tcad/bvmethods_sentaurus2018/`：ABA Poisson、ABA Coupled、外接
  电阻、Voltage-to-Current、Continuation 和 Transient。
- 在 `configs/templates/` 新增外接电阻与 Voltage-to-Current 两份 Vela 模板，
  并接入配置生成器。
- Continuation 和 Transient 仅归档 Sentaurus 参考脚本；Vela 尚未实现等价模式。

### 5. 测试和提交

- `dc_sweep` 测试通过。
- impact-ionization：52 个测试用例、734 条断言全部通过。
- 配置模板：14 个测试全部通过；两份新模板均可成功渲染。
- 当前本地 `main` 已包含：
  - `b23798c feat(simulation): add BV boundary control methods`；
  - `c2f85f6 feat(reference): add Sentaurus BVmethods decks and templates`。

## 当前风险与待办

1. Vela 与 Sentaurus 的运行范围、硬件和线程环境不同，现有速度倍数是系统级观测，
   不是严格同机基准。
2. 尚未插桩 Vela 内部阶段耗时，单次 Newton 慢的具体归因仍是待验证推断。
3. 外接电阻的嵌套标量根带来 7 次完整 DD 求解，是已确认的算法级开销。
4. IIC rank-1 的物理通道选择和 rank-3 的低场停止语义仍未完全闭合。
5. Sentaurus Continuation/Transient 已归档，但 Vela 对应能力尚未实现。

## 建议日报表述

今日完成 BVmethods NMOS 外接电阻和 Voltage-to-Current 两种击穿计算方法，
相对 Sentaurus 的击穿电压误差分别为 `0.252%` 和 `0.199%`，均通过 3% 验收；
同电流场量比较显示电势、载流子和准费米势高度一致，高场形状相关系数超过
`0.99`。性能初诊确认高场 Newton 更新数、边界外层重复 DD 求解和单次 Newton
成本是主要瓶颈。已完成 Sentaurus 六种 BVmethods 输入脚本及两种 Vela 可执行模板
归档，下一步将通过分阶段计时和受控基准定位残差装配、雪崩 Jacobian 与稀疏分解
的真实占比，再按测量结果优化。
