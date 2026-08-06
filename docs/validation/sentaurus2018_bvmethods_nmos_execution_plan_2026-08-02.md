# Sentaurus 2018 BVmethods NMOS 对比验证执行记录与计划

日期：2026-08-02

## 1. 当前结论

Sentaurus Applications Library 的 `GettingStarted/sdevice/BVmethods` 已在独立虚拟机目录完整运行。六个 Sentaurus Device 物理节点均成功，只有依赖图形环境的 Inspect 后处理节点失败；本地分析器已经替代 Inspect 并复现官方 BV 提取准则。

六种方法的结果分为两组：

| 方法 | BV 提取准则 | BV (V) |
|---|---|---:|
| ABA Poisson | `min(eIonIntegral,hIonIntegral)=1.05` | 5.305526 |
| ABA coupled / IIC | `q*Integral(AvalancheGeneration)=abs(Idrain)` | 6.377494 |
| 外接电阻 | `abs(Idrain)=1e-4 A/um` | 6.379792 |
| voltage-to-current | `abs(Idrain)=1e-4 A/um` | 6.383184 |
| continuation | `abs(Idrain)=1e-4 A/um` | 6.383727 |
| transient | `abs(Idrain)=1e-4 A/um` | 6.378835 |

IIC 与四种全耦合/延拓方法集中在 6.377--6.384 V，最大差约 6.23 mV；ABA 的 5.306 V 是预期的非自洽低估，不能直接作为全耦合 BV 金标准。

## 2. 数据与可复现入口

- 虚拟机原始独立运行目录：`/home/tcad/sentaurus_runs/vela_oracle/bvmethods_nmos_20260802_run01`
- 本地原始归档：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/bvmethods_run01_results.tgz`
- 解包结果：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw`
- 中性 TDR 导出：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/imported`
- 统一曲线与 BV 汇总：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/analysis`
- Vela 首版网格和 deck：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela`
- 统一分析命令：

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\analyze_sentaurus_bvmethods.py `
  --input-dir build-release\reference_tcad\bvmethods_sentaurus2018\run01\full_raw `
  --output-dir build-release\reference_tcad\bvmethods_sentaurus2018\run01\analysis
```

## 3. 官方脚本的物理模型基线

共同二维 NMOS 网格含 2719 个节点、5210 个三角形、4 个电极和 6 个材料区域：Si、SiO2 与 Nitride。Sentaurus 2018 Nitride 默认参数已从 `sdevice -P:Nitride` 核对：`epsilon=7.5`、`Eg0=5 eV`、`Chi0=1.9 eV`。

全耦合节点的主要模型为：

- Fermi 统计；
- Old Slotboom 有效本征浓度；
- DopingDep 迁移率；
- HighFieldSaturation，驱动力为 GradQuasiFermi；
- Enormal 表面迁移率；
- 掺杂相关 SRH；
- Band2Band(E2)；
- Avalanche(Eparallel)，并启用雪崩导数。

比较时必须逐层开启模型，不能把“几何/掺杂误差”和“输运/雪崩模型误差”混在同一条最终 BV 曲线中。

## 4. Vela 首轮启动结果

Sentaurus TDR 已成功转换成 Vela `mesh.json`、逐节点 `doping.csv` 和算例级 Nitride 材料文件。网格报告为：

- 5210 个三角形，退化单元为 0；
- 56 个负余切单元已采用非负重心盒回退；
- 最小角 0.1466 度，最大角 134.23 度；
- 最小边长为 `3.90625e-5` 个内部长度单位。

首轮求解尚未形成可比较的 Vela BV 曲线：

1. 常数迁移率、无复合、无雪崩的 Gummel 路径可以完成迭代，但 0 V 在 200 次内未收敛。
2. 加入 Masetti/SRH/OldSlotboom 后，0 V 连续性线性系统在高掺杂内部节点 1172、1176 出现非有限对角项。
3. 常数迁移率的 coupled-Newton + Poisson-block 路径完成网格和 Poisson 初值阶段，但第一次耦合步因 `carrier_invalid` 被线搜索拒绝；报告中的非正载流子节点主要来自纯绝缘区的 pinned carrier 状态。
4. Sentaurus 的 gate 使用 `Barrier=-0.55`。Vela deck 暂以 `metal_gate` 和 `flatband_voltage=-0.55 V` 表示；在完整 DD 运行前还需验证该边界在 Newton 路径中的载流子行屏蔽与 Sentaurus 势垒定义是否一致。

因此，当前状态是“参考金标准已建立、网格/字段已导入、Vela 首启动阻塞点已定位”，不是“Vela BV 已通过”。

## 5. 分阶段关键物理量验证计划

### P0：输入封印

目标：确保两个求解器看到相同器件，而不是相似器件。

- 节点坐标、三角形连接、区域和接触节点集合；
- 逐节点 Donor、Acceptor、净掺杂与总杂质；
- Si/SiO2/Nitride 的介电常数；
- gate 势垒/平带电压、source/drain/substrate 电压；
- 2-D 电流归一化宽度。

门槛：拓扑计数和接触节点集合完全一致；坐标与掺杂仅允许格式转换舍入误差。

### P1：0 V 平衡态与低场输运

先禁用 avalanche 和 Band2Band，逐项比较：

- electrostatic potential；
- electron/hole density；
- e/h quasi-Fermi potential；
- space charge；
- source/drain/substrate 的平衡接触状态；
- 载流子连续性残差与端电流守恒。

首轮建议门槛：半导体区势中位绝对误差小于 20 mV；主要载流子活跃区 log10 密度误差小于 0.3 decade；端电流不守恒小于 1%。

### P2：关断态无雪崩电场

在 drain = 1、2、4、5、6 V 保存共同偏压状态，比较：

- 电势与电场矢量；
- 最大电场及其坐标；
- drain 边缘、栅边缘和结区的场剖面；
- 电子/空穴电流密度矢量；
- SRH 与 B2B 生成/复合率。

首轮建议门槛：最大场相对差小于 10%，热点位置误差不超过一个局部网格单元；主要剖面中位相对差小于 10%。

### P3：ABA 与 IIC 算子级验证

在同一无雪崩状态上离线比较：

- `eAlphaAvalanche`、`hAlphaAvalanche`；
- `eIonIntegral`、`hIonIntegral`、`MeanIonIntegral`；
- `AvalancheGeneration`、电子/空穴雪崩生成率；
- `q * integral(AvalancheGeneration)`；
- `abs(Idrain)` 与积分雪崩电流之比；
- 雪崩热点支撑区域、结区占比和左右肩部占比。

ABA 只用于验证积分路径和趋势；IIC 的电流交点作为全耦合 BV 的主要低成本判据。

### P4：全耦合 BV 曲线

完成 0--7 V 自适应扫描后比较：

- InnerVoltage 与 OuterVoltage；
- drain 总电流、电子电流、空穴电流；
- 差分电导与电流跳变；
- 1e-4 A/um 阈值 BV；
- IIC、外接电阻、voltage-to-current、continuation、transient 五种结果的一致性。

验收建议：Vela BV 对 IIC/全耦合 Sentaurus 中位基准的误差先控制在 5%，再收紧到 2%；击穿前电流在有效动态范围内控制到 0.3 decade。

## 6. Vela 后续执行顺序

1. 修复/澄清绝缘节点 carrier validation：纯绝缘节点的 `n=p=0` 应视为 pinned 合法状态，不应触发 coupled-Newton 的全局 `carrier_invalid`。
2. 对节点 1172、1176 输出电子/空穴连续性矩阵构成，分离 Masetti、SRH、OldSlotboom 和 SG Bernoulli 项，找到非有限对角的首个产生项。
3. 用常数迁移率、无复合、无雪崩跑通 0 V，再按 `SRH -> OldSlotboom -> Masetti -> high-field -> Enormal prototype -> avalanche postprocess -> avalanche self-consistent` 顺序加模型。
4. 明确 `Barrier=-0.55` 到 Vela `flatband_voltage` 的符号和参考能级映射，并增加氧化层金属栅 Newton 边界测试。
5. 先生成 0--6 V 无雪崩状态做 P1--P3；算子级通过后才运行 0--7 V 自洽雪崩。

## 7. 本地开源算例候选及 Sentaurus 转换优先级

### 候选 A：Genius PN diode Avalanche（优先）

源文件：`D:/code-repo/Genius-TCAD-Open/examples/PN_Diode/Avalanche/avalanche.inp`

优点：纯二维 Si PN 二极管、脚本自建三角网格、解析掺杂、Ohmic 接触、明确的反向 DC 扫描（0 到 -48 V，再以 0.1 V 扫到 -60 V）和局部冲击电离。没有外部网格或预计算 restart 依赖，最适合作为第二个独立 BV 金标准。

转换计划：

- 用 Sentaurus SDE 重建 3 um x 6 um 的二维矩形和等价网格密度；
- donor `1e16 cm^-3`，顶部局部 acceptor `1e19 cm^-3`，按 Genius analytic profile 复现结深/横向扩散；
- 建立无雪崩正向 IV 和反向漏电基线；
- 建立 Van Overstraeten 与 Genius Local-II 两套 Sentaurus 变体，分别用于软件间模型差异和 Vela 已支持模型的对照；
- 提取 IV、BV、最大场、alpha、雪崩生成率和积分生成电流。

### 候选 B：Charon silicon diode breakdown（第二优先）

源文件：`D:/code-repo/tcad-charon/test/nightlyTests/impact_ionization/si_diode.brk.crowellsze.dd.inp`

优点：二维 Si 二极管、Crowell-Sze 雪崩、Fermi-Dirac、Arora 迁移率、SRH 和 -30 到 -34.6 V 延拓，物理量定义完整。

限制：依赖预计算 Exodus restart 和 optical-generation 文件；Sentaurus 转换前应先去除光生项并用解析平衡态替代 restart，随后再决定是否保留 Crowell-Sze 作为模型差异用例。

### 暂不优先

- DEVSIM 本地 examples 未检出带 avalanche/breakdown 的现成二维算例。
- PISCES 有多个 PN IV 教学输入和 `breakdown.f` 实现，但算例/模型组织较旧，作为第三方历史交叉检查优于首轮 Sentaurus 转换模板。

## 8. 下一里程碑

下一里程碑不是直接追逐 6.38 V，而是得到一个满足以下条件的 Vela 0 V 状态：混合材料和四电极边界正确、所有半导体载流子为正且有限、绝缘节点被合法屏蔽、端电流守恒。完成后再一次性导出 1/2/4/5/6 V 的无雪崩状态，进入字段和 IIC 算子级对比。
