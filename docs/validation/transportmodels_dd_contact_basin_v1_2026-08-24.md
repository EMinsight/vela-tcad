# TransportModels DD contact_basin 21点 Id-Vg 回归

## 结论

- 计算状态：完成，共21个固定比较点，另使用 126 个内部延续点。
- 主曲线验收：通过。
- 深关断数值验收：通过。
- 总体验收：通过。

## 计算方法说明

- 前 5 个固定点来自同方向的 `contact_basin` 单调扫描。
- 细化延续在 `Vg=-0.345 V` 附近因 `line_search_non_decrease` 停滞；后续 16 个固定点采用相同物理契约、相同偏压物理状态作为初值，在 `contact_basin` 坐标下独立重闭合。
- 21 个报告点均为真实求解点，未使用曲线插值；每点保留独立配置、最终状态和 KCL 证据。
- 因此，本结果可用于 21 点数值/物理回归，但不应解读为“一次不间断单调扫描已全程收敛”。

## 分区误差

| 区域 | 最大相对误差 | 最大对数误差 |
|---|---:|---:|
| 深关断 | 2.409% | 0.010339 dex |
| 过渡区 | 7.610% | 0.031852 dex |
| 导通区 | 2.930% | 0.012542 dex |

## 深关断前三点

| Vg (V) | Sentaurus Id (A/um) | Vela Id (A/um) | 相对误差 | 对数误差 (dex) | Id/abs(KCL) | 状态 |
|---:|---:|---:|---:|---:|---:|---|
| -1.00 | 1.634684e-15 | 1.664426e-15 | 1.819% | 0.007831 | 1264.915 | pass |
| -0.84 | 1.622317e-15 | 1.647763e-15 | 1.569% | 0.006759 | 1208.930 | pass |
| -0.68 | 1.738638e-15 | 1.780526e-15 | 2.409% | 0.010339 | 118.465 | pass |

## 固定证据

- 物理契约：`D:\code-repo\vela-tcad\configs\regression\transportmodels_dd_dg_sentaurus2022_v1.json`
- DD 数值契约：`D:\code-repo\vela-tcad\configs\regression\transportmodels_dd_contact_basin_v1.json`
- 运行配置：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_contact_basin_fixed_contract_v1_2026-08-24\config.json`
- 21点曲线：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_contact_basin_fixed_contract_v1_2026-08-24\dd_idvg_21_point.csv`
- 对齐结果：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_contact_basin_fixed_contract_v1_2026-08-24\dd_idvg_21_point_aligned.csv`
- 对比图：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_contact_basin_fixed_contract_v1_2026-08-24\dd_idvg_21_point_comparison.png`
