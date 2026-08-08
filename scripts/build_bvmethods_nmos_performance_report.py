#!/usr/bin/env python3
"""Build the portable BVmethods field, I-V, and performance report input."""

from __future__ import annotations

import argparse
import base64
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_REPORT_DIR = Path(
    "build-release/reference_tcad/bvmethods_sentaurus2018/run01/"
    "sentaurus_boundary_state_20260808/report_20260808"
)


def image_html(path: Path, alt: str, caption: str) -> str:
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        '<figure style="margin:0">'
        f'<img alt="{alt}" src="data:image/png;base64,{encoded}" '
        'style="display:block;width:100%;height:auto;border-radius:8px">'
        f'<figcaption style="margin-top:8px;color:#667085">{caption}</figcaption>'
        "</figure>"
    )


def source(source_id: str, label: str, path: str,
           query: dict[str, object] | None = None) -> dict[str, object]:
    result: dict[str, object] = {"id": source_id, "label": label, "path": path}
    if query is not None:
        result["query"] = query
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    report_dir = args.report_dir.resolve()
    output = (args.output or report_dir / "artifact.json").resolve()
    perf = json.loads((report_dir / "performance_summary.json").read_text(encoding="utf-8"))

    sent = perf["sentaurus"]
    vela_current = perf["vela_current_switch"]
    vela_resistor = perf["vela_external_resistor"]
    comparison = perf["comparisons"]
    fields = perf["field_summary"]

    sent_reported = sent["reported_step_total_s"]
    phase_rows = [
        {"phase": "RHS / residual and assembly", "seconds": sent["reported_rhs_s"],
         "share_percent": 100.0 * sent["reported_rhs_s"] / sent_reported},
        {"phase": "Linear solve", "seconds": sent["reported_linear_solve_s"],
         "share_percent": 100.0 * sent["reported_linear_solve_s"] / sent_reported},
        {"phase": "Jacobian", "seconds": sent["reported_jacobian_s"],
         "share_percent": 100.0 * sent["reported_jacobian_s"] / sent_reported},
        {"phase": "Other", "seconds": sent_reported - sent["reported_rhs_s"]
         - sent["reported_linear_solve_s"] - sent["reported_jacobian_s"],
         "share_percent": 100.0 * (sent_reported - sent["reported_rhs_s"]
         - sent["reported_linear_solve_s"] - sent["reported_jacobian_s"]) / sent_reported},
    ]
    runtime_rows = [
        {"run": "Sentaurus full run", "wallclock_minutes": sent["wallclock_s"] / 60.0,
         "newton_updates": sent["newton_updates"],
         "seconds_per_update": sent["wallclock_s_per_update"]},
        {"run": "Vela current switch", "wallclock_minutes": vela_current["wallclock_s"] / 60.0,
         "newton_updates": vela_current["newton_updates"],
         "seconds_per_update": vela_current["wallclock_s_per_update"]},
        {"run": "Vela external resistor", "wallclock_minutes": vela_resistor["wallclock_s"] / 60.0,
         "newton_updates": vela_resistor["newton_updates"],
         "seconds_per_update": vela_resistor["wallclock_s_per_update"]},
    ]
    summary_rows = [{
        "sentaurus_wallclock_minutes": sent["wallclock_s"] / 60.0,
        "vela_current_wallclock_minutes": vela_current["wallclock_s"] / 60.0,
        "vela_resistor_wallclock_minutes": vela_resistor["wallclock_s"] / 60.0,
        "vela_current_speed_ratio": comparison["vela_current_wallclock_over_sentaurus"],
        "potential_p95_millivolts": fields["potential_abs_error_p95_V"] * 1000.0,
    }]

    sql_dir = Path(__file__).resolve().parent / "sql"
    runtime_sql = (sql_dir / "bvmethods_runtime_summary.sql").read_text(encoding="utf-8").strip()
    phase_sql = (sql_dir / "bvmethods_sentaurus_phase_summary.sql").read_text(encoding="utf-8").strip()
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE runtime_summary "
        "(sort_order INTEGER, run TEXT, wallclock_minutes REAL, "
        "newton_updates INTEGER, seconds_per_update REAL)"
    )
    connection.executemany(
        "INSERT INTO runtime_summary VALUES (?, ?, ?, ?, ?)",
        [(index, row["run"], row["wallclock_minutes"], row["newton_updates"],
          row["seconds_per_update"]) for index, row in enumerate(runtime_rows)],
    )
    connection.execute(
        "CREATE TABLE sentaurus_phase_summary "
        "(sort_order INTEGER, phase TEXT, seconds REAL, share_percent REAL)"
    )
    connection.executemany(
        "INSERT INTO sentaurus_phase_summary VALUES (?, ?, ?, ?)",
        [(index, row["phase"], row["seconds"], row["share_percent"])
         for index, row in enumerate(phase_rows)],
    )
    runtime_rows = [dict(zip(
        [column[0] for column in connection.execute(runtime_sql).description], row
    )) for row in connection.execute(runtime_sql).fetchall()]
    phase_rows = [dict(zip(
        [column[0] for column in connection.execute(phase_sql).description], row
    )) for row in connection.execute(phase_sql).fetchall()]
    connection.close()

    rel = "build-release/reference_tcad/bvmethods_sentaurus2018/run01/"
    rel += "sentaurus_boundary_state_20260808/report_20260808/"
    query_time = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    sources = [
        source("field_compare", "同网格场量逐节点比较", rel + "field_node_comparison.csv"),
        source("iv_compare", "Sentaurus 与 Vela I-V 对齐数据", rel + "iv_curve_comparison.csv"),
        source("runtime_summary", "运行日志汇总查询",
               "scripts/sql/bvmethods_runtime_summary.sql", {
                   "engine": "sqlite",
                   "sql": runtime_sql,
                   "description": "从解析后的三项运行日志汇总表读取墙钟时间和 Newton 成本。",
                   "executed_at": query_time,
                   "tables_used": ["runtime_summary"],
               }),
        source("sentaurus_phases", "Sentaurus 阶段耗时查询",
               "scripts/sql/bvmethods_sentaurus_phase_summary.sql", {
                   "engine": "sqlite",
                   "sql": phase_sql,
                   "description": "从 Sentaurus 日志解析后的阶段汇总表读取耗时占比。",
                   "executed_at": query_time,
                   "tables_used": ["sentaurus_phase_summary"],
               }),
        source("performance_summary", "性能解析与派生指标", rel + "performance_summary.json"),
    ]
    generated_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()

    def figure_block(block_id: str, filename: str, alt: str, caption: str,
                     source_id: str) -> dict[str, str]:
        return {
            "id": block_id,
            "type": "html",
            "sourceId": source_id,
            "body": image_html(report_dir / filename, alt, caption),
        }

    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "BVmethods NMOS：Sentaurus–Vela 场量、I-V 与性能对比",
            "description": "在 Id=1e-4 A/μm 电流边界及已有外接电阻验收结果上进行的可复核运行分析。",
            "generatedAt": generated_at,
            "charts": [
                {"id": "runtime_chart", "title": "观测到的端到端运行时间",
                 "subtitle": "Sentaurus 包含完整预偏置；两项 Vela 仅包含检查点后的高场求解。",
                 "type": "bar", "dataset": "runtime", "sourceId": "runtime_summary",
                 "encodings": {
                     "x": {"field": "run", "type": "nominal", "label": "运行"},
                     "y": {"field": "wallclock_minutes", "type": "quantitative",
                           "label": "墙钟时间 (min)"},
                 }, "yAxisTitle": "墙钟时间 (min)", "layout": "full"},
                {"id": "sentaurus_phase_chart", "title": "Sentaurus 已报告核心步骤耗时构成",
                 "subtitle": "RHS/残差与装配占 58.0%，线性求解占 33.0%。",
                 "type": "bar", "dataset": "sentaurus_phases", "sourceId": "sentaurus_phases",
                 "encodings": {
                     "x": {"field": "phase", "type": "nominal", "label": "阶段"},
                     "y": {"field": "share_percent", "type": "quantitative", "label": "占比 (%)"},
                 }, "yAxisTitle": "已报告核心步骤占比 (%)", "layout": "full"},
            ],
            "sources": sources,
            "blocks": [
                {"id": "title", "type": "markdown",
                 "body": "# BVmethods NMOS：Sentaurus–Vela 场量、I-V 与性能对比"},
                {"id": "technical_summary", "type": "markdown",
                 "body": "## 技术摘要\n\n在 **Id = 1e-4 A/μm** 的同一电流边界上，Vela 击穿电压为 **6.395904 V**，Sentaurus 新运行结果为 **6.384112 V**，相对误差 **0.1847%**。电势分布和高场热点位置一致，电势逐节点绝对误差 P95 为 **15.81 mV**。性能方面，Vela 电流切换恢复段为 **17.21 min**，外接电阻恢复段为 **64.75 min**；两者分别需要 **279** 和 **1356** 次 Newton 更新。首要瓶颈是高场非线性收敛与外层负载线重复调用完整 DD 求解，其次才是单次 Newton 更新成本。"},
                {"id": "potential_text", "type": "markdown",
                 "body": "## 电势分布保持高度一致\n\n两套结果在漏结弯曲区、耗尽区和接触附近呈现相同的电势拓扑。差异主要集中于器件底部高梯度区域；逐节点绝对误差中位数为 **0.764 mV**，P95 为 **15.81 mV**。"},
                figure_block("potential_figure", "potential_distribution_comparison.png",
                             "电势分布对比", "同网格电势、Vela 电势及逐节点差值。", "field_compare"),
                {"id": "field_text", "type": "markdown",
                 "body": "## 高场热点形状一致，峰值需按同一离散定义比较\n\n两者均把高场热点定位在漏结底部及弯曲结附近。图中 Sentaurus 使用导入的节点电场向量幅值，Vela 使用三角形电势梯度面积加权重构到节点；对应峰值分别为 **2.843e8** 与 **2.769e8 V/m**。该节点重构值不能与此前的边投影峰值直接混用。"},
                figure_block("field_figure", "electric_field_distribution_comparison.png",
                             "电场分布对比", "相同对数色标下的电场幅值和逐节点场强比。", "field_compare"),
                {"id": "iv_text", "type": "markdown",
                 "body": "## I-V 曲线在击穿区保持同一拐点\n\n从低漏电到雪崩区，两条曲线的数量级、斜率变化和击穿拐点相符。在判据电流 **1e-4 A/μm** 上，电压差为 **11.79 mV**。"},
                figure_block("iv_figure", "iv_curve_comparison.png",
                             "电压电流曲线对比", "全反偏扫描及击穿区局部放大。", "iv_compare"),
                {"id": "runtime_text", "type": "markdown",
                 "body": "## Vela 时间集中在高场 Newton 求解\n\nSentaurus 完整运行耗时 **103.24 s**，含 0→6 V 预偏置和电流边界；Vela 的 **1032.37 s** 与 **3885.10 s** 仅是检查点恢复后的高场段，因此原始倍数只能作为当前系统级观测。即便如此，Vela 每次 Newton 更新为 **2.87–3.70 s**，Sentaurus 全流程有效均值为 **0.279 s**，差约 **10.3–13.3 倍**。"},
                {"id": "runtime_chart_block", "type": "chart", "chartId": "runtime_chart", "layout": "full"},
                figure_block("newton_figure", "vela_newton_runtime_accumulation.png",
                             "Vela Newton 时间累积", "运行时间几乎全部随新计算的 Newton 状态数累积。", "performance_summary"),
                {"id": "sentaurus_phase_text", "type": "markdown",
                 "body": "## Sentaurus 内部以残差/装配和线性求解为主\n\nSentaurus 日志可直接拆分的 77.75 s 核心步骤中，RHS/残差与装配为 **45.13 s（58.0%）**，线性求解为 **25.67 s（33.0%）**，Jacobian 为 **5.58 s（7.2%）**。这提供了优化参照，但不能替代 Vela 自身的阶段计时。"},
                {"id": "sentaurus_phase_chart_block", "type": "chart", "chartId": "sentaurus_phase_chart", "layout": "full"},
                {"id": "scope", "type": "markdown",
                 "body": "## 数据范围、定义与比较口径\n\n场量使用同一器件网格、同一 1e-4 A/μm 电流边界。I-V 数据来自新 Sentaurus PLT、Vela 预偏置 CSV 与电流切换 CSV。Sentaurus 在虚拟机 O-2018.06-SP2 上运行并报告 4 CPU、1 个装配/求解线程；Vela 为 Windows UCRT64 Release、Eigen SparseLU。硬件、运行时与时间区间不同，所以速度倍数不是受控基准测试。"},
                {"id": "method", "type": "markdown",
                 "body": "## 方法\n\nSentaurus TDR 先导入 Vela 同网格格式；电势直接逐节点比较，Vela 电场由单元内电势梯度重构并按面积加权到节点。I-V 统一为漏极总电流绝对值与漏极电压。性能脚本解析三份运行日志、Newton 历史和 Sentaurus 每步阶段计时，并交叉核对外层边界评估次数。"},
                {"id": "limitations", "type": "markdown",
                 "body": "## 限制与稳健性\n\nVela 尚未对残差、雪崩源/Jacobian 装配、行缩放、SparseLU factorize/solve、线搜索和诊断分别计时，因此无法给出 Vela 内部阶段的可信百分比。现有日志可确认检查点 I/O 不是主要耗时，且最终接受状态的线搜索通常只尝试一次；对装配与稀疏分解的排序仍属于代码证据支持的待验证推断。"},
                {"id": "recommendations", "type": "markdown",
                 "body": "## 推荐优化顺序\n\n1. 先加入低开销分阶段计时，建立单次 Newton 的火焰图式账本。\n2. 优先减少外层完整 DD 调用：保留括区、割线/切线预测，并把外接电阻/电流接触写入分块或 Schur 联立约束。\n3. 降低高场 Newton 更新数：改进连续性缩放、初值预测和转折区 continuation；保持物理参数不变。\n4. 再优化单次迭代：缓存不变装配项、减少临时稀疏矩阵复制，评估更合适的稀疏直接/迭代求解器及数值分解复用条件。\n5. 使用同机、同线程、同预偏置范围做受控复测。"},
                {"id": "questions", "type": "markdown",
                 "body": "## 后续问题\n\n完成插桩后需要回答：高场每次 Newton 中残差/雪崩装配与 SparseLU 各占多少？外接电阻联立约束能否把 7 次完整 DD 评估降到 1 次耦合求解？相同线程和同扫描范围下，Vela 的真实端到端差距是多少？"},
            ],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "summary": summary_rows,
                "runtime": runtime_rows,
                "sentaurus_phases": phase_rows,
            },
        },
        "sources": sources,
        "package_info": {
            "root": "report_20260808",
            "manifestPath": "artifact.json",
            "snapshotPath": "artifact.json",
        },
    }
    artifact["snapshot"]["datasets"]["summary"][0]["breakdown_voltage_error_percent"] = 0.184717

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2), encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
