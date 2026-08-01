#!/usr/bin/env python3
"""Move one DF-ISE vertex to create a controlled obtuse-mesh probe."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


VERTEX_BLOCK = re.compile(
    r"(?P<head>\n\s*Vertices\s*\((?P<count>\d+)\)\s*\{\s*\n)"
    r"(?P<body>.*?)"
    r"(?P<tail>\n\s*\})",
    re.DOTALL,
)


def triangle_angles(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    result = []
    for index in range(3):
        center = points[index]
        left = points[(index + 1) % 3]
        right = points[(index + 2) % 3]
        u = (left[0] - center[0], left[1] - center[1])
        v = (right[0] - center[0], right[1] - center[1])
        denom = math.hypot(*u) * math.hypot(*v)
        cosine = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / denom))
        result.append(math.degrees(math.acos(cosine)))
    return tuple(result)  # type: ignore[return-value]


def replace_vertex(
    text: str,
    vertex: int,
    x_um: float,
    y_um: float,
) -> tuple[str, tuple[float, float]]:
    match = VERTEX_BLOCK.search(text)
    if match is None:
        raise ValueError("DF-ISE Vertices block not found")
    lines = [line for line in match.group("body").splitlines() if line.strip()]
    count = int(match.group("count"))
    if len(lines) != count:
        raise ValueError(f"declared {count} vertices but parsed {len(lines)}")
    if vertex < 0 or vertex >= count:
        raise IndexError(f"vertex {vertex} outside [0, {count})")
    old = tuple(float(value) for value in lines[vertex].split()[:2])
    lines[vertex] = f" {x_um:.15e} {y_um:.15e}"
    body = "\n".join(lines)
    return (
        text[: match.start()]
        + match.group("head")
        + body
        + match.group("tail")
        + text[match.end() :],
        old,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--vertex", type=int, default=1)
    parser.add_argument("--x-um", type=float, default=0.05)
    parser.add_argument("--y-um", type=float, default=0.0)
    args = parser.parse_args()

    source = args.input.resolve()
    output = args.output.resolve()
    if source == output:
        raise ValueError("input and output must differ")
    text, old = replace_vertex(
        source.read_text(encoding="utf-8"),
        args.vertex,
        args.x_um,
        args.y_um,
    )
    output.write_text(text, encoding="utf-8")
    print(
        f"moved vertex {args.vertex}: ({old[0]:.15g}, {old[1]:.15g}) -> "
        f"({args.x_um:.15g}, {args.y_um:.15g})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
