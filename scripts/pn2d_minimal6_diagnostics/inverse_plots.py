"""Deterministic static figures for the Minimal6 physics inverse audit."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import textwrap

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("SOURCE_DATE_EPOCH", "0")

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


FIGURE_NAMES = (
    "potential_field", "qf_gradient", "current_density",
    "alpha_generation", "replacement_matrix",
)
BLUE, ORANGE, INK, GREY, LIGHT_GREY = (
    "#2563A6", "#C66A1B", "#263238", "#8A949B", "#D8DEE3",
)
PALETTE = (BLUE, ORANGE, INK, GREY)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        return hashlib.sha256(image.convert("RGBA").tobytes()).hexdigest()


def _text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, *, fill: str = INK) -> None:
    draw.text(xy, value, fill=fill, font=ImageFont.load_default())


def _render_lines(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    x = [float(value) for value in panel["x"]]
    all_y = [float(value) for series in panel["series"] for value in series["values"]]
    if not x or not all_y:
        _text(draw, (left + 10, top + 10), "No compatible finite samples", fill=GREY)
        return
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(0.0, min(all_y)), max(0.0, max(all_y))
    if ymax == ymin:
        ymax = ymin + 1.0
    for fraction in (0.25, 0.5, 0.75):
        y_grid = int(bottom - fraction * (bottom - top))
        draw.line((left, y_grid, right, y_grid), fill=LIGHT_GREY, width=1)
    for index, series in enumerate(panel["series"]):
        color = PALETTE[index % len(PALETTE)]
        points = []
        for xv, yv in zip(x, series["values"]):
            px = left if xmax == xmin else int(left + (xv - xmin) / (xmax - xmin) * (right - left))
            py = int(bottom - (float(yv) - ymin) / (ymax - ymin) * (bottom - top))
            points.append((px, py))
        for first, second in zip(points, points[1:]):
            if index % 2:
                midpoint = ((first[0] + second[0]) // 2, (first[1] + second[1]) // 2)
                draw.line((first, midpoint), fill=color, width=2)
            else:
                draw.line((first, second), fill=color, width=2)
        for px, py in points:
            radius = 3 + index % 2
            shape = (px - radius, py - radius, px + radius, py + radius)
            if index % 2:
                draw.rectangle(shape, fill="white", outline=color)
            else:
                draw.ellipse(shape, fill="white", outline=color)
        _text(draw, (right - 190, top + 4 + index * 14), series["label"], fill=color)


def _render_bars(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    values = [float(value) for value in panel["values"]]
    extent = max((abs(value) for value in values), default=1.0) or 1.0
    label_width = 185
    zero = left + label_width + (right - left - label_width) // 2
    draw.line((zero, top, zero, bottom), fill=INK, width=1)
    height = max(12, (bottom - top) // max(1, len(values)))
    for index, (label, value) in enumerate(zip(panel["labels"], values)):
        y = top + index * height + 3
        length = int((right - left - label_width) * 0.45 * value / extent)
        x0, x1 = sorted((zero, zero + length))
        draw.rectangle((x0, y, x1, y + height - 6), fill=BLUE if value >= 0 else "white", outline=BLUE)
        if value < 0:
            for hatch in range(x0, x1, 6):
                draw.line((hatch, y + height - 6, min(hatch + height, x1), y), fill=BLUE)
        _text(draw, (left, y), str(label))
        _text(draw, (x1 + 4, y), f"{value:.3g}")


def _render_panel(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    _text(draw, (left, top - 18), panel["title"])
    draw.line((left, top, left, bottom), fill=INK, width=1)
    draw.line((left, bottom, right, bottom), fill=INK, width=1)
    _text(draw, (left, bottom + 7), panel.get("x_label", ""), fill=GREY)
    _text(draw, (left, top + 3), panel.get("y_label", ""), fill=GREY)
    if panel.get("kind") == "bar":
        _render_bars(draw, panel, box)
    else:
        _render_lines(draw, panel, box)


def render_inverse_figures(out_dir: str | Path, figure_specs: dict[str, dict]) -> dict:
    """Render the fixed five PNG/PDF pairs and return their chart map."""
    if tuple(figure_specs) != FIGURE_NAMES:
        raise ValueError("figure specifications must use the fixed deterministic order")
    root = Path(out_dir)
    figure_dir = root / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    for name in FIGURE_NAMES:
        spec = figure_specs[name]
        panels = spec["panels"]
        width, height = 864, 456 + 210 * (len(panels) - 1)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)
        _text(draw, (45, 24), spec["title"])
        for line_index, line in enumerate(textwrap.wrap(spec["subtitle"], width=115)[:2]):
            _text(draw, (45, 43 + line_index * 13), line, fill=GREY)
        panel_top = 95
        panel_height = (height - panel_top - 55) // len(panels)
        for index, panel in enumerate(panels):
            top = panel_top + index * panel_height + 15
            _render_panel(draw, panel, (80, top, width - 45, top + panel_height - 55))
        png, pdf = figure_dir / f"{name}.png", figure_dir / f"{name}.pdf"
        image.save(png, format="PNG", optimize=False, compress_level=9, dpi=(120, 120))
        pdf_canvas = canvas.Canvas(str(pdf), pagesize=(width, height), invariant=1, pageCompression=1)
        pdf_canvas.setTitle(spec["title"])
        pdf_canvas.setAuthor("Vela TCAD")
        pdf_canvas.setSubject(spec["subtitle"])
        pdf_canvas.setCreator("Vela TCAD deterministic inverse audit")
        pdf_canvas.drawImage(ImageReader(image), 0, 0, width=width, height=height,
                             preserveAspectRatio=True, mask=None)
        pdf_canvas.showPage()
        pdf_canvas.save()
        contract = dict(spec["chart_contract"])
        contract["output_paths"] = [f"figures/{name}.png", f"figures/{name}.pdf"]
        entries.append({
            "name": name, "title": spec["title"], "subtitle": spec["subtitle"],
            "png_sha256": _sha256(png), "png_pixel_sha256": pixel_sha256(png),
            "pdf_sha256": _sha256(pdf), "chart_contract": contract,
        })
    return {
        "schema": "vela.pn2d_minimal6_inverse_figure_manifest.v1",
        "renderer": {"backend": "Agg-static fallback", "raster": "Pillow",
                     "pdf": "ReportLab invariant mode", "dpi": 120,
                     "font": "Pillow deterministic default bitmap"},
        "palette_roots": {"primary": BLUE, "secondary": ORANGE,
                          "neutrals": [INK, GREY, LIGHT_GREY]},
        "figures": entries,
    }


def write_figure_manifest(path: str | Path, manifest: dict) -> None:
    Path(path).write_text(json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8", newline="\n")
