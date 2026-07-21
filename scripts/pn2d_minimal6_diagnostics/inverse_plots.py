"""Deterministic Pillow/ReportLab figures for the Minimal6 inverse audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import textwrap

import PIL
from PIL import Image, ImageDraw, ImageFont
import reportlab
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


FIGURE_NAMES = (
    "potential_field", "qf_gradient", "current_density",
    "alpha_generation", "replacement_matrix",
)
BLUE, ORANGE, INK, GREY, LIGHT_GREY = (
    "#2563A6", "#C66A1B", "#263238", "#66737B", "#D8DEE3",
)
PALETTE = (BLUE, ORANGE, INK, GREY)
_FONT_ROOT = Path(reportlab.__file__).resolve().parent / "fonts"
FONT = ImageFont.truetype(str(_FONT_ROOT / "Vera.ttf"), 12)
FONT_BOLD = ImageFont.truetype(str(_FONT_ROOT / "VeraBd.ttf"), 13)
FONT_SMALL = ImageFont.truetype(str(_FONT_ROOT / "Vera.ttf"), 10)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        return hashlib.sha256(image.convert("RGBA").tobytes()).hexdigest()


def _text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, *,
          fill: str = INK, bold: bool = False, small: bool = False) -> None:
    font = FONT_BOLD if bold else FONT_SMALL if small else FONT
    draw.text(xy, value, fill=fill, font=font)


def _finite(value: object) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number == number and abs(number) != float("inf")


def _line_marker(draw: ImageDraw.ImageDraw, x: int, y: int, color: str, index: int) -> None:
    draw.line((x, y, x + 27, y), fill=color, width=2)
    if index % 2:
        draw.rectangle((x + 11, y - 4, x + 19, y + 4), fill="white", outline=color, width=1)
    else:
        draw.ellipse((x + 11, y - 4, x + 19, y + 4), fill="white", outline=color, width=1)


def _render_lines(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    x = [float(value) for value in panel["x"]]
    all_y = [float(value) for series in panel["series"] for value in series["values"]
             if _finite(value)]
    if not x or not all_y:
        _text(draw, (left + 12, top + 28), "No compatible finite samples", fill=GREY)
        return
    xmin, xmax = min(x), max(x)
    ymin, ymax = min(0.0, min(all_y)), max(0.0, max(all_y))
    if ymax == ymin:
        ymax = ymin + 1.0
    for fraction in (0.0, 0.5, 1.0):
        py = int(bottom - fraction * (bottom - top))
        draw.line((left, py, right, py), fill=LIGHT_GREY, width=1)
        value = ymin + fraction * (ymax - ymin)
        _text(draw, (4, py - 6), f"{value:.3g}", fill=GREY, small=True)
    tick_indices = sorted({0, len(x) // 2, len(x) - 1})
    for position in tick_indices:
        px = left if xmax == xmin else int(left + (x[position] - xmin) / (xmax - xmin) * (right - left))
        draw.line((px, bottom, px, bottom + 4), fill=INK, width=1)
        _text(draw, (px - 10, bottom + 6), f"{x[position]:g}", fill=GREY, small=True)
    for index, series in enumerate(panel["series"]):
        color = PALETTE[index % len(PALETTE)]
        points: list[tuple[int, int] | None] = []
        for xv, yv in zip(x, series["values"]):
            if not _finite(yv):
                points.append(None)
                continue
            px = left if xmax == xmin else int(left + (xv - xmin) / (xmax - xmin) * (right - left))
            py = int(bottom - (float(yv) - ymin) / (ymax - ymin) * (bottom - top))
            points.append((px, py))
        for first, second in zip(points, points[1:]):
            if first is None or second is None:
                continue
            if index % 2:
                midpoint = ((first[0] + second[0]) // 2, (first[1] + second[1]) // 2)
                draw.line((first, midpoint), fill=color, width=2)
            else:
                draw.line((first, second), fill=color, width=2)
        for point in points:
            if point is None:
                continue
            px, py = point
            shape = (px - 4, py - 4, px + 4, py + 4)
            if index % 2:
                draw.rectangle(shape, fill="white", outline=color)
            else:
                draw.ellipse(shape, fill="white", outline=color)
        legend_x, legend_y = right - 250, top + 8 + index * 19
        _line_marker(draw, legend_x, legend_y + 5, color, index)
        _text(draw, (legend_x + 35, legend_y), series["label"], fill=color, small=True)


def _render_bars(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    values = [float(value) for value in panel["values"]]
    extent = max((abs(value) for value in values), default=1.0) or 1.0
    label_width = 205
    zero = left + label_width
    draw.line((zero, top, zero, bottom), fill=INK, width=1)
    height = max(14, (bottom - top) // max(1, len(values)))
    for index, (label, value) in enumerate(zip(panel["labels"], values)):
        y = top + index * height + 3
        length = int((right - zero - 45) * value / extent)
        x0, x1 = sorted((zero, zero + length))
        draw.rectangle((x0, y, x1, y + height - 6), fill=BLUE if value >= 0 else "white", outline=BLUE)
        if value < 0:
            for hatch in range(x0, x1, 6):
                draw.line((hatch, y + height - 6, min(hatch + height, x1), y), fill=BLUE)
        _text(draw, (left, y), str(label), small=True)
        _text(draw, (x1 + 4, y), f"{value:.3g}", small=True)
    for fraction in (0.0, 0.5, 1.0):
        px = int(zero + fraction * (right - zero - 45))
        _text(draw, (px - 8, bottom + 6), f"{fraction * extent:.3g}", fill=GREY, small=True)


def _render_panel(draw: ImageDraw.ImageDraw, panel: dict, box: tuple[int, int, int, int]) -> None:
    left, top, right, bottom = box
    _text(draw, (left, top - 23), panel["title"], bold=True)
    draw.line((left, top, left, bottom), fill=INK, width=1)
    draw.line((left, bottom, right, bottom), fill=INK, width=1)
    _text(draw, (left + 190, bottom + 25), panel.get("x_label", ""), fill=GREY, small=True)
    _text(draw, (left + 4, top + 4), panel.get("y_label", ""), fill=GREY, small=True)
    (_render_bars if panel.get("kind") == "bar" else _render_lines)(draw, panel, box)


def render_inverse_figures(out_dir: str | Path, figure_specs: dict[str, dict]) -> dict:
    if tuple(figure_specs) != FIGURE_NAMES:
        raise ValueError("figure specifications must use the fixed deterministic order")
    figure_dir = Path(out_dir) / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    entries = []
    for name in FIGURE_NAMES:
        spec, panels = figure_specs[name], figure_specs[name]["panels"]
        width, height = 960, 520 + 240 * (len(panels) - 1)
        image = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(image)
        _text(draw, (45, 24), spec["title"], bold=True)
        for line_index, line in enumerate(textwrap.wrap(spec["subtitle"], width=118)[:2]):
            _text(draw, (45, 49 + line_index * 15), line, fill=GREY, small=True)
        panel_top = 105
        panel_height = (height - panel_top - 65) // len(panels)
        for index, panel in enumerate(panels):
            top = panel_top + index * panel_height + 18
            _render_panel(draw, panel, (105, top, width - 45, top + panel_height - 65))
        png, pdf = figure_dir / f"{name}.png", figure_dir / f"{name}.pdf"
        image.save(png, format="PNG", optimize=False, compress_level=9, dpi=(120, 120))
        output = canvas.Canvas(str(pdf), pagesize=(width, height), invariant=1, pageCompression=1)
        output.setTitle(spec["title"]); output.setAuthor("Vela TCAD")
        output.setSubject(spec["subtitle"]); output.setCreator("Vela TCAD deterministic inverse audit")
        output.drawImage(ImageReader(image), 0, 0, width=width, height=height, mask=None)
        output.showPage(); output.save()
        contract = dict(spec["chart_contract"])
        contract["output_paths"] = [f"figures/{name}.png", f"figures/{name}.pdf"]
        entries.append({
            "name": name, "title": spec["title"], "subtitle": spec["subtitle"],
            "png_sha256": _sha256(png), "png_pixel_sha256": pixel_sha256(png),
            "pdf_sha256": _sha256(pdf), "chart_contract": contract,
        })
    return {
        "schema": "vela.pn2d_minimal6_inverse_figure_manifest.v1",
        "renderer": {
            "backend": "Pillow PNG + ReportLab invariant PDF", "dpi": 120,
            "font": "ReportLab Vera TrueType", "pillow_version": PIL.__version__,
            "reportlab_version": reportlab.Version,
            "determinism_scope": "same runtime and library versions",
        },
        "palette_roots": {"primary": BLUE, "secondary": ORANGE,
                          "neutrals": [INK, GREY, LIGHT_GREY]},
        "figures": entries,
    }


def write_figure_manifest(path: str | Path, manifest: dict) -> None:
    Path(path).write_text(json.dumps(manifest, sort_keys=True, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8", newline="\n")
