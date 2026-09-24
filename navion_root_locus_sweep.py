"""Quét Vv, Vh, SM từ cấu hình Navion chuẩn và xuất sáu đồ thị trị riêng.

Chạy: python navion_root_locus_sweep.py
Phụ thuộc: numpy, matplotlib và navion_linear_models.py trong cùng thư mục.
"""

from dataclasses import replace
from itertools import permutations
from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
import numpy as np

from navion_linear_models import (
    LateralParameters,
    LongitudinalParameters,
    build_lateral_matrix,
    build_longitudinal_matrix,
    lateral_geometry,
    longitudinal_aerodynamics,
    longitudinal_geometry,
)


FRACTIONS = np.arange(11, dtype=float) * 0.05
OUT_DIR = Path(__file__).resolve().parent / "results" / "root_locus_sweep"
MODELS = ("longitudinal", "lateral_directional")
MODE_LABELS = {
    "longitudinal": ("Phugoid +", "Phugoid −", "Short-period +", "Short-period −"),
    "lateral_directional": ("Dutch roll +", "Dutch roll −", "Roll", "Spiral"),
}
TITLES = {"Vv": "$V_v$", "Vh": "$V_H$", "SM": "SM"}


def starting_roots(matrix, model):
    roots = np.linalg.eigvals(matrix)
    oscillatory = sorted((z for z in roots if z.imag > 1e-8), key=abs)
    if model == "longitudinal" and len(oscillatory) == 2:
        return np.array([oscillatory[0], oscillatory[0].conjugate(),
                         oscillatory[1], oscillatory[1].conjugate()])
    real = sorted((z for z in roots if abs(z.imag) <= 1e-8), key=lambda z: z.real)
    if model == "lateral_directional" and len(oscillatory) == 1 and len(real) == 2:
        return np.array([oscillatory[0], oscillatory[0].conjugate(), real[0], real[1]])
    raise ValueError(f"Dạng nghiệm baseline không khớp 4 mode dự kiến: {model}: {roots}")


def track_roots(previous, matrix):
    roots = np.linalg.eigvals(matrix)
    scale = np.maximum(1.0, np.abs(previous))
    order = min(permutations(range(4)),
                key=lambda perm: sum(abs(previous[j] - roots[i]) / scale[j]
                                     for j, i in enumerate(perm)))
    return roots[list(order)]


def cases():
    lon0, lat0 = LongitudinalParameters(), LateralParameters()
    geom_lon0 = longitudinal_geometry(lon0)
    geom_lat0 = lateral_geometry(lat0)
    sm0 = longitudinal_aerodynamics(lon0)["SM"]
    if sm0 <= 0:
        raise ValueError("SM chuẩn cần dương để áp dụng quét tăng theo phần trăm.")
    h_np = lon0.x_cg_c + sm0

    for parameter in ("Vv", "Vh", "SM"):
        for step, fraction in enumerate(FRACTIONS):
            if parameter == "Vv":
                lon = lon0
                lat = replace(lat0, S_v=lat0.S_v * (1 + fraction))
                nominal, initial = lateral_geometry(lat)["V_v"], geom_lat0["V_v"]
            elif parameter == "Vh":
                lon = replace(lon0, S_t=lon0.S_t * (1 + fraction))
                lat = lat0
                nominal, initial = longitudinal_geometry(lon)["V_H"], geom_lon0["V_H"]
            else:
                nominal, initial = sm0 * (1 + fraction), sm0
                cg_c = h_np - nominal
                lon = replace(lon0, x_cg_c=cg_c)
                lat = replace(lat0, x_cg_c=cg_c)
                actual = longitudinal_aerodynamics(lon)["SM"]
                if not np.isclose(actual, nominal, rtol=1e-10, atol=1e-12):
                    raise ValueError(f"SM thực không khớp SM đặt tại bước {step}.")

            if not np.isclose(nominal, initial * (1 + fraction), rtol=1e-10):
                raise ValueError(f"Bước quét {parameter} không đúng {fraction:.0%}.")
            yield parameter, step, nominal, lon, lat


def sweep():
    grouped = {(p, m): [] for p in TITLES for m in MODELS}
    for parameter, step, value, lon, lat in cases():
        for model, matrix in (
            ("longitudinal", build_longitudinal_matrix(lon)[0]),
            ("lateral_directional", build_lateral_matrix(lat)[0]),
        ):
            series = grouped[(parameter, model)]
            roots = starting_roots(matrix, model) if step == 0 else track_roots(series[-1][2], matrix)
            series.append((step, value, roots))
    return grouped


def plot_locus(parameter, model, records, output_dir):
    roots = np.stack([row[2] for row in records])
    unchanged = np.allclose(roots, roots[0], rtol=0, atol=1e-12)
    fig, ax = plt.subplots(figsize=(10.5, 7))
    colors = ("#0f766e", "#2563eb", "#e07a19", "#8b5cf6")
    for i, label in enumerate(MODE_LABELS[model]):
        branch = roots[:, i]
        ax.plot(branch.real, branch.imag, color=colors[i], lw=1.7,
                marker="o", markersize=3, label=label)
        ax.scatter(branch[0].real, branch[0].imag, s=66, marker="o",
                   facecolor="white", edgecolor=colors[i], lw=1.8, zorder=5)
        ax.scatter(branch[-1].real, branch[-1].imag, s=60, marker="x",
                   color=colors[i], lw=2, zorder=6)
    ax.axvline(0, color="black", lw=0.8, alpha=0.55)
    ax.axhline(0, color="black", lw=0.8, alpha=0.55)
    ax.grid(True, alpha=0.22)
    ax.set_xlabel(r"$\operatorname{Re}(\lambda)$ (s$^{-1}$)")
    ax.set_ylabel(r"$\operatorname{Im}(\lambda)$ (s$^{-1}$)")
    channel = "dọc" if model == "longitudinal" else "ngang-hướng"
    ax.set_title(f"Navion | {channel} | thay đổi {TITLES[parameter]} từ 0 đến +50%")
    ax.legend(loc="best", fontsize=9, ncol=2)
    ax.text(0.01, 0.01,
            "Không đổi trong ma trận của kênh này" if unchanged else
            "○ Cấu hình chuẩn (0%)     × Bước 10 (+50%)     Mỗi bước +5% so với giá trị chuẩn",
            transform=ax.transAxes, va="bottom", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="#dadada", alpha=0.90))
    ax.margins(x=0.12, y=0.16)
    fig.tight_layout()
    if not unchanged:
        index = 0 if model == "longitudinal" else 3
        placement = (0.57, 0.54, 0.33, 0.27) if model == "longitudinal" else (0.24, 0.54, 0.34, 0.27)
        detail = fig.add_axes(placement)
        branch = roots[:, index]
        detail.plot(branch.real, branch.imag, color=colors[index],
                    lw=1.5, marker="o", markersize=3)
        detail.scatter(branch[0].real, branch[0].imag, s=55, marker="o",
                       facecolor="white", edgecolor=colors[index], lw=1.5, zorder=5)
        detail.scatter(branch[-1].real, branch[-1].imag, s=60, marker="x",
                       color=colors[index], lw=1.8, zorder=6)
        if model == "lateral_directional":
            detail.axvline(0, color="black", lw=0.8, linestyle="--")
            detail.set_ylim(-max(0.002, np.ptp(branch.real) * 0.1),
                            max(0.002, np.ptp(branch.real) * 0.1))
        detail.set_title("Phóng to: " + MODE_LABELS[model][index], fontsize=9)
        detail.set_xlabel(r"$\operatorname{Re}(\lambda)$", fontsize=8)
        detail.set_ylabel(r"$\operatorname{Im}(\lambda)$", fontsize=8)
        detail.tick_params(labelsize=8)
        detail.grid(alpha=0.25)
        detail.margins(x=0.15, y=0.3)
    stem = f"root_locus_{parameter}_{model}"
    for extension in ("svg", "png"):
        fig.savefig(output_dir / f"{stem}.{extension}", dpi=180)
    plt.close(fig)
    return unchanged


def save_values(parameter, model, records, output_dir):
    path = output_dir / f"root_locus_{parameter}_{model}.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        columns = [f"lambda_{label.lower().replace(' ', '_')}_{part}"
                   for label in MODE_LABELS[model] for part in ("real_s-1", "imag_s-1")]
        writer.writerow(["step", "delta_from_baseline_pct", parameter, *columns])
        for step, value, roots in records:
            writer.writerow([step, FRACTIONS[step] * 100, value,
                             *(part for z in roots for part in (z.real, z.imag))])


def main(output_dir=OUT_DIR):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped = sweep()
    for parameter in TITLES:
        for model in MODELS:
            records = grouped[(parameter, model)]
            unchanged = plot_locus(parameter, model, records, output_dir)
            save_values(parameter, model, records, output_dir)
            print(f"{parameter:>2} {model:>19}: {records[0][1]:.9g} → "
                  f"{records[-1][1]:.9g}; không đổi={unchanged}")
    print("Đồ thị và CSV:", output_dir)
    return grouped


if __name__ == "__main__":
    main()
