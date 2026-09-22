"""Chạy python navion_reports.py để xuất đồ thị và hai báo cáo trị riêng."""
from dataclasses import asdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from navion_linear_models import simulate_models, spectrum_groups, OUTPUT_DIR

def plot_response(t, x, model, output_dir):
    longitudinal = model == "longitudinal"
    labels = ([r"$\Delta u$ (ft/s)", r"$\Delta w$ (ft/s)",
               r"$\Delta q$ (rad/s)", r"$\Delta\theta$ (rad)"] if longitudinal else
              [r"$\Delta\beta$ (rad)", r"$\Delta p$ (rad/s)",
               r"$\Delta r$ (rad/s)", r"$\Delta\phi$ (rad)"])
    name = "Chuyển động dọc" if longitudinal else "Chuyển động ngang-hướng"
    colors = ["#2563eb", "#16836b", "#dc6b25", "#8056b3"]
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)
    fig.suptitle(f"Navion | {name}\nĐáp ứng tự do 0–{t[-1]:g} s · RK4 · Δt = {t[1]-t[0]:g} s", fontsize=14)
    for i, ax in enumerate(axes):
        ax.plot(t, x[:, i], color=colors[i], lw=1.5)
        ax.set_ylabel(labels[i], fontsize=11)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
        ax.grid(True, alpha=0.25)
        ax.set_xlim(t[0], t[-1])
    axes[-1].set_xlabel("Thời gian (s)")
    note = "Hình học theo tài liệu; hệ số hiệu dụng hiệu chỉnh bằng dữ liệu Navion."
    fig.text(0.5, 0.018, note, ha="center", fontsize=9, color="#555555")
    fig.tight_layout(rect=(0, 0.035, 1, 0.935))
    for ext in ("svg", "png"):
        fig.savefig(output_dir / f"{model}_response_60s.{ext}", dpi=170)
    plt.close(fig)


NELSON_EIGENVALUES = {
    'longitudinal': {'Phugoid':complex(-0.0171,0.213), 'Short-period':complex(-2.5,2.59)},
    'lateral_directional': {'Dutch roll':complex(-0.487,2.335),'Roll':complex(-8.435,0),'Spiral':complex(-0.00877,0)}
}


NELSON_PRINTED_MATRIX = {
    'longitudinal':np.array([[-.045,.036,0,-32.2],[-.369,-2.02,176,0],
                             [.0019,-.0396,-2.948,0],[0,0,1,0]]),
    'lateral_directional':np.array([[-.254,0,-1,.182],[-16.02,-8.40,2.19,0],
                                    [4.488,-.350,-.760,0],[0,1,0,0]])
}


def fmt(z):
    if abs(z.imag)<1e-10:
        return f'{z.real:+.9f}'
    return f'{z.real:+.9f} {z.imag:+.9f}i'


def eigen_comparison(A, model):
    groups = spectrum_groups(A, model)
    printed = spectrum_groups(NELSON_PRINTED_MATRIX[model], model)
    lines = ["| Mode | Mô hình | Nelson công bố | Tính từ A in trong sách | Sai lệch |",
             "|---|---:|---:|---:|---:|"]
    for mode, z in groups.items():
        ref = NELSON_EIGENVALUES[model].get(mode)
        if ref is None:
            lines.append(f"| {mode} | {fmt(z)} | — | — | — |")
        else:
            error = 100 * abs(z - ref) / abs(ref)
            lines.append(f"| {mode} | {fmt(z)} | {fmt(ref)} | {fmt(printed[mode])} | {error:.4f}% |")
    return "\n".join(lines)

def parameter_table(p):
    lines=['| Tham số đầu vào | Giá trị |','|---|---:|']
    for key,v in asdict(p).items():
        lines.append(f'| `{key}` | {v} |')
    return '\n'.join(lines)


def calculated_table(c):
    lines=['| Đại lượng đã tính | Giá trị |','|---|---:|']
    for key,v in c.items():lines.append(f'| `{key}` | {v:.12g} |')
    return '\n'.join(lines)


def write_eigen_report(A,p,c,model,output_dir=OUTPUT_DIR):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    lon=model=='longitudinal';name='dọc' if lon else 'ngang-hướng'
    eigs=np.linalg.eigvals(A)
    text=[f'# Navion — trị riêng mô hình {name}, phiên bản đã hiệu chỉnh','',
          '## Phương pháp và nguồn','',
          'Hình học, trọng lượng, quán tính và đạo hàm tham khảo là dữ liệu độc lập. Các đại lượng dẫn xuất được tính bằng công thức. Hệ số hiệu dụng được suy từ đạo hàm Navion ở cấu hình chuẩn cố định, không tối ưu trị riêng. Khớp dữ liệu dùng để hiệu chỉnh không phải validation độc lập.',
          'Nguồn: Nelson (1998), 2nd ed.; NASA CR-96008, Tables X-B/X-C; Seckel & Morris (1971), Table I chỉ dùng bổ sung hình học.',
          f'Nelson: Example {"4.3, trang in 158 (PDF 169)" if lon else "5.3, trang in 199–200 (PDF 210–211)"}.',
          'Hai hệ 4×4 chạy riêng; RK4 Δt=0.01 s, t=0…60 s, 6001 mẫu; trạng thái đầu [0,0,0.1,0]. Trị riêng tính trực tiếp từ A, không suy từ RK4.',
          'Thứ tự trạng thái: '+('[Δu,Δw,Δq,Δθ], đơn vị ft/s, ft/s, rad/s, rad.' if lon else '[Δβ,Δp,Δr,Δφ], đơn vị rad, rad/s, rad/s, rad.'),
          '', '## Ma trận mới','', '```text',np.array2string(A,precision=10),'```','',
          '## So sánh trị riêng','',eigen_comparison(A,model),'',
          'Đơn vị λ: s⁻¹. Với cặp phức, bảng dùng nghiệm có phần ảo dương; nghiệm kia là liên hợp. Sai lệch chuẩn hóa = |λ_mới−λ_N|/|λ_N| × 100%.',
          '**Kết luận: '+('mọi trị riêng có phần thực âm; ổn định tiệm cận.' if np.all(eigs.real<0) else 'có trị riêng không âm; kiểm tra ổn định và cách nhận dạng mode.')+'**','',
          'Hai cột Nelson được giữ riêng: (1) nghiệm tác giả công bố; (2) eig của các phần tử A đã làm tròn in trong sách. Chúng không bằng nhau hoàn toàn. Không ép code khớp một ma trận làm tròn hoặc các con số nội bộ không nhất quán.',
          'Code dùng m=W/g, Q=ρu0²/2, CL0=W/(QS), khác các giá trị làm tròn m=85.4, Q=36.8, CL0=0.41. Các hệ số phụ thuộc CL0 thay đổi tương ứng; mode spiral rất nhạy với những chênh lệch nhỏ này.',
          '', '## Giả thiết và giới hạn','',
          ('Bỏ Zq và Z_w_dot theo Eq.(4.51); bỏ ảnh hưởng nén được theo cấu hình incompressible=True. eta_effective và eps_effective là tham số tương đương suy từ dữ liệu, không phải phép đo dòng tại đuôi. Điểm trung hòa tính từ Cm_alpha=0 của chính công thức mới với cánh tay đòn cập nhật theo CG.' if lon else
           'Tính đầy đủ Cy_p/Cy_r và Y_p_full/Y_r_full; đặt Yp=Yr=0 trong A khi neglect_Yp_Yr=True để dùng cùng xấp xỉ Nelson Example 5.3. Đặt False sẽ giữ hai số hạng này. K_beta, K_r và z_effective là các tham số tương đương suy từ bộ dữ liệu, không phải hình học đo. Chưa đủ dữ liệu để xác định riêng a_v, eta_v, sidewash; không gán tùy ý chúng.'),
          '', '## Tham số đầu vào','',parameter_table(p),'',
          '## Các đại lượng tính bằng công thức','',calculated_table(c),'']
    (output_dir/f'{model}_eigenvalues.md').write_text('\n'.join(text),encoding='utf-8')



def main(output_dir=OUTPUT_DIR):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    for model, result in simulate_models().items():
        plot_response(result["t"], result["x"], model, output_dir)
        write_eigen_report(result["A"], result["parameters"],
                           result["coefficients"], model, output_dir)
    print("Results:", output_dir)


if __name__ == "__main__":
    main()
