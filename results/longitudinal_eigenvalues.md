# Navion — trị riêng mô hình dọc, phiên bản đã hiệu chỉnh

## Phương pháp và nguồn

Hình học, trọng lượng, quán tính và đạo hàm tham khảo là dữ liệu độc lập. Các đại lượng dẫn xuất được tính bằng công thức. Hệ số hiệu dụng được suy từ đạo hàm Navion ở cấu hình chuẩn cố định, không tối ưu trị riêng. Khớp dữ liệu dùng để hiệu chỉnh không phải validation độc lập.
Nguồn: Nelson (1998), 2nd ed.; NASA CR-96008, Tables X-B/X-C; Seckel & Morris (1971), Table I chỉ dùng bổ sung hình học.
Nelson: Example 4.3, trang in 158 (PDF 169).
Hai hệ 4×4 chạy riêng; RK4 Δt=0.01 s, t=0…60 s, 6001 mẫu; trạng thái đầu [0,0,0.1,0]. Trị riêng tính trực tiếp từ A, không suy từ RK4.
Thứ tự trạng thái: [Δu,Δw,Δq,Δθ], đơn vị ft/s, ft/s, rad/s, rad.

## Ma trận mới

```text
[[-4.4805050582e-02  3.5698447894e-02  0.0000000000e+00 -3.2200000000e+01]
 [-3.6590909091e-01 -2.0241924240e+00  1.7600000000e+02 -0.0000000000e+00]
 [ 1.8899826936e-03 -3.9512090101e-02 -2.9857526359e+00  0.0000000000e+00]
 [ 0.0000000000e+00  0.0000000000e+00  1.0000000000e+00  0.0000000000e+00]]
```

## So sánh trị riêng

| Mode | Mô hình | Nelson công bố | Tính từ A in trong sách | Sai lệch |
|---|---:|---:|---:|---:|
| Phugoid | -0.017025513 +0.211967347i | -0.017100000 +0.213000000i | -0.017048749 +0.213544122i | 0.4845% |
| Short-period | -2.510349542 +2.591787828i | -2.500000000 +2.590000000i | -2.489451251 +2.597763769i | 0.2918% |

Đơn vị λ: s⁻¹. Với cặp phức, bảng dùng nghiệm có phần ảo dương; nghiệm kia là liên hợp. Sai lệch chuẩn hóa = |λ_mới−λ_N|/|λ_N| × 100%.
**Kết luận: mọi trị riêng có phần thực âm; ổn định tiệm cận.**

Hai cột Nelson được giữ riêng: (1) nghiệm tác giả công bố; (2) eig của các phần tử A đã làm tròn in trong sách. Chúng không bằng nhau hoàn toàn. Không ép code khớp một ma trận làm tròn hoặc các con số nội bộ không nhất quán.
Code dùng m=W/g, Q=ρu0²/2, CL0=W/(QS), khác các giá trị làm tròn m=85.4, Q=36.8, CL0=0.41. Các hệ số phụ thuộc CL0 thay đổi tương ứng; mode spiral rất nhạy với những chênh lệch nhỏ này.

## Giả thiết và giới hạn

Bỏ Zq và Z_w_dot theo Eq.(4.51); bỏ ảnh hưởng nén được theo cấu hình incompressible=True. eta_effective và eps_effective là tham số tương đương suy từ dữ liệu, không phải phép đo dòng tại đuôi. Điểm trung hòa tính từ Cm_alpha=0 của chính công thức mới với cánh tay đòn cập nhật theo CG.

## Tham số đầu vào

| Tham số đầu vào | Giá trị |
|---|---:|
| `S` | 184.0 |
| `b` | 33.4 |
| `c` | 5.7 |
| `S_t` | 43.0 |
| `AR_t` | 4.0 |
| `W` | 2750.0 |
| `g` | 32.2 |
| `rho` | 0.002378 |
| `u0` | 176.0 |
| `speed_of_sound` | 1117.0 |
| `I_y` | 3000.0 |
| `x_cg_c` | 0.295 |
| `x_ac_c` | 0.25 |
| `x_tail_ac` | 17.6815 |
| `gamma0_deg` | 0.0 |
| `a0_w_per_deg` | 0.097 |
| `a0_t_per_deg` | 0.1 |
| `incompressible` | True |
| `fuselage_dx` | (1.5, 1.5, 1.5, 1.5, 2.9, 2.9, 2.9, 2.9, 2.9) |
| `fuselage_width` | (3.0, 3.4, 3.8, 4.2, 3.8, 3.1, 2.3, 1.5, 0.8) |
| `fuselage_flow_gradient` | (1.2, 1.3, 1.4, 3.2, 0.06, 0.18, 0.31, 0.43, 0.55) |

## Các đại lượng tính bằng công thức

| Đại lượng đã tính | Giá trị |
|---|---:|
| `m` | 85.4037267081 |
| `Q` | 36.830464 |
| `QS` | 6776.805376 |
| `QSc` | 38627.7906432 |
| `AR` | 6.06282608696 |
| `b_t` | 13.1148770486 |
| `l_t` | 16 |
| `V_H` | 0.655987795576 |
| `a_w` | 4.30231969828 |
| `a_t` | 3.93529746444 |
| `eps_geometry` | 0.45175991323 |
| `fuselage_integral` | 169.84924 |
| `Cmf_geometry` | 0.25421472173 |
| `theta0` | 0 |
| `Mach` | 0.157564905998 |
| `CL0` | 0.4057959241 |
| `eta_effective` | 0.687243894881 |
| `eps_factor` | 0.968990366778 |
| `Cmf_factor` | 0.475564958452 |
| `CLalpha_residual` | -0.217679073278 |
| `e_effective` | 0.579239138994 |
| `CD_profile` | 0.0347635135135 |
| `eps_effective` | 0.437751004016 |
| `Cmf_effective` | 0.120895613577 |
| `CLalpha_geometry` | 4.65767907328 |
| `C_Lalpha` | 4.44 |
| `C_D0` | 0.0496891368835 |
| `C_Dalpha` | 0.326616231593 |
| `C_malpha` | -0.683 |
| `C_mq` | -9.96 |
| `C_malpha_dot` | -4.36 |
| `C_Du` | 0 |
| `C_mu` | 0 |
| `C_Lu_compressible` | 0.0103310596212 |
| `C_Lu` | 0 |
| `x_np_c` | 0.441639557869 |
| `SM` | 0.146639557869 |
| `C_malpha_geometry` | -0.967466573083 |
| `C_mq_geometry` | -14.4926714871 |
| `C_malpha_dot_geometry` | -6.54720801346 |
| `X_u` | -0.0448050505815 |
| `X_w` | 0.0356984478936 |
| `Z_u` | -0.365909090909 |
| `Z_w` | -2.02419242404 |
| `M_u` | 0 |
| `M_w` | -0.0499673882752 |
| `M_w_dot` | -0.0051651700944 |
| `M_q` | -2.07668269924 |
