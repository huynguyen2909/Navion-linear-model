# Navion — trị riêng mô hình ngang-hướng, phiên bản đã hiệu chỉnh

## Phương pháp và nguồn

Hình học, trọng lượng, quán tính và đạo hàm tham khảo là dữ liệu độc lập. Các đại lượng dẫn xuất được tính bằng công thức. Hệ số hiệu dụng được suy từ đạo hàm Navion ở cấu hình chuẩn cố định, không tối ưu trị riêng. Khớp dữ liệu dùng để hiệu chỉnh không phải validation độc lập.
Nguồn: Nelson (1998), 2nd ed.; NASA CR-96008, Tables X-B/X-C; Seckel & Morris (1971), Table I chỉ dùng bổ sung hình học.
Nelson: Example 5.3, trang in 199–200 (PDF 210–211).
Hai hệ 4×4 chạy riêng; RK4 Δt=0.01 s, t=0…60 s, 6001 mẫu; trạng thái đầu [0,0,0.1,0]. Trị riêng tính trực tiếp từ A, không suy từ RK4.
Thứ tự trạng thái: [Δβ,Δp,Δr,Δφ], đơn vị rad, rad/s, rad/s, rad.

## Ma trận mới

```text
[[ -0.2542814196   0.            -1.             0.1829545455]
 [-15.9823971062  -8.4022936875   2.171254767    0.          ]
 [  4.4948457504  -0.3462519649  -0.7605199151   0.          ]
 [  0.             1.             0.             0.          ]]
```

## So sánh trị riêng

| Mode | Mô hình | Nelson công bố | Tính từ A in trong sách | Sai lệch |
|---|---:|---:|---:|---:|
| Dutch roll | -0.486070361 +2.333148611i | -0.487000000 +2.335000000i | -0.486162486 +2.333575284i | 0.0869% |
| Roll | -8.435807446 | -8.435000000 | -8.432762053 | 0.0096% |
| Spiral | -0.009146854 | -0.008770000 | -0.008912975 | 4.2971% |

Đơn vị λ: s⁻¹. Với cặp phức, bảng dùng nghiệm có phần ảo dương; nghiệm kia là liên hợp. Sai lệch chuẩn hóa = |λ_mới−λ_N|/|λ_N| × 100%.
**Kết luận: mọi trị riêng có phần thực âm; ổn định tiệm cận.**

Hai cột Nelson được giữ riêng: (1) nghiệm tác giả công bố; (2) eig của các phần tử A đã làm tròn in trong sách. Chúng không bằng nhau hoàn toàn. Không ép code khớp một ma trận làm tròn hoặc các con số nội bộ không nhất quán.
Code dùng m=W/g, Q=ρu0²/2, CL0=W/(QS), khác các giá trị làm tròn m=85.4, Q=36.8, CL0=0.41. Các hệ số phụ thuộc CL0 thay đổi tương ứng; mode spiral rất nhạy với những chênh lệch nhỏ này.

## Giả thiết và giới hạn

Tính đầy đủ Cy_p/Cy_r và Y_p_full/Y_r_full; đặt Yp=Yr=0 trong A khi neglect_Yp_Yr=True để dùng cùng xấp xỉ Nelson Example 5.3. Đặt False sẽ giữ hai số hạng này. K_beta, K_r và z_effective là các tham số tương đương suy từ bộ dữ liệu, không phải hình học đo. Chưa đủ dữ liệu để xác định riêng a_v, eta_v, sidewash; không gán tùy ý chúng.

## Tham số đầu vào

| Tham số đầu vào | Giá trị |
|---|---:|
| `S` | 184.0 |
| `b` | 33.4 |
| `c` | 5.7 |
| `S_v` | 12.5 |
| `W` | 2750.0 |
| `g` | 32.2 |
| `rho` | 0.002378 |
| `u0` | 176.0 |
| `Ix` | 1048.0 |
| `Iz` | 3530.0 |
| `x_cg_c` | 0.295 |
| `x_fin_ac` | 17.6815 |
| `wing_sweep_deg` | 2.9961111111111114 |
| `wing_dihedral_deg` | 7.5 |
| `wing_taper` | 0.54 |
| `a0_w_per_deg` | 0.097 |
| `gamma0_deg` | 0.0 |
| `neglect_Yp_Yr` | True |

## Các đại lượng tính bằng công thức

| Đại lượng đã tính | Giá trị |
|---|---:|
| `m` | 85.4037267081 |
| `Q` | 36.830464 |
| `QS` | 6776.805376 |
| `QS_b` | 226345.299558 |
| `AR` | 6.06282608696 |
| `l_v` | 16 |
| `V_v` | 0.0325436084353 |
| `theta0` | 0 |
| `sweep` | 0.0522920036445 |
| `dihedral` | 0.1308996939 |
| `a_w` | 4.30231969828 |
| `CL0` | 0.4057959241 |
| `K_beta` | 8.30208 |
| `K_r` | 4.00904375 |
| `z_effective` | 0.576 |
| `Cnb_wf_reference` | -0.200079640719 |
| `dihedral_gain` | -0.491013622869 |
| `Clp_factor` | 0.672174222223 |
| `Cnp_factor` | 1.12195121951 |
| `Cy_beta` | -0.564 |
| `Cn_beta` | 0.0701 |
| `Cn_beta_wf` | -0.200079640719 |
| `Cy_p` | 0.0149124703571 |
| `Cy_r` | 0.2609375 |
| `Cl_beta` | -0.074 |
| `Cl_p` | -0.41 |
| `Clp_geometry` | -0.609960909604 |
| `Cn_p` | -0.0569104039896 |
| `Cl_r` | 0.105948981025 |
| `Cn_r` | -0.125 |
| `Y_beta` | -44.7535298445 |
| `Y_p_full` | 0.112279789006 |
| `Y_r_full` | 1.96466492421 |
| `Y_p` | 0 |
| `Y_r` | 0 |
| `L_beta` | -15.9823971062 |
| `L_p` | -8.40229368752 |
| `L_r` | 2.17125476699 |
| `N_beta` | 4.49484575044 |
| `N_p` | -0.346251964866 |
| `N_r` | -0.760519915059 |
