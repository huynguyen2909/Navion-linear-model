# Công thức tính mô hình tuyến tính Navion

*Trang dưới đây là số trang in của tài liệu. Nelson: Flight Stability and Automatic Control, 2nd ed., 1998. Dấu * chỉ đạo hàm tham chiếu; r chỉ cấu hình chuẩn. “Suy đại số/hiệu chỉnh trong code” không phải phương trình do Nelson công bố.*

## 1. Hình học và các đại lượng dẫn xuất

| Đại lượng | Công thức | Nguồn tham chiếu |
|---|---|---|
| Khối lượng | m=W/g | Định nghĩa trọng lượng; dữ liệu NASA CR-96008, tr.101, Table X-A (không đánh số phương trình). |
| Áp suất động | Q=ρu0²/2; QS=Q·S; QSc=QS·c; QSb=QS·b | Định nghĩa áp suất động; Nelson, tr.150, Table 4.2; tr.199, Table 5.1. |
| Diện tích cánh | S=184 ft² | Nelson, tr.57, Fig.2.16; NASA CR-96008, tr.100 (dữ liệu hình học). |
| Vị trí CG | x_cg=h·c; l_t=x_tail_ac−x_cg; l_v=x_fin_ac−x_cg | Quan hệ hình học trong code; Nelson, tr.57, Fig.2.16 (cánh tay đòn chuẩn). |
| Aspect ratio và sải đuôi | AR_w=b²/S; b_t=sqrt(AR_t·S_t), với AR_t là tỷ số hình học độc lập nguồn cho | Định nghĩa aspect ratio; Nelson, tr.57, Example 2.2; Seckel & Morris (1971), tr.4, Table I. |
| Hệ số thể tích | V_H=S_t l_t/(Sc); V_v=S_v l_v/(Sb) | Nelson, tr.56, Eq.(2.35) và tr.76, Eq.(2.78)–(2.79): định nghĩa V_H, V_v. |
| Lift-curve slope cánh/đuôi | a=a0/[1+a0/(πAR)], a0 đổi từ per-degree sang per-radian | Nelson, tr.57, Example 2.2, công thức không đánh số. |
| CL0 | CL0=W cosθ0/(QS), giả thiết lực nâng cân bằng thành phần trọng lượng | Cân bằng lực nâng trong code; Nelson, tr.57–58, Example 2.2; Table 3.3, tr.116. |
| CD0, CDalpha | Polar CD=CD_profile+CL²/(π e AR); CDalpha=2 CL CLalpha/(π e AR) | Nelson, tr.116, Table 3.3; suy đại số polar trong code. |
| Downwash | ε_geom=2a_w/(πAR); hiệu chỉnh bằng gain suy từ C_madot/C_mq chuẩn | Nelson, tr.58, Example 2.2, công thức downwash không đánh số; gain do code hiệu chỉnh. |
| Fuselage Cm_alpha | Cmf_eff=Gain_Multhopp·Cmf_geometry (tổng ở bảng 2) | Nelson, tr.53, Eq.(2.32); tr.61, Fig.2.17; gain do code hiệu chỉnh. |
| C_mq và C_madot | C_mq=−2a_t eta_eff V_H l_t/c; C_madot=C_mq ε_eff | Nelson, tr.113, Eq.(3.79); tr.115, Eq.(3.92); eta hiệu dụng do code hiệu chỉnh. |
| M_u | C_mu=M·dCm/dM; M_u=C_mu QSc/(I_yu0); vẫn bằng 0 do dữ liệu cho dCm/dM=0 | Nelson, tr.116, Table 3.3; tr.150, Table 4.2. |
| Cybeta | Cybeta=−K_beta S_v/S | Nelson, tr.121, Table 3.4; K_beta suy từ dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cnbeta | Cnbeta=Cnbeta_wf−Cybeta l_v/b | Nelson, tr.76, Eq.(2.79); tr.121, Table 3.4; phần wing-body được code hiệu chỉnh. |
| a_v, eta_v, sidewash | K_beta=a_v·eta_v·(1+sigma_beta); K_r=a_v·eta_v (hiệu dụng) | Suy đại số trong code từ Nelson, tr.121, Table 3.4 và NASA CR-96008, tr.102, Table X-C. |
| Clr và Cyr/Cnr | Cyr=2K_r(S_v/S)(l_v/b); Cnr=−(l_v/b)Cyr; Clr=CL0/4+(z_eff/b)Cyr | Nelson, tr.120–121, §3.6.5 và Table 3.4; z_eff suy trong code. |
| Cyp | Cyp=CL0(AR+cosΛ)tanΛ/(AR+4cosΛ) | Nelson, tr.121, Table 3.4. |
| Clp, Cnp | Clp=Gain_roll·Clp_geometry; Cnp=Gain_Cnp·(−CL0/8) | Nelson, tr.121, Table 3.4; gain suy từ NASA CR-96008, tr.102, Table X-C. |
| Yp, Yr trong A | Yp=Yr=0 khi neglect_Yp_Yr=True; ngược lại dùng giá trị đầy đủ | Nelson, tr.199, Example 5.3. |
| Điểm trung hòa | Cm_alpha(h_NP)=0; SM=h_NP−h | Suy nghiệm Eq.(2.35), Nelson tr.56, với l_t thay đổi theo CG; không dùng trực tiếp Eq.(2.36). |

## 2. Hiệu chỉnh và khí động dọc

| Đại lượng | Công thức | Nguồn tham chiếu |
|---|---|---|
| eta_eff | −C_mq* / [2 a_tr V_Hr l_tr/c_r] | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.113, Eq.(3.79). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| ε_eff,r | C_madot* / C_mq* | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.113/115, Eq.(3.79)/(3.92). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| Gain downwash | ε_eff,r / [2a_wr/(πAR_r)] | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.58, Example 2.2 (downwash, không đánh số). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| Cmf_eff,r | C_malpha* − a_wr(h_r−h_ac,r) + eta_eff V_Hr a_tr(1−ε_eff,r) | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.56, Eq.(2.35). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| Gain Multhopp | Cmf_eff,r / Cmf_geometry,r | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.53, Eq.(2.32); tr.61, Fig.2.17. Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| Phần dư CLalpha | CLalpha* − [a_wr+eta_eff(S_tr/S_r)a_tr(1−ε_eff,r)] | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.116, Table 3.3 (tổng lực nâng wing–tail). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| e_eff | 2 CL0* CLalpha* / (π AR_r CDalpha*) | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.116, Table 3.3 (đạo hàm lực cản). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| CD_profile | CD0* − CL0*²/(π e_eff AR_r) | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.116, Table 3.3 (quan hệ polar). Dữ liệu: NASA CR-96008, tr.101, Table X-B. |
| Cmf_geometry | $\frac{180/\pi}{36.5S\bar c}\sum_i\Delta x_i w_i^2(\partial\epsilon_u/\partial\alpha)_i$ | Nelson, tr.53, Eq.(2.32); tr.61, Fig.2.17. |
| C_malpha | $a_w(h-h_{ac})+C_{mf}-\eta_{eff}V_Ha_t(1-\epsilon_{eff})$ | Nelson, tr.56, Eq.(2.35); hệ số hiệu dụng theo code. |
| T | $\eta_{eff}(S_t/S)a_t(1-\epsilon_{eff})$ | Đặt biến trong code từ Nelson, tr.56, Eq.(2.35). |
| h_NP; SM | $h_{NP}=[a_wh_{ac}+Tx_{tail,ac}/c-C_{mf}]/(a_w+T);\quad SM=h_{NP}-h$ | Suy nghiệm Eq.(2.35), Nelson tr.56, với l_t thay đổi theo CG. |
| C_Lu | $M^2C_{L0}/(1-M^2)$; dùng 0 khi incompressible=True | Nelson, tr.116, Table 3.3; giản lược theo Example 4.3, tr.158. |
| C_Du; C_mu | $M\partial C_D/\partial M;\quad M\partial C_m/\partial M$ | Nelson, tr.116, Table 3.3. |

## 3. Hiệu chỉnh và khí động ngang-hướng

| Đại lượng | Công thức | Nguồn tham chiếu |
|---|---|---|
| K_beta | −Cybeta* S_r/S_vr; tương đương a_v·eta_v·(1+sigma_beta) theo quy ước diện tích đang dùng | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| K_r | −Cnr* / [2 V_vr l_vr/b_r]; tương đương a_v·eta_v hiệu dụng | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cybeta | −K_beta S_v/S | Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cnbeta_wf,r | Cnbeta* + Cybeta* l_vr/b_r | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cnbeta | Cnbeta_wf,r S_r b_r/(S b) − Cybeta l_v/b | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cyr | 2K_r(S_v/S)(l_v/b) | Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cnr | −2K_r V_v(l_v/b) | Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| z_eff | (Clr*−CL0*/4)b_r/Cyr_r | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Clr | CL0/4 + (z_eff/b)Cyr | Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Gain dihedral | [Clbeta*−(z_eff/b_r)Cybeta*]/Γ_r | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Clbeta | Gain_dihedral·Γ+(z_eff/b)Cybeta | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Clp_geometry | −a_w(1+3λ_t)/[12(1+λ_t)] | Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Gain roll damping | Clp*/Clp_geometry,r | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Clp | Gain_roll·Clp_geometry | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Gain Cnp | Cnp*/(−CL0*/8) | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cnp | Gain_Cnp·(−CL0/8) | Suy đại số/hiệu chỉnh trong code từ Nelson, tr.121, Table 3.4; dữ liệu NASA CR-96008, tr.102, Table X-C. |
| Cy_p | $C_{L0}\frac{AR+\cos\Lambda}{AR+4\cos\Lambda}\tan\Lambda$ | Nelson, tr.121, Table 3.4. |

## 4. Đạo hàm có thứ nguyên và ma trận dọc

| Đại lượng | Công thức | Nguồn tham chiếu |
|---|---|---|
| X_u; X_w | $X_u=-(C_{Du}+2C_{D0})QS/(mu_0)$; $X_w=-(C_{D\alpha}-C_{L0})QS/(mu_0)$ | Nelson, tr.150, Table 4.2. |
| Z_u; Z_w | $Z_u=-(C_{Lu}+2C_{L0})QS/(mu_0)$; $Z_w=-(C_{L\alpha}+C_{D0})QS/(mu_0)$ | Nelson, tr.150, Table 4.2. |
| M_u; M_w | $M_u=C_{mu}QS\bar c/(I_yu_0)$; $M_w=C_{m\alpha}QS\bar c/(I_yu_0)$ | Nelson, tr.150, Table 4.2. |
| M_w_dot; M_q | $M_{\dot w}=C_{m\dot\alpha}QS\bar c^2/(2I_yu_0^2)$; $M_q=C_{mq}QS\bar c^2/(2I_yu_0)$ | Nelson, tr.150, Table 4.2. |
| A_lon, hàng 1–2 | $[X_u,X_w,0,-g\cos\theta_0]$; $[Z_u,Z_w,u_0,-g\sin\theta_0]$ | Nelson, tr.149, Eq.(4.51). |
| A_lon, hàng 3–4 | $[M_u+M_{\dot w}Z_u,M_w+M_{\dot w}Z_w,M_q+M_{\dot w}u_0,-M_{\dot w}g\sin\theta_0]$; $[0,0,1,0]$ | Nelson, tr.149, Eq.(4.51). |

## 5. Đạo hàm có thứ nguyên và ma trận ngang-hướng

| Đại lượng | Công thức | Nguồn tham chiếu |
|---|---|---|
| Y_beta; Y_p; Y_r | $Y_\beta=QSC_{Y\beta}/m$; $Y_p=QSbC_{Yp}/(2mu_0)$; $Y_r=QSbC_{Yr}/(2mu_0)$ | Nelson, tr.199, Table 5.1; đổi v=u_0β theo tr.195, Eq.(5.35). |
| L_beta; N_beta | $L_\beta=QSbC_{l\beta}/I_x$; $N_\beta=QSbC_{n\beta}/I_z$ | Nelson, tr.199, Table 5.1. |
| L_p; L_r | $L_p=QSb^2C_{lp}/(2I_xu_0)$; $L_r=QSb^2C_{lr}/(2I_xu_0)$ | Nelson, tr.199, Table 5.1. |
| N_p; N_r | $N_p=QSb^2C_{np}/(2I_zu_0)$; $N_r=QSb^2C_{nr}/(2I_zu_0)$ | Nelson, tr.199, Table 5.1. |
| A_lat, hàng 1 | $[Y_\beta/u_0,Y_p/u_0,Y_r/u_0-1,g\cos\theta_0/u_0]$ | Nelson, tr.195, Eq.(5.35); tr.199, Example 5.3: mặc định bỏ Y_p,Y_r. |
| A_lat, hàng 2–4 | $[L_\beta,L_p,L_r,0]$; $[N_\beta,N_p,N_r,0]$; $[0,1,\tan\theta_0,0]$ | Nelson, tr.195, Eq.(5.35), với I_xz=0; tr.199, Example 5.3. |

