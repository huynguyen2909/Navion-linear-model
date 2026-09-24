# Quỹ đạo trị riêng Navion: Vv, Vh và SM

Chạy `python navion_root_locus_sweep.py` ở thư mục gốc để tạo lại sáu đồ thị SVG/PNG và sáu bảng CSV trong thư mục này. Script dùng `numpy` và `matplotlib`, tính trực tiếp trị riêng của hai ma trận 4×4 từ `navion_linear_models.py`.

## Quy ước quét

| Tham số | Giá trị chuẩn | Bước 10 | Cách thực hiện |
|---|---:|---:|---|
| `Vv` | 0.0325436084 | 0.0488154127 | Tăng `S_v` và tính lại `V_v=S_v l_v/(S b)`; giữ hình học khác cố định. |
| `Vh` | 0.655987796 | 0.983981693 | Tăng `S_t` và tính lại `V_H=S_t l_t/(S c)`; giữ `AR_t` cố định. |
| `SM` | 0.146639558 | 0.219959337 | Tăng `SM=h_NP−h_CG` bằng cách dịch CG về phía trước; giữ vị trí đuôi trong hệ tọa độ máy bay. |

Giá trị bước `k=0,…,10` bằng giá trị chuẩn nhân `1+0.05k`. Ở mỗi bước chỉ thay một tham số, điều kiện bay và các tham số hiệu chỉnh tham chiếu giữ nguyên. Vòng tròn rỗng chỉ nghiệm ở bước 0; dấu × chỉ bước 10; ô phóng to hiển thị mode gần gốc.

## Kết quả cần chú ý

| Quét | Kênh dọc | Kênh ngang-hướng |
|---|---|---|
| `Vv` | Hai mode không đổi vì ma trận dọc không phụ thuộc `S_v`. | Dutch roll dịch chuyển, nghiệm spiral đổi từ −0.00914685 sang +0.01438639 s⁻¹; lần đầu dương ở bước 2 (+10%). |
| `Vh` | Short-period dịch chuyển rõ; phugoid thay đổi nhẹ. | Ba mode không đổi vì ma trận ngang-hướng không phụ thuộc `S_t`. |
| `SM` | Phugoid và short-period thay đổi. | Dutch roll và spiral thay đổi; roll gần như không đổi. |

Nghiệm spiral dương ở phép quét `Vv` là **kết quả của mô hình hiệu dụng hiện tại**. Các gain khí động được xác định tại cấu hình gốc; dải tăng đến +50% là khảo sát độ nhạy, chưa được xác thực bằng dữ liệu của những cấu hình đuôi mới. Hai đồ thị không đổi là đặc điểm của hai mô hình giảm bậc độc lập.
