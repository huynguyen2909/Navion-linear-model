# Đáp ứng tự do 60 giây: mô hình Python và ma trận Nelson

Chạy `python navion_reports.py` để tạo lại hai hình PNG/SVG. Mỗi hình có một khung chứa bốn biến trạng thái và hai bộ đáp ứng.

- Python: nét liền, marker đặc. Nelson: nét đứt, marker rỗng.
- Màu xanh dương/vuông: Δu hoặc Δβ; xanh lá/tròn: Δw hoặc Δp; cam/tam giác: Δq hoặc Δr; tím/ngôi sao: Δθ hoặc Δφ.
- Khung 8 × 8 inch; tiêu đề 16 pt; nhãn trục 14 pt; số trục 12 pt; legend và ghi chú 11 pt.
- Mô phỏng đáp ứng tự do ẋ = Ax: 0–60 s, RK4, Δt = 0,01 s, 6001 mẫu; cùng x(0) = [0, 0, 0.1, 0]ᵀ. Kênh dọc có Δq(0) = 0,1 rad/s; kênh ngang–hướng có Δr(0) = 0,1 rad/s.
- Các đường dùng giá trị và đơn vị nguyên gốc, không chuẩn hóa. Dọc: [Δu, Δw, Δq, Δθ], đơn vị [ft/s, ft/s, rad/s, rad]. Ngang–hướng: [Δβ, Δp, Δr, Δφ], đơn vị [rad, rad/s, rad/s, rad]. Vì chung trục tung, các biến có biên độ nhỏ hơn có thể nằm sát trục 0.

## Nguồn Nelson

Robert C. Nelson, *Flight Stability and Automatic Control*, 2nd ed., 1998: ma trận chuyển động dọc ở trang in 158 (PDF trang 169); ma trận ngang–hướng ở trang in 199 (PDF trang 210).

`nelson_reference.py` lưu đúng các phần tử A đã làm tròn in trong sách. Đường Nelson được tính lại từ các ma trận đó với điều kiện đầu nêu trên, không phải đường số hóa từ hình trong sách. Không dựng đáp ứng chỉ từ các trị riêng do Nelson công bố; riêng trị riêng không xác định đầy đủ đáp ứng các biến trạng thái. Ma trận này cũng được dùng cho cột “Tính từ A in trong sách” trong các báo cáo trị riêng.

Đây là đối chiếu giữa hai mô hình tuyến tính ở cấu hình chuẩn; không phải xác thực độc lập bằng dữ liệu bay.
