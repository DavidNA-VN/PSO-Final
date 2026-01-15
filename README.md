## Cell-free ISAC Beamforming (MATLAB)

Dự án mô phỏng hệ **cell-free ISAC MIMO** (Integrated Sensing and Communication) với nhiều Access Point (AP) phối hợp vừa phục vụ liên lạc cho các UE vừa “sensing” mục tiêu.

Repo này chạy mô phỏng theo luồng trong [run.m](run.m) → [simulation.m](simulation.m), và xuất kết quả ra thư mục [output/](output).

### Nội dung chính trong code

Trong [simulation.m](simulation.m):

- Sinh hình học (vị trí AP/UE/Target) bằng [utils/generate_positions.m](utils/generate_positions.m)
- Sinh kênh LOS cho liên lạc bằng [channel/LOS_channel.m](channel/LOS_channel.m)
- Tạo beam cho sensing theo hướng mục tiêu (beamsteering) và tạo bản **nulling** để giảm nhiễu lên UE
- Quét theo tỉ lệ công suất cho liên lạc $\rho$ (biến `params.P_comm_ratio`)
- So sánh các phương án:
	- **NS+RZF**: Sensing Null-Space + Communication Regularized ZF
	- **JSC+PSO-LD**: Giữ beam liên lạc cố định (RZF), tối ưu sensing bằng **PSO low-dimension** (tối ưu trọng số phức theo từng AP) trong [optimization/opt_sens_PSO_lowdim.m](optimization/opt_sens_PSO_lowdim.m)

Các metric được tính trong [utils/compute_metrics.m](utils/compute_metrics.m):

- `min_SINR`: SINR nhỏ nhất trong các UE
- `SSNR`: sensing SNR của mục tiêu

### Yêu cầu

- MATLAB (khuyến nghị R2021a+)
- (Tuỳ chọn) Parallel Computing Toolbox: code dùng `parfor` trong [simulation.m](simulation.m)
- (Tuỳ chọn) CVX: repo có các hàm tối ưu (SOCP/SDP) trong [optimization/](optimization), nhưng **luồng mặc định trong [run.m](run.m) hiện không bắt buộc CVX**. Nếu bạn không có CVX, có thể comment phần addpath CVX trong [run.m](run.m).

### Cách chạy nhanh

1) Mở MATLAB tại thư mục gốc của repo.

2) (Tuỳ chọn) chỉnh tham số trong [sim_params.m](sim_params.m):

- `params.N_t`: số antenna mỗi AP
- `params.M_t`: số AP
- `params.U`: số UE
- `params.P_comm_ratio`: danh sách $\rho$ để sweep
- `params.repetitions`: số lần lặp Monte Carlo
- `params.geo.*`: thiết lập hình học

3) Chạy script [run.m](run.m).

Script sẽ:

- Chạy mô phỏng theo sweep $\rho$ và lưu [output/power_data.mat](output/power_data.mat)
- Chạy mô phỏng theo khoảng cách tối thiểu UE–Target và lưu nhiều file `output/dist_data*.mat`
- Gọi các script plot trong [plots/](plots)

### Kết quả (figures có sẵn)

- Power ratio vs performance: ![power-vs-perf](plots/power_vs_perf.png)
- Min distance vs performance: ![dist-vs-perf](plots/dist_vs_perf.png)

### Ảnh minh hoạ trong thư mục png/

Các ảnh trong [png/](png) (ảnh chụp/ghi lại kết quả khi chạy) đã được include để tham khảo:

- ![png-1](png/92c0aca2-847e-4e38-8dd2-f247ac3cb11f.png)
- ![png-2](png/dd40365a-b851-43f3-ba01-3ff6f795428c.png)
- ![png-3](png/Screenshot%202026-01-15%20225735.png)
- ![png-4](png/Screenshot%202026-01-15%20225744.png)
- ![png-5](png/Screenshot%202026-01-15%20230826.png)
- ![png-6](png/Screenshot%202026-01-15%20233215.png)
- ![png-7](png/Screenshot%202026-01-15%20233225.png)

### Ghi chú về dữ liệu đầu ra

- Dữ liệu được lưu dạng `results` trong các file `.mat` ở [output/](output)
- Các script vẽ hình đọc dữ liệu bằng [plots/load_results.m](plots/load_results.m)

### Tham khảo

Mã nguồn/ý tưởng mô phỏng bám theo bài: U. Demirhan và A. Alkhateeb, “Cell-free ISAC MIMO systems: Joint sensing and communication beamforming,” arXiv:2301.11328.

### License

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0): https://creativecommons.org/licenses/by-nc-sa/4.0/
