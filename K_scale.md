# Phương pháp tính toán hệ số K_scale (Giải tích lý thuyết FOG)

Để tính toán hệ số chuyển đổi $K_{\text{scale}}$ (từ dòng điện ngõ ra $I_{\text{pd}}$ tính bằng $\text{mA}$ về tốc độ góc $\Omega$ tính bằng $\text{°/h}$) theo phương pháp giải tích lý thuyết, chúng ta thiết lập chuỗi phương trình liên kết từ đặc tính hình học cuộn sợi, thông số quang học của linh kiện (SLD, MIOC) cho đến đặc tính của Photodetector.

## Các bước xây dựng công thức lý thuyết

### 1. Bước 1: Từ tốc độ góc sang độ dịch pha Sagnac ($\Omega \rightarrow \Delta\phi_S$)
Tốc độ góc đầu vào $\Omega$ (đơn vị $\text{°/h}$) cần được đổi sang đơn vị tiêu chuẩn quốc tế $\text{rad/s}$:
$$\Omega_{\text{rad/s}} = \Omega_{\text{deg/h}} \times \left( \frac{\pi}{180} \right) \times \left( \frac{1}{3600} \right)$$

Độ dịch pha quang học Sagnac $\Delta\phi_S$ (rad) sinh ra do hiệu ứng quay của trái đất/bàn xoay:
$$\Delta\phi_S = \frac{2\pi L D}{\lambda c} \cdot \Omega_{\text{rad/s}} = \frac{2\pi L D}{\lambda c} \cdot \Omega_{\text{deg/h}} \cdot \left( \frac{\pi}{180 \times 3600} \right)$$

Trong đó:
- **$L$**: Chiều dài sợi quang trong cuộn ($\text{m}$).
- **$D$**: Đường kính trung bình của cuộn sợi quang ($\text{m}$).
- **$\lambda$**: Bước sóng trung tâm của nguồn quang SLD ($\text{m}$).
- **$c$**: Vận tốc ánh sáng trong chân không ($\approx 3 \times 10^8 \text{ m/s}$).

---

### 2. Bước 2: Từ pha Sagnac sang biến thiên công suất quang ($\Delta\phi_S \rightarrow \Delta P_{\text{opt}}$)
Trong hệ thống giao thoa kế vòng kín (hoặc tại điểm làm việc vuông góc - quadrature point của hệ thống), sự biến thiên công suất quang học thu được tại cổng giao thoa tỉ lệ tuyến tính với độ dịch pha:
$$\Delta P_{\text{opt}} = P_0 \cdot V \cdot \Delta\phi_S \quad (\text{mW})$$

Trong đó:
- **$P_0$**: Công suất quang đầu ra của nguồn SLD ($\text{mW}$).
- **$V$**: Hệ số tương phản giao thoa (Fringe Visibility, thường dao động trong khoảng $0.8 \div 0.98$).

---

### 3. Bước 3: Từ công suất quang sang dòng điện Photodetector ($I_{\text{pd}}$)
Thông qua hệ số đáp ứng quang điện (Responsivity $R$ của Photodiode):
$$I_{\text{pd}} = R \cdot \Delta P_{\text{opt}} \quad (\text{mA})$$

Trong đó:
- **$R$**: Độ nhạy của Photodetector (đơn vị $\text{mA/mW}$ hoặc $\text{A/W}$).

---

## Công thức tổng hợp tính $K_{\text{scale}}$

Thay thế các phương trình trên vào nhau, ta có mối quan hệ trực tiếp giữa dòng điện ngõ ra $I_{\text{pd}}$ ($\text{mA}$) và tốc độ góc $\Omega$ ($\text{°/h}$):
$$I_{\text{pd}} = \left[ \frac{2\pi^2 \cdot L \cdot D \cdot P_0 \cdot V \cdot R}{180 \times 3600 \cdot \lambda c} \right] \cdot \Omega_{\text{deg/h}}$$

Từ đó, suy ra hệ số lý thuyết $K_{\text{scale}}$ (để chuyển đổi ngược từ $\text{mA}$ về $\text{°/h}$):
$$K_{\text{scale}} = \frac{\Omega_{\text{deg/h}}}{I_{\text{pd}}} = \frac{180 \times 3600 \cdot \lambda c}{2\pi^2 \cdot L \cdot D \cdot P_0 \cdot V \cdot R} \quad \left( \frac{\text{°/h}}{\text{mA}} \right)$$
