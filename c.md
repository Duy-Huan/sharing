# Module: Superluminescent Diode (`sld.m`)

## 1. Tổng quan
Module `sld` mô phỏng nguồn phát quang siêu phát xạ (Superluminescent Diode - SLD), đóng vai trò là nguồn quang băng thông rộng, độ kết hợp thấp (low coherence) cung cấp cho hệ thống quang học FOG. Module tích hợp mô hình dòng điện - công suất (P-I curve), hiệu ứng dịch chuyển bước sóng theo nhiệt độ và dòng điện hoạt động, mô phỏng gợn sóng phổ (spectral ripple), tiếng ồn cường độ tương đối (RIN), phân bố phổ Gaussian, phân chia công suất theo độ phân cực (DoP), và mô hình nhiễu pha Gauss-Markov (AR1) cho hai nhánh TE và TM.

## 2. Giải thích các tham số (`properties`)

### Nhóm tham số cấu hình quang - điện chính:
* `LambdaAxis`: Trục bước sóng suy giảm/phân bố phổ (mảng $1\times M$).
* `Centrewavelength`: Bước sóng trung tâm danh định của nguồn SLD (nm).
* `Bandwidth`: Băng thông quang học (FWHM) của nguồn (nm).
* `ThresholdCurrent`: Dòng điện ngưỡng phát quang của diode (mA).
* `SlopeEfficiency`: Hệ số hiệu suất dốc của đặc tuyến P-I (mW/mA).
* `TempCoeff`: Hệ số dịch chuyển bước sóng theo nhiệt độ chip (nm/$^\circ$C).
* `DoP`: Độ phân cực (Degree of Polarization) trong khoảng $[0; 1]$.
* `Ripple_dB`: Biên độ gợn sóng phổ (rms gain ripple) (dB).
* `RipplePeriod`: Chu kỳ gợn sóng phổ (nm).
* `RIN_dB`: Tiếng ồn cường độ tương đối (Relative Intensity Noise) (dB/Hz).
* `Fs`: Tần số lấy mẫu mô phỏng (Hz).
* `CurrentCoeff`: Hệ số dịch chuyển bước sóng theo dòng điện điều khiển (nm/mA).
* `NominalCurrent`: Dòng điện định mức của SLD (mA).

### Nhóm tham số nhiễu pha Gauss-Markov:
* `Alpha`: Hệ số tương quan AR(1) ($0 < \text{Alpha} < 1$, mặc định $0.985$).
* `SigmaNoise`: Biên độ kích thích nhiễu trắng (mặc định $0.2$).

### Nhóm tham số riêng tư (`private`):
* `PrevPhaseTE` / `PrevPhaseTM`: Lưu trữ trạng thái pha của chu kỳ trước cho nhánh TE và TM.

---

## 3. Các hiệu ứng và công thức chính

* **Đặc tuyến dòng - công suất (P-I Curve)**:
  Tổng công suất quang phát ra được tính dựa trên dòng điện ngõ vào so với dòng điện ngưỡng:
  $$\text{total\_power} = \begin{cases} \text{SlopeEfficiency} \cdot (\text{current} - \text{ThresholdCurrent}), & \text{if } \text{current} > \text{ThresholdCurrent} \\ 0, & \text{otherwise} \end{cases}$$

* **Dịch chuyển bước sóng theo dòng điện và nhiệt độ ($(\lambda_c)$)**:
  $$\Delta\lambda_{\text{current}} = \text{CurrentCoeff} \cdot (\text{current} - \text{NominalCurrent})$$
  $$\Delta\lambda_{\text{temp}} = \text{TempCoeff} \cdot (\text{temp} - 25)$$
  $$\lambda_c = \text{Centrewavelength} + \Delta\lambda_{\text{temp}} + \Delta\lambda_{\text{current}}$$

* **Mô hình gợn sóng phổ (Spectral Ripple)**:
  $$\text{RippleAmplitude} = \frac{10^{\text{Ripple\_dB} / 10} - 1}{10^{\text{Ripple\_dB} / 10} + 1}$$
  $$\text{ripple\_modulation} = 1 + \text{RippleAmplitude} \cdot \cos\left(\frac{2\pi}{\text{RipplePeriod}} \cdot \lambda_c\right)$$
  $$\text{base\_power} = \max(\text{total\_power} \cdot \text{ripple\_modulation}, 0)$$

* **Tiếng ồn cường độ tương đối (RIN)**:
  $$\text{rin\_linear} = 10^{\text{RIN\_dB} / 10}$$
  $$\text{RINlinear} = \sqrt{\text{rin\_linear} \cdot \frac{F_s}{2}}$$
  $$\text{effective\_power} = \text{base\_power} \cdot \left(1 + \text{RINlinear} \cdot \text{rand}\right)$$

* **Phân bố phổ Gaussian và Phân chia DoP**:
  Độ lệch chuẩn phổ $\sigma$ và phổ bao envelope:
  $$\sigma = \frac{\text{Bandwidth}}{2\sqrt{2\ln(2)}}$$
  $$\text{spectral\_envelope} = \exp\left(-\frac{(\Lambda_{\text{axis}} - \lambda_c)^2}{2\sigma^2}\right)$$
  Công suất phân bố trên hai trục phân cực TE và TM:
  $$\text{total\_power\_te} = \frac{1 + \text{DoP}}{2} \cdot \text{effective\_power}$$
  $$\text{total\_power\_tm} = \frac{1 - \text{DoP}}{2} \cdot \text{effective\_power}$$

* **Mô hình nhiễu pha Gauss-Markov (AR1 Process)**:
  $$\text{current\_phase} = \text{Alpha} \cdot \text{PrevPhase} + \sqrt{1 - \text{Alpha}^2} \cdot \sum \left(\text{SigmaNoise} \cdot \text{randn}\right)$$

* **Định dạng cấu trúc tín hiệu ngõ ra (`sig`)**:
  $$\text{sig.TE} = \sqrt{\max(\text{total\_power\_te}, 0)} \cdot e^{j \cdot \text{current\_phase\_te}}$$
  $$\text{sig.TM} = \sqrt{\max(\text{total\_power\_tm}, 0)} \cdot e^{j \cdot \text{current\_phase\_tm}}$$
  $$\text{sig.centrewavelength} = \frac{\sum (\Lambda_{\text{axis}} \cdot \text{spectral\_envelope})}{\sum \text{spectral\_envelope}}$$
