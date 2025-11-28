# ERS - Extreme Response Style Correction Methods in PISA Data

[**Türkçe**](#türkçe) | [**English**](#english)

---

## Türkçe

### 📊 Proje Başlığı

**Aşırı Yanıt Tarzı (ERS) Düzeltme Yöntemlerinin PISA Verileri Üzerindeki Etkisi: Simülasyon ve Uygulama**

### 🎯 Amaç

Bu araştırma projesi, PISA gibi büyük ölçekli uluslararası değerlendirmelerde sıklıkla karşılaşılan Aşırı Yanıt Tarzı'nın (Extreme Response Style - ERS) veri kalitesi ve ülke sıralamaları üzerindeki etkisini incelemektedir.

Çalışmanın temel amacı, daha basit ve yaygın kullanılan ERS düzeltme yöntemlerinin (Z-skor standardizasyonu, ERS indeksi ile kovaryans kontrolü) etkinliğini, karmaşık IRT tabanlı yöntemlerle (ML-MNRM) karşılaştırarak, araştırmacılara pratik ve erişilebilir bir metodolojik kılavuz sunmaktır.

### 📋 Metodoloji

Çalışma iki ana aşamadan oluşmaktadır:

#### 1. Simülasyon Çalışması

**Amaç:** Farklı ERS düzeylerinin ve düzeltme yöntemlerinin ülke ortalamaları üzerindeki yanlılık (bias) ve hata (RMSE) düzeylerini karşılaştırmak.

**Koşullar:**
- ERS Düzeyleri: Düşük, Orta, Yüksek
- Gruplar: İki ülke/grup (ERS farklı vs. ERS eşit)
- Tekrarlama: 1000 iterasyon

**Değerlendirme:**
- Yanlılık (Bias)
- Kök Ortalama Kare Hata (RMSE)

#### 2. PISA Gerçek Veri Analizi

**Veri:** PISA 2018/2022 öğrenci anketi (sosyo-duygusal/motivasyonel ölçekler)

**Düzeltme Yöntemleri:**
- **Model 0:** Ham puanlar (kontrol grubu)
- **Model 1:** Z-skor standardizasyonu
- **Model 2:** Kovaryans kontrolü (Greenleaf ERS İndeksi)
- **Model 3:** ML-MNRM (Multilevel Multidimensional Nominal Response Model)

**Karşılaştırma:**
- Ülke ortalamalarındaki değişim
- Ülke sıralamalarındaki değişim
- Model performans karşılaştırması

### 🔑 Ana Bulgular (Beklenen)

1. Basit düzeltme yöntemlerinin (Z-skor) ülke sıralamaları üzerindeki etkisi
2. Karmaşık IRT modellerine göre basit yöntemlerin performansı
3. Farklı ERS düzeylerinin ülke karşılaştırmalarına etkisi
4. Pratik öneriler: Hangi yöntem ne zaman kullanılmalı?

### 📁 Proje Yapısı

```
ERS/
├── data/               # Veri dosyaları
│   ├── raw/            # Ham PISA verileri
│   ├── processed/      # İşlenmiş veriler
│   └── simulated/      # Simülasyon verileri
├── R/                  # R kaynak kodları
│   ├── simulation/     # Simülasyon çalışması
│   ├── analysis/       # PISA veri analizi
│   ├── functions/      # Yardımcı fonksiyonlar
│   └── visualization/  # Görselleştirme
├── output/             # Çıktılar
│   ├── figures/        # Grafikler
│   ├── tables/         # Tablolar
│   └── reports/        # Raporlar
└── docs/               # Dokümantasyon
```

### 🛠️ Kurulum

```r
# Gerekli paketleri yükleyin
install.packages(c("tidyverse", "mirt", "TAM", "lavaan", "psych", "here"))

# Proje dizinine gidin
setwd("path/to/ERS")

# (Opsiyonel) renv ile paket yönetimi
renv::restore()
```

### 🚀 Kullanım

```r
# Simülasyon çalışmasını çalıştır
source("R/simulation/01_data_generation.R")
source("R/simulation/02_correction_methods.R")
source("R/simulation/03_evaluation.R")

# PISA analizi çalıştır
source("R/analysis/01_data_preparation.R")
# ... diğer analiz scriptleri
```

### 📚 Kaynaklar

1. Ulitzsch, E., Lüdtke, O., & Robitzsch, A. (2023). The Role of Response Style Adjustments in Cross‐Country Comparisons—A Case Study Using Data from the PISA 2015 Questionnaire. *Educational Measurement: Issues and Practice, 42*(4), 101-112.

2. Lu, Y., & Bolt, D. M. (2015). Examining the attitude-achievement paradox in PISA using a multilevel multidimensional IRT model for extreme response style. *Large-scale Assessments in Education, 3*(1), 1-19.

3. Schoenmakers, M., Tijmstra, J., Vermunt, J., & Bolsinova, M. (2023). Correcting for Extreme Response Style: Model Choice Matters. *Educational and Psychological Measurement, 84*(1), 145-170.

4. Greenleaf, E. A. (1992). Improving rating scale measures by controlling for extreme response style. *Journal of Marketing Research, 29*(2), 176-185.

### 👥 İletişim

[İletişim bilgileri eklenecek]

### 📄 Lisans

[Lisans bilgisi eklenecek]

---

## English

### 📊 Project Title

**The Effects of Extreme Response Style (ERS) Correction Methods on PISA Data: Simulation and Application**

### 🎯 Purpose

This research project examines the impact of Extreme Response Style (ERS), frequently encountered in large-scale international assessments like PISA, on data quality and country rankings.

The main objective is to compare the effectiveness of simpler and more commonly used ERS correction methods (Z-score standardization, covariate control with ERS index) with complex IRT-based methods (ML-MNRM), providing researchers with a practical and accessible methodological guide.

### 📋 Methodology

The study consists of two main phases:

#### 1. Simulation Study

**Objective:** Compare bias and RMSE levels of different ERS levels and correction methods on country means.

**Conditions:**
- ERS Levels: Low, Medium, High
- Groups: Two countries/groups (different ERS vs. equal ERS)
- Replications: 1000 iterations

**Evaluation:**
- Bias
- Root Mean Square Error (RMSE)

#### 2. PISA Real Data Analysis

**Data:** PISA 2018/2022 student questionnaire (socio-emotional/motivational scales)

**Correction Methods:**
- **Model 0:** Raw scores (control group)
- **Model 1:** Z-score standardization
- **Model 2:** Covariate control (Greenleaf ERS Index)
- **Model 3:** ML-MNRM (Multilevel Multidimensional Nominal Response Model)

**Comparison:**
- Changes in country means
- Changes in country rankings
- Model performance comparison

### 🔑 Main Findings (Expected)

1. Effects of simple correction methods (Z-score) on country rankings
2. Performance of simple methods compared to complex IRT models
3. Impact of different ERS levels on country comparisons
4. Practical recommendations: Which method should be used when?

### 📁 Project Structure

```
ERS/
├── data/               # Data files
│   ├── raw/            # Raw PISA data
│   ├── processed/      # Processed data
│   └── simulated/      # Simulation data
├── R/                  # R source code
│   ├── simulation/     # Simulation study
│   ├── analysis/       # PISA data analysis
│   ├── functions/      # Utility functions
│   └── visualization/  # Visualization
├── output/             # Outputs
│   ├── figures/        # Figures
│   ├── tables/         # Tables
│   └── reports/        # Reports
└── docs/               # Documentation
```

### 🛠️ Installation

```r
# Install required packages
install.packages(c("tidyverse", "mirt", "TAM", "lavaan", "psych", "here"))

# Navigate to project directory
setwd("path/to/ERS")

# (Optional) Package management with renv
renv::restore()
```

### 🚀 Usage

```r
# Run simulation study
source("R/simulation/01_data_generation.R")
source("R/simulation/02_correction_methods.R")
source("R/simulation/03_evaluation.R")

# Run PISA analysis
source("R/analysis/01_data_preparation.R")
# ... other analysis scripts
```

### 📚 References

1. Ulitzsch, E., Lüdtke, O., & Robitzsch, A. (2023). The Role of Response Style Adjustments in Cross‐Country Comparisons—A Case Study Using Data from the PISA 2015 Questionnaire. *Educational Measurement: Issues and Practice, 42*(4), 101-112.

2. Lu, Y., & Bolt, D. M. (2015). Examining the attitude-achievement paradox in PISA using a multilevel multidimensional IRT model for extreme response style. *Large-scale Assessments in Education, 3*(1), 1-19.

3. Schoenmakers, M., Tijmstra, J., Vermunt, J., & Bolsinova, M. (2023). Correcting for Extreme Response Style: Model Choice Matters. *Educational and Psychological Measurement, 84*(1), 145-170.

4. Greenleaf, E. A. (1992). Improving rating scale measures by controlling for extreme response style. *Journal of Marketing Research, 29*(2), 176-185.

### 👥 Contact

[Contact information to be added]

### 📄 License

[License information to be added]

---

**Note:** For AI assistants working on this project, please refer to `CLAUDE.md` for detailed guidelines on code structure, conventions, and workflows.
