# MNRM vs IRTree: Teknik Karşılaştırma Kılavuzu

## 📖 Genel Bakış

Bu döküman, ERS araştırma projesinde MNRM (Multidimensional Nominal Response Model) ve IRTree (Item Response Tree) yöntemlerinin nasıl doğru bir şekilde karşılaştırıldığını açıklar.

**Referanslar:**
- Falk & Cai (2016) - MNRM modeli
- Böckenholt (2012) - IRTree modeli
- Plieninger (2017) - IRTree implementasyonu
- Schoenmakers et al. (2023) - MNRM ve IRTree karşılaştırması

---

## 🎯 Temel Farklar

### MNRM (Model 3)

**Yaklaşım:** Çok boyutlu nominal yanıt modeli

**Parametrizasyon:**
- İki boyut: θ (içerik), γ (ERS)
- s-matrisi ile kategori puanlaması:
  ```
  s_content: [0, 1, 2, 3]  (doğrusal puanlama)
  s_ERS:     [1, 0, 0, 1]  (uç kategoriler = 1)
  ```
- Alpha = 1.5 (her iki boyut için ayırt edicilik)
- m değerleri: [-0.5, 0.5] arası madde zorlukları

**Kategori Olasılıkları (Softmax):**
```
P(kategori k) = exp(θ'α_k + d_k) / Σ exp(θ'α_j + d_j)
```

**Avantajlar:**
- Tüm kategorileri eşzamanlı modelleme
- Kategori ayrımcılığı doğrudan kalibre edilir
- Hesaplama açısından daha verimli

**Dezavantajlar:**
- s-matrisinin spesifikasyonu subjektif olabilir
- Kategori yapısı hakkında güçlü varsayımlar


### IRTree (Model 4)

**Yaklaşım:** Hiyerarşik karar ağacı modeli

**Parametrizasyon:**
- Aynı iki boyut: θ (içerik), γ (ERS)
- Üç karar düğümü (4-nokta Likert için):
  ```
  Node 1: Katılma/Katılmama kararı (θ ile ilişkili)
  Node 2: Katılmama → Uç/Orta seçimi (γ ile ilişkili)
  Node 3: Katılma → Orta/Uç seçimi (γ ile ilişkili)
  ```

**Pseudo-Item Dönüşümü (4-nokta):**
```
Orijinal → (Node1, Node2, Node3)
    1    → (0,     0,     NA)   [Kesinlikle Katılmıyorum]
    2    → (0,     1,     NA)   [Katılmıyorum]
    3    → (1,     NA,    0)    [Katılıyorum]
    4    → (1,     NA,    1)    [Kesinlikle Katılıyorum]
```

**Düğüm Olasılıkları:**
```
Node 1: P(katılma) = logit^-1(a₁θ + d₁)
Node 2: P(uç | katılmama) = logit^-1(a₁θ - a₂γ + d₂)
Node 3: P(uç | katılma) = logit^-1(a₁θ + a₂γ + d₃)
```

**Avantajlar:**
- Yanıt sürecini psikolojik açıdan daha doğru modeller
- Farklı aşamalarda farklı parametreler
- Teorik olarak daha yorumlanabilir

**Dezavantajlar:**
- Hesaplama açısından daha karmaşık
- Daha fazla parametre → daha fazla belirsizlik
- Küçük örneklemlerde yakınsama sorunları

---

## 💻 İmplementasyon Detayları

### MNRM İmplementasyonu

**Dosya:** `R/simulation/00_mnrm_falk_generator.R`

```r
# s-matrix oluşturma
s_content <- matrix(0:(categories - 1), nrow = 1)
s_ers <- matrix(rep(0, categories), nrow = 1)
s_ers[1, 1] <- 1              # Birinci kategori = uç
s_ers[1, categories] <- 1     # Son kategori = uç
s_matrix <- rbind(s_content, s_ers)

# Alpha-ağırlıklı s-matrix
alpha_matrix <- diag(alpha) %*% s_matrix  # (2×2) %*% (2×K)

# Kategori olasılıkları (her kişi için)
linear_pred <- theta_matrix[person, ] %*% alpha_matrix + intercept
prob_matrix[person, ] <- exp(linear_pred) / sum(exp(linear_pred))
```

**Önemli Noktalar:**
- `diag(alpha)` kullanımı kritik (vektörle çarpım değil!)
- Intercept hesaplaması kategoriye göre kümülatif
- 4-nokta ve 5-nokta Likert otomatik olarak desteklenir


### IRTree İmplementasyonu

**Dosya:** `R/functions/irtree_pseudo_item_converter.R`

```r
# Pseudo-item dönüşümü (4-nokta)
convert_to_irtree_pseudoitems_4point <- function(data, n_items) {

  for (i in 1:n_items) {
    orig_resp <- data[, i]

    # Node 1: Katılma (1) vs Katılmama (0)
    node1 <- ifelse(orig_resp <= 2, 0, 1)

    # Node 2: Katılmama durumunda, uç (0) vs orta (1)
    node2 <- rep(NA, n_persons)
    node2[orig_resp == 1] <- 0  # Kesinlikle katılmıyorum = uç
    node2[orig_resp == 2] <- 1  # Katılmıyorum = orta

    # Node 3: Katılma durumunda, orta (0) vs uç (1)
    node3 <- rep(NA, n_persons)
    node3[orig_resp == 3] <- 0  # Katılıyorum = orta
    node3[orig_resp == 4] <- 1  # Kesinlikle katılıyorum = uç

    tree_data[, (i-1)*3 + 1] <- node1
    tree_data[, (i-1)*3 + 2] <- node2
    tree_data[, (i-1)*3 + 3] <- node3
  }
}
```

**Custom Node Fonksiyonları:**

**Node 2** (Katılmama → Uç/Orta):
```r
P.Node2TreePL <- function(par, Theta, ncat) {
  a1 <- par[1]  # İçerik ayırt ediciliği
  a2 <- par[2]  # ERS ayırt ediciliği
  d <- par[3]   # Eşik

  # NEGATİF loading on ERS (yüksek ERS → uç katılmama)
  P_extreme <- exp((a1*Theta[,1] - a2*Theta[,2] + d)) /
               (1 + exp((a1*Theta[,1] - a2*Theta[,2] + d)))

  cbind(1 - P_extreme, P_extreme)
}
```

**Node 3** (Katılma → Orta/Uç):
```r
P.Node3TreePL <- function(par, Theta, ncat) {
  a1 <- par[1]
  a2 <- par[2]
  d <- par[3]

  # POZİTİF loading on ERS (yüksek ERS → uç katılma)
  P_extreme <- exp((a1*Theta[,1] + a2*Theta[,2] + d)) /
               (1 + exp((a1*Theta[,1] + a2*Theta[,2] + d)))

  cbind(1 - P_extreme, P_extreme)
}
```

**Model Fitting:**
```r
fit_proper_irtree_model <- function(data, group_var, n_items, n_categories) {

  # 1. Pseudo-item dönüşümü
  tree_data <- convert_to_irtree_pseudoitems_4point(data, n_items)

  # 2. Model spesifikasyonu
  model_syntax <- paste0(
    "Content = ", paste(seq(1, n_items*3, by=3), collapse=","), "\n",
    "ERS = ", paste(c(seq(2, n_items*3, by=3),
                     seq(3, n_items*3, by=3)), collapse=","), "\n",
    "FREE = (GROUP, COV_21)"
  )

  # 3. Item tipleri (her madde için 3 düğüm)
  item_types <- rep(c("2PL", "Node2TreePL", "Node3TreePL"), n_items)

  # 4. Model fit
  irtree_fit <- multipleGroup(
    data = tree_data,
    model = model_syntax,
    group = group_var,
    itemtype = item_types,
    customItems = list(
      Node2TreePL = irtree_nodes$Node2,
      Node3TreePL = irtree_nodes$Node3
    ),
    method = "EM",
    invariance = c("free_mean", "free_var", "slopes", "intercepts"),
    technical = list(NCYCLES = 2000)
  )
}
```

---

## 🔬 Karşılaştırma Kriterleri

### 1. Bias (Yanlılık)

**Tanım:** Tahminlerin gerçek değerlerden sistematik sapması

```r
bias = mean(θ_estimated - θ_true)
```

**Yorumlama:**
- Bias ≈ 0: İdeal (yansız tahmin)
- Bias > 0: Pozitif yanlılık (olduğundan yüksek tahmin)
- Bias < 0: Negatif yanlılık (olduğundan düşük tahmin)


### 2. Variance Bias

**Tanım:** Tahminlerin varyansının gerçek varyansdan farkı

```r
variance_bias = var(θ_estimated) - var(θ_true)
```

**Yorumlama:**
- Variance bias > 0: Varyansı olduğundan büyük tahmin ediyor
- Variance bias < 0: Varyansı küçümseyor (shrinkage)


### 3. RMSE (Root Mean Square Error)

**Tanım:** Tahmin hatalarının karekök ortalaması

```r
rmse = sqrt(mean((θ_estimated - θ_true)^2))
```

**Yorumlama:**
- Bias ve varyans hatalarını birleştirir
- Küçük RMSE = daha iyi tahmin
- RMSE² = Bias² + Variance


### 4. ERS Sınıflandırma Doğruluğu

**Tanım:** Bireyleri yüksek/düşük ERS olarak doğru sınıflandırma oranı

```r
# Eşik: Median ERS
high_ers_true <- true_ers > median(true_ers)
high_ers_pred <- estimated_ers > median(estimated_ers)

accuracy = mean(high_ers_true == high_ers_pred)
sensitivity = TP / (TP + FN)  # Yüksek ERS'lileri yakalama
specificity = TN / (TN + FP)  # Düşük ERS'lileri doğru belirleme
```

---

## 📊 Beklenen Sonuçlar

### Koşullara Göre Performans

**Yüksek ERS Farkı Koşulları (high_ers_diff):**
- **MNRM:** Daha iyi bias performansı beklenir
- **IRTree:** Daha iyi ERS sınıflandırması beklenir
- **Sebep:** IRTree, ERS'yi daha detaylı modeller

**Küçük Örneklemler (N=100):**
- **MNRM:** Daha stabil yakınsama
- **IRTree:** Yakınsama sorunları olabilir
- **Sebep:** IRTree daha fazla parametre tahmin eder

**Uzun Ölçekler (30 madde):**
- **Her iki model:** İyi performans
- **IRTree:** Hesaplama süresi 3× daha uzun
- **Sebep:** 30 madde × 3 düğüm = 90 pseudo-item

**5-nokta Likert:**
- **MNRM:** Sorunsuz genişleme
- **IRTree:** Orta kategori (3) ataması kritik
- **Not:** Bizim implementasyonda 3 = katılmama olarak kodlandı

---

## ⚠️ Olası Sorunlar ve Çözümler

### Sorun 1: IRTree NA Sonuçları

**Sebep:**
- Model yakınsamadı
- Küçük örneklem + karmaşık model
- Pseudo-item dönüşümünde hata

**Çözüm:**
```r
# Debug scripti çalıştır
source("R/simulation/05_debug_na_values.R")

# IRTree convergence kontrolü
if (!irtree_result$convergence) {
  cat("IRTree yakınsamadı, MNRM kullanılıyor\n")
}

# Alternatif: NCYCLES artır
technical = list(NCYCLES = 5000, BURNIN = 2000)
```

### Sorun 2: Z-score NA Sonuçları

**Sebep:**
- Bazı bireyler tüm maddelere aynı yanıtı verdi (SD = 0)

**Çözüm:**
```r
# correction_methods.R'de zaten düzeltildi:
person_sd[person_sd == 0 | is.na(person_sd)] <- 1
```
**Not:** SD=0 olduğunda z-skorlar 0 olur (doğru davranış)

### Sorun 3: MNRM Matrix Dimension Hatası

**Sebep:**
- `alpha %*% t(s_matrix)` yerine `diag(alpha) %*% s_matrix` kullanılmalı

**Çözüm:**
```r
# YANLIŞGeçersiz
alpha_matrix <- alpha %*% t(s_matrix)  # HATALI

# DOĞRU
alpha_matrix <- diag(alpha) %*% s_matrix  # ✓
```

---

## 📈 Değerlendirme Workflow'u

### Adım 1: Veri Üretimi
```r
source("R/simulation/01_data_generation.R")
# Çıktı: data/simulated/all_conditions_full_factorial.rds
```

### Adım 2: Düzeltme Yöntemlerini Uygula
```r
source("R/simulation/02_correction_methods.R")
# Çıktı: output/tables/simulation_country_means_aggregated.csv
```

### Adım 3: Değerlendirme
```r
source("R/simulation/03_evaluation.R")
# Çıktılar:
#   - output/tables/simulation_evaluation_metrics_full.csv
#   - output/tables/simulation_factorial_summary.csv
#   - output/figures/*.png
```

### Adım 4: ERS Sınıflandırma Doğruluğu
```r
source("R/simulation/04_ers_classification_accuracy.R")
# Çıktılar:
#   - output/tables/ers_classification_summary.csv
#   - Sensitivity, Specificity, Accuracy, F1-score
```

---

## 📚 Literatür Karşılaştırması

### Schoenmakers et al. (2023) - Bizim İmplementasyon

| Özellik | Schoenmakers et al. | Bizim Çalışma |
|---------|---------------------|---------------|
| Veri Üretimi | MNRM-Falk | ✓ Aynı (MNRM-Falk) |
| IRTree Kalibrasyon | MNRM'den türetme | ✓ 00_irtree_calibration.R |
| Alpha | 1.5 | ✓ 1.5 |
| m değerleri | Var | ✓ seq(-0.5, 0.5) |
| Koşullar | 8 (4 ERS × 2 eşik) | 72 (4 ERS × 3 N × 3 items × 2 Likert) |
| Replikasyon | 500 | 100 (daha hızlı) |
| Değerlendirme | Bias, Variance Bias | ✓ Aynı + Classification Accuracy |
| Pseudo-item | ✓ Doğru IRTree | ✓ Doğru IRTree |

**Sonuç:** Bizim implementasyon Schoenmakers et al. (2023) metodolojisini takip ediyor ve genişletiyor.

---

## 🎯 Ana Çıkarımlar

1. **MNRM ve IRTree farklı modeller:** Aynı θ'yı farklı şekillerde tahmin ederler

2. **IRTree daha spesifik:** ERS'yi karar ağacıyla modelleyerek daha detaylı bilgi verir

3. **MNRM daha robust:** Küçük örneklemlerde daha az yakınsama sorunu

4. **Koşula göre seçim:** En iyi yöntem araştırma hedefine ve veri yapısına bağlı

5. **İki yöntem birlikte:** MNRM + IRTree karşılaştırması en zengin bilgiyi verir

---

## ✅ Kontrol Listesi

Simülasyon sonuçlarınız geçerli mi? Kontrol edin:

- [ ] MNRM s-matrix doğru (s_ERS = [1,0,0,1])?
- [ ] IRTree pseudo-item dönüşümü doğru (NA'ler doğru yerlerde)?
- [ ] Alpha = 1.5 her iki boyut için?
- [ ] m değerleri [-0.5, 0.5] arasında?
- [ ] 4-nokta eşikler: [-1, 0, 1]?
- [ ] 5-nokta eşikler: [-1.5, -0.5, 0.5, 1.5]?
- [ ] Bias ≈ 0 ideal koşulda (equal_ers)?
- [ ] RMSE high_ers_diff > low_ers?
- [ ] IRTree convergence ≥ %80?

---

**Son Güncelleme:** 2025-11-30
**Versiyon:** 1.0.0
