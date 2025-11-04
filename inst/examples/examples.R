# =========================================================
# File: examples/examples.R
# Demonstrasi Lengkap Penggunaan Paket dimvR
# =========================================================

# ---------------------------------------------------------
# 0. Persiapan
# ---------------------------------------------------------
library(dimvR)

set.seed(123)

# =========================================================
# 1. Simulasi Dataset dengan Missing Values
# =========================================================

# Buat dataset sederhana dengan 3 variabel numerik
X <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100) * 2 + rnorm(100, sd = 0.2),  # agak berkorelasi dengan x1
  x3 = rnorm(100)
)

# Introduksi missing values (10% tiap kolom)
for (j in seq_along(X)) {
  X[sample.int(nrow(X), 10), j] <- NA
}

# Cek proporsi missing
colSums(is.na(X))
#> x1 x2 x3 
#> 10 10 10

# =========================================================
# 2. Melatih DIMV Imputer
# =========================================================
# Parameter:
# - lambda = 0.1 → kekuatan regularisasi ridge
# - feature_select = FALSE → tidak melakukan seleksi fitur
# =========================================================

imp <- dimv_train(X, lambda = 0.1)

# Tampilkan ringkasan metadata imputer
print(imp)
# Output:
# R-native DIMV imputer (ridge-conditional)
# lambda = 0.1 | adaptive = FALSE | feature_select = FALSE | iterations = ...

# =========================================================
# 3. Diagnostik Konvergensi (opsional)
# =========================================================
# Mengevaluasi stabilitas iterasi dan varians residual
# =========================================================

X_imp_once <- dimv_impute_new(imp, X)
diag <- dimv_diagnostics(imp, X, X_imp_once)
print(diag)

# Alternatif tampilan ringkas:
# print_dimv_diag(imp)

# =========================================================
# 4. Imputasi Deterministik (Single Imputation)
# =========================================================
# Menghasilkan dataset lengkap tanpa elemen acak tambahan
# =========================================================

X_imputed <- dimv_impute_new(imp, X)
summary(X_imputed)

# =========================================================
# 5. Multiple Imputation (Rubin-Ready)
# =========================================================
# m = 5 menghasilkan 5 versi dataset berbeda,
# dengan penambahan residual Gaussian untuk meniru ketidakpastian.
# =========================================================

imputations <- dimv_impute_multiple(imp, X, m = 5, seed = 2024)

# Jumlah dataset hasil imputasi
length(imputations)
#> [1] 5

# Lihat ringkasan dataset pertama
summary(imputations[[1]])

# =========================================================
# 6. Contoh Analisis Setelah Multiple Imputation
# =========================================================
# Misal kita melakukan regresi sederhana pada setiap dataset
# =========================================================

models <- lapply(imputations, function(df) {
  lm(x1 ~ x2 + x3, data = df)
})

# Lihat ringkasan model pertama
summary(models[[1]])

# =========================================================
# 7. Visualisasi Relevansi Fitur (opsional)
# =========================================================
# Dapat digunakan untuk melihat fitur mana yang paling
# berpengaruh dalam struktur kovarians data.
# =========================================================

scores <- feature_select_score(X)
print(scores)

plot_feature_selection(scores)

# =========================================================
# 8. Contoh Diagnostik Adaptif Lambda (opsional)
# =========================================================
# Memberikan rekomendasi tingkat regularisasi berdasarkan
# median korelasi antar fitur.
# =========================================================

recommended_lambda <- adaptive_lambda(X)
cat("Rekomendasi lambda adaptif:", recommended_lambda, "\n")

# =========================================================
# 9. Catatan Akhir
# =========================================================
# Fungsi-fungsi tambahan yang tersedia:
# - dimv_convergence_diag() : ringkasan konvergensi imputer
# - generate_report() : menghasilkan laporan HTML otomatis
# - run_full_pipeline() : eksperimen lengkap dengan SHAP
# =========================================================

# Contoh (tidak dijalankan di sini):
# results <- run_full_pipeline(X, y = X$x1)
# generate_report(results, dataset_name = "Simulated Data")
