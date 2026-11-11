# ECG and echocardiography

R code for the paper. / 心電図・心エコー研究の解析コード\
↓\
[Machine-Reported Electrocardiographic Right Axis Deviation or Right Ventricular Hypertrophy and Echocardiographic Pulmonary Hypertension Phenotypes and Right-Heart Abnormalities](https://pmc.ncbi.nlm.nih.gov/articles/PMC13358884/)

## Setup

In R

```r
install.packages(c("data.table", "dplyr", "stringr", "lubridate", "janitor", "sandwich", "lmtest", "binom", "R.utils"))
```

MIMIC input paths [ecg.R](ecg.R) · [echo.R](echo.R) · [clinical.R](clinical.R)

## Run

Terminal (folder with the R files) / Rファイルがあるフォルダで実行

```bash
Rscript 1_prepare.R && Rscript 2_analyze.R && Rscript 3_tables.R
```

CSV output `data_processed/`, `results/`, `tables/`

Without the manual review CSV, the script stops after S5. / 手動判定CSVがない場合はS5まで出力

[Supporting Information (S6, S7)](https://pmc.ncbi.nlm.nih.gov/articles/PMC13358884/#_adsm93_)
