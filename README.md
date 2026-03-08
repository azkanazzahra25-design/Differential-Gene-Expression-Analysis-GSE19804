# Differential Gene Expression Analysis of Lung Cancer Dataset (GSE19804)

## Project Overview

Project ini bertujuan untuk mengidentifikasi gen yang mengalami perubahan ekspresi signifikan antara jaringan paru-paru normal dan kanker paru-paru menggunakan analisis transcriptomics.

Dataset yang digunakan berasal dari Gene Expression Omnibus (GEO) dengan ID **GSE19804**.

## Repository Structure

```
data/        : Dataset yang digunakan dalam analisis
scripts/     : Script R untuk analisis differential expression
outputs/     : Hasil analisis (visualisasi dan tabel DEG)
README.md    : Dokumentasi proyek
laporan_singkat.md : Laporan singkat analisis
```

## Methods

Analisis dilakukan menggunakan bahasa pemrograman **R** dengan paket:

* GEOquery
* limma
* EnhancedVolcano
* pheatmap

Langkah analisis:

1. Mengunduh dataset dari GEO
2. Normalisasi dan preprocessing data
3. Differential expression analysis menggunakan limma
4. Visualisasi hasil menggunakan volcano plot
5. Visualisasi pola ekspresi menggunakan heatmap

## Results

Analisis menghasilkan sejumlah gen yang menunjukkan perubahan ekspresi signifikan antara kondisi normal dan kanker paru-paru.

Visualisasi hasil meliputi:

* Volcano plot untuk melihat distribusi gen signifikan
* Heatmap untuk melihat pola ekspresi gen utama

Hasil lengkap analisis tersedia pada file:
**outputs/DEG_results_GSE19804.csv**

## Azka Nabilah Azzahra

Project ini dibuat sebagai bagian dari tugas capstone program bioinformatics.
