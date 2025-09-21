#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw
cd data/raw

echo "[*] Downloading GSE147821 raw tarball from GEO..."
wget -O GSE147821_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147821&format=file"

echo "[*] Extracting..."
tar -xvf GSE147821_RAW.tar

echo "[*] Done. Raw .h5 matrices are in data/raw/"