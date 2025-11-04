# 🧬 DNA Sequence Analyzer  
*A desktop GUI application for motif search, ORF detection, sequence comparison, and alignment.*

**Author:** Artemida Chadzinikolau  
**Version:** 1.2  

DNA Sequence Analyzer is a Tkinter-based GUI tool built for educational and research purposes.  
It lets you load DNA sequences (from FASTA, TXT, or NCBI), analyze motifs, detect open reading frames (ORFs), perform pairwise alignment, visualize motif distribution, generate heatmaps, and export results (CSV, PDF, FASTA).

---

## ✨ Features

### ✅ Sequence Loading
- Load sequences from **FASTA** or **plain text** files  
- Fetch sequences directly from **NCBI (Entrez API)**  
- Automatic computation of:
  - Sequence length  
  - GC content  
  - Sequence statistics table  

### ✅ Motif Analysis
- Search motifs (e.g., `ATG, CGCG, CATA`)  
- Adjustable segmentation for motif density  
- Highlight motifs inside the sequence text  
- Motif distribution plot (bar + scatter visualization)  

### ✅ Heatmap Visualization
- Compare motif frequencies across multiple sequences  
- Automatic heatmap sizing based on motif count  
- Interactive window (zoom/pan)  

### ✅ ORF Detection
- Detect ORFs on both strands and all three frames  
- Minimum AA length filter (default: 10 aa)  
- Show AA sequences and genomic coordinates  
- Export ORFs to **FASTA**  

### ✅ Pairwise Alignment
- Global alignment using **Biopython PairwiseAligner**  
- Automatic warnings for long sequences  
- Clean alignment display  

### ✅ Export Options
- **CSV** (motif statistics)  
- **FASTA** (detected ORFs)  
- **PDF report** including:
  - Sequence stats  
  - Motif table  
  - Motif distribution plot  
  - Motif heatmap (if multiple sequences loaded)

---

## Installation

```
git clone https://github.com/yourusername/dna-sequence-analyzer.git
cd dna-sequence-analyzer
```

## Requirements

```
biopython
pandas
matplotlib
seaborn
fpdf2
```

## Run the App

```
python3 dna_sequence_analyzer.py
```
