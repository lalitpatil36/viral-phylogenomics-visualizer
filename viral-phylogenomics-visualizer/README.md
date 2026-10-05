# Viral Phylogenomics & Sequence Identity Visualizer

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python: 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Antigravity: Skill](https://img.shields.io/badge/Antigravity-Custom%20Skill-00E5FF.svg)](https://deepmind.google)
[![Format: Multi--Export](https://img.shields.io/badge/export-PNG%20%7C%20TIFF%20%7C%20PDF%20%7C%20EMF-success.svg)](#export-formats)

A publication-grade, automated phylogenomic analysis and identity visualization pipeline for **ANY** nucleotide or amino acid sequence dataset (viral, microbial, plant, or animal).

Designed for research manuscripts submitting to high-impact journals (*Nature Microbiology*, *Phytopathology*, *Virology*, *Journal of General Virology*).

---

## 📖 Complete Word Manual Included

A complete, camera-ready user manual is provided in Microsoft Word format:
- **File**: [`Viral_Phylogenomics_Visualizer_Manual.docx`](./Viral_Phylogenomics_Visualizer_Manual.docx)
- **Formatting Specifications**: Times New Roman, 12 pt, 1.5 line spacing, justified paragraph alignment.
- **Includes**: Comprehensive protocol, execution steps, schematic layouts, troubleshooting matrices, and citation standards.

---

## 🌟 Visual Portfolio Architecture

### 1. Pure Lower-Triangular Stepped SDT Identity Matrix
Replicates the authentic **Sequence Demarcation Tool (SDT v1.2, Muhire et al., 2014)** interface:
- **Clean Stepped Geometry**: Lower triangle active, upper triangle clean white.
- **Official Rainbow Gradient**: Deep Navy (low) $\rightarrow$ Cyan $\rightarrow$ Grass Green $\rightarrow$ Vivid Yellow $\rightarrow$ Crimson (100%).
- **Integrated Colorbar**: Sleek colorbar embedded directly into the upper-right corner.
- **Flush Demarcations**: Dashed diagonal clade separation lines with labels flush against staircase cell borders.

```
       Taxon A  Taxon B  Taxon C  Taxon D
Taxon A [ 100 ]
Taxon B [  98 ][ 100 ]
Taxon C [  92 ][  91 ][ 100 ]
Taxon D [  85 ][  84 ][  86 ][ 100 ]
```

---

### 2. Midpoint-Rooted Linear Maximum-Likelihood (ML) Phylogeny
- Inferred using IQ-TREE with ultrafast bootstrap approximation ($n=1,000$).
- Midpoint-rooted with ladderized branch topologies.
- Calibrated branch scale bar (substitutions/site).
- Visual isolate badges (`★ Isolate`) and reference indicators (`● Reference`).

---

### 3. Polar Circular Cladogram with Isolated Legends
- Polar circular projection ($15^\circ$ to $345^\circ$).
- Concentric metadata rings (Host, Clade/Genogroup, Provenance).
- **Dedicated Non-Overlapping Legend Container**: Anchored cleanly in the lower canvas quadrant ($0.05 \le y \le 0.17$), completely preventing collisions with tree branches or radial leaf labels.

---

### 4. Combined Dual-Panel Publication Figures
- Side-by-side composite linking linear evolutionary phylogeny (Panel a) with the stepped SDT identity matrix (Panel b).

---

## 🚀 Quick Start

### Installation

Clone the repository and install the dependencies:
```bash
git clone https://github.com/lalitpatil36/viral-phylogenomics-visualizer.git
cd viral-phylogenomics-visualizer
pip install -r requirements.txt
```

### Requirements
- Python 3.9+
- `biopython`
- `matplotlib`
- `numpy`
- `pillow`
- `python-docx`
- `aspose-words` (for vector EMF generation)

---

## 💻 Command Line Usage

### Scenario 1: Single Gene / Dataset Analysis
Generate all individual figures (Stepped SDT, SDT Heatmap, Linear ML Tree, Circular Cladogram) and the Combined Dual-Panel plate:
```bash
python3 viral_phylogenomics_visualizer.py \
  --fasta path/to/sequences.fasta \
  --gene-name TargetGene \
  --output-dir ./output
```

### Scenario 2: Paired / Dual-Gene Comparative Pipeline
When analyzing two related genes:
```bash
python3 viral_phylogenomics_visualizer.py \
  --fasta1 path/to/gene1_sequences.fasta \
  --name1 Gene1 \
  --fasta2 path/to/gene2_sequences.fasta \
  --name2 Gene2 \
  --output-dir ./output
```

### Scenario 3: Converting Figures to Vector EMF for Microsoft Office
```bash
python3 -c "
import aspose.words as aw, glob, os
out_dir = './output'
for png in glob.glob(os.path.join(out_dir, '*.png')):
    emf = png.replace('.png', '.emf')
    doc = aw.Document()
    builder = aw.DocumentBuilder(doc)
    builder.insert_image(png)
    doc.save(emf)
    print(f'Converted {emf}')
"
```

### Scenario 4: Exporting 960 DPI Print TIFFs (macOS)
```bash
sips -s format tiff -s dpiHeight 960.0 -s dpiWidth 960.0 \
  output/Figure_TargetGene_SDT_Identity_Matrix.png \
  --out output/Figure_TargetGene_SDT_Identity_Matrix_960dpi.tiff
```

---

## 🧬 Supported Sequence Formats

The pipeline accommodates any sequence header style:
- **Clean Isolate IDs (Accession-Free)**: `>Taxon_A`, `>Isolate_1`, `>Ref_Seq` (renders crisp isolate badges without printing empty `[]` brackets).
- **NCBI Formatted Headers**: `>Isolate_1 [MW123456]` (automatically parses accession numbers).
- **Reference Recognition**: Automatically highlights reference genomes (`● Reference`) distinct from study isolates (`★ Isolate`).

---

## 🛠️ Fault-Tolerant Architecture
- **Apple Silicon Native**: Prioritizes ARM64 native binaries (`iqtree-2.4.0-macOS-arm`) to bypass Rosetta x86 execution errors.
- **Pure-Python MSA Fallback**: Runs progressive pairwise alignment using `Bio.Align.PairwiseAligner` when external MAFFT binaries are unavailable.
- **Polar Legend Separation**: Custom polar coordinate geometry eliminates legend overlapping.

---

## 📚 Citation

When utilizing this pipeline in published work, please cite:
1. **SDT**: Muhire, B. M., Varsani, A., & Martin, D. P. (2014). SDT: a computer program for classifying and dividing sequences based on pairwise sequence identity. *PLoS ONE*, 9(9), e108277.
2. **IQ-TREE 2**: Minh, B. Q., et al. (2020). IQ-TREE 2: New models and efficient methods for phylogenomic inference in the genomic era. *Molecular Biology and Evolution*, 37(5), 1530–1534.
3. **MAFFT**: Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7. *Molecular Biology and Evolution*, 30(4), 772–780.
4. **BioPython**: Cock, P. J., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423.

---

## 📄 License
This project is open-source and licensed under the [MIT License](LICENSE).
