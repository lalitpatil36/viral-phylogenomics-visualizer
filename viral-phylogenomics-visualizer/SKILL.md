---
name: viral-phylogenomics-visualizer
description: >-
  Generalized, publication-grade phylogenomics and sequence identity visualization pipeline for ANY nucleotide or amino acid sequence dataset. Generates: (1) Pure lower-triangular stepped SDT identity matrices (SDT v1.2 style), (2) All-against-all SDT heatmaps, (3) Midpoint-rooted Linear Maximum-Likelihood (ML) phylogenetic trees with study isolate badges and bootstrap supports, (4) Polar circular cladograms with concentric metadata rings (Host, Region, Clade) and dedicated non-overlapping legends, and (5) Combined dual-panel figures. Automatically exports in PNG, 300/960 DPI TIFF, PDF, and vector EMF formats.
---

# Generalized Viral Phylogenomics & Identity Visualization Skill

This skill defines the standardized, publication-grade phylogenomic analysis and visualization pipeline for **ANY** sequence dataset (e.g., *Citrus tristeza virus*, *Begomovirus*, *Potyvirus*, *Tospovirus*, plant/animal viruses, bacteria, or eukaryotic genes).

It supports sequences in standard FASTA format, with or without NCBI GenBank accession numbers, and seamlessly distinguishes study isolates from reference sequences.

---

## Key Visual Outputs

For each target gene/dataset (e.g. `p23`, `CP`, or any generic marker `GeneX`), the pipeline produces:

| Figure Name | Description | Publication Benchmark |
| :--- | :--- | :--- |
| **`Figure_[Gene]_SDT_Identity_Matrix`** | **Pure lower-triangular stepped SDT matrix** with white upper triangle, rainbow colorbar in upper right, dashed diagonal clade lines, and taxon labels flush against staircase cells. | Matches original Sequence Demarcation Tool (SDT v1.2, Muhire et al.) |
| **`Figure_[Gene]_SDT_Heatmap`** | **Full all-against-all pairwise identity matrix heatmap** ordered by phylogenetic leaf order. | Comprehensive identity overview |
| **`Figure_[Gene]_Linear_Phylogeny`** | **Midpoint-rooted, ladderized ML tree** with bootstrap supports ($\ge 50\%$), branch scale bar, dotted leader lines, study isolate badges (`★ Isolate`), and reference indicators (`● Reference`). | Standard linear phylogram |
| **`Figure_[Gene]_Circular_Phylogeny`** | **Polar circular cladogram** ($15^\circ$ to $345^\circ$) with clade sectors, concentric metadata rings, and an **isolated non-overlapping legend canvas** at the base. | Multi-dimensional metadata visualization |
| **`Combined_[Gene]_Phylogeny_and_Identity_Matrix`** | **Dual-panel composite figure**: Panel (a) Linear ML tree with clade brackets + Panel (b) Stepped SDT identity matrix. | High-impact dual journal figure |

---

## File Export Standard

Every figure is automatically exported in four industry-standard formats:

1. **PNG (`.png`)**: High-resolution raster (300 DPI) for web viewing and presentation slides.
2. **TIFF (`.tiff`)**: LZW-compressed 300 DPI raster, alongside **960 DPI ultra-high-resolution print TIFF** (`_960dpi.tiff`) generated via `sips` for top-tier journal submissions.
3. **PDF (`.pdf`)**: Vector graphic file with embedded vector fonts, ideal for Adobe Illustrator or inkscape touch-ups.
4. **EMF (`.emf`)**: Windows Enhanced Metafile vector format generated via `aspose.words` for loss-free scaling and editing in Microsoft Word, PowerPoint, and Excel.

---

## Pipeline Location & Architecture

The master automated pipeline script is located in this directory:
```bash
/Users/lalit/lp/MAPI/nwang/viral_phylogenomics_visualizer.py
```

### Universal Compatibility & Fallbacks:
- **Architecture Aware**: Seamlessly detects Apple Silicon ARM64 native binaries (`iqtree-2.4.0-macOS-arm/bin/iqtree2`), system x86_64 binaries, or falls back to fast BioPython TreeConstruction if external binaries are absent.
- **Multiple Sequence Alignment (MSA)**: Executes `mafft --auto` when available, with an automatic pure-Python progressive alignment fallback (`Bio.Align.PairwiseAligner`) to guarantee error-free execution on any environment.
- **Accession-Free Support**: Accepts clean sequence IDs (`D1`, `D3`, `Ref. Seq.`) as well as full NCBI headers (`Isolate_D1 [MW123456]`). Never outputs empty brackets `[]` or dummy accessions.
- **Non-Overlapping Polar Layout**: Circular cladograms use a fixed dedicated legend band in the bottom quadrant ($0.05 \le y \le 0.17$), preventing collisions between tree branches, outer text labels, and legend boxes.

---

## Command Line Usage

### Scenario 1: Single Gene / Dataset Analysis
Generate all individual figures (Stepped SDT, SDT Heatmap, Linear ML Tree, Circular Cladogram) and Combined Dual-Panel for a single sequence file:
```bash
python3 /Users/lalit/lp/MAPI/nwang/viral_phylogenomics_visualizer.py \
  --fasta /path/to/sequences.fasta \
  --gene-name p23 \
  --output-dir /Users/lalit/lp/MAPI/nwang
```

### Scenario 2: Paired / Dual Gene Comparison
When analyzing two related genes (e.g. Coat Protein and p23 suppressor, or Gene A and Gene B):
```bash
python3 /Users/lalit/lp/MAPI/nwang/viral_phylogenomics_visualizer.py \
  --fasta1 /path/to/gene1_sequences.fasta \
  --name1 CTV_CP \
  --fasta2 /path/to/gene2_sequences.fasta \
  --name2 CTV_p23 \
  --output-dir /Users/lalit/lp/MAPI/nwang
```

### Scenario 3: Generating EMF Vector Files
To convert any generated `.png` or `.pdf` figure in the output directory into vector EMF:
```bash
python3 -c "
import aspose.words as aw, glob, os
out_dir = '/Users/lalit/lp/MAPI/nwang'
for png in glob.glob(os.path.join(out_dir, '*.png')):
    emf = png.replace('.png', '.emf')
    doc = aw.Document()
    builder = aw.DocumentBuilder(doc)
    shape = builder.insert_image(png)
    doc.save(emf)
    print(f'Created {emf}')
"
```

### Scenario 4: Generating 960 DPI TIFF for Print
To create journal-ready 960 DPI TIFFs:
```bash
sips -s format tiff -s dpiHeight 960.0 -s dpiWidth 960.0 \
  Figure_CTV_p23_SDT_Identity_Matrix.png \
  --out Figure_CTV_p23_SDT_Identity_Matrix_960dpi.tiff
```

---

## Methodological Guidelines

1. **Pairwise Identity Calculation**:
   - Matches standard SDT v1.2 metric:
     $$\text{Identity} = \frac{\text{Identical Sites}}{\text{Alignment Length} - \text{Shared Internal Gaps}} \times 100$$
2. **Stepped Matrix Geometry**:
   - The matrix renders strictly lower-triangular cells $(row \ge col)$.
   - Cell annotations display integer percentages (e.g., `100`, `98`, `95`).
   - For datasets $> 25$ taxa, numeric cell text is automatically suppressed to preserve clean aesthetics, while color gradients remain prominent.
3. **Phylogeny Ordering**:
   - Both SDT matrix and heatmap order taxa according to phylogenetic tree leaf traversal, ensuring biologically related isolates cluster into distinct blocks along the diagonal.
4. **Provenance & Neutral Descriptors**:
   - Avoid hardcoding geographic assumptions. When geographic origin is unverified or spans mixed regions, use phylogenetic/biological descriptors (e.g., `Intermediate Sub-cluster`, `Divergent Genotype`, `Major Clade A/B`).
