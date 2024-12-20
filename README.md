# Protein-Variant-Phenotype Study of NBAS using AlphaFold in the Aspect of SOPH Syndrome

This repository contains the supplementary materials and code for the article "Protein-variant-phenotype study of NBAS using AlphaFold in the aspect of SOPH syndrome".

## Contents

1. [Heatmap Data for Binding Sites](#heatmap-data-for-binding-sites)
2. [Pymol Scripting and Mapping Predicted Binding Sites](#pymol-scripting-and-mapping-predicted-binding-sites)
3. [PAE Plots and Mapping Binding Sites](#pae-plots-and-mapping-binding-sites)
4. [Mapping Phenotyping Traits onto PAE Plot](#mapping-phenotyping-traits-onto-pae-plot)
5. [Supplementary Materials](#supplementary-materials)

---

## Heatmap Data for Binding Sites

This section contains the scripts and data to generate heatmaps for the binding sites of the proteins using seaborn, matplotlib, and other libraries.

**Libraries Used:**
- seaborn
- matplotlib.pyplot
- matplotlib.colors.LinearSegmentedColormap
- pandas
- numpy
- json

**Files:**
- `Heatmap_AF3.ipynb`
- `Heatmap_AF.ipynb`

**Example Output:**

Heatmap of USE1/NBAS binding site.

<div align="center">
  <img src="Heatmap example.png" alt="Heatmap Example" width="600">
</div>

---

## Pymol Scripting and Mapping Predicted Binding Sites

This section includes the Pymol scripts and notebooks for mapping the predicted binding sites onto the NBAS protein.

**Files:**
- `figure 5 mol AF3.pse`
- `figure 5 mol.pse`
- `Figure 5.ipynb`
- `Figure 5-AF3.ipynb`

**Example Output:**

Binding sites are colored as follows: green - USE1, purple - Rab18, yellow - ZW10 (cyan stands for UPF3B putative binding site).

<div align="center">
  <img src="C-term NBAS.png" alt="C-term NBAS" width="600">
</div>

---

## PAE Plots and Mapping Binding Sites

This section covers the scripts for generating PAE plots and mapping the binding sites onto them.

**Libraries Used:**
- jsonlite
- ggnewscale
- scales
- reshape2
- ggplot2

**Files:**
- `Figure 6 batch bs comb AF3.R`
- `Figure 6 batch bs comb.R`
- `Figure 6.ipynb`
- `Figure 6-AF3.ipynb`

**Example Output:**

Binding sites are colored as follows: blue - USE1, purple - Rab18, yellow - ZW10.

<div align="center">
  <img src="Figures data/PAE/test_bs/AF3-BS5.png" alt="PAE Plot with Binding Sites" width="600">
</div>

---

## Mapping Phenotyping Traits onto PAE Plot

This section details the code for mapping phenotyping traits onto the PAE plot.

**Libraries Used:**
- jsonlite
- ggnewscale
- scales
- reshape2
- ggplot2

**Files:**
- `Figure 7 layers.R`
- `Figure 7 layers.ipynb`
- `Figure 7 pheno w-o sg.csv`

**Example Outputs:**

Subdomain affected by SOPH variant.

<div align="center">
  <img src="Figures data/PAE/1914_count-104.png" alt="Phenotyping Trait Mapping 1" width="600">
</div>

Protein-Variant-Phenotype (PVP) map of Osteogenesis Imperfecta.

<div align="center">
  <img src="Figures data/PAE_pheno/filter 5/PAE_OI.png" alt="Phenotyping Trait Mapping 2" width="600">
</div>

---

## Supplementary Materials

The supplementary materials include phenotype and variants data, and details about AlphaFold runs performed for the article.

**Files:**
- `Supplementary table 1.xlsx`
- `Supplementary table 2.ods`

All data is too large to be loaded into GitHub and is available by request. Please contact me through my [personal GitHub page](https://github.com/LyonyaZhozhikov).

<div id="header" align="center">

#### Published in:
 
##### 
<a href="https://onlinelibrary.wiley.com/journal/10970134"><img src="https://shields.io/badge/PROTEINS:%20Structure%20Function%20and%20Bioinformatics-Wiley-green?logo=#E9711C&PROTEINS:%20Structure%20Function%20and%20Bioinformatics=Wiley" alt="ARTICLE"/></a>

##### 
<a href="https://onlinelibrary.wiley.com/doi/10.1002/prot.26764"><img src="https://shields.io/badge/Full%20Access-Article-green?logo=#E9711C&Full%20Access=Article" alt="FullAccess"/></a>

<table align="center">
  <tr>
    <td align="left" style="width: 50%;">
      <a href="https://www.scimagojr.com/journalsearch.php?q=14283&amp;tip=sid&amp;exact=no" title="SCImago Journal &amp; Country Rank">
        <img src="https://www.scimagojr.com/journal_img.php?id=14283" alt="SCImago Journal &amp; Country Rank" />
      </a>
    </td>
    <td align="right" style="width: 50%;">
      <a href="https://wos-journal.info/journalid/14272" title="WOS-Journal.info">
        <img src="https://wos-journal.info/journalide/14272" width="320" height="120" alt="WOS-Journal.info" />
      </a>
    </td>
  </tr>
</table>

</div>

Please, also check my ResearchGate page for full texts and more: <a href="https://www.researchgate.net/profile/Leonid-Zhozhikov-3"><img src="https://shields.io/badge/ResearchGate-Page-green?logo=#E9711C&ResearchGate=Page" alt="Page"/></a>
