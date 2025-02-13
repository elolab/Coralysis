# Coralysis

<br>

## Overview

<br>

Coralysis is an R package featuring a multi-level integration algorithm for sensitive integration, reference-mapping, and cell state identification in single-cell data, described in the paper _["Coralysis enables sensitive identification of imbalanced cell types and states in single-cell data via multi-level integration"](https://doi.org/10.1101/2025.02.07.637023)_.

<br>

<p style="text-align:center;"> <img src="man/figures/Coralysis_flowchart.png" width="90%" /></p>
<font size="2"> <b>Coralysis integration flowchart.</b> (A) An input of heterogeneous scRNA-seq datasets are overclustered batch wise into a training set modelled through the Iterative Clustering Projection (ICP) algorithm in order to predict the cell cluster probabilities and obtain an integrated embedding. Adaptations to the original ICP algorithm (Smolander et al., 2021): (B) batch wise cluster assignment at start, dependent on the cell distribution across Principal Component 1 (median as cutoff); (C) training cells selected from batch k nearest neighbours of the cell with the highest probability for every batch per cluster; and, (D) upon ICP clustering convergence, each cluster is further divided into two for the next clustering round, dependent on the batch wise cluster probability distribution (median as cutoff). (E) Multi-level integration is achieved through multiple divisive clustering rounds, blending the batch effect and highlighting the biological signal incrementally. Shapes represent cell types and colours batches. 
</font>


<br>

---

<br>

## Installation

<br>

The latest version of `Coralysis` can be downloaded from GitHub using the [`devtools`](https://devtools.r-lib.org/) R package.

```R
devtools::install_github("elolab/Coralysis")
```

<br>

---

<br>

## Citation

<br>

>Sousa A, Smolander J, Junttila S, Elo L (2025). **Coralysis enables sensitive identification of imbalanced cell types and states in single-cell data via multi-level integration.** _bioRxiv_. [doi:10.1101/2025.02.07.637023](https://doi.org/10.1101/2025.02.07.637023)

<br>

---

<br>

## Contact information

<br>

If you have questions related to `Coralysis`, please contact us [here](https://github.com/elolab/Coralysis/issues). 

<br>

---

<br>

## References

<br>

1. Johannes Smolander, Sini Junttila, Mikko S Venäläinen, Laura L Elo. "ILoReg: a tool for high-resolution cell population identification from single-cell RNA-seq data". Bioinformatics, Volume 37, Issue 8, 15 April 2021, Pages 1107–1114, [https://doi.org/10.1093/bioinformatics/btaa919](https://doi.org/10.1093/bioinformatics/btaa919).

<br>

<br>

<br>
