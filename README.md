# 6110 Assignment2: Global Transcriptome Changes Associated with Biofilm Development in Wine Yeast

## Introduction
Biofilm formation is a conserved microbial survival strategy that enables communities of cells to persist under environmental conditions that would otherwise be lethal to suspended individuals. In *Saccharomyces cerevisiae*, this capacity is demonstrated by the formation of a surface biofilm known as a velum, which develops at the air-liquid interface during the biological aging of sherry-type wines. This structure is produced by specialized flor strains that have acquired the physiological capacity to transition from fermentative to oxidative metabolism following the depletion of fermentable carbon sources, allowing sustained growth on ethanol and other non-fermentable substrates in an environment characterized by high ethanol concentrations, low pH, and severe nutrient limitation (Legras et al., 2016; Moreno-García et al., 2018).
Velum development proceeds through a series of morphologically and physiologically distinct stages, reflecting progressive adaptation to the wine surface environment. The transition from early to mature biofilm involves coordinated changes across multiple cellular systems, including restructuring of the cell wall composition to increase surface hydrophobicity, reprogramming of central carbon metabolism toward oxidative pathways, activation of stress tolerance mechanisms, and regulation of cell adhesion to maintain biofilm cohesion (Zara et al., 2005; Mardanov et al., 2020). Flor strains are distinguished from non-biofilm-forming laboratory strains in part by specific genetic modifications to the FLO11 promoter region that constitutively elevate its transcription, representing a key evolutionary adaptation to the wine aging environment (Fidalgo et al., 2006; Eldarov et al., 2018).

Despite the biological and industrial significance of velum formation, the genome-wide transcriptional dynamics underlying its staged development remain incompletely resolved. Mardanov et al. (2020) generated the first RNA-seq transcriptome dataset spanning three stages of velum development in the industrial flor strain I-329, identifying broad transcriptional changes in metabolic, structural, and stress-response gene categories. However, that study employed an older differential expression methodology and reference-guided alignment approach that predate significant advances in both quantification accuracy and statistical modeling. Reanalysis of this dataset using contemporary tools offers the opportunity to obtain more accurate transcript abundance estimates and more sensitive detection of differentially expressed genes.

The present study reanalyzes the Mardanov et al. (2020) dataset using a modern RNA-seq analysis pipeline. Transcript-level quantification was performed using Salmon, which employs a quasi-mapping algorithm with built-in correction for GC-content and sequence-specific biases, offering superior computational efficiency and quantification accuracy relative to traditional alignment-based approaches (Patro et al., 2017; Soneson et al., 2015). Differential expression analysis was conducted using DESeq2, which implements a negative binomial generalized linear model framework with empirical Bayes shrinkage of dispersion estimates, providing robust and well-calibrated statistical inference appropriate for experiments with limited biological replication (Love et al., 2014). To capture both stage-specific transcriptional changes and expression trajectories spanning the full developmental time course, pairwise contrasts were complemented by likelihood ratio testing against an intercept-only reduced model. Functional interpretation of differentially expressed gene sets was performed using clusterProfiler, enabling Gene Ontology over-representation analysis, KEGG pathway enrichment, and gene set enrichment analysis within a unified analytical framework (Wu et al., 2021). Together, this reanalysis aims to provide a statistically rigorous characterization of the transcriptional program governing velum maturation in *S. cerevisiae*.
## Methods

## Results
<img width="2000" height="1428" alt="image" src="https://github.com/user-attachments/assets/57e6ff3e-e227-4fce-929c-c077da30cb52" />

<img width="2000" height="1500" alt="image" src="https://github.com/user-attachments/assets/eb8658ad-853f-43d2-b8a6-f91c036f4bc1" />

<img width="2000" height="997" alt="image" src="https://github.com/user-attachments/assets/b773a55c-4e17-4607-927f-540bae1c13fb" />

<img width="2000" height="1428" alt="image" src="https://github.com/user-attachments/assets/dae4191f-735a-49c3-bd1d-daf86a771b0e" />

<img width="2000" height="1428" alt="image" src="https://github.com/user-attachments/assets/25b43085-5f77-4926-bd4c-90b226979384" />

<img width="2000" height="1428" alt="image" src="https://github.com/user-attachments/assets/4a45362e-065a-46f9-bf65-670c6d0be93f" />

<img width="1259" height="983" alt="image" src="https://github.com/user-attachments/assets/9abd8e02-234a-4c8b-b7eb-8186102ab8ee" />

<img width="2000" height="2500" alt="image" src="https://github.com/user-attachments/assets/b9a236f5-1ec2-4342-bea8-9abed8b380e3" />

## Discussion

## Conclusion
