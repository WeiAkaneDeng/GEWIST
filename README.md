# GEWIST
Prioritizing SNPs using Gene Environment Wide Interaction Search Threshold (GEWIST)

It is challenging to search for gene-environment interactions in a genome-wide setting because of the low statistical power and the heavy computational burden involved due to the large number of genetic variants. 

Pare et al. (2010) proposed a novel method - variance prioritization (VP) to prioritize single nucleotide polymorphisms (SNPs) by exploiting the interaction effects on the variance of quantitative traits.The prioritization is achieved by comparing the variance of a quantitative trait conditioned on three possible genotypes using Levene's test (Levene, 1960) for variance inequality. The variance prioritization procedure consists of two steps: 1) Select SNPs with Levene's test p-value lower than their individual optimal variance prioritization thresholds (eta_0); 2) Test the selected SNPs against all other SNPs (i.e. gene-gene) or environmental covariates (i.e. gene-environment) using a linear regression for interactions while correcting for eta_0*M tests, where M is the number of total SNPs tested).

To reduce both the stringent significance threshold due to multile hypothesis and the heavy computational burden, we introduced a fast algorithm, Gene Environment Wide Interaction Search Threshold (GEWIST; Deng and Pare 2011), to efficiently and accurately determine the optimal variance p-value threshold for individual SNPs, i.e. each SNP will be chosen according to different significant threshold based on the minor allele frequency and the expected effect sizes. The GEWIST package provides functions to prioritize SNPs using the algorithms described in Deng and Pare (2001).

The original submission can be found here: http://bioconductor.org/packages/release/bioc/html/GEWIST.html

The next release will include a likelihood ratio test for variance heteogeneity test that is expected to be more powerful than Levene's test, as well as functions that faciliate a meta-analysis of Levene's test across multiple studies.
