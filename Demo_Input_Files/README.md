# TerraTools demo input files

Small, purpose-built inputs for all 21 tools. Each folder holds `input/` (the files to load), `expected_output/` (what the tool produced when we ran it, so you can compare)
and `README.txt` (the exact settings to enter, the run time and the numbers you should see).

| # | Tool | Run time* | Headline result |
|---|------|-----------|-----------------|
| 01 | Phenotypic Data Analysis | 5-30 s per analysis | Yield heritability 0.49 (plot) / 0.74 (entry mean); Wet_2024 yields most |
| 02 | Mendelian Segregation | 10 s | 37 of 334 markers distorted after FDR correction (95 uncorrected), 297 kept; a second, dense simulated RIL map (4,499 markers) shows three clear distortion regions (chr3, chr6, chr9) |
| 03 | Linkage Mapping | 25 s | 345 markers read, 334 mapped, 11 faulty markers caught; 12 chromosomes |
| 04 | QTL Mapping | 5 s per run | CIM finds the 2 planted Yield QTL (chr1, chr5); the control trait gives none |
| 05 | QTL-Seq | 10 s | One QTL on chr2, 8.8-15.2 Mb, peak 11.9 Mb (planted: 12 Mb) |
| 06 | Diversity Analysis | 10 s | Shannon 0.36 -> 1.09 over 5 varieties; clusters equal the varieties |
| 07 | GWAS | 12 s | OLS 31 / MLM 4 / FarmCPU 2 significant SNPs; FarmCPU: chr1 and chr5 |
| 08 | ET-GWAS | 10 s | Two MTAs: the planted common SNP (chr3) and the planted rare SNP (chr8) |
| 09 | PyBSA-Seq | 50 s | A single significant region on chr2 (about 4.7-15 Mb) |
| 10 | Genes Extraction | 2 s | 22 genes in the QTL regions (QTL table in cM plus a marker position file or VCF) |
| 11 | Gene Prioritization | about 3 min (needs internet) | 6 known yield genes score 49-62, 2 controls score 0 |
| 12 | Clustering & Genomic Selection | 5 s | 5 clusters; 140 train / 36 test |
| 13 | Additive Dominance | 10 s | 176 x 176 additive and dominance matrices |
| 14 | BreedSight | 1-2 min | Real data of the Rice Diversity Panel 1 (seed width, 377 accessions, source in its README): BreedSight R2 0.57 / r 0.76 beats GBLUP 0.41, RBLUP 0.38 and Lasso 0.38 |
| 15 | Optimal Cross Selection | 12 s | 20 crosses, 27 parents, usefulness predicted for 1,770 crosses; far lower coancestry than truncation selection |
| 16 | Heterosis | 90 s | Cross-validation r 0.63; top hybrid BSP2-121 x BSP2-80 |
| 17 | Genotype Viewer | 5 s | 12 chromosomes drawn for 12 rice RILs, 308 markers; same picture from the CSV and the VCF |
| 18 | Favourable Allele Analysis | 1 s + 5 s | 327 / 327 trait-gene markers matched in a 1,027-marker rice 4K file (24 lines); balloon heatmap of 12 genes |
| 19 | Genofile Converter | 1 s | The start of the chain: HapMap, VCF or raw allele files of a genotyping company become one coded 0/1/2 table (used by tools 02 and 03) (345 markers read, 334 coded, 11 left out) |
| 20 | Recombination Analysis | 5 s | 405 crossovers in 80 simulated lines (489 planted; 83 % found at a 2 % block); mean donor 15.3 %; 37 intragenic |
| 21 | In silico BSA | 5 s | 3 regions in 300 simulated F2 plants, one of them fully over-dominant (found only by Gw) |

\* Measured on the development PC with the settings in each README (some use fewer simulation replicates, epochs or folds than the defaults to stay well under 5 minutes).

**What is real and what is simulated.** The genotypes in tools 02-04, 07, 12, 13, 15-17 and 19 come from 176 real rice RILs; tool 14 (BreedSight) uses real genotypes and a real trait of the public Rice Diversity Panel; the traits (and all data of tools 01, 05, 06, 08, 09, 20, 21), bulks, haplotype blocks and panels are simulated with known answers
(each README says what was planted), so you can see that the tool finds it. Gene Prioritization uses real rice genes and live public databases, so its scores can drift slightly.
Results of the neural-network tools (BreedSight, Heterosis) vary a little between runs.
