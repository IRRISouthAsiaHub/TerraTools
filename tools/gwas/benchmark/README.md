# Benchmarking the GWAS models against established programs

mrMLM and FASTmrMLM are not benchmarked here because the app runs the mrMLM R package itself (installed in the bundled R). None of the other reference programs (GAPIT, GEMMA, EMMAX, rMVP) was installed where these scripts were written, so **no
agreement with any of them has been measured**. The scripts below do the comparison once you have the programs. Until you have
run them, the models in the app are described as "GEMMA-style", "EMMAX-style" and so on, and no equivalence is claimed.

## 1. Write identical inputs

```bash
python -m tools.gwas.benchmark.make_inputs --scenario causal_3 --seed 7 --pcs 3 --out benchmark_run
```

This simulates rice-like data (three subpopulations, siblings, missing calls), runs this program's OLS, MLM, EMMAX, GEMMA, FarmCPU
and BLINK on it, and writes the genotypes, map, phenotype, relationship matrix and principal components in plain CSV so every
program uses the same samples, markers and covariates.

## 2. Run the reference programs

| Program | Install | Run |
|---|---|---|
| GAPIT (GLM, MLM, FarmCPU, Blink) | `source("http://zzlab.net/GAPIT/gapit_functions.txt")` in R, or https://github.com/jiabowang/GAPIT | `Rscript tools/gwas/benchmark/run_references.R benchmark_run` |
| rMVP | `install.packages("rMVP")` | same R script (prints the steps for your version) |
| GEMMA | binary from https://github.com/genetics-statistics/GEMMA/releases | commands below |
| EMMAX | https://genome.sph.umich.edu/wiki/EMMAX | commands below |

GEMMA (centred relationship matrix, Wald test, `-c` covariates = intercept + the PCs in `pcs.csv`):

```bash
# BIMBAM genotype file from numeric_genotypes.csv: one row per marker  SNP,A,T,g1,g2,...  (python snippet in gemma_external.write_inputs)
gemma -g geno.txt -p pheno.txt -gk 1 -o k
gemma -g geno.txt -p pheno.txt -c covariates.txt -k output/k.cXX.txt -lmm 1 -o gemma_wald
```

`tools/gwas/backend/gemma_external.py` writes those files and does the comparison automatically when `gemma` is on the PATH or the
"GEMMA executable" option is set.

EMMAX (p-values from `--emmax` with the same kinship, variance components estimated once):

```bash
emmax-kin -v -d 10 emmax_geno
emmax -v -d 10 -t emmax_geno -p pheno.tfam -k emmax_geno.hBN.kinf -c covariates.txt -o emmax_out
```

Save each reference result as a CSV with a marker column and a p-value column.

## 3. Compare

```bash
python -m tools.gwas.benchmark.compare benchmark_run --ours MLM --reference "GAPIT MLM=benchmark_run/gapit_MLM.csv:SNP:P"
python -m tools.gwas.benchmark.compare benchmark_run --ours GEMMA --reference "GEMMA Wald=benchmark_run/gemma_wald.csv:SNP:P"
```

Tolerances for the same model on identical inputs: Spearman correlation of -log10 p >= 0.99, largest difference in -log10 p <= 0.5,
overlap of the top 1% of markers >= 0.8. A reference that has not been run is reported as "not run".

Expect differences that are not errors: GAPIT compresses the relationship matrix and re-estimates variance components differently
for MLM; FarmCPU and BLINK differ in the order they select markers; mrMLM uses its own selection criterion. Report those
differences instead of tuning this program to match.

## Simulation study (no outside program needed)

```bash
python -m tools.gwas.validation.run_validation --reps 20 --out validation_results
```

See `tools/gwas/validation/run_validation.py` for the quantities it reports.
