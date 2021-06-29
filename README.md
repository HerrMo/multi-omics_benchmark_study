# Code to the multi-omics benchmark study by Herrmann et al.

This repo allows to
- download the (preprocessed) data used in the study 
- reproduce the results (table, figures etc.) presented in the paper 
- rerun the entire benchmark experiment

**If you use the code or data please cite:** 

Moritz Herrmann, Philipp Probst, Roman Hornung, Vindi Jurinovic, Anne-Laure Boulesteix, Large-scale benchmark study of survival prediction methods using multi-omics data, Briefings in Bioinformatics, Volume 22, Issue 3, May 2021, bbaa167, https://doi.org/10.1093/bib/bbaa167


### To download the data: 
- The preprocessed data (described in the study) is available via [OpenML](https://www.openml.org/)
- The OpenML dataset IDs can be found in `data/datset_ids.txt` or `data/datset_ids.RData`
- Note that the datasets had to be split into two to three parts in order to be uploaded to OpenML
- **R users** can use the code in `R/bench_experiment.R` (lines [44-81](https://github.com/HerrMo/multi-omics_benchmark_study/blob/0b58b56a7ef80905d812f0a2644f9ae549394363/R/bench_experiment.R#L46-L81)) to directly download the data (and convert it to `mlr` tasks)

### To reprocude the results (in R):
- to only reproduce the tables, figures etc. displayed in the paper without rerunning the benchmark experiments use `reproduce_table-and-figures.Rmd`
- to rerun the full experiments (this takes several days or weeks, depending on the available resources) use `R/bench_experiment.R`
  - see the instructions in `R/packages.R`!
  - make sure the required packages are installed 
  - make sure to use correct package versions via [checkpoint](https://github.com/HerrMo/multi-omics_benchmark_study/blob/0b58b56a7ef80905d812f0a2644f9ae549394363/R/bench_experiment.R#L13-L14)
  - not all packages are covered by checkpoint, this is specifically relevant for `mlr` (s. `R/packages.R`)!
- to merge the benchmark results use `R/merge_bmr_results.R`

Note, `mlr` has deprecated (https://github.com/mlr-org/mlr) in the meantime. There is now the new framework `mlr3` (https://mlr3.mlr-org.com/).

