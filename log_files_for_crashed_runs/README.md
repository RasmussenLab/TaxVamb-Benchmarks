# Log files for crashed runs

For benchmarking for figure 3 in the paper the following runs crashed internally.

**Note on Snakemake resource fields:** In the log files, the relevant resources used by this pipeline are `mem_gb` and `walltime`. The other fields (`mem_mb`, `mem_mib`, `disk_mb`, `disk_mib`) are Snakemake defaults and are not used by this pipeline.

## SemiBin
SemiBin crashes in the Vaginal and Salvia datasets orginally. These errors were fixed.
- **Vaginal**: SemiBin crashes because one sample produces no bins, causing the `--write-pre-reclustering-bins` flag to fail with `FileNotFoundError: output_recluster_bins`. Version 2.2.0 (https://github.com/BigDataBiology/SemiBin/releases/tag/v2.2.0) is required as it includes a fix for "Do not fail if no bins are produced" (#170 & #173), but this fix alone does not solve the crash when `--write-pre-reclustering-bins` is also passed. Removing the `--write-pre-reclustering-bins` flag fixes the issue. See `semibin_vaginal_v2.2.0.log` for the log file.
- **Salvia**: SemiBin crashes with `ValueError: Input X contains NaN` in sklearn's NearestNeighbors using Semibin v.2.1.0. See https://github.com/BigDataBiology/SemiBin/issues/211 and https://github.com/BigDataBiology/SemiBin/issues/201. Upgrading SemiBin fixed this issue. See `semibin_salvia.log` for the log file for the crashed run.
- **Resources**: mem_gb=500, walltime=288:00:00 (v2.1.0 Vaginal, Salvia) / walltime=1088:00:00 (v2.2.0 Vaginal), GPU=1, 64 threads.

## COMEBin (multi-sample)
COMEBin (v1.0.3) crashes in the Human Gut (IBS) and Forest Soil datasets.
In both datasets the error is described in the following GitHub issue: https://github.com/ziyewang/COMEBin/issues/17.
The log files for these runs can be found in: `comebin_multisample_human_gut_ibs.log` and `comebin_multisample_forest_soil.log`. Upgrading comebin to version v.1.0.3 did not fix the issue (newest version as of Feb 12, 2026)
Since multi-sample COMEBin crashed on these datasets, they were run in single-sample mode instead. This is equivalent to assigning each read pair and their corresponding contig file to a different sample name in the `sample` column of the config files for our benchmarking pipeline. Additionally, for completeness, single-sample COMEBin was also run as an additional benchmark.
- **Resources**: Human Gut IBS: mem_gb=500, walltime=1088:00:00, GPU=1. Forest Soil: mem_gb=250, walltime=288:00:00, no GPU.

## COMEBin (single-sample)
COMEBin (single-sample) crashes on some of the samples in the Bee Hives (11 samples), Human Gut Antibiotics (1 sample), Human Gut IBS (2 samples), and Vaginal (7 samples) datasets.
COMEBin's internal `markerCmd` (`test_getmarker_2quarter.pl`) fails, This only happens to some samples in the single-sample mode. 
- **Resources**: mem_gb=500, walltime=1088:00:00, 64 threads, no GPU.
- **Log files**: `comebin_singlesample/bee_hives/`, `comebin_singlesample/human_gut_antibiotics/`, `comebin_singlesample/human_gut_ibs/`, `comebin_singlesample/vaginal/`.

## GUNC
GUNC crashes on Bee Hives.
GUNC's Diamond step finds no genes mapped to the reference database, producing no output files ("No diamond output files").
This was resolved by implementing the fix described in https://github.com/grp-bork/gunc/issues/42.
- **Resources**: mem_gb=250, walltime=20-00:00:00, 64 threads.
- **Log files**: `gunc_bee_hives_comebin_multisample.log`, `gunc_bee_hives_comebin_singlesample/`, `gunc_bee_hives_semibin.log`, `gunc_bee_hives_comebin_singlesample`.
