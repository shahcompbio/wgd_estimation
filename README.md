# Calc WGD pipeline
A simple pipeline for running estimating WGD status with Isabl consensus somatic maf and ReMixT post-processing results.

## Install
Cloning the directory is enough. However, you need to provide an accurate Isabl sample ID. See example in the next section.
```
git clone git@github.com:shahcompbio/wgd_estimation.git
```

## Example
At least running the following code will result in a successful run with `${out_dir}/results/${sample_id}/${sample_id}.wgd.tsv` as an output.
```
out_dir=/path/to/out/dir  # you need write privilege
sample_id=DAH123a  # matching Isabl sample ID

bash run_snakemake.sh $out_dir $sample_id
```

## Settings
1. Create an `$out_dir` for running `run_snakemake.sh`.
2. Make sure that the `$sample_id` input for `run_snakemake.sh` is exact with the Isabl sample ID included in the metadata file.
3. For first-time runs I recommend running `run_snakemake.sh` with `--dry-run`.
