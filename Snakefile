include: "rules/common.smk"
include: "rules/calc_wgd.smk"

SAMPLES = runinfo.sample_ids

rule all:
    input:
        expand(os.path.join(config['results_dir'], '{sample_id}.wgd.tsv'), sample_id=SAMPLES),

