
rule all:
    input:
        os.path.join(config['results_dir'], f'{config["sample_id"]}.wgd.tsv'),

include: "rules/common.smk"
include: "rules/process_variants.smk"
include: "rules/calc_wgd.smk"
