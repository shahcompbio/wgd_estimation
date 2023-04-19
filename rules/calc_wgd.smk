import os
import yaml
import pandas as pd

def extract_purity_and_ploidy(sample_id):
    paths_path = '/juno/work/shah/users/chois7/metacohort/paths.WGS-REMIXT-POSTPROCESS.tsv'
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==sample_id]
    metas = paths[paths['result_type']=='meta'][['isabl_sample_id', 'result_filepath']]
    meta_paths = dict(zip(metas['isabl_sample_id'], metas['result_filepath']))
    assert len(meta_paths) == 1, meta_paths
    sample2pp = {}
    for sample_id, meta_path in meta_paths.items():
        meta = yaml.load(open(meta_path))
        purity = sum(meta['mix'][1:])
        ploidy = meta['ploidy']
    return purity, ploidy

rule parse_purity_and_ploidy:
    input:
        vcf=os.path.join(config['results_dir'], '{sample_id}.vcf') # not used
    output:
        os.path.join(config['results_dir'], '{sample_id}.purity,ploidy.txt')
    run:
        # get purity, ploidy
        print("No problem before extract_purity")
        purity, ploidy = extract_purity_and_ploidy(wildcards.sample_id)

        with open(output[0], 'w') as outfile:
            outfile.write(f'{purity},{ploidy}\n')

def get_ploidy(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return ploidy

def get_purity(purity_and_ploidy):
    purity, ploidy = open(purity_and_ploidy, 'r').read().strip().split(',')
    return purity

rule calc_wgd:
    input:
        vcf=os.path.join(config['results_dir'], '{sample_id}.vcf'),
        cn=os.path.join(config['results_dir'], '{sample_id}.cn.tsv'),
        purity_and_ploidy=os.path.join(config['results_dir'], 
                                       '{sample_id}.purity,ploidy.txt')
    output:
        wgd=os.path.join(config['results_dir'], '{sample_id}.wgd.tsv'),
    params:
        purity = lambda wildcards, input: get_purity(input.purity_and_ploidy),
        ploidy = lambda wildcards, input: get_ploidy(input.purity_and_ploidy)
    log:
        os.path.join(config['log_dir'], '{sample_id}.pdf.log')
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/calc_wgd.R '
        '{input.vcf} {input.cn} '
        '{params.purity} ' # clonal_freq in MutationTimeR
        '{params.ploidy} ' # ploidy in MutationTimeR
        '{output.wgd} ' # wgd estimation
        '&> {log}'
