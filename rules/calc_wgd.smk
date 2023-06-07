import os
import yaml
import pandas as pd

def get_remixtpp_input_paths(wildcards):
    return {
        'remixt_cn': runinfo.paths[(wildcards.sample_id, 'WGS-REMIXT-POSTPROCESS', 'remixt_cn')],
    }

rule proc_remixtpp_for_mt:
    input:
        unpack(get_remixtpp_input_paths)
        #remixtpp = config['remixtpp']
    output:
        os.path.join(config['results_dir'], '{sample_id}.cn.tsv')
    singularity: 
        "docker://soymintc/clickpdvcf:0.0.2"
    shell:
        'python scripts/proc_remixtpp_for_mt.py '
        '--in_pp {input.remixt_cn} --out_cn {output} '

def extract_purity_and_ploidy(sample_id):
    paths_path = config['metadata']
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==sample_id]
    metas = paths[paths['result_type']=='meta'][['isabl_sample_id', 'result_filepath']]
    meta_paths = dict(zip(metas['isabl_sample_id'], metas['result_filepath']))
    assert len(meta_paths) == 1, meta_paths
    sample2pp = {}
    for sample_id, meta_path in meta_paths.items():
        meta = yaml.load(open(meta_path), Loader=yaml.Loader)
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

def get_ploidy(sample_id):
    purity, ploidy = extract_purity_and_ploidy(sample_id)
    return ploidy

def get_purity(sample_id):
    purity, ploidy = extract_purity_and_ploidy(sample_id)
    return purity

rule calc_wgd:
    input:
        cn=os.path.join(config['results_dir'], '{sample_id}.cn.tsv'),
    output:
        wgd=os.path.join(config['results_dir'], '{sample_id}.wgd.tsv'),
    params:
        purity = lambda wildcards: get_purity(wildcards.sample_id),
        ploidy = lambda wildcards: get_ploidy(wildcards.sample_id),
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/calc_wgd.R '
        '{input.cn} '
        '{params.purity} ' # clonal_freq in MutationTimeR
        '{params.ploidy} ' # ploidy in MutationTimeR
        '{output.wgd} ' # wgd estimation
