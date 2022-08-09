rule mutationtimer:
    input:
        vcf=os.path.join(config['results_dir'], '{aliquot_id}.vcf'),
        cn=os.path.join(config['results_dir'], '{aliquot_id}.cn.tsv')
    output:
        pdf=os.path.join(config['results_dir'], '{aliquot_id}.pdf'),
        rdata=os.path.join(config['results_dir'], '{aliquot_id}.RData')
    log:
        os.path.join(config['log_dir'], '{aliquot_id}.pdf.log')
    params:
        clonal_freq = 0.54 # TODO: soft code
    singularity: 
        "docker://soymintc/mutationtimer:latest"
    shell:
        'Rscript scripts/run_mutationtimer.R '
        '{input.vcf} {input.cn} {params.clonal_freq} '
        '{output.pdf} {output.rdata} '
        '&> {log}'
