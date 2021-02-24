configfile: "conf/config.yml"
SAMPLES = config["samples"]
CORES = config["resources"]["cores"]
RAM = config["resources"]["ram"]

rule all:
    input:
        'data/reference/gencode.vM25.transcripts.fa.gz',
        'data/reference/gencode.vM25.annotation.gtf',
        'intermediate/salmon_index',
        expand('results/tables/salmon/{sample}/quant.sf.gz', sample=SAMPLES),
        expand('results/tables/salmon/{sample}/quant.genes.sf.gz', sample=SAMPLES)

rule download_transcriptome:
    output:
        'data/reference/gencode.vM25.transcripts.fa.gz'
    log:
        "logs/wget_transcriptome_fasta.log"
    threads: 1
    params:
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz',
        dir = 'data/reference'
    shell:
        '''
        export ftp_proxy=http://192.168.100.1:80
        wget -v -o {log} -P {params.dir} {params.url}
        '''

rule download_GTF:
    output:
        'data/reference/gencode.vM25.annotation.gtf'
    log:
        "logs/wget_GTF_fasta.log"
    threads: 1
    params:
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz',
        dir = 'data/reference'
    shell:
        '''
        export ftp_proxy=http://192.168.100.1:80
        wget -v -o {log} -P {params.dir} {params.url}
        gunzip {output}.gz
        '''

rule salmon_index:
    input:
        'data/reference/gencode.vM25.transcripts.fa.gz'
    output:
        dir = directory('intermediate/salmon_index'),
    threads:
        CORES
    log:
        'logs/salmon/salmon_index.log'
    shell:
        '''
        salmon index \
            -p {threads} \
            --gencode \
            -t {input} \
            -i {output.dir} \
            2> {log}
        '''

rule salmon_quant:
    input:
        ref = 'intermediate/salmon_index',
        gtf = 'data/reference/gencode.vM25.annotation.gtf',
        R1 = 'data/fastq/{sample}_1.fastq.gz',
        R2 = 'data/fastq/{sample}_2.fastq.gz'
    output:
        quant = 'results/tables/salmon/{sample}/quant.sf',
        quant_genes = 'results/tables/salmon/{sample}/quant.genes.sf'
    params:
        out_dir = 'results/tables/salmon/{sample}'
    threads:
        CORES
    log:
        'logs/salmon/{sample}_quant.log'
    shell:
        '''
        salmon quant \
            -p {threads} \
            -l A \
            -i {input.ref} \
            -g {input.gtf}\
            -1 {input.R1} \
            -2 {input.R2} \
            --numBootstraps 100 \
            --seqBias \
            --gcBias \
            --validateMappings \
            -o {params.out_dir} 2> {log}
        '''

rule sf_gzip:
    input:
        quant = 'results/tables/salmon/{sample}/quant.sf',
        quant_genes = 'results/tables/salmon/{sample}/quant.genes.sf'
    output:
        quant = 'results/tables/salmon/{sample}/quant.sf.gz',
        quant_genes = 'results/tables/salmon/{sample}/quant.genes.sf.gz'
    shell:
        '''
        gzip {input.quant} && gzip {input.quant_genes}
        '''
