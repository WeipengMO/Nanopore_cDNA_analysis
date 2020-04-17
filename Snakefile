'''
A pipline for calling polya tail length
You can run like this:
snakemake -j 6 -c "qsub -N {rulename} -l nodes=1:ppn={threads} -l mem=32g -j oe -l walltime=900:00:00"
'''

import os
import glob

configfile: 'config.yml'


SAMPLE_NAME = [os.path.basename(path).replace('.fastq.gz', '') for path in glob.glob('basecalled_data/*.fastq.gz')]
STRAND = ['positive', 'negative']
    
rule all:
    input:
        #expand('aligned_data/{sample_name}.adjust.mm2.sorted.bam', sample_name=SAMPLE_NAME),
        expand('transitional/{sample_name}.total_coverage.tsv', sample_name=SAMPLE_NAME),
        #expand('transitional/{sample_name}.{strand}.genelist.bed', sample_name=SAMPLE_NAME, strand=STRAND)
        expand('transitional/{sample_name}.full_genelist.bed', sample_name=SAMPLE_NAME)
        

# qcat only accept ungzip file
rule decompress:
    input:
        'basecalled_data/{sample_name}.fastq.gz'
    output:
        temp('aligned_data/{sample_name}.fastq')
    threads: 20
    shell:
        '''
        pigz -p {threads} -cd {input} > {output}
        '''


# Remove adapter from reads.
# Search for adapters in the whole read
rule qcat_trim:
    input:
        'aligned_data/{sample_name}.fastq'
    output:
        'trim_data/{sample_name}.qcat.fastq.gz'
    threads: 1
    shell:
        '''
        qcat -f {input} -o {output} -t {threads} --trim --detect-middle
        '''
        
# Adjust read strandedness and polya filter
rule adjust_strandedness:
    input:
        fastq='trim_data/{sample_name}.qcat.fastq.gz',
        polya_file=config['polya_data']
    output:
        'trim_data/{sample_name}.adjust.fastq.gz'
    threads: 1
    shell:
        '''
        python script/adjust_strandedness.py -i {input.fastq} -o {output} -p {input.polya_file}
        '''


rule minimap2_mapping:
    input:
        reads='trim_data/{sample_name}.adjust.fastq.gz',
        ref=config['reference']
    output:
        temp_bam=temp('aligned_data/{sample_name}.adjust.tmp.bam'),
        bam='aligned_data/{sample_name}.adjust.mm2.sorted.bam'
    threads: 64
    shell:
        '''
        minimap2 -t {threads} -ax splice -ub -k14 -G 10000 \
            --secondary=no --cs {input.ref} {input.reads} \
            | samtools view -Sb - > {output.temp_bam}
        samtools sort -@ {threads} -o {output.bam} {output.temp_bam}
        samtools index {output.bam}
        '''


# add polya length and gene_id to bam tags
rule add_tags:
    input:
        infile='aligned_data/{sample_name}.adjust.mm2.sorted.bam',
        polya_file=config['polya_data'],
        gff_file=config['gff_file'],
    output:
        'aligned_data/{sample_name}.tagged.mm2.sorted.bam'
    threads: 1
    shell:
        '''
        python script/add_tag_to_bam.py -i {input.infile} -o {output} -p {input.polya_file} -g {input.gff_file}
        samtools index {output}
        '''
        

rule get_full_length_transcripts:
    input:
        infile='aligned_data/{sample_name}.tagged.mm2.sorted.bam',
        first_exon_path = config['first_exon_path']
    output:
        full_len='full_length_transcripts/{sample_name}.full_length.sorted.bam',
        non_full_len='full_length_transcripts/{sample_name}.non_full_length.sorted.bam'
    threads: 1
    shell:
        '''
        python script/get_full_length_transcripts.py -i {input.infile} --full_len {output.full_len} --non_full_len {output.non_full_len} --first_exon_path {input.first_exon_path}
        samtools index {output.full_len}
        samtools index {output.non_full_len}
        '''
        
rule bam_to_bed:
    input:
        'aligned_data/{sample_name}.tagged.mm2.sorted.bam'
    output:
        bed = 'transitional/{sample_name}.tagged.mm2.sorted.bed',
        split_bed = 'transitional/{sample_name}.tagged.mm2.sorted.split.bed'
    threads: 1
    shell:
        '''
        bedtools bamtobed -i {input} > {output.bed}
        bedtools bamtobed -split -i {input} > {output.split_bed}
        '''


rule separate_bed_file:
    input:
        'transitional/{sample_name}.tagged.mm2.sorted.bed'
    output:
        positive = 'transitional/{sample_name}.positive.bed',
        negative = 'transitional/{sample_name}.negative.bed'
    threads: 1
    run:
        with open(output.positive, "w") as positive, open(output.negative, "w") as negative, open(input[0], 'r') as infile:
            for line in infile:
                parts=line.rstrip().split('\t')
                if len(parts) >= 6 and parts[5] == "+":
                    positive.write(line)
                elif len(parts) >= 6 and parts[5] == "-":
                    negative.write(line)


rule bed_to_fasta:
    input:
        bed = 'transitional/{sample_name}.{strand}.bed',
        reference = config['reference']
    output:
        temp('transitional/{sample_name}.{strand}.fasta')
    threads: 1
    shell:
        '''
        bedtools getfasta -fi {input.reference} -bed {input.bed} -fo {output}
        '''


rule remapping:
    input:
        reads='transitional/{sample_name}.{strand}.fasta',
        ref=config['reference']
    output:
        temp_bam=temp('transitional/{sample_name}.{strand}.tmp.bam'),
        bam='transitional/{sample_name}.{strand}.mm2.sorted.bam'
    threads: 64
    shell:
        '''
        minimap2 -t {threads} -ax splice -ub -k14 -G 10000 \
            --secondary=no --cs {input.ref} {input.reads} \
            | samtools view -Sb - > {output.temp_bam}
        samtools sort -@ {threads} -o {output.bam} {output.temp_bam}
        samtools index {output.bam}
        '''

# 得到基因组单碱基覆盖度
rule genomecov:
    input:
        'transitional/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'transitional/{sample_name}.{strand}_coverage.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -d -ibam {input} > {output}
        '''

# 将正负链覆盖度文件合并
rule combine_coverage:
    input:
        positive='transitional/{sample_name}.positive_coverage.tsv',
        negative='transitional/{sample_name}.negative_coverage.tsv'
    output: 
        'transitional/{sample_name}.total_coverage.tsv'
    threads: 1
    shell:
        '''
        python script/combine_coverage.py --positive {input.positive} --negative {input.negative} --output {output}
        '''

rule find_coverage_breaks:
    input:
        'transitional/{sample_name}.{strand}_coverage.tsv'
    output:
        'transitional/{sample_name}.{strand}_genelist_incomplete.bed'
    params:
        min_coverage_threshold = 6
    threads: 1
    shell:
        '''
        python script/find_coverage_breaks.py --infile {input} --outfile {output} --min_coverage_threshold {params.min_coverage_threshold}
        '''
        

rule genome_cov_3dash:
    input:
        'transitional/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'transitional/{sample_name}.{strand}.3dash_positions.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -dz -3 -ibam {input} > {output}
        '''


rule genome_cov_5dash:
    input:
        'transitional/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'transitional/{sample_name}.{strand}.5dash_positions.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -dz -5 -ibam {input} > {output}
        '''

rule insert_dash_cuts:
    input:
        in3dash = 'transitional/{sample_name}.{strand}.3dash_positions.tsv',
        in5dash = 'transitional/{sample_name}.{strand}.5dash_positions.tsv',
        chromosome_data = 'transitional/{sample_name}.{strand}_genelist_incomplete.bed'
    output:
        'transitional/{sample_name}.{strand}.genelist.bed'
    threads: 1
    shell:
        '''
        python script/insert_dash_cuts.py \
            --in3dash_file {input.in3dash} \
            --in5dash_file {input.in5dash} \
            --outfile {output} \
            --chromosome_data {input.chromosome_data}
        '''
        

rule combine_chr_lists:
    input:
        positive = 'transitional/{sample_name}.positive.genelist.bed',
        negative = 'transitional/{sample_name}.negative.genelist.bed'
    output:
        'transitional/{sample_name}.full_genelist.bed'
    threads: 1
    shell:
        '''
        python script/combine_chr_lists.py --positive {input.positive} --negative {input.negative} --output {output}
        '''
    