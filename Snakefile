'''
A pipline for nanopore cdna
You can run like this:
snakemake -j 6 -c "qsub -N {rulename} -l nodes=1:ppn={threads} -l mem=32g -j oe -l walltime=900:00:00"
'''

import os
import glob

configfile: 'config.yml'


SAMPLE_NAME = [os.path.basename(path).replace('.fastq.gz', '') for path in glob.glob('basecalled_data/*.fastq.gz')]
STRAND = ['positive', 'negative']
TYPE = ['spliced', 'intronless']


rule all:
    input:
        # full_length transcript:
        expand('full_length_transcripts/{sample_name}.full_length.sorted.bam', sample_name=SAMPLE_NAME),
        expand('full_length_transcripts/{sample_name}.non_full_length.sorted.bam', sample_name=SAMPLE_NAME),
        # gene list
        expand('gene_list/{sample_name}.full_genelist.bed', sample_name=SAMPLE_NAME),
        # isoform
        expand('splicing_isoforms/{sample_name}.full_length_isoform.bed', sample_name=SAMPLE_NAME)


#######################################
# Pre-process
#######################################

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
        ref=config['genome']
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

#######################################
# Get gene list
#######################################

rule bam_to_bed:
    input:
        'aligned_data/{sample_name}.tagged.mm2.sorted.bam'
    output:
        'gene_list/{sample_name}.tagged.mm2.sorted.bed',
    threads: 1
    shell:
        '''
        bedtools bamtobed -i {input} > {output}
        '''

# separating bedfile between positive and negative reads
rule separate_bed_file:
    input:
        'gene_list/{sample_name}.tagged.mm2.sorted.bed'
    output:
        positive = 'gene_list/{sample_name}.positive.bed',
        negative = 'gene_list/{sample_name}.negative.bed'
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
        bed = 'gene_list/{sample_name}.{strand}.bed',
        reference = config['genome']
    output:
        temp('gene_list/{sample_name}.{strand}.fasta')
    threads: 1
    shell:
        '''
        bedtools getfasta -s -fi {input.reference} -bed {input.bed} -fo {output}
        '''


rule remapping:
    input:
        reads='gene_list/{sample_name}.{strand}.fasta',
        ref=config['genome']
    output:
        temp_bam=temp('gene_list/{sample_name}.{strand}.tmp.bam'),
        bam='gene_list/{sample_name}.{strand}.mm2.sorted.bam'
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
        'gene_list/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'gene_list/{sample_name}.{strand}_coverage.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -d -ibam {input} > {output}
        '''

# 将正负链覆盖度文件合并
rule combine_coverage:
    input:
        positive='gene_list/{sample_name}.positive_coverage.tsv',
        negative='gene_list/{sample_name}.negative_coverage.tsv'
    output: 
        'gene_list/{sample_name}.total_coverage.tsv'
    threads: 1
    shell:
        '''
        python script/combine_coverage.py --positive {input.positive} --negative {input.negative} --output {output}
        '''

# 找连续的overlap的reads区间
rule find_coverage_breaks:
    input:
        'gene_list/{sample_name}.{strand}_coverage.tsv'
    output:
        'gene_list/{sample_name}.{strand}_genelist_incomplete.bed'
    params:
        min_coverage_threshold = 6
    threads: 1
    shell:
        '''
        python script/find_coverage_breaks.py --infile {input} --outfile {output} --min_coverage_threshold {params.min_coverage_threshold}
        '''
        

rule genome_cov_3dash:
    input:
        'gene_list/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'gene_list/{sample_name}.{strand}.3dash_positions.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -dz -3 -ibam {input} > {output}
        '''


rule genome_cov_5dash:
    input:
        'gene_list/{sample_name}.{strand}.mm2.sorted.bam'
    output:
        'gene_list/{sample_name}.{strand}.5dash_positions.tsv'
    threads: 1
    shell:
        '''
        bedtools genomecov -dz -5 -ibam {input} > {output}
        '''


# Intersect genes start and end positions from the coverage analysis and cuts to refine our analysis
rule insert_dash_cuts:
    input:
        in3dash = 'gene_list/{sample_name}.{strand}.3dash_positions.tsv',
        in5dash = 'gene_list/{sample_name}.{strand}.5dash_positions.tsv',
        chromosome_data = 'gene_list/{sample_name}.{strand}_genelist_incomplete.bed'
    output:
        'gene_list/{sample_name}.{strand}.genelist.bed'
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
        positive = 'gene_list/{sample_name}.positive.genelist.bed',
        negative = 'gene_list/{sample_name}.negative.genelist.bed'
    output:
        'gene_list/{sample_name}.full_genelist.bed'
    threads: 1
    shell:
        '''
        python script/combine_chr_lists.py --positive {input.positive} --negative {input.negative} --outfile {output}
        '''


#######################################
# Find splice sites
#######################################

rule bam_to_bed_split:
    input:
        'aligned_data/{sample_name}.tagged.mm2.sorted.bam'
    output:
        'splicing_isoforms/{sample_name}.tagged.mm2.sorted.split.bed'
    threads: 1
    shell:
        '''
        bedtools bamtobed -split -i {input} > {output}
        '''

# Finding splicing isoforms
rule find_splice_sites:
    input:
        'splicing_isoforms/{sample_name}.tagged.mm2.sorted.split.bed'
    output:
        'splicing_isoforms/{sample_name}.spliced.tsv'
    threads: 1
    shell:
        '''
        python script/find_splice_sites.py --infile {input} --outfile {output}
        '''


rule find_intronless_transcripts:
    input:
        'splicing_isoforms/{sample_name}.tagged.mm2.sorted.split.bed'
    output:
        'splicing_isoforms/{sample_name}.intronless.tsv'
    threads: 1
    shell:
        '''
        python script/find_intronless_transcripts.py --infile {input} --outfile {output}
        '''


# Translating to bed format for the final results
rule convert_isoforms_to_bed:
    input:
        'splicing_isoforms/{sample_name}.{type}.tsv'
    output:
        'splicing_isoforms/{sample_name}.{type}.bed'
    threads: 1
    shell:
        '''
        python script/convert_isoforms_to_bed.py --infile {input} --outfile {output}
        '''


rule get_full_length_isoform:
    input:
        'splicing_isoforms/{sample_name}.{type}.bed'
    output:
        'splicing_isoforms/{sample_name}.{type}.full_length.bed'
    threads: 1
    shell:
        '''
        bedtools intersect -a {input} -b supplementary_data/representative_exon/representative_gene_first_exon.bed -wa -s > {output}
        '''

rule merge_bedfile:
    input:
        spliced = 'splicing_isoforms/{sample_name}.spliced.full_length.bed',
        intronless = 'splicing_isoforms/{sample_name}.intronless.full_length.bed'
    output:
        'splicing_isoforms/{sample_name}.full_length_isoform.bed'
    threads: 1
    shell:
        '''
        cat {input.spliced} {input.intronless} | bedtools sort -i - > {output}
        '''