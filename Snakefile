
import os
import pandas as pd
import numpy as np
import json


chrom_regex = re.compile('(chr)?([1-9][0-9]?|[XY])')


def get_max_reads(sample, sampleinfos):
    infos = pd.read_table(sampleinfos).set_index("id", drop=False)
    return 20*float(infos.loc[sample].depth)


def genome_depth(contigs_info):
    depth = []
    for contig in contigs_info:
        if chrom_regex.match(contig['name']):
            depth.append(float(contig['depth']))
    return np.mean(depth)


def sample_info(idxdepth_file):
    with open(idxdepth_file) as json_fin:
        loaded_json = json.load(json_fin)
        read_length = loaded_json["read_length"]
        bam_path = loaded_json["bam_path"]
        read_depth = genome_depth(loaded_json["contigs"])
    return (bam_path, read_depth, read_length)


samples = pd.read_table(config["sample_file"]).set_index("sample", drop=False)

reference = os.path.abspath(config['reference'])
inputvcf = os.path.abspath(config['inputvcf'])

workdir: config['workdir']

localrules: sampleinfo, variants, merge

wildcard_constraints:
    sample="[A-Za-z0-9\-]+"

rule all:
    input:
        "genotypes.vcf.gz"

rule depth:
    input:
        bam = lambda wildcards: samples.loc[wildcards.sample].bam,
        reference = reference
    output:
        "{sample}/depth.json"
    params:
        reference = reference
    threads:
        8
    shell:
        "idxdepth --threads {threads} -b {input.bam} -r {params.reference} "
        " > {output}"

rule sampleinfo:
    input:
        "{sample}/depth.json"
    output:
        "{sample}/sample.txt"
    run:
        with open(output[0], "w") as fout:
            fout.write("id\tpath\tdepth\tread length\n")
            bam_path, read_depth, read_length = sample_info(input[0])
            fout.write("%s\t%s\t%3.2f\t%d\n" %
                       (wildcards.sample, bam_path, read_depth, read_length))

rule variants:
    input:
        inputvcf
    output:
        temp("{sample}/input.vcf")
    shell:
        "bcftools view -s {wildcards.sample} {input} > {output}"


rule paragraph:
    message: """ --- trigger the paragraph genotyping --- """
    input:
        inputvcf = "{sample}/input.vcf",
        sampleinfo = "{sample}/sample.txt",
        reference = reference
    output:
        "{sample}/genotypes.vcf.gz"
    run:
        max_reads = get_max_reads(wildcards.sample, input.sampleinfo)
        cmd = "multigrmpy.py "
        cmd += " -i %s" % input.inputvcf
        cmd += " -m %s" % input.sampleinfo
        cmd += " -r %s" % input.reference
        cmd += " -M %d" % max_reads
        cmd += " -o %s" % wildcards.sample
        shell(cmd)
        shell("tabix %s" % output[0])


rule merge:
    input:
        expand("{sample}/genotypes.vcf.gz", sample=samples.index)
    output:
        "genotypes.vcf.gz"
    params:
        num_sample=lambda w: len(samples.index)
    shell:
        """
        if [ "{params.num_sample}" != 1 ];
        then
            bcftools merge {input} -Oz -o {output}
        else
            cp {input} {output}
        fi
        """
