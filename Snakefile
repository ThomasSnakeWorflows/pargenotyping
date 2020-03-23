
import os
import numpy as np
import json
import subprocess
from collections import defaultdict

chrom_regex = re.compile('(chr)?([1-9][0-9]?|[XY])')

DEFAULT_THREADS = 2


def get_threads(rule, default=DEFAULT_THREADS):
    cluster_config = snakemake.workflow.cluster_config
    if rule in cluster_config and "threads" in cluster_config[rule]:
        return cluster_config[rule]["threads"]
    if "default" in cluster_config and "threads" in cluster_config["default"]:
        return cluster_config["default"]["threads"]
    return default


def get_max_reads(sample, sampleinfos):
    with open(sampleinfos, 'r') as file_in:
        headers = next(file_in, None).rstrip().split("\t")
        samples = defaultdict()
        for row in file_in:
            fields = row.rstrip().split("\t")
            samples[fields[0]] = dict(zip(headers, fields))
    return 20*float(samples[sample]['depth'])


def genome_depth(contigs_info):
    depth = []
    for contig in contigs_info:
        if chrom_regex.match(contig['name']):
            depth.append(float(contig['depth']))
    return np.mean(depth)


def sample_info(raw_sampleinfo_file, samplename):
    with open(raw_sampleinfo_file, 'r') as file_in:
        first_line = next(file_in, None).rstrip()[2:]
        headers = first_line.split("\t")
        samples = defaultdict()
        for row in file_in:
            fields = row.rstrip().split("\t")
            sample_d = dict(zip(headers, fields))
            samples[sample_d['sampleName']] = sample_d
    sample_info =  samples[samplename]
    return float(sample_info['meanCoverage'])


def read_samples(sample_file):
    with open(sample_file, 'r') as file_in:
        headers = next(file_in, None)
        samples = defaultdict()
        for row in file_in:
            fields = row.rstrip().split("\t")
            samples[fields[0]] = fields[1]
    return samples


samples = read_samples(config['sample_file'])
reference = os.path.abspath(config['reference'])
inputvcf = os.path.abspath(config['inputvcf'])

workdir: config['workdir']

localrules: sampleinfo, variants, merge, depth, allinfo

wildcard_constraints:
    sample="[A-Za-z0-9\-]+"


rule all:
    "samples.tar.gz",
    "genotypes.vcf.gz"


rule tarsamples:
    input:
        expand("{sample}.tar.gz", sample=samples.keys())
    output:
        "samples.tar.gz"
    shell:
        """
        mv -t samples {input}
        tar zcvf {output} samples
        rm -fr samples
        """


rule cleanparagraph:
    input:
        "{sample}/genotypes.vcf.gz",
        "genotypes.vcf.gz"
    output:
        tar="{sample}.tar.gz"
    shell:
        """
        tar cvzf {wildcards.sample}.tar.gz {wildcards.sample}
        rm -fr {wildcards.sample}
        """



rule allinfo:
    input:
        expand("{sample}/sampleinfo.txt", sample=samples.keys())
    output:
        "allsampleinfo.txt"
    shell:
        "cat {input} > {output}"

rule depth:
    input:
        bam = lambda wildcards: samples[wildcards.sample]
    output:
        "{sample}/depth.txt"
    params:
        region = config['statregion']
    threads:
        2
    shell:
        """
         sambamba depth region -L {params.region} {input.bam} -o {output}
        """

rule sampleinfo:
    input:
        bam = lambda wildcards: samples[wildcards.sample],
        depth = "{sample}/depth.txt"
    output:
        "{sample}/sampleinfo.txt"
    params:
        region = config['statregion']
    threads:
        1
    run:
        cmd = "samtools stats  "
        cmd += "%s " % input.bam
        cmd += "%s " % params.region
        cmd += " | grep ^FRL | cut -f 2-"
        result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        lines = result.stdout.decode("utf-8").split("\n")
        length, counts = [], []
        for line in lines:
            if line:
                l, c = line.split("\t")
                length.append(int(l))
                counts.append(int(c))
        read_length = np.average(length, weights=counts)
        with open(output[0], "w") as fout:
                fout.write("id\tpath\tdepth\tread length\n")
                read_depth = sample_info(input.depth, wildcards.sample)
                fout.write("%s\t%s\t%3.2f\t%d\n" %
                           (wildcards.sample, input.bam, read_depth, read_length))


rule variants:
    input:
        inputvcf
    output:
        "variants.vcf"
    shell:
        "bcftools view --drop-genotypes {input} -Ov -o {output}"


rule paragraph:
    message: """ --- trigger the paragraph genotyping --- """
    input:
        inputvcf = "variants.vcf",
        sampleinfo = "{sample}/sampleinfo.txt",
        reference = reference
    output:
        "{sample}/genotypes.vcf.gz"
    threads:
        get_threads("paragraph")
    run:
        max_reads = get_max_reads(wildcards.sample, input.sampleinfo)
        cmd = "multigrmpy.py "
        cmd += " -t %d" % threads
        cmd += " -i %s" % input.inputvcf
        cmd += " -m %s" % input.sampleinfo
        cmd += " -r %s" % input.reference
        cmd += " -M %d" % max_reads
        cmd += " -o %s" % wildcards.sample
        print(cmd)
        shell(cmd)
        shell("tabix %s" % output[0])


rule merge:
    input:
        expand("{sample}/genotypes.vcf.gz", sample=samples.keys())
    output:
        "genotypes.vcf.gz"
    params:
        num_sample=lambda w: len(samples.keys())
    shell:
        """
        if [ "{params.num_sample}" != 1 ];
        then
            bcftools merge {input} -Oz -o {output}
        else
            cp {input} {output}
        fi
        """
