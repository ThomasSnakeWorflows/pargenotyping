# pargenotyping
Workflow for genotyping using paragraph

```bash
module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load bioinfo/snakemake-4.8.0
paragraph_bin = "/home/faraut/dynawork/SeqOccin/giab/softwares/paragraph/bin"
export PATH=$paragraph_bin:$PATH
```

Running the CNVPipeline

```bash
snakemake --jobs 30 --cluster-config cluster.yaml --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" -p -n
```
