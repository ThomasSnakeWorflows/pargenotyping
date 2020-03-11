# pargenotyping
Workflow for genotyping using paragraph

```bash
module load bioinfo/samtools-1.9
module load bioinfo/bcftools-1.9
module load system/Python-3.6.3
python3 -m venv parenv
source parenv/bin/activate
pip install -r requirements.txt

paragraph_bin="/home/faraut/dynawork/SeqOccin/giab/softwares/paragraph/bin"
export PATH=$paragraph_bin:$PATH
```

Running the CNVPipeline

```bash
snakemake --configfile config.yaml \
          --cluster-config cluster.yaml \
          --drmaa " --mem-per-cpu={cluster.mem}000 --mincpus={threads} --time={cluster.time} -J {cluster.name} -N 1=1" \
          --jobs 30  \
          -p -n
```
