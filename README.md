# MetalNet 

MetalNet is designed for detecting metal binding sites with coevolution data. And metal binding sites and metal type will be offered if predicted.

# preparation of coevolution profile

Could be downloaded from GREMLIN server.

(reference) 
jackhmmer -E 1E-20 -N 8 -A $gene.sto $gene.fasta $sequence_database
sto2fas.pl $gene.sto > $gene.fas
hhfilter -i $gene.fas -id 90 -cov 75 -o $gene.msa
gremlin_cpp -i $gene.msa -o $gene.pair

# running script

```bash
trim_coevolution.sh $gene
python predict.py $gene
```
