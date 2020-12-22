# MetalNet 

MetalNet is designed for detecting metal binding sites with coevolution data. And metal binding sites and metal type will be offered if predicted.

## make prediction
what you need:

1. msa file  (*$gene.msa*)

2. coevolution profile  (*$gene.pair*)

`#`**$gene**: input file name, such as **8f4a248217074d25a8bf3c71391e0bc9**.msa

```bash
cd demo/input
bash ../../trim_coevolution.sh 8f4a248217074d25a8bf3c71391e0bc9
python ../../predict.py 8f4a248217074d25a8bf3c71391e0bc9
```
## result

*$gene.cnn*: coevolved CHED pairs with predicted possibility

*$gene.dat*: prediction result after graph-based filter 

*$gene.gv.eps*: visualization of predicted network cluster 

## preparation of coevolution profile

Could be downloaded from GREMLIN server or calculated locally.

(reference instruction:)

jackhmmer -E 1E-20  -N 8  -A *$gene.sto*   *$gene.fasta*   {sequence_database}

sto2fas.pl  *$gene.sto*  >  *$gene.fas*

hhfilter  -i  *$gene.fas*  -id 90  -cov 75  -o *$gene.msa*

gremlin_cpp  -i *$gene.msa*  -o *$gene.pair*



