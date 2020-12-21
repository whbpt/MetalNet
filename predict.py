#!/usr/bin/env python
# coding: utf-8

import keras
import sys,os
import h5py
from keras.models import load_model
import pickle
import numpy as np
import networkx as nx
from graphviz import Graph as DOTGraph
from collections import Counter

model=load_model('database/dense_model_balanced8314.h5')


alphabet='CHDENSTKGQYLAVRIMFWP-'

states = len(alphabet)
a2n = {}
for a,n in zip(alphabet,range(states)):
    a2n[a] = n

def aa2num(aa):
    if aa in a2n: return a2n[aa]
    else: return a2n['-']
    

def parse_aln2seqmtx(filename,limit=-1):
    sequence = []
    f = open(filename, "r")
    for line in f.readlines():
        if line[0]=='>':
            continue
        line = line.strip()
        sequence.append(line)
    f.close()
    return np.array(sequence)

def frequency_matrix(seqs_mtx,i,j):
    fq_mtx=np.zeros((21,21))
    i=int(i)-1
    j=int(j)-1
    for r in range(seqs_mtx.shape[0]):   
        iaa=seqs_mtx[r][i]
        jaa=seqs_mtx[r][j]
        fq_mtx[aa2num(iaa),aa2num(jaa)]+=1
    return fq_mtx

def get_frequency_mtx(msa_file,i,j):
    fas_path=msa_file
    m=parse_aln2seqmtx(fas_path)
    fq_mtx=frequency_matrix(m,i,j)
    fq_mtx=fq_mtx/np.sum(fq_mtx)
    return fq_mtx



def main(gene,msa_file,contact_file):
    if os.path.getsize(msa_file) ==0:
        print( gene,'blank msafile')
        return
    with open(msa_file,'r') as f:
        if f.readline()[0]!= '>':
            print(gene,'no msa')
            return
    CHED_pair=['HH','HC','CH','CC','HD','DD','DH','CD','DC','EE','EC','CE','DE','ED','HE','EH']
    coevolution=[]
    state=True
    contact_handle=open(contact_file,'r')
    for line in contact_handle.readlines():
        line=line.strip()
        if len(line)==0:
            state=False
            break
        i_id=line.split()[-2][1:]
        j_id=line.split()[-1][1:]
        i_res=line.split()[-2][0]
        j_res=line.split()[-1][0]
        
        pair=i_res+j_res
        if pair in CHED_pair:
            coevolution.append((i_id+'_'+i_res,j_id+'_'+j_res))
    if state==False:
        print( gene,'break confile')
        return
    print(coevolution)
#########################predict the frequency matrix of each pair
    G=nx.Graph()
    output_file='%s.cnn'%gene 
    output_handle=open(output_file,'w')
    for pair in coevolution:
        i_AA=pair[0]
        j_AA=pair[1]
        i=pair[0].split('_')[0]
        j=pair[1].split('_')[0]
        data=get_frequency_mtx(msa_file,i,j)
        x_prediction = data.reshape(1,441)  
        weighted_prediction = model.predict(x_prediction,batch_size=1)
        pred = weighted_prediction[:,1]
        if pred > 0.5:
            G.add_edge(i_AA,j_AA)
        output_handle.write(i_AA+"\t"+j_AA+"\t"+str(pred)+"\n")        
    output_handle.close()
    nx.write_gpickle(G,'%s.gpickle'%gene )




def color_def(elem):
    if elem=="C":
        color="#FFED97"
    elif elem=="H":
        color="#66B3FF"
    elif elem=="E":
        color="#FFA6FF"
    elif elem =="D":
        color="#D3A4FF"
    return color

def find_ring(subG):
    g=nx.Graph(subG)
    while ([d<2 for (n,d) in list(g.degree())].count(True)>0):
        for (n,d) in list(g.degree()):
            if d<2:
                g.remove_node(n)
    return g

def resi_comparison(G1,G2):
    return G1['residue']==G2['residue']
def gap_comparison(G1,G2):
    return G1['gap_label']==G2['gap_label']




def motif_scan(G1,stdout):
    f=open('database/motif_bank.pkl','rb')
    motif_bank=pickle.load(f)
    f.close()
    result=[]
    for subG_metal in list(motif_bank.keys()):
        graph_list=motif_bank[subG_metal]
        for count,G2 in enumerate(graph_list):
            GM = nx.isomorphism.GraphMatcher(G1,G2,node_match=resi_comparison, edge_match=gap_comparison)
            if GM.is_isomorphic():  
                result.append(subG_metal)
                print(subG_metal,count,file=stdout)
                print(GM.mapping,file=stdout)
    result=dict(Counter(result))
    for k in list(result.keys()):
        ratio="%.2f"%(result[k]/len(motif_bank[k]))
        print(k,result[k],len(motif_bank[k]),ratio,file=stdout)




def option(gene):
    G=nx.read_gpickle('%s.gpickle'%gene)
    dot = DOTGraph(comment=gene,format='eps')
    GT_export=open("%s.dat"%gene,'w')
    for c in nx.connected_components(G):
        nodeset=G.subgraph(c).nodes()
        if len(nodeset)<3:          
            continue
        subG=nx.Graph(G.subgraph(c))
        subG=find_ring(subG)
        pair_list = list(subG.edges)
        print(gene,pair_list,file=GT_export)
        if len(pair_list) >=2:      
            for pair in pair_list:        
                i_tag,iAA=pair[0].split('_')
                j_tag,jAA=pair[1].split('_')
                i_id=i_tag+iAA            
                j_id=j_tag+jAA
                attrs={pair[0]:{'residue':iAA,'id':i_id},pair[1]:{'residue':jAA,'id':j_id}}
                nx.set_node_attributes(subG,attrs)
                gap=abs(int(i_tag)-int(j_tag))
                attrs={(pair[0],pair[1]):{'gap':gap,'gap_label': (gap<3)}}
                nx.set_edge_attributes(subG,attrs) 
                ##draw graph
                dot.node(str(i_id),style="radial",fillcolor=color_def(iAA))
                dot.node(str(j_id),style="radial",fillcolor=color_def(jAA))
                if gap < 3:
                    penwidth= '2'
                else:
                    penwidth= '1'
                dot.edge(str(i_id),str(j_id),penwidth=penwidth)
            motif_scan(subG,GT_export)
    dot.render('%s.gv'%gene,view=False)
    GT_export.close()





gene=sys.argv[1]
msa_file='%s.msa'%gene
contact_file='%s.copair'%gene
main(gene,msa_file,contact_file)
if os.path.exists("%s.gpickle"%gene):
    if os.path.getsize("%s.gpickle"%gene) >0:
        option(gene)

