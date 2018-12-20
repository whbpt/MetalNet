#!/usr/bin/env python
from __future__ import print_function
import os,sys
import pickle
import keras
from keras import backend as K
from keras.models import load_model
import numpy as np
import requests
from networkx import nx
from utils import *
from collections import Counter
################
import warnings
warnings.filterwarnings("ignore")
################
try:
    input_ID = str(sys.argv[1])
except:
    print("please input ID number")
    sys.exit()
coevolution_dict=dict()
contact_handle = requests.get('https://gremlin2.bakerlab.org/preds_cst.php?db=SUB&id='+input_ID)
for line in contact_handle.text.split("\n")[1:]:
	line_list = line.split()
	if len(line_list) == 9:
		i_id = line_list[0]
		j_id = line_list[1]
		i_res = line_list[2]
		j_res = line_list[3]
		pair = i_res + j_res
		if pair in CHED_pair:
			coevolution_dict[(i_id,j_id)]=np.zeros((21,21))
seq_num=0
# function for parsing FASTA file
def parse_fasta(lines):
	'''function to parse fasta'''
	header = []
	sequence = []
	for line in lines:
		if len(line) > 0:
			if line[0] == ">":
				header.append(line[1:])
				sequence.append([])
			else:
				sequence[-1].append(line)
	sequence = [''.join(seq) for seq in sequence]
	return header, sequence

msa_handle = requests.get('http://openseq.org/sub_fasta.php?id='+input_ID)
msa_lines = msa_handle.text.split("\n")
headers, sequences = parse_fasta(msa_lines)

origin_seq = sequences[0]
for sequence in sequences:
	sequence=list(sequence)
	seq_num+=1
	for pair in coevolution_dict.keys():
		i_id=pair[0]
		j_id=pair[1]
		i_res=sequence[int(i_id)-1]
		j_res=sequence[int(j_id)-1]
		index_i=AA_dict[i_res]
		index_j=AA_dict[j_res]
		coevolution_dict[pair][index_i,index_j]+=1
seq_len=len(sequence)
#########################predict the frequency matrix of each pair
states = 21
model = load_model('database/contact_model.h5')
Protein=PROTEIN()
for pair in coevolution_dict.keys():
	data1=coevolution_dict[pair]/seq_num
	data1=transform_matrix(data1)
	x_prediction = data1.reshape(1, states, states, 1)
	classes = model.predict(x_prediction,batch_size=1)
	i_id=pair[0]
	j_id=pair[1]
	i_tag=i_id+"_"+origin_seq[int(i_id)-1]
	j_tag=j_id+"_"+origin_seq[int(j_id)-1]
	for i,num in enumerate(classes.argmax(axis=1)):
		if num==1:
			# print(i_tag,j_tag,num,classes[i][num])
			Protein.metal.append((i_tag,j_tag))
			if i_tag not in Protein.site:
				Protein.site.append(i_tag)
			if j_tag not in Protein.site:
				Protein.site.append(j_tag)
			Protein.prob_dict[(i_tag,j_tag)]=classes[i][num]


##################################################
# to load the 2D bank information
###############################################
bench_file="database/benchmark.pkl"
if os.path.exists(bench_file):
	f2=file(bench_file,'rb')
	bench_r_dict=pickle.load(f2)
	f2.close()
else:
    print("benchmark database file missing")

#########################################
#to separate each subgraph as output 
# and to predict the metal type
#########################################
flag=0
dot = Graph(comment='The Round Table',format='png')
new_matrix=list2mtx(Protein.site,Protein.metal)
for (i,mtx) in enumerate(subgraph(new_matrix)):
    mtx=find_ring(mtx)
    pair_list=mtx2list(Protein.site,mtx)
    mtx_site=node_count(pair_list)
    if len(mtx_site)<=3:
        continue
    if find_ring(mtx,False):
    	temp_str=str()
    	for line in pair_list:
            for a in line:
    	        temp_str=temp_str+" "+a
    	        temp_str=temp_str+";"
    	#### set a graph G for comparison with 2D bank
    	G1=nx.Graph()
    	for pair in pair_list:
    	    i_tag=pair[0]
    	    j_tag=pair[1]
    	    i_id=i_tag.split("_")[0]+i_tag.split("_")[1]
    	    j_id=j_tag.split("_")[0]+j_tag.split("_")[1]
    	    i=int(i_tag[0:-2])
    	    j=int(j_tag[0:-2])
    	    if fabs(i-j) <=3:
    	    	gap_type='N'
    	    else:
    	    	gap_type='F'
    	    #####trying to read the probability with wrong order
    	    try:
    	    	prob=Protein.prob_dict[pair]				
    	    except:
    	    	prob=Protein.prob_dict[(j_tag,i_tag)]
    	    ##### adding nodes and edges in the G graph
    	    G1.add_node(i_id,residue=i_tag.split("_")[1])
    	    G1.add_node(j_id,residue=j_tag.split("_")[1])
    	    G1.add_edge(i_id,j_id,gap_type=gap_type)
            dot.node(i_id,style="radial",fillcolor=color_def(i_tag.split("_")[1]))
            dot.node(j_id,style="radial",fillcolor=color_def(j_tag.split("_")[1]))
            dot.edge(i_id,j_id,penwidth="2",color=color_digit(float(prob)*2-1))
    	    flag=1
        pred_metal=list()
    	#########transfer the metal type from the same graph
    	for bench_G in bench_r_dict.keys():
    	    GM=nx.isomorphism.GraphMatcher(bench_G,G1,node_match=nx.isomorphism.categorical_node_match(['residue','residue','residue','residue'],['C','D','E','H']),edge_match=nx.isomorphism.categorical_edge_match(['gap_type','gap_type'],['N','F']))
            if G1.number_of_nodes() == bench_G.number_of_nodes():
                if GM.is_isomorphic() or GM.subgraph_is_isomorphic():
    	            pred_metal=bench_r_dict[bench_G]
    	            break
    	result = []
    	if len(pred_metal)>0:
            result.append(str(pred_metal[0]))
            dot.edge(str(pred_metal[0]),j_id,penwidth="0",color="white")
    	else:
    	    result.append(str(pred_metal))
if flag >0:
    dot.render(input_ID, view=False)
