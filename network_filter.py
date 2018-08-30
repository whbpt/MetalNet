#!/usr/bin/env python
import os
import sys
import pickle
import math
import numpy as np 
from subgraph import *
from graphviz import Graph
class PROTEIN:
    def __init__(self):
        self.site=[]
        self.metal_both=[]
        self.prob_dict={}
def color_def(residue):
	elem=residue.split("_")[1]
	if elem=="C":
		color="#FFED97"
	elif elem=="H":
		color="#66B3FF"
	elif elem=="E":
		color="#FFA6FF"
	elif elem =="D":
		color="#D3A4FF"
	return color
def color_digit(num):
	if num>0.9:
		color="#000000"
	elif num>0.8:
		color="#3C3C3C"
	elif num>0.7:
		color="#5B5B5B"
	elif num >0.6:
		color="#7B7B7B"
	elif num >0.5:
		color="#9D9D9D"
#	elif num >0.4:
#		color="#BEBEBE"
#	elif num >0.3:
#		color="#E0E0E0"
	else:
		color="#BEBEBE"
	return color
def node_count(pair_list):
	mtx_site=[]
	for line in pair_list:
		i_tag=line[0]
		j_tag=line[1]
		if i_tag not in mtx_site:
			mtx_site.append(i_tag)
		if j_tag not in mtx_site:
			mtx_site.append(j_tag)
	return mtx_site
input_ID=sys.argv[1]
cnn_file="test_data/output_"+input_ID+".cnn"
nf_file="test_data/output_"+input_ID+".nf"#network_filter
render_file="test_data/output_"+input_ID+".render.gv"#network_filter
Protein=PROTEIN()
for line in open(cnn_file,'r'):
	line=line.strip()
	pair_i=line.split()[0]
	pair_j=line.split()[1]
	prob=line.split()[2]
	Protein.metal_both.append((pair_i,pair_j))
	if pair_i not in Protein.site:
		Protein.site.append(pair_i)
	if pair_j not in Protein.site:
		Protein.site.append(pair_j)
	Protein.prob_dict[(pair_i,pair_j)]=prob
nf_handle=open(nf_file,'w')
dot = Graph(comment='The Round Table',format='png')
flag=0
temp_pair=list()
new_matrix=list2mtx(Protein.site,Protein.metal_both)
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
	#if True:
		for pair in pair_list:
			i_tag=pair[0]
			j_tag=pair[1]
			i_id=i_tag.split("_")[0]+i_tag.split("_")[1]
			j_id=j_tag.split("_")[0]+j_tag.split("_")[1]
			try:
				prob=Protein.prob_dict[pair]				
			except:
				prob=Protein.prob_dict[(j_tag,i_tag)]				
			dot.node(i_id,style="radial",fillcolor=color_def(i_tag))
			dot.node(j_id,style="radial",fillcolor=color_def(j_tag))
			dot.edge(i_id,j_id,penwidth="2",color=color_digit(float(prob)*2-1))
			nf_handle.write(i_id+"\t"+j_id+"\t"+prob+"\n")
			flag=1
			temp_pair.append(pair)
if flag>0:
	dot.render(render_file, view=False)
nf_handle.close()
