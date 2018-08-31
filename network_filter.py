#!/usr/bin/env python
import os
import sys
from graphviz import Graph
import numpy as np
#data=[[0,1,1,0],[1,0,1,0],[1,1,0,1],[0,0,1,0]]
def find_ring(data,mod_mtx_bool=True):
	connectmatrix=np.array(data)
	dim=connectmatrix.shape
	numnode=dim[1]
	Node=np.array(range(0,numnode))
	pair_dict={}
	connect_dict={}
	for node in Node:
    		currentline=connectmatrix[node,:]
    		for currentnode in range(0,numnode):
        		if currentline[currentnode] == 1:
            			pair_dict.setdefault(node,[]).append(currentnode)
	connect=int()
	for node in range(0,numnode):
		if node in pair_dict.keys():
			connect=len(pair_dict[node])
		else:
			connect=0
		connect_dict.setdefault(node,connect)
	#print connect_dict
	while 1 in connect_dict.values():
		for n in range(0,numnode):
			if connect_dict[n]==1:
				connect_dict[n]=0               
				for connect_node in pair_dict[n]:
					connect_dict[connect_node]=connect_dict[connect_node]-1
					pair_dict[n].remove(connect_node)
					pair_dict[connect_node].remove(n)                    
	ring_node=[]
	for i in connect_dict.keys():
		if connect_dict[i]>0:
			ring_node.append(i)
	#print ring_node
	ring_matrix=np.zeros([numnode,numnode])
	ring_submatrix=connectmatrix[ring_node]
	ring_submatrix=ring_submatrix[:,ring_node]
#print ring_submatrix
	if mod_mtx_bool==True:
		for i in range(0,len(ring_node)):
			for j in range(0,len(ring_node)):
				ring_matrix[ring_node[i],ring_node[j]]=ring_submatrix[i,j]
		return ring_matrix
	else:
		if len(ring_node)>0:
			return True
		else:
			return False

def subgraph(data):
	connectmatrix=np.array(data)
	dim=connectmatrix.shape
	numnode=dim[1]
	Node=np.array(range(0,numnode))
	Branches=[]
	quence=[]
	neighborNode=[]
	while len(np.where(Node >-1)[0])>0:
		quence=np.where(Node>-1)[0][0:1]
		subField=[]
		while len(quence)!=0:
			currentNode=quence[0]
			quence=np.delete(quence,0)
			subField.append(currentNode)
			Node[int(currentNode)]=-1
			currentline=connectmatrix[currentNode,:]
			neighborNode=np.where(currentline>0)[0]
			for k in neighborNode:
				if Node[k]!=-1 :
					quence=np.append(quence,k)
					Node[k]=-1
		#for i in range(numnode-len(subField)):
		#	subField.append(-1)
		Branches.append(subField)
	Branch_mtx=[]
	for line in Branches:
		sub_mtx=np.zeros((numnode,numnode))
		for i in line:
			for j in line:
				sub_mtx[i,j]=connectmatrix[i,j]
		Branch_mtx.append(sub_mtx)
	return Branch_mtx	
def list2mtx(elem_list,adj_list):
	length=len(elem_list)
	new_matrix=np.zeros((length,length))
	for line in adj_list:
		i=line[0]
		j=line[1]
		new_matrix[elem_list.index(i),elem_list.index(j)]=1
		new_matrix[elem_list.index(j),elem_list.index(i)]=1
	return new_matrix
def mtx2list(elem_list,adj_mtx):
	elem_index=np.where(adj_mtx>0)
	adj_list=list()
	for i in range(0,len(elem_index[0])):
		elem_i=elem_index[0][i]
		elem_j=elem_index[1][i]
		if elem_i<elem_j:
			adj_list.append((elem_list[elem_i],elem_list[elem_j]))
	return adj_list
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

class PROTEIN:
    def __init__(self):
        self.site=[]
        self.metal_both=[]
        self.prob_dict={}

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
