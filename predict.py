#!/usr/bin/env python
from __future__ import print_function
import keras
from keras import backend as K
from keras.models import load_model
import numpy as np
from graphviz import Graph

################ transfrom the matrix into a specific order
################, which related with metal binding orders
def transform_matrix(matrix):
        origin_seq='ARNDCQEGHILKMFPSTWYV-'
        new_seq='CHDENSTKGQYLAVRIMFWP-'
        new_matrix=np.zeros((21,21))
        for (i,line) in enumerate(matrix):
                for (j,elem) in enumerate(line):
                        new_i=new_seq.index(origin_seq[i])
                        new_j=new_seq.index(origin_seq[j])
                        new_matrix[new_i,new_j]=matrix[i,j]
        return new_matrix


############find the closed ring from the adjacent matrix
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
	ring_matrix=np.zeros([numnode,numnode])
	ring_submatrix=connectmatrix[ring_node]
	ring_submatrix=ring_submatrix[:,ring_node]

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

###########finding subgraphs from the whole predicted graph

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


########count nodes
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

#####some funky function that transform adjacent list to matrix
def list2mtx(elem_list,adj_list):
	length=len(elem_list)
	new_matrix=np.zeros((length,length))
	for line in adj_list:
		i=line[0]
		j=line[1]
		new_matrix[elem_list.index(i),elem_list.index(j)]=1
		new_matrix[elem_list.index(j),elem_list.index(i)]=1
	return new_matrix
 
#####some funky function that transform adjacent matrix to list
def mtx2list(elem_list,adj_mtx):
	elem_index=np.where(adj_mtx>0)
	adj_list=list()
	for i in range(0,len(elem_index[0])):
		elem_i=elem_index[0][i]
		elem_j=elem_index[1][i]
		if elem_i<elem_j:
			adj_list.append((elem_list[elem_i],elem_list[elem_j]))
	return adj_list

#####some helper function to transform the residue type into diffrent colors
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

#####some helper function to transform the probablity into diffrent colors
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
	else:
		color="#BEBEBE"
	return color

######define a class that can store the results
class PROTEIN:
    def __init__(self):
        self.site=[]
        self.metal_both=[]
        self.prob_dict={}
###########################download the msa file and contact map file from website
msa_file="examples/input_1535637622.msa"
contact_file="examples/input_1535637622.con"
############################get the coevolution pair from .con file
########a list than contains the combination of all the CHEDs.
CHED_pair=['HH','HC','CH','CC','HD','DD','DH','CD','DC','EE','EC','CE','DE','ED','HE','EH']
##########transform the amino acid into number,the special amino acids are all treated with 20
AA_dict={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'-':20,'X':20,'U':20,'Z':20,'B':20,'J':20,'O':20}

coevolution_dict=dict()
contact_handle=open(contact_file,'r')
contact_handle.readline()#discard the header
for line in contact_handle:
	line=line.strip()
	i_id=line.split()[0]
	j_id=line.split()[1]
	i_res=line.split()[2].split("_")[1]
	j_res=line.split()[3].split("_")[1]
	pair=i_res+j_res
	if pair in CHED_pair:
		coevolution_dict[(i_id,j_id)]=np.zeros((21,21))
seq_num=0
msa_handle=open(msa_file,'r').readlines()
msa_handle="".join(msa_handle)#.split(">")
origin_seq="".join(msa_handle.split(">")[1].split("\n")[1:])
for line in msa_handle.split(">")[1:]:
	sequence="".join(line.split("\n")[1:])
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
model = load_model('models/contact_model.h5')
Protein=PROTEIN()
for pair in coevolution_dict.keys():
	data1=coevolution_dict[pair]/seq_num
	data1=transform_matrix(data1)
	x_prediction = data1.reshape(1, states, states, 1)
	classes =model.predict(x_prediction,batch_size=1)
	i_id=pair[0]
	j_id=pair[1]
	i_tag=i_id+"_"+origin_seq[int(i_id)-1]
	j_tag=j_id+"_"+origin_seq[int(j_id)-1]
	for i,num in enumerate(classes.argmax(axis=1)):
		if num==1:
			print(i_tag,j_tag,num,classes[i][num])
			Protein.metal_both.append((i_tag,j_tag))
			if i_tag not in Protein.site:
				Protein.site.append(i_tag)
			if j_tag not in Protein.site:
				Protein.site.append(j_tag)
			Protein.prob_dict[(i_tag,j_tag)]=classes[i][num]

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
			flag=1
			temp_pair.append(pair)
if flag>0:
	dot.render("test_out", view=False)
