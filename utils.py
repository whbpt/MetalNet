import numpy as np
from graphviz import Graph
from math import fabs
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
######define a class that can store the results
class PROTEIN:
    def __init__(self):
        self.site=[]
        self.metal=[]
        self.prob_dict={}
	self.pair=[]
	self.subgraphs=[]
########a list than contains the combination of all the CHEDs.
CHED_pair=['HH','HC','CH','CC','HD','DD','DH','CD','DC','EE','EC','CE','DE','ED','HE','EH']
##########transform the amino acid into number,the special amino acids are all treated with 20
AA_dict={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'-':20,'X':20,'U':20,'Z':20,'B':20,'J':20,'O':20}
metal_list=["Calcium","Cobalt","Copper","Iron","Iron-sulfur","Magnesium","Manganese","Molybdenum","Nickel","Potassium","Sodium","Zinc"]
