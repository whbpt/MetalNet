import numpy as np
from find_ring import *
data=[[0,1,1,0],[1,0,1,0],[1,1,0,1],[0,0,1,0]]
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

def degree_count(elem_list,pair_list):
	mtx=list2mtx(elem_list,pair_list)
	numnode=mtx.shape[1]
	connect_dict={}
	for node in range(0,numnode):
		connect=int(sum(mtx[node,:]))
		connect_dict[node]=connect
	return connect_dict	
#def find_ring(data,mod_mtx_bool=True,full_sub=False):
#	connectmatrix=np.array(data)
#	none_zero=np.where(sum(connectmatrix[0:]>0))[0]
#	new_mtx=np.zeros((len(none_zero),len(none_zero)))
#	for (i,ii) in enumerate(none_zero):
#		for (j,jj) in enumerate(none_zero):
#			new_mtx[i,j]=connectmatrix[ii,jj]
#	connectmatrix=new_mtx
#	dim=connectmatrix.shape
#	numnode=dim[1]
#	Node=np.array(range(0,numnode))
#	pair_list=list()
#	connect_dict={}
#	for node in Node:
#   		currentline=connectmatrix[node,:]
#   		for currentnode in range(0,node):
#   			if currentline[currentnode] == 1:
#				pair_list.append((int(currentnode),int(node)))
#	connect=int()
#	elem_list=list()
#	for line in pair_list:
#		if line[0] not in elem_list:
#			elem_list.append(line[0])
#		if line[1] not in elem_list:
#			elem_list.append(line[1])
#	connect_dict=degree_count(elem_list,pair_list)
#	while int(1) in connect_dict.values():
#		for n in range(0,numnode):
#			if connect_dict[n]==1:
#				temp=list()
#				for pair in pair_list:
#					if int(n) not in pair:
#						temp.append(pair)
#				pair_list=temp
#				connect_dict=degree_count(elem_list,pair_list)
#	ring_node=[]
#	for i in connect_dict.keys():
#		if connect_dict[i]>0:
#			ring_node.append(i)
##	print ring_node
#	ring_matrix=np.zeros([numnode,numnode])
#	ring_submatrix=connectmatrix[ring_node]
#	ring_submatrix=ring_submatrix[:,ring_node]
##print ring_submatrix
#	if mod_mtx_bool==True:
#		for i in range(0,len(ring_node)):
#			for j in range(0,len(ring_node)):
#				ring_matrix[ring_node[i],ring_node[j]]=ring_submatrix[i,j]
#		return ring_matrix
#	else:
#		if len(ring_node)>0:
#			return True 
#		else:
#			return False 
##print find_ring(data,True)
