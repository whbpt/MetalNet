import numpy as np
data=[[0,1,1,0],[1,0,1,0],[1,1,0,0],[0,0,0,0]]
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

