from __future__ import print_function
import pickle
import numpy as np
import random,math
class PROTEIN:
	def __init__(self):
		self.site=[]
		self.metal=[]
		self.contact_map=[]
		self.metal_pair=[]
		self.none_pair=[]
		self.CHED_pair=[]
		self.prob=[]
		self.site_dict={}
		self.fulseq=str()
CHED_pair=['HH','HC','CH','CC','HD','DD','DH','CD','DC','EE','EC','CE','DE','ED','HE','EH']
metal_list=["Calcium","Cobalt","Copper","Iron","Iron-sulfur","Magnesium","Manganese","Molybdenum","Nickel","Potassium","Sodium","Zinc"]
hard_metal=["Calcium","Magnesium","Potassium","Sodium","Molybdenum"]
soft_metal=["Cobalt","Copper","Iron","Iron-sulfur","Manganese","Nickel","Zinc","Other"]
def transform_matrix(matrix):
	origin_seq='ARNDCQEGHILKMFPSTWYV-'
	new_seq='CHDENSTKGQYLAVRIMFWP-'
	new_seq2='CHDENQSTYWPGLAVIMFKR-'
	new_matrix=np.zeros((21,21))
	for (i,line) in enumerate(matrix):
		for (j,elem) in enumerate(line):
			new_i=new_seq.index(origin_seq[i])
			new_j=new_seq.index(origin_seq[j])
			new_matrix[new_i,new_j]=matrix[i,j]
	return new_matrix#[0:dimentions,0:dimentions]
def loaddata():
#################load metal info
########### using protein name to classify into three groups, train:validation:test=4:1:1
	classify_dict=dict()
	metal_label=1
	positive_list=open("data/negative_list",'r').readlines()
	positive_dict=dict()
	for pair_file in positive_list:
		pair_file=pair_file.strip()
		gene=pair_file.split('_')[0]
		positive_dict.setdefault(gene,[])
		positive_dict[gene].append(pair_file)
	for (i,gene) in enumerate(positive_dict.keys()):
		if i%3!=0:
			classify_dict[gene]="train"
		else:
			if i%6!=0:
				classify_dict[gene]="validation"
			else:
				classify_dict[gene]="test"
######## put positive samples and negative samples into train validation and test set
	train=[]
	validation=[]
	test=[]
	metal_label=1
	positive_list=open("data/positive_list",'r').readlines()
	for (i,line) in enumerate(positive_list): 
		gene_file=line.strip()
		try:
			f1=open("data/positive/"+gene_file,'rb')
		except:
			print(gene_file)
			raise
		data1=pickle.load(f1)
		data1=transform_matrix(data1)
		gene=gene_file.split("_")[0]
		if gene in classify_dict.keys():
			if classify_dict[gene]=="train":
				train.append((data1,metal_label))
			elif classify_dict[gene]=="validation":
				validation.append((data1,metal_label))
			elif classify_dict[gene]=="test":
				test.append((data1,metal_label))
	positive_samples=len(positive_list)
	negative_list=open("data/negative_list",'r').readlines()
	random.shuffle(negative_list)
	for (i,line) in enumerate(negative_list):
		if i < positive_samples:
			gene_file=line.strip()
			f1=open("data/negative/"+gene_file,'rb')
			data1=pickle.load(f1)
			data1=transform_matrix(data1)
			gene=gene_file.split("_")[0]
			if gene in classify_dict.keys():
				if classify_dict[gene]=="train":
					train.append((data1,0))
				elif classify_dict[gene]=="validation":
					validation.append((data1,0))
				elif classify_dict[gene]=="test":
					test.append((data1,0))
	random.shuffle(train)
	random.shuffle(validation)
	random.shuffle(test)
	x_train=[]
	y_train=[]
	x_validation=[]
	y_validation=[]
	x_test=[]
	y_test=[]
	for line in train:
		x_train.append(line[0])
		y_train.append(line[1])
	for line in validation:
		x_validation.append(line[0])
		y_validation.append(line[1])
	for line in test:
		x_test.append(line[0])
		y_test.append(line[1])
	x_train=np.array(x_train)
	y_train=np.array(y_train)
	x_validation=np.array(x_validation)
	y_validation=np.array(y_validation)
	x_test=np.array(x_test)
	y_test=np.array(y_test)
	return (x_train,y_train,x_validation,y_validation,x_test,y_test)	
