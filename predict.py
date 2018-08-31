#!/usr/bin/env python
from __future__ import print_function
#import urllib
import keras
from keras import backend as K
import os
from keras.models import load_model
import sys
import numpy as np
from loaddata import transform_matrix
###################################################
#the GREMLIN software has not been released yet, so please calculate the coevolution result on the online version
#website: gremlin.bakerlab.org/submit.php
#after finish this calculation, please use the ID number as the script input;
#for example, gremlin.bakerlab.org/submit.php&id=1535631343, the ID number is 1535631343
#run the script, python prepare.py 1535631343
input_ID=sys.argv[1]
###########################download the msa file and contact map file from website
msa_file="test_data/input_"+input_ID+".msa"
contact_file="test_data/input_"+input_ID+".con"
msa_website="http://gremlin.bakerlab.org/sub_fasta.php?id="+input_ID
#contact_website="http://gremlin.bakerlab.org/sub_txt.php?db=SUB&id="+input_ID
contact_website="http://gremlin.bakerlab.org/sub_txt.php?id="+input_ID
os.system("wget "+msa_website+" -O "+msa_file)
os.system("wget "+contact_website+" -O "+contact_file)
if "ERROR" in open(contact_file,'r').readlines():
	print("the GREMLIN Result is being analyzed, please wait~~~~~~~~~~~~~~~~~~~~~~~~~")
	sys.exit(1)
#filehandle=open(contact_file,'w')
#filehandle.write(newfile)
#filehandle.close()
############################get the coevolution pair from .con file
CHED_pair=['HH','HC','CH','CC','HD','DD','DH','CD','DC','EE','EC','CE','DE','ED','HE','EH']
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
img_rows, img_cols = 21, 21
model = load_model('models/contact_model.h5')
output_file="test_data/output_"+input_ID+".cnn"
output_handle=open(output_file,'w')
for pair in coevolution_dict.keys():
	data1=coevolution_dict[pair]/seq_num
	data1=transform_matrix(data1)
	x_prediction = data1.reshape(1, img_rows, img_cols, 1)
	classes =model.predict(x_prediction,batch_size=1)
	i_id=pair[0]
	j_id=pair[1]
	i_tag=i_id+"_"+origin_seq[int(i_id)-1]
	j_tag=j_id+"_"+origin_seq[int(j_id)-1]
	for i,num in enumerate(classes.argmax(axis=1)):
		if num==1:
			print(i_tag,j_tag,num,classes[i][num])
			output_handle.write(i_tag+"\t"+j_tag+"\t"+str(classes[i][num])+"\n")
output_handle.close()
