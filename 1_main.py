#!/usr/bin/env python
from __future__ import print_function
import sys
import argparse
import keras
import functools
from itertools import product
from keras import regularizers
from keras.utils import plot_model
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv3D, MaxPooling2D,Conv2D
from keras import backend as K
import os, random
import numpy as np
import h5py
from sklearn.datasets import make_classification
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix,average_precision_score
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
import math
import pickle
from loaddata import *
def w_categorical_crossentropy(y_true, y_pred, weights):
	nb_cl = len(weights)
	final_mask = K.zeros_like(y_pred[:, 0])
	y_pred_max = K.max(y_pred, axis=1)
	y_pred_max = K.reshape(y_pred_max, (K.shape(y_pred)[0], 1))
	y_pred_max_mat = K.cast(K.equal(y_pred, y_pred_max), K.floatx())
	for c_p, c_t in product(range(nb_cl), range(nb_cl)):
		final_mask += (weights[c_t, c_p] * y_pred_max_mat[:, c_p] * y_true[:, c_t])
	return K.categorical_crossentropy(y_pred, y_true) * final_mask
result=0
for a in range(0,10):
	num_classes=2
	(x_train,y_train,x_validation,y_validation,x_test,y_test)=loaddata()
	img_rows, img_cols = 21, 21
	############start nn
	batch_size = 32
	epochs = 50
	class_factor=float((len(y_test)-np.sum(y_test)))/float(np.sum(y_test))
	print("the class_factor is %0.1f" % (class_factor))
	x_train = x_train.reshape(x_train.shape[0],img_rows, img_cols, 1)
	x_validation = x_validation.reshape(x_validation.shape[0],img_rows, img_cols, 1)
	x_test = x_test.reshape(x_test.shape[0], img_rows, img_cols, 1)
	input_shape = (img_rows, img_cols, 1)
	print('x_train shape:', x_train.shape)
	print(x_train.shape[0], 'train samples')
	print(x_test.shape[0], 'test samples')
	print(x_validation.shape[0], 'validation samples')
	print(x_train.shape[0],x_train.shape[1],x_train.shape[2])
	# convert class vectors to binary class matrices
	y_train = keras.utils.to_categorical(y_train, num_classes)
	y_validation = keras.utils.to_categorical(y_validation, num_classes)
	y_test = keras.utils.to_categorical(y_test, num_classes)
	model = Sequential()
	model.add(Conv2D(16, kernel_size=(3, 3),
	                 activation='relu',
					 padding='same',
	                 input_shape=input_shape))
	model.add(Conv2D(16, (3, 3), activation='relu'))
	model.add(MaxPooling2D(pool_size=(2, 2)))
	model.add(Dropout(0.25))
	model.add(Flatten())
	model.add(Dense(16, activation='relu',
				activity_regularizer=regularizers.l1(0.0001),))
	model.add(Dropout(0.25))
	model.add(Dense(num_classes, activation='softmax'))
	model.summary()
	model.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=keras.optimizers.Adadelta(),
	              metrics=['accuracy'])
	model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
	          verbose=1, validation_data=(x_validation, y_validation))
	class_weight= {	0:1.,
					1:class_factor,}
	tbCAllBack=keras.callbacks.TensorBoard(log_dir='./logs/run_d',histogram_freq=0,write_graph=True,write_images=True)
	score = model.evaluate(x_validation, y_validation, verbose=0)
	print('validation loss:', score[0])
	print('validation accuracy:', score[1])
	print('validation ?:', score)
	accuracy=[]
	F1_score=[]
	Recall=[]
	Precision=[]
	weighted_prediction=model.predict(x_test)
	for i in range(0,num_classes):
		l=y_test[:,i]
		p=np.array([round(x[i]) for x in weighted_prediction])
		accuracy.append(accuracy_score(l,p))
		F1_score.append(f1_score(l,p))
		Recall.append(recall_score(l,p))
		Precision.append(precision_score(l,p))
	print('Accuracy','F1 score','Recall','Precision')
	for i in range(1,num_classes):
		print("flag\t%0.3f\t%0.3f\t%0.3f\t%0.3f" % (accuracy[i],F1_score[i],Recall[i],Precision[i]))
if result <F1_score[i]:
	result = F1_score[i]
	model.save('models/contact_model.h5')
