# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 18:05:30 2018

@author: xiang
"""

print(__doc__)

from IPython import get_ipython
get_ipython().magic('reset -sf')
get_ipython().magic('clear')

from sklearn.pipeline import Pipeline
from sklearn import neighbors, preprocessing
from sklearn import datasets
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier, export_graphviz
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score, ShuffleSplit
from time import time

from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.datasets import load_digits
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.mixture import GaussianMixture
from sklearn import mixture
import itertools
from scipy import linalg
import matplotlib as mpl
from sklearn.decomposition import FastICA
from sklearn import random_projection
import TensorFlow_DNN
import tensorflow as tf
from sklearn.decomposition import PCA
from sklearn.manifold import LocallyLinearEmbedding
from sklearn.manifold import MDS

filename = 'Letter.csv'
data = np.genfromtxt(filename, delimiter=',', dtype='str')
sample_data = data[:,1:].astype(np.float32)
num_of_features = sample_data.shape[1];
sample_data = scale(sample_data)

sample_target = data[:,0]
#Transform the labels
le = preprocessing.LabelEncoder()
le.fit(sample_target.ravel())
sample_target = le.transform(sample_target.ravel()).reshape(sample_target.shape) #.astype(np.float32)

labels = np.unique(sample_target)
num_labels = labels.size

for index, item in enumerate(["LLE"]): #enumerate(["PCA", "ICA", "Random", "LLE"]):
    
    reduced_num_components = 13
    
    if item == "PCA":
        estimator = PCA(n_components = reduced_num_components)
        sample_data_transformed = estimator.fit_transform(sample_data)
        estimator.explained_variance_ratio_
        
        dims = []
        variance_ratios = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]
        for ratio in variance_ratios:
            pca2 = PCA(n_components = ratio)
            sample_data_transformed2 = pca2.fit_transform(sample_data)
            dims.append(sample_data_transformed2.shape[1])
            
        plt.figure()
        plt.plot(dims, variance_ratios, label='asdfds')
        plt.legend(["Ratio-dims"])
        plt.xlabel("Dimensions")
        plt.ylabel("Explained variance")
            
    elif item == "ICA":        
        estimator = FastICA(n_components=reduced_num_components,random_state=0)
        sample_data_transformed = estimator.fit_transform(sample_data)
    elif item == "Random":
        estimator = random_projection.GaussianRandomProjection(n_components = reduced_num_components)
        sample_data_transformed = estimator.fit_transform(sample_data)
    elif item == "LLE":
        #estimator = LocallyLinearEmbedding(n_components=reduced_num_components, n_neighbors=10)
        #sample_data_transformed = estimator.fit_transform(sample_data)  
        estimator = MDS(n_components=2)
        sample_data_transformed = estimator.fit_transform(sample_data)

    max_score = 0.0
    
    #**********************Train the deep neural network!!!**********************
    for activation in [tf.nn.relu, tf.nn.tanh, tf.nn.sigmoid, tf.nn.leaky_relu, tf.nn.elu]:
        for learning_rate in [0.01, 0.05, 0.1]:
            for initiation in [tf.contrib.layers.variance_scaling_initializer(), tf.contrib.layers.xavier_initializer()]:
                
                return_dum = TensorFlow_DNN.TensorFlow_DNN(sample_data = sample_data_transformed, sample_target = sample_target, initiation = initiation,  
                                                   Percentage_testing = 0.2, random_key = 31, learning_rate = learning_rate, n_epochs = 5, batch_size = 1000,
                                                   activation = activation)
                
                validation_score = return_dum[2]
    
                print("Returned val accuracy: ", initiation, activation, learning_rate, validation_score)
                
                if ( validation_score > max_score ):
                    max_score = validation_score
                    activation_chosen = activation
                    learning_rate_chosen = learning_rate
                    initiation_chosen = initiation
                    
    print("The ultimate one: %s, %f, %s, %s, %s" %(item, max_score, activation_chosen, learning_rate_chosen, initiation_chosen))
                
    #**********************Make predictions!!!**********************
    X = tf.placeholder(tf.float32, shape=(None, reduced_num_components), name = "X")
    y = tf.placeholder(tf.int64, shape=(None), name = "y")
    n_hidden1 = 64
    n_hidden2 = 64
    n_outputs = 26
    
    tf.get_variable_scope().reuse_variables()
    
    with tf.name_scope("dnn"):
        hidden1 = tf.layers.dense(X, n_hidden1, name="hidden1", activation=tf.nn.relu)
        hidden2 = tf.layers.dense(hidden1, n_hidden2, name="hidden2", activation=tf.nn.relu)
        logits = tf.layers.dense(hidden2, n_outputs, name="outputs")
    
    with tf.Session() as sess:
        saver = tf.train.Saver()
        saver.restore(sess, return_dum[1])
        X_new_scaled = sample_data_transformed[:20]
        X_new_scaled = return_dum[0].transform(X_new_scaled)
        Z = logits.eval(feed_dict={X:X_new_scaled})
        y_pred = np.argmax(Z, axis = 1)