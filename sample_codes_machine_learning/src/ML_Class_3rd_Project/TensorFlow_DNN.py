print(__doc__)
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 15:46:34 2018

@author: xiang
"""

import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
#from tensorflow.keras.callbacks import TensorBoard
from datetime import datetime

#np.set_printoptions(threshold=np.nan) #print the full arrays

def TensorFlow_DNN(sample_data, sample_target, initiation, Percentage_testing = 0.2, random_key = 31, learning_rate = 0.1, n_epochs = 100, batch_size = 50, 
                   activation = tf.nn.relu):
    
    #use train_test_split in sklearn to split the train and test data
    X_train, X_test, y_train, y_test = train_test_split(sample_data, sample_target, test_size= Percentage_testing, random_state=random_key)
    
    #Scale the training data
    scaler = StandardScaler()
    scaler.fit(X_train)
    X_train = scaler.transform(X_train)
    #Scale the testing data
    X_test = scaler.transform(X_test)
    
    tf.reset_default_graph()

    sample_size = X_train.shape[0]
    n_inputs = X_train.shape[1]
    n_hidden1 = 64
    n_hidden2 = 64
    n_outputs = 26

    X = tf.placeholder(tf.float32, shape=(None, n_inputs), name = "X")
    y = tf.placeholder(tf.int64, shape=(None), name = "y")

    with tf.name_scope("dnn"):
        hidden1 = tf.layers.dense(X, n_hidden1, name="hidden1", activation=activation, kernel_initializer=initiation)
        hidden2 = tf.layers.dense(hidden1, n_hidden2, name="hidden2", activation=activation)
        logits = tf.layers.dense(hidden2, n_outputs, name="outputs")

    with tf.name_scope("loss"):
        #For classification:
        xentropy = tf.nn.sparse_softmax_cross_entropy_with_logits(labels=y, logits=logits)
        loss = tf.reduce_mean(xentropy, name = "loss")
        #For regression:
        #loss = tf.losses.mean_squared_error(y, y_predict)

    with tf.name_scope("train"):
        optimizer = tf.train.GradientDescentOptimizer(learning_rate)
        training_op = optimizer.minimize(loss)
    
    with tf.name_scope("eval"):
        correct = tf.nn.in_top_k(logits, y, 1)
        accuracy = tf.reduce_mean(tf.cast(correct, tf.float32))
   
    init = tf.global_variables_initializer()
    saver = tf.train.Saver()

    with tf.Session() as sess:
        init.run()
        for epoch in range(n_epochs):
            for iteration in range(sample_size // batch_size):
                X_batch = X_train[iteration * batch_size : (iteration + 1) * batch_size]
                y_batch = y_train[iteration * batch_size : (iteration + 1) * batch_size].ravel()
                sess.run(training_op, feed_dict={X:X_batch, y:y_batch})
            
            acc_train = accuracy.eval(feed_dict={X:X_train, y: y_train})
            #dum = logits.eval(feed_dict={X:X_train, y: y_train})
            #print("Test logits: ", dum, dum.shape)
            acc_val = accuracy.eval(feed_dict={X:X_test, y: y_test})
            
            #print(epoch, "Train accuracy: ", acc_train, "Val accuracy: ", acc_val)
        
        save_path = saver.save(sess, "./my_model_final.ckpt")
        
        return scaler, save_path, acc_val










    