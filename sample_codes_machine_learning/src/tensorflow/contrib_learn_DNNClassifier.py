# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:09:07 2018

@author: xiang
"""

import tensorflow as tf
from tensorflow import keras
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.metrics import accuracy_score

print(tf.__version__)

#Load the data
boston_housing = keras.datasets.boston_housing
(train_x, train_y), (test_x, test_y) = boston_housing.load_data()

#Scale the data
scaler = StandardScaler()
scaler.fit(train_x)
scaled_housing_data_plus_bias = scaler.transform(train_x)

#Cast the data set from float64 to float32
train_x = np.float32(train_x) #tf.cast(train_x, tf.float32)

#Get the number of features
feature_columns = tf.contrib.learn.infer_real_valued_columns_from_input(train_x) #list(range(train_x.shape[1])) #range(train_x.get_shape().as_list()[1])

#Get the number of classes
n_classes = train_y.shape[0]

#Build the DNNClassifier
dnn_clf = tf.contrib.learn.DNNClassifier(hidden_units=[300,100], n_classes=n_classes, feature_columns = feature_columns)
dnn_clf.fit(x = train_x, y=train_y, batch_size = 50, steps = 40000)

y_pred = list(dnn_clf.predict(test_x))

accuracy_score(test_y, y_pred)