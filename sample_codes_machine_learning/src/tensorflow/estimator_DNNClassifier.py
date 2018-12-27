# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:09:07 2018

@author: xiang
"""

import tensorflow as tf
from tensorflow import keras
#from sklearn.preprocessing import StandardScaler
import numpy as np
#from sklearn.metrics import accuracy_score

print(tf.__version__)

#Load the data
fashion_mnist = keras.datasets.fashion_mnist
(train_x, train_y), (test_x, test_y) = fashion_mnist.load_data()

#Scale the data
#scaler = StandardScaler()
#scaler.fit(train_x)
#scaled_housing_data_plus_bias = scaler.transform(train_x)

#Cast the data set from float64 to float32
#train_x = np.float32(train_x) #tf.cast(train_x, tf.float32)

#Get the number of features
#feature_columns = tf.contrib.learn.infer_real_valued_columns_from_input(train_x) #list(range(train_x.shape[1])) #range(train_x.get_shape().as_list()[1])

#Get the number of classes
n_classes = 10 #train_y.shape[0]

#Build the DNNClassifier
# Specify feature
feature_columns = [tf.feature_column.numeric_column("x", shape=[28, 28])]

classifier = tf.estimator.DNNClassifier(
    feature_columns=feature_columns,
    hidden_units=[256, 32],
    optimizer=tf.train.AdamOptimizer(1e-4),
    n_classes=10,
    dropout=0.1,
    model_dir="./tmp/mnist_model"
)

# Define the training inputs
train_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"x": train_x},
    y=np.int32(train_y),
    num_epochs=None,
    batch_size=50,
    shuffle=True
)

classifier.train(input_fn=train_input_fn, steps=1000)

# Define the test inputs
test_input_fn = tf.estimator.inputs.numpy_input_fn(
    x={"x": test_x},
    y=np.int32(test_y),
    num_epochs=1,
    shuffle=False
)

# Evaluate accuracy
accuracy_score = classifier.evaluate(input_fn=test_input_fn)["accuracy"]
print("\nTest Accuracy: {0:f}%\n".format(accuracy_score*100))