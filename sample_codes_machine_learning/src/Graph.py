# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 08:02:21 2018

@author: xiang
"""

import tensorflow as tf

graph = tf.Graph()

returned = graph is tf.get_default_graph()

print("Test graph1: ", returned)

x1 = tf.Variable(3)

print("Test graph2: ", x1.graph is tf.get_default_graph())



