# -*- coding: utf-8 -*-
"""
Created on Sat Dec 22 08:10:21 2018

@author: xiang
"""

import tensorflow as tf

A = tf.placeholder(tf.float32, shape=(None, 3))

B = A + 3

with tf.Session() as sess:
    B_val_1 = B.eval(feed_dict={A:[[1,2,3]]})
    B_val_2 = B.eval(feed_dict={A:[[4,5,6],[7,8,9]]})

print("B_val_1: \n", B_val_1)

print("B_val_1: \n", B_val_2)