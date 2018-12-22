# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 07:50:38 2018

@author: xiang
"""

import tensorflow as tf

x = tf.Variable(2, name = 'x')

f = x

saver = tf.train.Saver()

with tf.Session() as sess:
    saver.restore(sess, "./temp/temp.ckpt")
    print("Evaluate f the third time: ", f.eval())

    