# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 07:50:38 2018

@author: xiang
"""

import tensorflow as tf

x = tf.Variable(3, name = 'x')

y = tf.Variable(4, name = 'y')

f = x*x*y + y + 2
    
init = tf.global_variables_initializer()

#Can save specific nodes, e.g. tf.train.Saver({"weights": theta})
saver = tf.train.Saver()

with tf.Session() as sess:
    init.run()
    print("Evaluate f the third time: ", f.eval()) 
    save_path = saver.save(sess, "./temp/temp.ckpt")

    