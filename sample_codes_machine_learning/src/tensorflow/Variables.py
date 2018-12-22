# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 07:50:38 2018

@author: xiang
"""

import tensorflow as tf

x = tf.Variable(3, name = 'x')

y = tf.Variable(4, name = 'y')

f = x*x*y + y + 2

#Method 1 of using tf.Session()
sess = tf.Session()

sess.run(x.initializer)

sess.run(y.initializer)

result = sess.run(f)

print(result)

sess.close()

#Method 2 of using tf.Session()
with tf.Session() as sess:
    x.initializer.run()
    y.initializer.run()
    print("Evaluate f again: ", f.eval())
    
#Method 3 of using global_variables_initializer()
    
init = tf.global_variables_initializer()

with tf.Session() as sess:
    init.run()
    print("Evaluate f the third time: ", f.eval())
    