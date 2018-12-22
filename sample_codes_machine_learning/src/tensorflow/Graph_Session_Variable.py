#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 19:31:06 2018

@author: xiang
"""
#
import tensorflow as tf

x = tf.Variable(3, name = "x")
y = tf.Variable(4, name = "y")
f = x*x*y + y + 2

sess = tf.Session()

sess.run(x.initializer)
sess.run(y.initializer)
result = sess.run(f)
print(result)
sess.close()

with tf.Session() as sess:
    x.initializer.run()
    y.initializer.run()
    result = f.eval()
    
init = tf.global_variables_initializer()

with tf.Session() as sess:
    init.run()
    result2 = f.eval()
    
#Interactive Session
sess2 = tf.InteractiveSession()

init.run()

result3 = f.eval()

print(result3)

sess2.close()

x1 = tf.Variable(1)

bool_return = x1.graph is tf.get_default_graph()

tf.reset_default_graph()

bool_return2 = x1.graph is tf.get_default_graph()

w = tf.constant(3)

x = w + 2

y = x + 5

z = x + 3

with tf.Session() as sess:
    print(y.eval())
    print(z.eval())

#evaluate y and z simultaneously to avoid duplicate variable evaluations
with tf.Session() as sess:
    y_val, z_val = sess.run([y,z])
    print(y_val)
    print(z_val)    
    
    
    
    
    
    
    
