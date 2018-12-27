# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 10:03:21 2018

@author: xiang
"""

def gen_model_info(Crack_type, Crack_Center, Crack_Radius, Axis, angle):
    file_new = "Crack_info.inp"
    wr_file = open(file_new, "wb")
    wr_file.write("Crack_Type: " + Crack_type + "\n")
    wr_file.write("Crack_Center: " + str(Crack_Center[0]) + ", " + str(Crack_Center[1]) + ", " + str(Crack_Center[2]) + "\n")
    wr_file.write("Crack_Radius: " + str(Crack_Radius) + "\n")
    wr_file.write("Rotation_Axis: " + str(Axis) + "\n")
    wr_file.write("Rotation_Angle: " + str(angle) + "\n")