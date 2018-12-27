import abaqus_gen
import Input_Written
import abaqus_crack
import subprocess
import gen_model_info
from subprocess import check_output
import time
import process_crack_file
from shutil import copyfile
import os
import datetime

#seed_sizes = [0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02]
#angles = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]

seed_sizes = [0.1] #[0.1, 0.05]
angles = [45.0] #[45.0, 60.0]
Axis = "+X"

Crack_type = "Penny"
Crack_Center = [1.0, 1.0, 1.0]

for i0 in seed_sizes:
    
    seed_size = i0
    
    abaqus_gen.abaqus_gen(seed_size = seed_size)
    
    for i in angles:
                
    	theta = i
    
    	abaqus_crack.generate_crack(theta = theta)

        #Write the model information to a file
        gen_model_info.gen_model_info(Crack_type = Crack_type, Crack_Center = Crack_Center, Crack_Radius = 0.1, 
                                      Axis = Axis, angle = i)
        
    	process_crack_file.process_crack_file()
    
    	#time.sleep(10)
    
    	a = check_output(["abaqus", "-j", "Job-1_new", "int", "standard_parallel=all", "cpus=3"], shell=True).decode()
    	
    	try:
    		os.mkdir("Pre_Results")
    	except OSError:
    		print("The folder has been created.\n")
    	
    	src = "./Outputs/Job-1_new_Crack_1.FatigueOut"
    	
    	now = datetime.datetime.now()
    	
    	new_file = "Job-1_new_Crack_1_" + str(now.year) + "_" + str(now.month) \
    			+ "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) \
    			+ "_" + str(now.second) + "_" + str(now.microsecond) + ".FatigueOut"
    	
    	dst = "./Pre_Results/" + new_file
    	
    	copyfile(src, dst)
        
        src2 = "./Outputs/Job-1_new_Crack_1.Quality"
        
    	new_file2 = "Job-1_new_Crack_1_" + str(now.year) + "_" + str(now.month) \
    			+ "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) \
    			+ "_" + str(now.second) + "_" + str(now.microsecond) + ".Quality"
    	
    	dst2 = "./Pre_Results/" + new_file2
    	
    	copyfile(src2, dst2)
        
	 
#a = check_output(["chdir"], shell = True)

#import sys
#import subprocess

#theproc = subprocess.Popen([sys.executable, "abaqus -j Job-1_new int standard_parallel=all cpus=3"])
#theproc.communicate()

#subprocess.check_call([sys.executable, "abaqus -j Job-1_new int standard_parallel=all cpus=3"])

#subprocess.call("abaqus -j Job-1_new int standard_parallel=all cpus=3") #, shell=True)

#import os

#p = os.popen('abaqus -j Job-1_new int standard_parallel=all cpus=3')
