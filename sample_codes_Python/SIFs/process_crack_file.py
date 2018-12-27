def process_crack_file():

	file = "Cracks_Ori.inp"

	op_file = open(file, "r")

	string = op_file.read()

	strings = string.split("\n")

	dum = 0
	dum1 = 0
	for index, i in enumerate(strings):
		if i.upper().find("*NODE") >= 0 and i.upper().find("**NODE") == -1 and dum == 0:
			strings[index] = "*NODE, Crack_1"
			dum = dum + 1
		elif i.upper().find("*ELEMENT") >= 0 and i.upper().find("**ELEMENT") == -1 and \
			 i.upper().find("TYPE=S3") >=0 and dum1 == 0:
			
			strings[index] = "*Element, type=S3, Crack_1"
			dum1 = dum1 + 1
			
		elif i.upper().find("*INSTANCE") >= 0 and i.upper().find("**INSTANCE") == -1:
			strings[index] = strings[index] + ", Crack_1"
		
		#print("Test i: ", i)

	op_file.close()

	file_new = "Cracks.inp"
	wr_file = open(file_new, "w")

	for index, i in enumerate(strings):
		wr_file.write(i+"\n")

	wr_file.close()


