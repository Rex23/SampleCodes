file = "Job-1.inp"

op_file = open(file, "r")

string = op_file.read()

strings = string.split("\n")

for index, i in enumerate(strings):
	if i.upper().find("*ELEMENT") >= 0 and i.upper().find("**ELEMENT") == -1 and i.upper().find("C3D8RH") >= 0:
		pos = i.upper().find("TYPE")
		strings[index] = i.replace("C3D8RH", "U1381, ELSET = XFEM_SET")
		strings.insert(index, "*userelement, nodes=8, type=U1381, properties=8, coordinates=3\n1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15")
	elif i.upper().find("*MATERIAL") >= 0 and i.upper().find("**MATERIAL") == -1 and i.upper().find("MATERIAL-XFEM") >= 0:
		print(len(strings))
		xfem_material_str = strings[index + 2]
		#print("Test idfs: ", xfem_material_str)
	
	#*Material, name=Material-XFEM

op_file.close()

file_new = "Job-1_new.inp"
wr_file = open(file_new, "wb")

dum = 0
dum2 = -1
for index, i in enumerate(strings):
	if i.upper().find("*END PART") >= 0 and i.upper().find("**END PART") == -1 and dum == 0:
		wr_file.write("*uelproperty, elset=XFEM_SET\n")
		wr_file.write(xfem_material_str+"\n")
		dum = dum + 1
		wr_file.write(i+"\n")
	elif i.upper().find("*SOLID SECTION") >= 0 and i.upper().find("**SOLID SECTION") == -1 and i.upper().find("MATERIAL-XFEM") >= 0:
		wr_file.write("**" + i +"\n")
		dum2 = index
	elif index == dum2 + 1:
		wr_file.write("**" + i + "\n")
	else:
		wr_file.write(i+"\n")
		
wr_file.close()