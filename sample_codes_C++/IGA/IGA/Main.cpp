#include <map>
#include <vector>
#include <iostream>
#include <string>
#include "Mesh.h"
#include "IGA_Object.h"

int main(int argc, char * argv[])
{
	IGA_Object Instance_IGA;

	Instance_IGA.ReadMesh("Mesh_2D.inp");
	
	Instance_IGA.ProcessMesh();

	Instance_IGA.ObtainKnotInformation();


	system("pause");
	return 1;
}