#include "IGA_Object.h"

IGA_Object::~IGA_Object()
{
	delete Knots;
}

int IGA_Object::ReadMesh(const std::string& FileName)
{
	std::string tmp2;
	std::vector< int > row;
	std::vector< double > row2;
	std::vector< std::string > row_string;
	std::vector< std::string > Terms;

	std::ifstream infile(FileName);

	if (infile.fail())
	{
		std::cout << "Error in opening the IGA input file.\n" << std::endl;
		system("pause");
		exit(-3);
	}

	std::string tmp;
	while (getline(infile, tmp)){
		if (tmp[0] == '*' && tmp[1] != '*'){
			Terms.clear();
			DivideTerms(tmp, &Terms);

			if (MatchString(Terms[0].c_str(), "*NODE") && !(Terms[0][0] == '*' && Terms[0][1] == '*') && !MatchString(Terms[0].c_str(), "*NODEOUTPUT")){
				while (true){
					row2.clear();
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if (tmp2[0] == '*' || tmp2.find_first_not_of(' ') == std::string::npos){
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else{
						ReadLine(tmp2, row2);
						std::vector<double> row_Temp(row2.begin() + 1, row2.end());
						//Nodes[(int)row2[0]] = row_Temp;
						Original_Mesh.Assign_Node((int)row2[0], row_Temp);
					}
				}
			}
			else if (MatchString(Terms[0].c_str(), "*ELEMENT")){
				while (true){
					row.clear();
					std::string tmp2;
					std::streamoff previous = infile.tellg();
					getline(infile, tmp2);
					if (tmp2[0] == '*' || tmp2.find_first_not_of(' ') == std::string::npos){
						infile.seekg(previous, std::ios::beg);
						break;
					}
					else{
						ReadLine(tmp2, row);
						std::vector<int> row_Temp(row.begin() + 1, row.end());

						if (row_Temp.size() != 4 && row_Temp.size() != 5) {
							std::cout << "Error: The element " << (int) row[0] <<" has "<< row_Temp.size() <<" edges which can cause trouble.\n";
							system("pause");
							exit(-4);
						}

						Original_Mesh.Assign_Element((int)row[0], row_Temp);
					}
				}
			}
		} //if (tmp[0] == '*' && tmp[1] != '*'){
	}

	infile.close();

	return 1;
}

int IGA_Object::ProcessMesh()
{
	Original_Mesh.IdentifyElementEdges();
	Original_Mesh.ObtainNodesInverseConnectivities();
	Original_Mesh.ObtainNodesAdjacentNodes();

	Original_Mesh.ObtainEdgeIndexes();
	Original_Mesh.ObtainElementsEdgeFlags();
	Original_Mesh.Obtain_TJunction_ExtraOrdinaryNodes();

	std::set <int> Elements;
	for (size_t m = 0; m < 10; m++){
		Original_Mesh.Shared_Elements_of_Edge(m+1, Elements);
	}

	int Current_Edge_Index;
	int Next_Edge_Index;
	Original_Mesh.ObtainNextEdge_of_Node_Edge(13, 5, Current_Edge_Index, Next_Edge_Index);

	std::cout << Current_Edge_Index << ", " << Next_Edge_Index << std::endl;

	//Original_Mesh.PrintNodesInformation();
	//Original_Mesh.PrintElementsInformation();

	Original_Mesh.Paraview_Visualizations();

	//////////int k = 3;
	//////////std::vector <double> Knot_Vector{ 0,0,0,1,1,3,3,3 };
	//////////std::vector <std::vector <std::vector <double> > > Coefficients;
	//////////std::vector <int> Initial_Intervals;

	//////////Cox_de_Boor_Recursion(k, Knot_Vector, Initial_Intervals, Coefficients);

	//////////PrintSplines(2, Knot_Vector, Initial_Intervals, Coefficients);

	//////////std::string FileName = "BSplines.vtk";
	//////////Paraview_Visualization(FileName, 1000, Knot_Vector, Initial_Intervals, Coefficients);

	return 1;
}

int IGA_Object::ObtainKnotInformation()
{
	Knots = new Class_Knots(Original_Mesh);
	Knots -> ObtainKnotIntervals_on_edges();

	return 1;
}