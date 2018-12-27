#include "Mesh.h"

int Class_Mesh::New_Edges_and_Nodes_Based_on_T_Junction(const int& T_Junction_Node_Index, const std::map <int, std::set <int> >& nodes_inverse_connectivities,
	const std::map <int, std::vector <int> >& elements_map, std::map <int, int> Nodes_Flags, std::map <int, int>& TJunction_Nodes_Ele)
{
	//for (size_t m = 0; m < TJunction_Nodes.size(); m++) {
	int Node_Index = T_Junction_Node_Index; // TJunction_Nodes[m];
	std::set <int> Shared_Elements = nodes_inverse_connectivities.find(Node_Index)->second;

	for (std::set <int>::const_iterator iter1 = Shared_Elements.begin(); iter1 != Shared_Elements.end(); iter1++) {

		int Ele_Index = *iter1;

		std::vector <int> Connectivities = elements_map.find(Ele_Index)->second;

		if (Connectivities.size() != 4) {
			TJunction_Nodes_Ele[Node_Index] = Ele_Index;
			break;
		}
	}

	//}

	return 1;
}

int Class_Mesh::New_Edges_and_Nodes_Based_on_Super_T_Junction(const int& Super_T_Junction_Node_Index, const std::map <int, std::set <int> >& nodes_inverse_connectivities,
	const std::map <int, std::vector <int> >& elements_map, std::map <int, int> Nodes_Flags, std::map <int, int>& TJunction_Nodes_Ele)
{
	//for (size_t m = 0; m < Super_TJunction_Nodes.size(); m++) {
	int Node_Index = Super_T_Junction_Node_Index; // Super_TJunction_Nodes[m];
	std::set <int> Shared_Elements = nodes_inverse_connectivities.find(Node_Index)->second;

	for (std::set <int>::const_iterator iter1 = Shared_Elements.begin(); iter1 != Shared_Elements.end(); iter1++) {

		int Ele_Index = *iter1;

		std::vector <int> Connectivities = elements_map.find(Ele_Index)->second;

		if (Connectivities.size() != 4) {

			int Temp = 0;
			for (size_t mm = 0; mm < Connectivities.size(); mm++) {
					
				int Node_Temp = Connectivities[mm];

				int Type = Nodes_Flags.find(Node_Temp)->second;

				if (Type == 1 || Type == 5) {
					Temp++;
				}
			}

			if (Temp == 1) {
				TJunction_Nodes_Ele[Node_Index] = Ele_Index;
				break;
			}

		}
	}

	//}

	return 1;
}

void Class_Mesh::Paraview_Visualizations()
{
	Paraview_Mesh_Visualization("Original_Mesh.vtk", nodes_map, elements_map, Nodes_Flags);
	Paraview_Mesh_Visualization_Edges("Original_Mesh_Edges.vtk", nodes_map, edges_nodes);
}

int Class_Mesh::Paraview_Mesh_Visualization(const std::string& FileName, const std::map <int, std::vector <double> >& nodes_map, 
	                                        const std::map <int, std::vector <int> >& elements_map, const std::map <int, int>& Nodes_Flags)
{
	//It is assumed that the nodes and elements are continuously indexed

	int NumberofElements = elements_map.size();
	int NumberofNodes = nodes_map.size();

	std::ofstream FileOutput;
	FileOutput.open(FileName);

	if (FileOutput.fail() != true) {
		FileOutput << "# vtk DataFile Version 4.0" << "\r\n";
		FileOutput << "XFEM Outputs" << "\r\n";
		FileOutput << "ASCII" << std::endl;
		FileOutput << "DATASET UNSTRUCTURED_GRID" << std::endl;
		FileOutput << std::endl;
		FileOutput << "POINTS" << " " << NumberofNodes << " " << "float" << std::endl;
		FileOutput << std::fixed;
		FileOutput << std::setprecision(8);

		for (std::map <int, std::vector <double> >::const_iterator iter1 = nodes_map.begin(); iter1 != nodes_map.end(); iter1++){
			FileOutput << std::setprecision(15) << iter1 ->second[0] << " " << iter1->second[1] << " " << 0.0 << std::endl;
		}

		FileOutput << std::endl;

		//*****************************
		int Num_Ele_Data = 0;
		for (std::map <int, std::vector <int> >::const_iterator iter1 = elements_map.begin(); iter1 != elements_map.end(); iter1++) {
			Num_Ele_Data += (iter1->second).size() + 1;
		}
		//*****************************

		FileOutput << "CELLS" << " " << NumberofElements << " " << Num_Ele_Data << std::endl;

		for (std::map <int, std::vector <int> >::const_iterator iter1 = elements_map.begin(); iter1 != elements_map.end(); iter1++){
			FileOutput << (iter1->second).size();
			for (size_t m = 0; m < (iter1->second).size(); m++) {
				FileOutput << " " << (iter1->second)[m] - 1;
			}
			FileOutput << std::endl;
		}

		FileOutput << std::endl;
		FileOutput << "CELL_TYPES" << " " << NumberofElements << std::endl;
		for (int i = 0; i < NumberofElements; i++) {
			FileOutput << "7" << std::endl;
		}

		FileOutput << std::endl;
		FileOutput << "POINT_DATA " << Nodes_Flags.size() << std::endl;
		FileOutput << std::endl;
		FileOutput << "SCALARS " << "Node_Flags float 1"<< std::endl;
		FileOutput << "LOOKUP_TABLE default" << std::endl;
		for (std::map <int, int>::const_iterator iter1 = Nodes_Flags.begin(); iter1 != Nodes_Flags.end(); iter1++) {
			FileOutput << iter1->second << std::endl;
		}
	}
	else {
		//fprintf(stdout, "Failed in opening the file for SurfaceMeshVisualization.\n");
		std::cout << "Failed in opening the file for " << FileName << std::endl;
		system("pause");
		exit(-54);
	}

	return 1;
}

int Class_Mesh::Paraview_Mesh_Visualization_Edges(const std::string& FileName, const std::map <int, std::vector <double> >& nodes_map,
	const std::map < int, std::set<int> >& edges_nodes)
{
	int NumberofElements = edges_nodes.size();
	int NumberofNodes = nodes_map.size();

	std::ofstream FileOutput;
	FileOutput.open(FileName);

	if (FileOutput.fail() != true) {
		FileOutput << "# vtk DataFile Version 4.0" << "\r\n";
		FileOutput << "XFEM Outputs" << "\r\n";
		FileOutput << "ASCII" << std::endl;
		FileOutput << "DATASET UNSTRUCTURED_GRID" << std::endl;
		FileOutput << std::endl;
		FileOutput << "POINTS" << " " << NumberofNodes << " " << "float" << std::endl;
		FileOutput << std::fixed;
		FileOutput << std::setprecision(8);

		for (std::map <int, std::vector <double> >::const_iterator iter1 = nodes_map.begin(); iter1 != nodes_map.end(); iter1++) {
			FileOutput << std::setprecision(15) << iter1->second[0] << " " << iter1->second[1] << " " << 0.0 << std::endl;
		}

		FileOutput << std::endl;

		FileOutput << "CELLS" << " " << NumberofElements << " " << NumberofElements *3 << std::endl;

		for (std::map < int, std::set<int> >::const_iterator iter1 = edges_nodes.begin(); iter1 != edges_nodes.end(); iter1++) {
			FileOutput << (iter1->second).size();
			for (std::set<int>::const_iterator iter2 = (iter1->second).begin(); iter2 != (iter1->second).end(); iter2++) {
				FileOutput << " " << *iter2 - 1;
			}
			FileOutput << std::endl;
		}

		FileOutput << std::endl;
		FileOutput << "CELL_TYPES" << " " << NumberofElements << std::endl;
		for (int i = 0; i < NumberofElements; i++) {
			FileOutput << "3" << std::endl;
		}

	}
	else {
		std::cout << "Failed in opening the file for " << FileName << std::endl;
		system("pause");
		exit(-55);
	}

	return 1;
}

int Class_Mesh::ObtainNextEdge_of_Node_Edge(const int& Node_Index, const int& Adjacent_Node_Index, int& Current_Edge_Index, int& Next_Edge_Index) //The adjacent node shall not be an extraordinary node
{
	Next_Edge_Index = 0;

	std::set <int> Edge_Set = { Node_Index, Adjacent_Node_Index };

	Current_Edge_Index = edges_map.find(Edge_Set)->second;

	std::set <int> Shared_Elements;
	Shared_Elements_of_Edge(Current_Edge_Index, Shared_Elements);

	int Node_Type = Nodes_Flags.find(Adjacent_Node_Index)->second;

	//Condition 1: Normal, Boundary Edge, Boundary Corner
	//Condition 2: T-Junction
	//Condition 3: Super T-Junction
	//Note: For conditions 2 and 3, the prerequisite is that the elements at most have 5 edges, i.e. at most one super T-Junction per element

	//*********************************************************************************
	if (Node_Type != 1 && Node_Type != 2 && Node_Type != 5){

		std::set <int> Adjacent_Nodes2 = nodes_adjacent_nodes.find(Adjacent_Node_Index)->second;

		for (std::set <int>::const_iterator iter1 = Adjacent_Nodes2.begin(); iter1 != Adjacent_Nodes2.end(); iter1++){

			int node_index_temp = *iter1;

			if (node_index_temp != Node_Index){

				std::set <int> node_set_temp = { node_index_temp, Adjacent_Node_Index };

				int edge_index_temp = edges_map.find(node_set_temp)->second;

				std::set <int> Shared_Elements2;
				Shared_Elements_of_Edge(edge_index_temp, Shared_Elements2);

				std::vector <int> Inter_Temp(max_intersection_num);

				std::vector <int>::iterator it = std::set_intersection(Shared_Elements.begin(), Shared_Elements.end(), Shared_Elements2.begin(), Shared_Elements2.end(), Inter_Temp.begin());

				Inter_Temp.resize(it - Inter_Temp.begin());

				if (Inter_Temp.size() == 0){
					Next_Edge_Index = edge_index_temp;
					break;
				}
			}
		}

		if (Next_Edge_Index == 0) {
			return 0;
		}

	}

	//*********************************************************************************
	//From T-Junction to its adjacent node
	else if (Node_Type == 1){ //T Junction
		
		int Flag = 0;

		std::set <int> Adjacent_Nodes2 = nodes_adjacent_nodes.find(Adjacent_Node_Index)->second;

		for (std::set <int>::const_iterator iter1 = Adjacent_Nodes2.begin(); iter1 != Adjacent_Nodes2.end(); iter1++){
			int node_index_temp = *iter1;

			if (node_index_temp != Node_Index){
				
				std::set <int> node_set_temp = { node_index_temp, Adjacent_Node_Index };

				int edge_index_temp = edges_map.find(node_set_temp)->second;

				std::set <int> Shared_Elements2;
				Shared_Elements_of_Edge(edge_index_temp, Shared_Elements2);

				for (std::set <int>::const_iterator iter1 = Shared_Elements2.begin(); iter1 != Shared_Elements2.end(); iter1++){

					int Ele_Index = *iter1;

					if ((elements_map.find(Ele_Index)->second).size() != 4){
						Next_Edge_Index = edge_index_temp;
						//std::cout << "NExt Edge " << Next_Edge_Index << std::endl;
						Flag = 1;
						goto flag1;
					}

				}

			}
		}

	flag1:;

		if (Flag == 0){
			std::cout << "Error: The next edge is not found.\n";
			system("pause");
			exit(-2);
		}

	}

	//*********************************************************************************
	//From super T-Junction to its adjacent node (cannot go to the adjacent T-Junction) 
	else if (Node_Type == 5){ //Super T Junction
		int Node_Type1 = Nodes_Flags.find(Node_Index)->second;

		if (Node_Type1 == 1){
			std::cout << "Did not find!" << std::endl;
		}
		else{

			int Flag = 0;

			std::set <int> Adjacent_Nodes2 = nodes_adjacent_nodes.find(Adjacent_Node_Index)->second;

			for (std::set <int>::const_iterator iter1 = Adjacent_Nodes2.begin(); iter1 != Adjacent_Nodes2.end(); iter1++){
				int node_index_temp = *iter1;

				int Node_Type_temp = Nodes_Flags.find(node_index_temp)->second;

				if (node_index_temp != Node_Index && Node_Type_temp != 1){

					std::set <int> node_set_temp = { node_index_temp, Adjacent_Node_Index };

					int edge_index_temp = edges_map.find(node_set_temp)->second;

					std::set <int> Shared_Elements2;
					Shared_Elements_of_Edge(edge_index_temp, Shared_Elements2);

					for (std::set <int>::const_iterator iter1 = Shared_Elements2.begin(); iter1 != Shared_Elements2.end(); iter1++){

						int Ele_Index = *iter1;

						if ((elements_map.find(Ele_Index)->second).size() != 4){
							Next_Edge_Index = edge_index_temp;
							//std::cout << "Next Edge " << Next_Edge_Index << std::endl;
							Flag = 1;
							goto flag2;
						}

					}

				}
			}

		flag2:;

			if (Flag == 0){
				return 0;
			}
		}

	}

	return 1;
}

int Class_Mesh::Shared_Elements_of_Edge(const int& Edge_Index, std::set<int>& Elements)
{
	std::set <int> shared_two_nodes = edges_nodes.find(Edge_Index) -> second;

	if (shared_two_nodes.size() != 2){
		std::cout << "Error in identifying the shared nodes of the edge!\n";
		exit(-1);
	}

	int first_node = *shared_two_nodes.begin();
	int second_node = *std::next(shared_two_nodes.begin(),1);
	
	std::set <int> set1 = nodes_inverse_connectivities.find(first_node)->second;
	std::set <int> set2 = nodes_inverse_connectivities.find(second_node)->second;

	std::vector <int> Intersection(max_intersection_num);
	std::vector <int>::iterator it;
	it = std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), Intersection.begin());
	Intersection.resize(it - Intersection.begin());

	Elements.clear();

	for (size_t m = 0; m < Intersection.size(); m++){
		Elements.insert(Intersection[m]);
	}

	//std::cout << first_node << ", " << second_node << ", " << Elements.size() << "\n";

	return 1;
}

int Class_Mesh::ObtainEdgeIndexes()
{
	//std::map < std::set<int>, int > edges_map;
	//std::map < int, int > nodes_adjacent_nodes_edge_index;
	//std::map <int, std::set <int> > nodes_adjacent_nodes;

	int Temp = 0;

	for (std::map <int, std::set <int> >::const_iterator iter1 = nodes_adjacent_nodes.begin();
		 iter1 != nodes_adjacent_nodes.end(); iter1++){
		int node_index = iter1->first;
		std::set<int> adjacent_nodes = iter1->second;
		for (std::set<int>::const_iterator iter2 = adjacent_nodes.begin(); iter2 != adjacent_nodes.end(); iter2++){
			int node_index2 = *iter2;

			std::set<int> set_temp = { node_index, node_index2 };

			if (edges_map.find(set_temp) == edges_map.end()){
				Temp++;
				edges_map[set_temp] = Temp;
				edges_nodes[Temp] = set_temp;

				//std::cout << "Test the edge_map: " << node_index << ", " << node_index2 << ", " << Temp << std::endl;
			}
		}
	}

	return 1;
}

int Class_Mesh::Assign_Node(const int& Node_Index, const std::vector <double>& Coords)
{
	nodes_map[Node_Index] = Coords;
	return 1;
}

int Class_Mesh::Assign_Element(const int& Element_Index, const std::vector <int>& Connectivities)
{
	elements_map[Element_Index] = Connectivities;
	return 1;
}

void Class_Mesh::PrintNodesInformation()
{
	PrintScreen_Map("Nodes", nodes_map, false);
	PrintScreen_MapSet("Nodes_Inverse_Connectivities", nodes_inverse_connectivities, false);
	PrintScreen_MapIntInt("Nodes_Flags", Nodes_Flags, false);
	PrintScreen_MapSet("Nodes_Adjacent_Nodes", nodes_adjacent_nodes, true);
}

void Class_Mesh::PrintElementsInformation()
{
	PrintScreen_Map("Elements", elements_map, false);
	PrintScreen_MapVector("Test Edges", elements_edges, false);
}

int Class_Mesh::IdentifyElementEdges()
{
	//Save the edges
	for (std::map <int, std::vector <int> >::const_iterator iter1 = elements_map.begin();
		iter1 != elements_map.end(); iter1++){

		std::vector <int> Connec_Temp = iter1->second;

		for (size_t m = 0; m < Connec_Temp.size(); m++){
			if (m < Connec_Temp.size() - 1){
				int P1 = Connec_Temp[m];
				int P2 = Connec_Temp[m + 1];
				elements_edges[iter1->first].push_back(std::make_pair(P1, P2));
			}
			else{
				int P1 = Connec_Temp[m];
				int P2 = Connec_Temp[0];
				elements_edges[iter1->first].push_back(std::make_pair(P1, P2));
			}
		}

	}

	//Assign flags to the edges

	return 1;
}

int Class_Mesh::ObtainNodesInverseConnectivities()
{
	for (std::map <int, std::vector <int> >::const_iterator iter1 = elements_map.begin();
		iter1 != elements_map.end(); iter1++){

		std::vector <int> Connec = iter1->second;

		for (size_t mm = 0; mm < Connec.size(); mm++){
			int Node_Index = Connec[mm];
			nodes_inverse_connectivities[Node_Index].insert(iter1->first);
		}
	}
	
	return 1;
}

int Class_Mesh::ObtainNodesAdjacentNodes()
{
	std::set<int> Temp;

	//nodes_adjacent_nodes
	for (std::map <int, std::set <int> >::const_iterator iter1 = nodes_inverse_connectivities.begin();
		iter1 != nodes_inverse_connectivities.end(); iter1++){
		
		int Node_Index = iter1->first;
		std::set<int> Eles = iter1->second;

		Temp.clear();

		//for (size_t mm = 0; mm < Eles.size(); mm++){
		for (std::set<int>::const_iterator iter2 = Eles.begin(); iter2 != Eles.end(); iter2++){
			int Ele_Index = *iter2;

			std::vector <int> Connec = elements_map.find(Ele_Index)->second;

			for (size_t mm1 = 0; mm1 < Connec.size(); mm1++){
				if (Connec[mm1] == Node_Index){
					if (mm1 == 0){
						Temp.insert(Connec.back());
						Temp.insert(Connec[mm1+1]);
					}
					else if (mm1 == Connec.size() - 1){
						Temp.insert(Connec.front());
						Temp.insert(Connec[mm1-1]);
					}
					else{
						Temp.insert(Connec[mm1 + 1]);
						Temp.insert(Connec[mm1 - 1]);
					}
					break;
				}
			}

		}

		nodes_adjacent_nodes[Node_Index] = Temp;
	}

	return 1;
}

int Class_Mesh::ObtainElementsEdgeFlags()
{
	//std::map <int, std::vector <bool> > elements_edges_flags; //first: ele index; second: true: on the boundary, false: not on the boundary
	//std::map < std::set <int>, int > edges_map; //first: pair of two nodes; second: edge index (starting from 1)
	//std::map < int, std::set<int> > edges_nodes; //first: edge index (starting from 1); second: pair of two nodes
	//std::map < int, bool > edges_flags; //first: edge index (starting from 1); second: true: on the boundary, false: not on the boundary

	//std::map <int, std::vector <int> > elements_edges_indexes; //first: ele index; second: edge index (starting from 1)
	//std::map <int, std::vector < std::pair <int, int> > > elements_edges; //first: ele index; second: ele edge with two nodes

	for (std::map <int, std::vector < std::pair <int, int> > >::const_iterator iter1 = elements_edges.begin(); iter1 != elements_edges.end(); iter1++){

		int ele_index = iter1->first;

		std::vector < std::pair <int, int> > edges = iter1 -> second;

		//std::cout << "DSFDS" << edges.size() << std::endl;

		std::vector <int> Edge_Indexes;

		for (std::vector <std::pair <int, int>>::const_iterator iter2 = edges.begin(); iter2 != edges.end(); iter2++){

			std::pair <int, int> node_pair = *iter2;

			std::set <int> node_set = { node_pair.first, node_pair.second };

			int Edge_Index = edges_map.find(node_set)->second;

			Edge_Indexes.push_back(Edge_Index);

			//std::cout << ele_index << ", " << Edge_Index << std::endl;
		}

		elements_edges_indexes[ele_index] = Edge_Indexes;

	}

	//std::map <int, std::vector < std::pair <int, int> > > elements_edges;
	//std::map <int, std::vector <bool> > elements_edges_flags;

	for (std::map <int, std::vector < std::pair <int, int> > >::const_iterator iter1 = elements_edges.begin();
		iter1 != elements_edges.end(); iter1++){

		int Ele_Index = iter1->first;

		std::vector < std::pair <int, int> > Edges = iter1->second;

		for (size_t m = 0; m < Edges.size(); m++){
			std::pair <int, int> Node_Pair = Edges[m];

			int N1 = Node_Pair.first;
			int N2 = Node_Pair.second;

			std::set<int> set1 = nodes_inverse_connectivities.find(N1)->second;
			std::set<int> set2 = nodes_inverse_connectivities.find(N2)->second;


			//////////std::vector<int>::iterator it;

			//////////std::sort(first, first + 5);     //  5 10 15 20 25
			//////////std::sort(second, second + 5);   // 10 20 30 40 50

			//////////it = std::set_intersection(first, first + 5, second, second + 5, v.begin());
			//////////// 10 20 0  0  0  0  0  0  0  0
			//////////v.resize(it - v.begin());

			std::vector <int> Intersection(max_intersection_num);
			std::vector <int>::iterator it;
			it = std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), Intersection.begin());
			Intersection.resize(it - Intersection.begin());

			int Edge_Index = (elements_edges_indexes.find(Ele_Index)->second)[m];

			if (Intersection.size() == 1){
				elements_edges_flags[Ele_Index].push_back(true);
				edges_flags[Edge_Index] = true;
			}
			else{
				elements_edges_flags[Ele_Index].push_back(false);
				edges_flags[Edge_Index] = false;
			}

			//PrintScreen_Vector("Intersection", Intersection, 1);
		}
	}

	//std::map < int, bool > edges_flags; //first: edge index (starting from 1); second: true: on the boundary, false: not on the boundary
	//std::map <int, std::vector <bool> > elements_edges_flags; //first: ele index; second: true: on the boundary, false: not on the boundary

	return 1;
}

int Class_Mesh::Obtain_TJunction_ExtraOrdinaryNodes()
{
	for (std::map <int, std::set <int> >::const_iterator iter1 = nodes_inverse_connectivities.begin(); iter1 != nodes_inverse_connectivities.end(); iter1++){
		int node_index = iter1->first;

		if (iter1->second.size() == 3){
			
			std::set <int> Elements = iter1->second;

			int Temp = 0;
			for (std::set <int>::const_iterator iter1 = Elements.begin(); iter1 != Elements.end(); iter1++){
				int Ele_Index = *iter1;

				size_t Ele_Connec_Size = (elements_map.find(Ele_Index)->second).size();

				if (Ele_Connec_Size != 4){
					Temp++;
				}
			}

			if (Temp > 1){
				Nodes_Flags[node_index] = 5;
				Super_TJunction_Nodes.push_back(node_index);
			}
			else{
				Nodes_Flags[node_index] = 1;
				TJunction_Nodes.push_back(node_index);
			}
		}
		else if (iter1->second.size() > 4){
			Extraordinary_Nodes.push_back(node_index);
			Nodes_Flags[node_index] = 2;
		}
		else if (iter1->second.size() == 2){
			BoundaryEdge_Nodes.push_back(node_index);
			Nodes_Flags[node_index] = 3;
		}
		else if (iter1->second.size() == 1){
			Corner_Nodes.push_back(node_index);
			Nodes_Flags[node_index] = 4;
		}
		else{
			Nodes_Flags[node_index] = 0;
		}
	}

	//std::cout << Super_TJunction_Nodes.size() << std::endl;

	return 1;
}