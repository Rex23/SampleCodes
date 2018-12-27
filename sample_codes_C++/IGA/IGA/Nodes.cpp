//////////#include "Nodes.h"
//////////
//////////int Class_Nodes::ObtainNodesInverseConnectivities()
//////////{
//////////	for (std::map <int, std::vector <int> >::const_iterator iter1 = Class_Elements::elements_map.begin();
//////////		iter1 != Class_Elements::elements_map.end(); iter1++){
//////////
//////////		std::cout << "Test 1" << std::endl;
//////////
//////////	}
//////////	
//////////	return 1;
//////////}
//////////
//////////int Class_Nodes::Assign(const int& Node_Index, const std::vector <double>& Coords)
//////////{
//////////	nodes_map[Node_Index] = Coords;
//////////	return 1;
//////////}
//////////
//////////void Class_Nodes::PrintNodesInformation()
//////////{
//////////	PrintScreen_Map("Nodes", nodes_map, false);
//////////}
//////////
