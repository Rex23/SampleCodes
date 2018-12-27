//////////#include <map>
//////////#include <vector>
//////////#include "Utilities.h"
//////////#include "Elements.h"
//////////#pragma once
//////////
//////////class Class_Nodes : virtual New_Class_Utilities::Class_Utilities
//////////{
//////////public:
//////////	//Nodes(const int& Node_Index, const std::vector <double>& Nodes_Coords){ nodes_map[Node_Index] = Nodes_Coords; };
//////////	Class_Nodes(){};
//////////	int Assign(const int& Node_Index, const std::vector <double>& Coords);
//////////	virtual ~Class_Nodes(){};
//////////
//////////	void PrintNodesInformation();
//////////
//////////	std::vector <double>& operator[](const int& Node_Index){ return nodes_map[Node_Index]; };
//////////
//////////	int ObtainNodesInverseConnectivities();
//////////
//////////public:
//////////	std::map <int, std::vector <int> > nodes_inverse_connectivities;
//////////
//////////private:
//////////	std::map <int, std::vector <double> > nodes_map;
//////////	//std::vector <int> m_list;
//////////};