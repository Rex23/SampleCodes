//////////#include "Elements.h"
//////////
//////////std::map <int, std::vector <int> > Class_Elements::elements_map;
//////////
//////////int Class_Elements::Assign(const int& Element_Index, const std::vector <int>& Connectivities)
//////////{
//////////	Class_Elements::elements_map[Element_Index] = Connectivities;
//////////	return 1;
//////////}
//////////
//////////void Class_Elements::PrintElementsInformation()
//////////{
//////////	PrintScreen_Map("Elements", Class_Elements::elements_map, false);
//////////	PrintScreen_MapVector("Test Edges", elements_edges, 1);
//////////}
//////////
//////////int Class_Elements::IdentifyElementEdges()
//////////{
//////////	//Save the edges
//////////	for (std::map <int, std::vector <int> >::const_iterator iter1 = Class_Elements::elements_map.begin();
//////////		iter1 != Class_Elements::elements_map.end(); iter1++){
//////////
//////////		vector <int> Connec_Temp = iter1->second;
//////////
//////////		for (size_t m = 0; m < Connec_Temp.size() - 1; m++){
//////////			int P1 = Connec_Temp[m];
//////////			int P2 = Connec_Temp[m + 1];
//////////			elements_edges[iter1 -> first].push_back( std::make_pair(P1, P2) );
//////////		}
//////////
//////////	}
//////////
//////////	//Assign flags to the edges
//////////
//////////	return 1;
//////////}