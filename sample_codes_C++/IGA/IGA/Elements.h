//////////#include <vector>
//////////#include <map>
//////////#include "Utilities.h"
//////////#pragma once
//////////
//////////class Class_Elements : virtual New_Class_Utilities::Class_Utilities
//////////{
//////////public:
//////////	Class_Elements(){};
//////////	int Assign(const int& Element_Index, const std::vector <int>& Connectivities);
//////////	virtual ~Class_Elements(){};
//////////
//////////	void PrintElementsInformation();
//////////	std::vector <int>& operator[](const int& Ele_Index){ return elements_map[Ele_Index]; };
//////////
//////////	int IdentifyElementEdges();
//////////
//////////public:
//////////	static std::map <int, std::vector <int> > elements_map;
//////////
//////////private:
//////////	std::map <int, std::vector < std::pair <int, int> > > elements_edges;
//////////	std::map <int, std::vector <bool> > elements_edges_flags;
//////////};