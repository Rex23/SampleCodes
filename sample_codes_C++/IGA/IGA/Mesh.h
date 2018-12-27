#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include "Utilities.h"
#pragma once

class Class_Mesh : virtual New_Class_Utilities::Class_Utilities
{
public:
	Class_Mesh(){};
	virtual ~Class_Mesh() {};

	int Assign_Node(const int& Node_Index, const std::vector <double>& Coords);
	int Assign_Element(const int& Element_Index, const std::vector <int>& Connectivities);

	int IdentifyElementEdges();
	int ObtainNodesInverseConnectivities();
	int ObtainNodesAdjacentNodes();
	int ObtainElementsEdgeFlags();
	int Obtain_TJunction_ExtraOrdinaryNodes();
	int ObtainEdgeIndexes();
	int ObtainNextEdge_of_Node_Edge(const int& Node_Index, const int& Adjacent_Node_Index, int& Current_Edge_Index, int& Next_Edge_Index); //Current and Next Edge Indexes are from 1 //The adjacent node shall not be an extraordinary node
	int Shared_Elements_of_Edge(const int& Edge_Index, std::set<int>& Elements);

	int New_Edges_and_Nodes_Based_on_T_Junction(const int& T_Junction_Node_Index, const std::map <int, std::set <int> >& nodes_inverse_connectivities,
		const std::map <int, std::vector <int> >& elements_map, std::map <int, int> Nodes_Flags, std::map <int, int>& TJunction_Nodes_Ele);

	int New_Edges_and_Nodes_Based_on_Super_T_Junction(const int& Super_T_Junction_Node_Index, const std::map <int, std::set <int> >& nodes_inverse_connectivities,
		const std::map <int, std::vector <int> >& elements_map, std::map <int, int> Nodes_Flags, std::map <int, int>& TJunction_Nodes_Ele);

	int ObtainTJunction_Edges(const int& Node_Index, ...)
	
	void PrintNodesInformation();
	void PrintElementsInformation();

	void Paraview_Visualizations();

	int Paraview_Mesh_Visualization(const std::string& FileName, const std::map <int, std::vector <double> >& nodes_map, 
		const std::map <int, std::vector <int> >& elements_map, const std::map <int, int>& Nodes_Flags);

	int Paraview_Mesh_Visualization_Edges(const std::string& FileName, const std::map <int, std::vector <double> >& nodes_map,
		const std::map < int, std::set<int> >& edges_nodes);

protected:
	std::map <int, std::vector <double> > nodes_map;
	std::map <int, std::set <int> > nodes_inverse_connectivities;
	std::map <int, std::set <int> > nodes_adjacent_nodes;

	std::map <int, std::vector <int> > elements_map; //first: ele index; second: element connectivities
	std::map <int, std::vector < std::pair <int, int> > > elements_edges; //first: ele index; second: ele edge with two nodes
	std::map <int, std::vector <bool> > elements_edges_flags; //first: ele index; second: true: on the boundary, false: not on the boundary

	std::map < std::set <int>, int > edges_map; //first: pair of two nodes; second: edge index (starting from 1)
	std::map < int, std::set<int> > edges_nodes; //first: edge index (starting from 1); second: pair of two nodes
	
	std::map <int, std::vector <int> > elements_edges_indexes; //first: ele index; second: edge index (starting from 1)
	std::map < int, bool > edges_flags; //first: edge index (starting from 1); second: true: on the boundary, false: not on the boundary

	std::vector <int> Super_TJunction_Nodes;
	std::vector <int> TJunction_Nodes;
	std::map <int, int> TJunction_Nodes_Ele; //First: Index in 'TJunction_Nodes', Second: The corresponding element
	std::vector <int> Extraordinary_Nodes;
	std::vector <int> BoundaryEdge_Nodes;
	std::vector <int> Corner_Nodes;
	std::map <int, int> Nodes_Flags; //0: Normal; 1: T Junction; 2: Extraordinary; 3: Boundary edge; 4: Boundary corner; 5: Super T Junction

	const int max_intersection_num = 50;
};