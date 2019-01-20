/*!
Used to save local parameters for local stiffness assembly.
*/


#pragma once
#include <vector>
#include <string>
#include "C3D8.h"

struct struct_Local_Para
{
	struct_Local_Para() { };
	virtual ~struct_Local_Para() { };

	//UEL:
	std::vector <double> Rhs;
	std::vector <std::vector <double> > Amatrx;
	std::vector <double> Svars;
	double Energy;
	int Ndofel;
	int Nrhs;
	int Nsvars;
	std::vector <double> Props;
	int Nprops;
	std::vector <std::vector <double> > Coords;
	int Mcrd;
	int Nnode;
	std::vector <std::vector <double> > U;
	std::vector <double> Du;
	std::vector <double> V;
	std::vector <double> A;
	int Jtype;
	double Time;
	double Dtime;
	int Kstep;
	int Kinc;
	int Jelem;
	std::vector <double> Params;
	int Ndload;
	int Jdltyp;
	double Adlmag;
	double Predef;
	int Npredf;
	int Lflags;
	int Mlvarx;
	double Ddlmag;
	int Mdload;
	double Pnewdt;
	std::vector <int> Jprops;
	int Njprop;
	double Period;

	std::string Ele_Type;
	std::string Material_Name;
	std::vector <std::vector <double> > CC;
	std::vector <double> Sig;
	std::string Step_Name;
};
