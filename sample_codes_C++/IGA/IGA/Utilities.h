#define _CRT_SECURE_NO_DEPRECATE
#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include <fstream>
#include <iomanip>
#pragma once

namespace New_Class_Utilities{
	
	using namespace std;

	class Class_Utilities
	{
	public:

		int Cox_de_Boor_Recursion(const int& k, const std::vector <double>& Knot_Vector, std::vector <int>& Initial_Intervals, std::vector <std::vector <std::vector <double> > >& Coefficients);

		int Cox_de_Boor_Recursion_Single(const int& k, const std::vector <std::vector <std::vector <double> > >& Coefficients, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, 
			                             std::vector <std::vector <std::vector <double> > >& Coefficients_New);

		int PrintSplines(const int& Mode, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients);

		int Paraview_Visualization(const std::string& FileName, const int& Intervals, const std::vector <double>& Knot_Vector, const std::vector <int>& Initial_Intervals, const std::vector <std::vector <std::vector <double> > >& Coefficients);

		void DivideTerms(string tmp, vector<string>* Terms);

		int MatchString(const char* str, const char* key2);

		char* StringToUpper(char* str);

		int no(char* param);

		template <typename Type1>
		int PrintScreen_Vector(const string& TestName, const vector<Type1>& Data, bool Paused)
		{
			for (size_t m1 = 0; m1 < Data.size(); m1++){
				cout << TestName << ", " << m1 + 1 << ", " << Data[m1] << endl;
			}

			if (Paused){
				system("pause");
			}

			return 1;
		}

		template <typename Type1, typename Type2>
		int PrintScreen_Map(const std::string& TestName, const std::map<Type1, Type2>& Data, bool Paused)
		{
			for (std::map<Type1, Type2>::const_iterator iter1 = Data.begin(); iter1 != Data.end(); iter1++){

				if (std::is_same<Type1, int>::value && std::is_same<Type2, std::vector<double>>::value){

					Type1 First = iter1->first;
					Type2 Second = iter1->second;

					std::cout.precision(6);
					std::cout << std::scientific;
					std::cout << TestName << ": " << First << ", ";
					for (size_t m = 0; m < Second.size(); m++){
						std::cout << Second[m] << ", ";
					}
					std::cout << std::endl;
				}
				else if (std::is_same<Type1, int>::value && std::is_same<Type2, std::vector<int>>::value){

					Type1 First = iter1->first;
					Type2 Second = iter1->second;

					std::cout << TestName << ": " << First << ", ";
					for (size_t m = 0; m < Second.size(); m++){
						std::cout << Second[m] << ", ";
					}
					std::cout << std::endl;
				}

			}

			if (Paused){
				system("pause");
			}

			return 1;
		}

		template <typename Type1, typename Type2>
		int PrintScreen_MapIntInt(const std::string& TestName, const std::map<Type1, Type2>& Data, bool Paused)
		{
			for (std::map<Type1, Type2>::const_iterator iter1 = Data.begin(); iter1 != Data.end(); iter1++){


				if (std::is_same<Type1, int>::value && std::is_same<Type2, int>::value){

					Type1 First = iter1->first;
					Type2 Second = iter1->second;

					std::cout << TestName << ": " << First << ", " << Second << std::endl;
				}

			}

			if (Paused){
				system("pause");
			}

			return 1;
		}

		template <typename Type1, typename Type2, typename Type3>
		int PrintScreen_MapVector(const std::string& TestName, const std::map<Type1, std::vector <std::pair<Type2, Type3> > >& Data, bool Paused)
		{
			for (std::map<Type1, std::vector <std::pair<Type2, Type3> > >::const_iterator iter1 = Data.begin(); iter1 != Data.end(); iter1++){

				if (std::is_same<Type1, int>::value && std::is_same<Type2, int>::value && std::is_same<Type3, int>::value){

					Type1 First = iter1->first;

					//std::cout.precision(6);
					//std::cout << std::scientific;
					std::cout << TestName << ": " << std::endl; 
					std::cout << "Index " << First << ": " << std::endl;
					for (size_t mm = 0; mm < (iter1->second).size(); mm++){
						std::cout << (iter1->second)[mm].first << ", " << (iter1->second)[mm].second << std::endl;
					}
					std::cout << std::endl;
				}
			}

			if (Paused){
				system("pause");
			}

			return 1;
		}

		template <typename Type1, typename Type2>
		int PrintScreen_MapSet(const std::string& TestName, const std::map<Type1, std::set <Type2>>& Data, bool Paused)
		{
			for (std::map<Type1, std::set <Type2>>::const_iterator iter1 = Data.begin(); iter1 != Data.end(); iter1++){

				if (std::is_same<Type1, int>::value && std::is_same<Type2, int>::value){

					Type1 First = iter1->first;
					std::cout << TestName << ": " << std::endl;
					std::cout << "Index " << First << ": " << std::endl;
					std::set <Type2> the_set = iter1->second;
					for (std::set<Type2>::const_iterator iter2 = the_set.begin(); iter2 != the_set.end(); iter2++){
						std::cout << *iter2 << std::endl;
					}
					std::cout << std::endl;
				}
			}

			if (Paused){
				system("pause");
			}

			return 1;
		}

		template <typename ValueType>
		int ReadLine(std::string tmp, std::vector <ValueType>& row) //Read a line of file separated by comma
		{
			//printf("Test tmp, %s\n", tmp.c_str());
			row.clear();
			int pre_pos = 0;
			int RowIndex = 0;
			//int Pre;
			//bool Pre_Bool = false;
			for (size_t i = 0; i < tmp.size(); i++)
			{
				string tmp_str;
				if (tmp[i] == ',')
				{
					stringstream r;
					tmp_str = tmp.substr(pre_pos, i - pre_pos);
					r << tmp_str;
					ValueType tt;
					r >> tt;
		
					row.push_back(tt);
					pre_pos = i + 1;
					//Pre = i + 1;
					//Pre_Bool = true;
				}
				else if (i == tmp.size() - 1)
				{
					//if ( Pre_Bool == true && Pre == i + 1 ){
					//	break;
					//}
		
					stringstream r;
					tmp_str = tmp.substr(pre_pos, i + 1 - pre_pos);
		
					size_t f1 = tmp_str.find_first_of("0", 0);
					size_t f2 = tmp_str.find_first_of("1", 0);
					size_t f3 = tmp_str.find_first_of("2", 0);
					size_t f4 = tmp_str.find_first_of("3", 0);
					size_t f5 = tmp_str.find_first_of("4", 0);
					size_t f6 = tmp_str.find_first_of("5", 0);
					size_t f7 = tmp_str.find_first_of("6", 0);
					size_t f8 = tmp_str.find_first_of("7", 0);
					size_t f9 = tmp_str.find_first_of("8", 0);
					size_t f10 = tmp_str.find_first_of("9", 0);
		
					//cout << "Test tmp_str" << tmp_str << endl;
		
					if (
						(f1 != string::npos) ||
						(f2 != string::npos) ||
						(f3 != string::npos) ||
						(f4 != string::npos) ||
						(f5 != string::npos) ||
						(f6 != string::npos) ||
						(f7 != string::npos) ||
						(f8 != string::npos) ||
						(f9 != string::npos) ||
						(f10 != string::npos)
						)
					{
						r << tmp_str;
						ValueType tt;
						r >> tt;
						row.push_back(tt);
						break;
					}
					else{
						break;
					}
				}
			}
			return 1;
		}

	};

} //namespace New_Class_Utilities{