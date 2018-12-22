//Author: Ren, Rex

//Readin a file like data

#include "utilities.h"

int main(int argc, char** argv)
{
	string a_string;
	
	vector <string> strings_read;

	while(getline(cin, a_string)){
		strings_read.push_back(a_string);
	}

	printvector("Test: ", strings_read);

	cout << a_string << "\n";

	return 1;
}
