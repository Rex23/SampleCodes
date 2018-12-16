#include<string.h>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <map>
#include <typeinfo>
#include <vector>

using namespace std;

int main(int argc, char** argv)
{
	char abc[256] = "test";

	cout << "Test abc.size(): " << strlen(abc) << endl;
	
	cout << sizeof(abc) << endl;

	string dd = "sacd";
	string dd2;

	transform(dd.begin(), dd.end(), dd.begin(), ::toupper);

	cout << "Test dd2: " << dd << endl;

	char abcde[] = "sadfdasfasdfadsf";

	cout << sizeof(abcde) << endl;
	
	stringstream ss;

	string a("a");
	string b("b");
	string c("c");

	ss << a << " " << b << c;
	
	string d;

	ss >> d;

	cout << ss.str() << endl;

	cout << d << endl;
	
	string temp = "abc";
	string temp2 = "d";
	
	temp.insert(temp.end(), temp2[0]);

	cout << temp << endl;

	map <int, int> mm = {{1,1},{2,2},{3,3}};
	
	cout << mm.size() << endl;
	cout << mm.count(1) << endl;

	mm.insert({4,4});

	cout << mm.count(4) << endl;

	char dum2[256] = "dsfdsfdfsd";
	
	dum2[0] = 'x';

	cout << typeid(dum2).name() << endl;

	//char*& dum3 = dum2;

	//cout << sizeof(dum3) << endl; 
	
	cout << "test max: " << max(1,2) << endl;

	int a_num = 12345;

	string dumm = to_string(a_num);

	vector <int> digits;
	
	cout << "Test string: " << dumm << endl;

	for (auto i : dumm){
		digits.push_back( (int) i - 48 );
	}

	for (auto m : digits){
		cout << m << endl;
	}


	return 1;
}
