#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <iostream>
#include <sstream>
#include <string>
#include <set>
#include <list>
#include <unordered_set>
#include <unordered_set>
#include <typeinfo>
#include <string.h>
#include <bitset>

using namespace std;

template <typename T>
void printvector(const string& a_string, const vector <T>& a_vector)
{
	for (auto it = a_vector.begin(); it != a_vector.end(); it++){
		cout << a_string << ", " << it - a_vector.begin() << ", " << *it << endl;
	}
}

template <typename T>

class Node_n_ary
{
public:
	T val;

	vector <Node_n_ary*> children;

	Node_n_ary(){}

	Node_n_ary(T _val) {
        	val = _val;
        }
};
