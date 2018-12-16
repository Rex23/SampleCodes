#include <string>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
	string dum = "aabbbcccc";
	int max_length = 0;

	for (auto i = 0; i < dum.size(); i++){
		
		char dum_char = dum[i];
		
		int len = 1;

		for (auto j = i + 1; j < dum.size(); j++){
			cout << "Test1: " << dum[j] << ", " << dum_char << endl;	
			if ( dum[j] == dum_char ){
				len++;	
			}
			else{
				cout << "Break!" << endl;
				//i = j;
				break;
			}
		}

		max_length = max(max_length, len);
		//cout << "Test " << i << ", " << max_length << endl;
	}

	cout << "max_length: " << max_length << endl;

	return 1;
}
