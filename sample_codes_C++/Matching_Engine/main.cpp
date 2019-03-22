#include "Matching_Engine.h"
#include "Book.h"
#include "string.h"

void run_testing_cases(const unsigned& start, const unsigned& end)
{
	for (auto m = start; m <= end; m++)
	{
		string filename_in = "./Case" + to_string(m) + "/in.txt";
		string filename_out = "./Case" + to_string(m) + "/out.txt";
		string filename_out_err = "./Case" + to_string(m) + "/out_err.txt";

		ifstream instream(filename_in);
		ofstream ostream(filename_out);
		ofstream ostream_err(filename_out_err);

		unique_ptr <Base_Matching_Engine> Matching_Engine_Instance(new Matching_Engine(instream, ostream, ostream_err));

		instream.close();
		ostream.close();
		ostream_err.close();
	}
}

int main(int argc, char* argv[])
{
	//Run testing cases 1 to 20
	if (argc == 2 && strcmp(argv[1],"RUN_TESTS") == 0) run_testing_cases(1, 20);

	unique_ptr <Base_Matching_Engine> Matching_Engine_Instance(new Matching_Engine(cin, cout, cerr));

	return 0;
}