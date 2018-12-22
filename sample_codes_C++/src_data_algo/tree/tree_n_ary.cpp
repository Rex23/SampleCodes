//Author: Ren, Rex

//Postorder traversal

#include "../template/utilities.h"

void generate_tree_n_ary()
{
	Node_n_ary <string> * root = new Node_n_ary <string> ();
	
	root -> children.resize(3);

	root -> children[0] = new Node_n_ary <string>("B");

	root -> children[1] = new Node_n_ary <string>("C");

	root -> children[2] = new Node_n_ary <string>("D");

	root -> children[0] -> children.resize(2);

	root -> children[0] -> children[0] = new Node_n_ary <string>("E");
	root -> children[0] -> children[1] = new Node_n_ary <string>("F");

	root -> children[0] -> children[1] -> children.resize(1);

	root -> children[0] -> children[1] -> children[0] = new Node_n_ary <string>("K");

	root -> children[2] -> children.resize(4);

	root -> children[2] -> children[0] = new Node_n_ary <string>("G");
	root -> children[2] -> children[1] = new Node_n_ary <string>("H");
	root -> children[2] -> children[2] = new Node_n_ary <string>('I');
	root -> children[2] -> children[3] = new Node_n_ary <string>('J');
}

int main(int argc, char** argv)
{
	
	return 1;
}
