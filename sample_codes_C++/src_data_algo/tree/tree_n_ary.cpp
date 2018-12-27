//Author: Ren, Rex

//Postorder traversal for n-ary tree

#include "../template/utilities.h"

//                A
//             / |   \
//            B  C    D
//           / \   / | \ \
//          E  F  G  H  I J
//            /
//            K

Node_n_ary <string>* generate_tree_n_ary()
{
	Node_n_ary <string> * root = new Node_n_ary <string> ();

	root -> val = "A";	
	
	root -> children.resize(3);

	root -> children[0] = new Node_n_ary <string>();

	root -> children[0] -> val = "B";
	
	root -> children[1] = new Node_n_ary <string>();

	root -> children[1] -> val = "C";

	root -> children[2] = new Node_n_ary <string>();

	root -> children[2] -> val = "D";

	root -> children[0] -> children.resize(2);

	root -> children[0] -> children[0] = new Node_n_ary <string>();

	root -> children[0] -> children[0] -> val = "E";

	root -> children[0] -> children[1] = new Node_n_ary <string>();

	root -> children[0] -> children[1] -> val = "F";

	root -> children[0] -> children[1] -> children.resize(1);

	root -> children[0] -> children[1] -> children[0] = new Node_n_ary <string>();

	root -> children[0] -> children[1] -> children[0] -> val = "K";

	root -> children[2] -> children.resize(4);

	root -> children[2] -> children[0] = new Node_n_ary <string>();
	
	root -> children[2] -> children[0] -> val = "G";
	
	root -> children[2] -> children[1] = new Node_n_ary <string>();

	root -> children[2] -> children[1] -> val = "H";

	root -> children[2] -> children[2] = new Node_n_ary <string>();

	root -> children[2] -> children[2] -> val = "I";

	root -> children[2] -> children[3] = new Node_n_ary <string>();

	root -> children[2] -> children[3] -> val = "J";

	return root;
}

void postorder_traversal(Node_n_ary <string>* root, vector <string>& vals)
{
	if (root == NULL)
		return;
	
	for (auto& m : root -> children){
		if ( m != NULL){
			postorder_traversal(m, vals);
		}
	}

	vals.push_back(root -> val);

	return;
}

int main(int argc, char** argv)
{
	Node_n_ary <string> * root = generate_tree_n_ary();

	vector <string> vals;

	postorder_traversal(root, vals);

	printvector("Test nodal val: ", vals);

	return 1;
}
