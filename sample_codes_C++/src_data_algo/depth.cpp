//Depth first search implemented by stack
#include <stack>
#include <vector>
#include <iostream>

using namespace std;

struct tree
{
	int val;
	tree* left;
	tree* right;
	tree(int x) : val(x), left(NULL), right(NULL){}
};

void generate_tree_structure(tree*& root)
{
	//[1 2 3 6 7 4 5]
	
	root = new tree(1);

	root -> left = new tree(2);
	
	root -> left -> left = new tree(6);

	root -> left -> right = new tree(7);
	
	root -> right = new tree(3);

	root -> right -> left = new tree(4);

	root -> right -> right = new tree(5);
}

int main(int argc, char* argv[])
{
	tree* root;

	generate_tree_structure(root);
	
	stack <tree*> a_stack;

	vector <int> values;
	
	if (root != NULL){
		a_stack.push(root);
		//values.push_back(root -> val);
	}

	while (! a_stack.empty()){
		tree* node = a_stack.top();
		
		values.push_back( node -> val );

		a_stack.pop();
		
		if ( node -> right != NULL ){
			a_stack.push(node -> right);
			//values.push_back(node -> right -> val);
		}

		if ( node -> left != NULL ){
			a_stack.push(node -> left);
			//values.push_back(node -> left -> val);
		}
	}

	for (auto mm=0; mm < values.size(); mm++){
		cout << "Test values: " << values[mm] << endl;
	}

	return 1;
}
