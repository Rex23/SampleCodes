#include <iostream> 
#include <queue>
#include <stack>
using namespace std;

struct TreeNode
{
	TreeNode* left;
	TreeNode* right;
	int val;
	TreeNode(int x) : left(NULL), right(NULL), val(x) {}
};

class Solution
{
	int max_length = 0;

	int depthOfBinaryTree(TreeNode* root)
	{
		if (root == nullptr)
			return 0;

		auto depth_left = depthOfBinaryTree(root->left);
		auto depth_right = depthOfBinaryTree(root->right);

		max_length = max(max_length, depth_left + depth_right);

		return max(depth_left, depth_right) + 1;
	}

public:
	//Recursive method
	void TreeMethod_Recursive_Search(TreeNode* root) {

		if (!root) {
			return;
		}

		if (root->left)
			TreeMethod_Recursive_Search(root->left);

		//Do something
		cout << root->val << endl;

		if (root->right)
			TreeMethod_Recursive_Search(root->right);

	}

	void TreeMethod_InOrder_Iterative(TreeNode* root)
	{
		stack <TreeNode*> a_stack;

		while (root != nullptr || !a_stack.empty())
		{
			while (root != nullptr)
			{
				a_stack.push(root);
				root = root->left;
			}
			
			root = a_stack.top();
			a_stack.pop();
			cout << root->val << endl;
			root = root -> right;
		}
	}

	void TreeMethod_BFS_Queue(TreeNode* root)
	{
		queue <TreeNode*> a_queue;

		if (root != nullptr)
		{
			a_queue.push(root);
		}

		while (!a_queue.empty())
		{
			TreeNode * a_root = a_queue.front();

			a_queue.pop();

			if (a_root->left != nullptr)
			{
				a_queue.push(a_root->left);
			}

			if (a_root->right != nullptr)
			{
				a_queue.push(a_root->right);
			}

			cout << a_root->val << endl;
		}
	}

	void TreeMethod_DFS_Stack(TreeNode* root)
	{
		stack <TreeNode*> a_stack;

		if (root != nullptr)
		{
			a_stack.push(root);
		}

		while (!a_stack.empty())
		{
			TreeNode * a_root = a_stack.top();

			a_stack.pop();

			if (a_root->right != nullptr)
			{
				a_stack.push(a_root->right);
			}

			if (a_root->left != nullptr)
			{
				a_stack.push(a_root->left);

			}

			cout << a_root->val << endl;
		}
	}

	int TreeMethod_Maximum_Depth(TreeNode* root)
	{
		if (root == nullptr)
			return 0;

		return max(TreeMethod_Maximum_Depth(root->left), TreeMethod_Maximum_Depth(root->right)) + 1;
	}

	TreeNode* TreeMethod_Clone(TreeNode* root) {
		TreeNode* new_root = new TreeNode(root ->val);
		new_root->left = (root->left) ? TreeMethod_Clone(root->left) : nullptr;
		new_root->right = (root->right) ? TreeMethod_Clone(root->right) : nullptr;
		return new_root;
	}

	void TreeMethod_Clone2(TreeNode* original, TreeNode* duplicate) {
		duplicate->val = original->val;
		//	cout<<duplicate->val<<endl;
		if (original->left) {
			TreeNode* newLeft = new TreeNode(original->left->val);
			//newLeft->parent = duplicate;
			duplicate->left = newLeft;
			//		printf("in left copying %d %d %d\n",duplicate->val,duplicate->left,duplicate->parent);
			TreeMethod_Clone2(original->left, duplicate->left);

		}
		if (original->right) {
			TreeNode* newRight = new TreeNode(original->right->val);
			//newRight->parent = duplicate;
			duplicate->right = newRight;
			//		printf("in right copying %d %d %d\n",duplicate->val,duplicate->right,duplicate->parent);
			TreeMethod_Clone2(original->right, duplicate->right);
		}
	}

	TreeNode* TreeMethod_Merge_Two_Trees(TreeNode* root1, TreeNode* root2)
	{
		if (root1 == nullptr && root2 == nullptr) return nullptr;

		int val = (root1 == nullptr ? 0 : root1 -> val) + (root2 == nullptr ? 0 : root2 -> val);
		
		TreeNode* newNode = new TreeNode(val);

		newNode -> left = TreeMethod_Merge_Two_Trees(root1 == nullptr ? nullptr : root1 -> left, root2 == nullptr ? nullptr : root2 -> left);
		newNode -> right = TreeMethod_Merge_Two_Trees(root1 == nullptr ? nullptr : root1->right, root2 == nullptr ? nullptr : root2->right);

		return newNode;
	}

	//The diameter of a binary tree is the length of the longest path between any two nodes in a tree.This path may or may not pass through the root.
	int diameterOfBinaryTree(TreeNode* root)
	{
		depthOfBinaryTree(root);

		return max_length;
	}
};

void single_tree()
{
	//     1 
	//   /  \
	//  2    3
	// / \  / \
        //4  5 6   7
	//[1,2,3,4,5,6,7]

	cout << "     1\n";
	cout << "   /  \\\n";
	cout << "  2    3\n";
	cout << " / \\  / \\\n";
	cout << "4  5 6   7\n";

	TreeNode* root = new TreeNode(1);

	root->left = new TreeNode(2);
	root->right = new TreeNode(3);
	root->left->left = new TreeNode(4);
	root->left->right = new TreeNode(5);
	root->right->left = new TreeNode(6);
	root->right->right = new TreeNode(7);

	Solution sol;
	//Pre Order:  1 2 4 5 3 6 7 
	//In Order:   4 2 5 1 6 3 7 
	//Post Order: 4 5 2 6 7 3 1 
	//BFS:        1 2 3 4 5 6 7

	cout << "In order:\n";
	sol.TreeMethod_Recursive_Search(root);
	cout << "In order (second approach):\n";
	sol.TreeMethod_InOrder_Iterative(root);
	
	cout << "BFS_QUEUE:\n";
	sol.TreeMethod_BFS_Queue(root);
	
	cout << "DFS_STACK:\n";
	sol.TreeMethod_DFS_Stack(root);

	cout << "Maximum Depth of Tree:\n";
	cout << sol.TreeMethod_Maximum_Depth(root) << endl;

	cout << "Diameter of the Tree:\n";
	cout << sol.diameterOfBinaryTree(root) << endl;

	cout << "Clone the tree: (Method 1)\n";
	TreeNode* root_clone = sol.TreeMethod_Clone(root);
	cout << "Maximum Depth of the Clone Tree:\n";
	cout << sol.TreeMethod_Maximum_Depth(root_clone) << endl;
	cout << "In order of the Clone Tree:\n";
	sol.TreeMethod_Recursive_Search(root_clone);

	cout << "Clone the tree: (Method 2)\n";
	TreeNode* duplicate = new TreeNode(root->val);
	sol.TreeMethod_Clone2(root, duplicate);
	cout << "Maximum Depth of the Clone Tree:\n";
	cout << sol.TreeMethod_Maximum_Depth(duplicate) << endl;
	cout << "In order of the Clone Tree:\n";
	sol.TreeMethod_Recursive_Search(duplicate);

}

void two_trees()
{
	//     1 
	//   /  \
	//  2    3
	// / \  / \
    //4  5 6   7
	//[1,2,3,4,5,6,7]
	cout << "     1\n";
	cout << "   /  \\\n";
	cout << "  2    3\n";
	cout << " / \\  / \\\n";
	cout << "4  5 6   7\n";

	TreeNode* root1 = new TreeNode(1);
	root1->left = new TreeNode(2);
	root1->right = new TreeNode(3);
	root1->left->left = new TreeNode(4);
	root1->left->right = new TreeNode(5);
	root1->right->left = new TreeNode(6);
	root1->right->right = new TreeNode(7);

	//    10 
	//   /   \
	//  20    30
	// / \    / \
    //40  50 60  70
	//[10,20,30,40,50,60,70]
	cout << "    10\n";
	cout << "   /  \\\n";
	cout << "  20    30\n";
	cout << " / \\   / \\\n";
	cout << "40  50 60   70\n";

	TreeNode* root2 = new TreeNode(10);
	root2->left = new TreeNode(20);
	root2->right = new TreeNode(30);
	root2->left->left = new TreeNode(40);
	root2->left->right = new TreeNode(50);
	root2->right->left = new TreeNode(60);
	root2->right->right = new TreeNode(70);

	Solution sol;
	TreeNode* merged_tree = sol.TreeMethod_Merge_Two_Trees(root1, root2);
	cout << "BFS_QUEUE:\n";
	sol.TreeMethod_BFS_Queue(merged_tree);
}

void test_tree()
{
	single_tree();
	two_trees();
}

int main()
{
	test_tree();
	system("pause");
	return 0;
}
