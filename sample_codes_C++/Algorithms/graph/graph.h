#pragma once
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>

using namespace std;

class Node {
public:
	int val;
	vector<Node*> neighbors;

	Node() {}
	Node(int _val) {}
	Node(int _val, vector<Node*> _neighbors) {
		val = _val;
		neighbors = _neighbors;
	}
};

unordered_set <Node*> visited;

void GraphMethod_DFS(Node* root)
{
	if (root == nullptr) return;

	visited.insert(root);

	auto neighbors = root->neighbors;

	for (auto m : neighbors)
	{
		if (visited.find(m) == visited.end()) GraphMethod_DFS(m);
	}

	cout << "value: " << root->val << endl;
}

bool GraphMethod_DFS_Find_A_Target(Node* root, int target) //Find a target within the graph by using DFS
{
	if (root == nullptr) return false;

	if (root->val == target) return true;

	visited.insert(root);

	auto neighbors = root->neighbors;

	for (auto m : neighbors)
	{
		if (visited.find(m) == visited.end())
		{
			if (GraphMethod_DFS_Find_A_Target(m, target)) return true; //If any neighbor can have the target, return true. Note "return GraphMethod_DFS_Find_A_Target(m, target)" is wrong.
		}
	}

	cout << "value: " << root->val << endl;

	return false;
}

void GraphMethod_BFS(Node* root)
{
	unordered_set <Node*> visited;
	visited.insert(root);

	if (root == nullptr)
		return;

	queue <Node*> a_queue;
	
	a_queue.push(root);

	while (!a_queue.empty())
	{
		Node* a_root = a_queue.front();
		a_queue.pop();

		for (auto iter = a_root->neighbors.begin(); iter != a_root->neighbors.end(); iter++)
		{
			if (visited.find(*iter) != visited.end())
			{
				continue;
			}
			else {
				a_queue.push(*iter);
				visited.insert(*iter);
			}
		}

		cout << "value: " << a_root->val << endl;
	}
}

bool GraphMethod_BFS_Find_A_Target(Node* root, int target)
{
	unordered_set <Node*> visited;
	visited.insert(root);

	if (root == nullptr)
		return false;

	queue <Node*> a_queue;

	if (root->val == target) return true;

	a_queue.push(root);

	while (!a_queue.empty())
	{
		Node* a_root = a_queue.front();
		a_queue.pop();

		for (auto iter = a_root->neighbors.begin(); iter != a_root->neighbors.end(); iter++)
		{
			if (visited.find(*iter) != visited.end())
			{
				continue;
			}
			else {
				a_queue.push(*iter);
				visited.insert(*iter);

				if ((*iter)->val == target) return true;
			}
		}

		cout << "value: " << a_root->val << endl;
	}

	return false;
}

void Initialize_Graph()
{
	//  2       9
	// / \     / \
	//1 - 4 - 6 - 7 - 5
	// \ /     \ 
	//  3       8

	Node *n1 = new Node();
	Node *n2 = new Node();
	Node *n3 = new Node();
	Node *n4 = new Node();
	Node *n5 = new Node();
	Node *n6 = new Node();
	Node *n7 = new Node();
	Node *n8 = new Node();
	Node *n9 = new Node();

	n1->val = 1;
	n2->val = 2;
	n3->val = 3;
	n4->val = 4;
	n5->val = 5;
	n6->val = 6;
	n7->val = 7;
	n8->val = 8;
	n9->val = 9;

	n1->neighbors = { n2,n3,n4 };
	n2->neighbors = { n1,n4 };
	n3->neighbors = { n1,n4 };
	n4->neighbors = { n1,n2,n3,n6 };
	n5->neighbors = { n7 };
	n6->neighbors = { n4,n8,n9,n7 };
	n7->neighbors = { n5,n6,n9 };
	n8->neighbors = { n6 };
	n9->neighbors = { n6,n7 };

	cout << "Graph BFS: " << endl;
	GraphMethod_BFS(n1);

	cout << "Graph DFS: " << endl;
	GraphMethod_DFS(n1);

	visited.clear();
	cout << "Graph DFS to find a target: " << endl;
	cout << boolalpha << GraphMethod_DFS_Find_A_Target(n1, 4) << endl;

	cout << "Graph BFS to find a target: " << endl;
	cout << boolalpha << GraphMethod_BFS_Find_A_Target(n1, 4) << endl;
}

void test_graph()
{
	Initialize_Graph();
}
