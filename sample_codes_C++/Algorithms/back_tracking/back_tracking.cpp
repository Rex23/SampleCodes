/*
Leetcode: 46. Permutations
Given a collection of distinct integers, return all possible permutations.

Example:

Input: [1, 2, 3]
	Output :
	[
		[1, 2, 3],
		[1, 3, 2],
		[2, 1, 3],
		[2, 3, 1],
		[3, 1, 2],
		[3, 2, 1]
	]
*/
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void helper(std::vector <std::vector <int> >& combinations, std::vector <int> solution, std::vector <int> numbers)
{
	//method 1:
	//if (numbers.size() == 0)
	//{
	//	combinations.push_back(solution);
	//	return;
	//}

	//for (int i = 0; i < numbers.size(); ++i)
	//{
	//	solution.push_back(numbers[i]);
	//	std::vector <int> num_copy(numbers);
	//	num_copy.erase(std::remove(num_copy.begin(), num_copy.end(), numbers[i]), num_copy.end());
	//	helper(combinations, solution, num_copy);
	//	solution.pop_back();
	//}

	//method 2:
	if (solution.size() == numbers.size())
	{
		combinations.push_back(solution);
		return;
	}
	for (int i = 0; i < numbers.size(); ++i)
	{
		if (find(solution.begin(), solution.end(), numbers[i]) != solution.end()) continue;
		solution.push_back(numbers[i]);
		helper(combinations, solution, numbers);
		solution.pop_back();
	}
}

int main()
{
	std::vector <int> numbers{ 1,2,3 };
	std::vector <std::vector <int> > combinations;
	std::vector <int> solution;

	helper(combinations, solution, numbers);

	for (int i = 0; i < combinations.size(); ++i)
	{
		std::cout << "Test com:" << std::endl;
		for (int j = 0; j < combinations[i].size(); ++j)
		{
			std::cout << "test: " << combinations[i][j] << std::endl;
		}
	}

	system("pause");

	return 0;
}
