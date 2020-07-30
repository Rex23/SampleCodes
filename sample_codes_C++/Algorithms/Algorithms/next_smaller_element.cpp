#include <vector>
#include <iostream>

using namespace std;

int main()
{
    
    vector<int> A{1,2,3,2,1};
    
    vector<int> stack;
    for (int i = 0; i < A.size(); ++i) {
        while (stack.size() && A[stack.back()] >= A[i]) {
            std::cout << "Test the pair: " << A[stack.back()] << ", " << A[i] << std::endl;
            stack.pop_back();
        }
        stack.push_back(i);
    }
    
    return 0;
}
