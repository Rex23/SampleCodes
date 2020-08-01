#pragma once
#include <list>

 //definition for singly-linked list.
 struct ListNode 
 {
      int val;
      ListNode *next;
      ListNode() : val(0), next(nullptr) {}
      ListNode(int x) : val(x), next(nullptr) {}
      ListNode(int x, ListNode *next) : val(x), next(next) {}
 };

//loop through the list from the head
while (head != nullptr)
{
    binary.push_back(head->val);
    head = head->next;
}

//problem 1290:
//convert a list of binary numbers to a decimal number

//method 1: bit manipulation
class Solution {
public:
    int getDecimalValue(ListNode* head) {
        int ret = 0;
        while(head)
        {
            ret <<= 1;
            ret |= head->val;
            head = head->next;
        }
        return ret;
    }
};

//method 2: (python) (Read the binary number from MSB to LSB?)
class Solution:
    def getDecimalValue(self, head: ListNode) -> int:
        answer = 0
        while head: 
            answer = 2*answer + head.val 
            head = head.next 
        return answer 

