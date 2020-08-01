#pragma once
#include <list>

 //Definition for singly-linked list.
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


