//Problem 345: reverse vowels of a string
//Write a function that takes a string as input and reverse only the vowels of a string.

//Example 1:

//Input: "hello"
//Output: "holle"
//Example 2:

//Input: "leetcode"
//Output: "leotcede"
//Note:
//The vowels does not include the letter "y".

//Solution 1:
class Solution {
public:
    string reverseVowels(string s) {
        int i = 0, j = s.size() - 1;
        while (i < j) {
            i = s.find_first_of("aeiouAEIOU", i);
            j = s.find_last_of("aeiouAEIOU", j);
            if (i < j) {
                swap(s[i++], s[j--]);
            }
        }
        return s;
    }
};

//Solution 2:
public:
    string reverseVowels(string s) {
        int i = 0, j = s.length()-1;
        while(i<j) {
            while(i<j && !isVowel(s[i])) i++;
            while(i<j && !isVowel(s[j])) j--;
            if(i<j) swap(s[i++], s[j--]);
        }
        return s;
    }
private:
    bool isVowel(char c) const{
        return c=='a' || c=='e' || c=='i' || c=='o' || c=='u' || 
            c=='A' || c=='E' || c=='I' || c=='O' || c=='U';
    }
