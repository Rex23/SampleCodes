/*convert a number into its left to right digits*/
#include <iostream>
#include <deque>

using namespace std;

std::deque <int> convert_int_to_digits(int number)
{
    std::deque <int> digits;
    
    if (number == 0) return {0};
    
    while (number > 0)
    {
        digits.push_front( number % 10 );    
        number = number / 10;
    }
    
    return digits;
}

int main()
{
   int number = 123789;
   
   std::deque <int> digits = convert_int_to_digits(number);
   
   for (const auto& a_digit : digits)
   {
       std::cout << "Test digit: " << a_digit << std::endl;
   }
   
   return 0;
}
