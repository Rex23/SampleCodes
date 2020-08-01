//https://www.geeksforgeeks.org/trapping-rain-water/
//algorithm:
//generate two arrays a and b.
//from left to right, find the maximum number and give the number into a. 
//from right to left, find the maximum number and give the number into b.
//for any element i, the accumaluated water is: min(a[i], b[i]) - array_itself[i]
//lesson learned: an array can be generated to obtain the maximum number ever found by going through either from left or from right.

#include <iostream>
#include <vector>

using namespace std;

int main()
{
   std::vector <int> mountains{0, 1, 0, 2, 1, 0, 1, 3, 2, 1, 2, 1}; //{3,0,2,0,4};
   std::vector <int> left_to_right_max(mountains.size());
   std::vector <int> right_to_left_max(mountains.size());
   int max_num = 0;
   for (size_t i = 0; i < left_to_right_max.size(); ++i)
   {
       max_num = std::max(max_num, mountains[i]);
       left_to_right_max[i] = max_num;
   }

   max_num = 0;
   for (int i = right_to_left_max.size() - 1; i >= 0; i--)
   {
       max_num = std::max(max_num, mountains[i]);
       right_to_left_max[i] = max_num;
   }

   int accumulated_water = 0;
   
   for (size_t i = 0; i < mountains.size(); i++)
   {
       double water = std::min(left_to_right_max[i], right_to_left_max[i]) - mountains[i];
       if ( water >= 0)
        accumulated_water += water; 
   }

   std::cout << "total accumulated water: " << accumulated_water << std::endl;
   
   return 0;
}

