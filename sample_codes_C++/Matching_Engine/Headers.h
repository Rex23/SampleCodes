#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <unordered_map>
#include <list>
#include <tuple>
#include <memory>

using namespace std;

typedef double real;
typedef int32_t integer;
typedef uint32_t uinteger;

//**********If the quantity of the order is larger than 2^31 - 1, define the macro "LARGE_QUANTITY" either here or on compile line.**********
//#define LARGE_QUANTITY
#ifdef LARGE_QUANTITY
typedef uint64_t ulinteger;
typedef int64_t linteger;
#else
typedef uint32_t ulinteger;
typedef int32_t linteger;
#endif
