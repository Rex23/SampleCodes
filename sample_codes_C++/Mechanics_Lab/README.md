The program is for obtaining a 3D solid structure's failure response under extern loadings by using the finite element approach. It is written in C++ using my spare time. I will continue to add advanced features of the package in terms of both physics, data structure, and testing.

To build thie package, type:

make -f Makefile

To clean the package, type:

make -f Makefile clean

To run google test, type:

make -f Makefile google

To clean the google test objects, type:

make -f Makefile cleangoogle

To run the executable, go to ./lab and type:

./LAB Test.in
