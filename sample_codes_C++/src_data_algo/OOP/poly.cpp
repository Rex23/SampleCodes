//Author: Ren, Rex
//Test polymorphism:
//1. if a function's argument is the reference of the base class, it can take the derived class as the argument.
//2. if a function's argument is the refernece of the derived class, it cannot take the base class as the argument.
//3. members of a class cannot use the keywork virtual. If the derived class is the argument, the member used is in the base class.

#include "utilities.h"

class A
{
	public:
	int val = 3;
	virtual void message(){cout << "Test message." << endl;}
};

class B : public A
{
	public:
	int val = 4;
	virtual void message(){cout << "Test message2." << endl;}
};

int a_fun(A& dum)
{
	cout << "a_fun" << endl;
	
	cout << dum.val << endl;
	
	dum.message();
	return 1;
}

int main(int argc, char** argv)
{
	A A_instance;
	
	B B_instance;

	a_fun(A_instance);

	return 1;
}

