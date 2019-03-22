#pragma once
#include "Headers.h"

using namespace std;

class book
{
public:
	map <real, list < pair<uinteger, ulinteger> >>* get_resting_queue_of_sell() { return &resting_sell; }
	map <real, list < pair<uinteger, ulinteger> >>* get_resting_queue_of_buy() { return &resting_buy; }
	unordered_map <uinteger, tuple <uinteger, real, list < pair<uinteger, ulinteger> >::iterator > >* get_resting_order_IDs(){return &resting_order_ID;}

private:
	//Resting Sell Queue
	map <real, list < pair<uinteger, ulinteger> >> resting_sell; //first: price, second: first: order ID, second: quantity
	//Resting Buy Queue
	map <real, list < pair<uinteger, ulinteger> >> resting_buy; //first: price, second: first: order ID, second: quantity
	//Resting Order ID
	unordered_map <uinteger, tuple <uinteger, real, list < pair<uinteger, ulinteger> >::iterator > > resting_order_ID; //first: order ID, second: first: side, second: price, third: iterator of the list
};
