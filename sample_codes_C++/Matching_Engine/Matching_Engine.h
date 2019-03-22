#pragma once
#include "Headers.h"
#include "Book.h"

class Base_Matching_Engine
{
public:
	virtual void print_bad_message(const unsigned type, ostream& o_stream) const = 0;
	virtual void process_order(string& tmp, istream& in_stream, ostream& o_stream, ostream& o_stream_err) = 0;
	virtual bool separate_order_items(string& tmp, tuple<uinteger, uinteger, uinteger, ulinteger, real>& items, ostream& o_stream) = 0;
};

class Matching_Engine : public Base_Matching_Engine
{
public:
	Matching_Engine() = default;

	virtual ~Matching_Engine() = default;

	Matching_Engine(istream& in_stream, ostream& o_stream, ostream& o_stream_err, book* order_book_in = nullptr) : in_stream(in_stream), o_stream(o_stream), o_stream_err(o_stream_err)
	{
		//Either load the pre-existing order book or construct it here:
		order_book = order_book_in == nullptr ? unique_ptr <book>(new book()) : unique_ptr <book>(order_book_in);

		resting_sell = order_book -> get_resting_queue_of_sell();
		resting_buy = order_book->get_resting_queue_of_buy();
		resting_order_ID = order_book->get_resting_order_IDs();

		string tmp;

		while (getline(in_stream, tmp)) //read input from istream
		{
			process_order(tmp, in_stream, o_stream, o_stream_err);
		}
	}

private:
	istream& in_stream;
	ostream& o_stream;
	ostream& o_stream_err;

	void process_order(string& tmp, istream& in_stream, ostream& o_stream, ostream& o_stream_err);
	bool separate_order_items(string& tmp, tuple<uinteger, uinteger, uinteger, ulinteger, real>& new_order, ostream& o_stream);
	void insert_order_into_book(tuple<uinteger, uinteger, uinteger, ulinteger, real>& new_order);

	void make_a_deal(const uinteger& msgtype, const uinteger& orderid, const uinteger& side, ulinteger& quantity, const real& price); //try to make a deal
	void book_keeping(const uinteger& side, const real& price, const uinteger& orderid, const ulinteger& quantity); //rest the order
	void cancel_order(const uinteger& orderid);

	void print_bad_message(const unsigned type, ostream& o_stream) const;

	void print_trade_event(ostream& o_stream, const ulinteger& quantity, const real& price) const
	{
		o_stream << 2 << "," << quantity << "," << price << endl;
	}

	void print_orderfullyfilled(ostream& o_stream, const uinteger& orderid) const
	{
		o_stream << 3 << "," << orderid << endl;
	}

	void print_orderpartiallyfilled(ostream& o_stream, const uinteger& orderid, const ulinteger& quantity) const
	{
		o_stream << 4 << "," << orderid << "," << quantity << endl;
	}

	unique_ptr<book> order_book;

	//Resting Sell Queue
	map <real, list < pair<uinteger, ulinteger> >>* resting_sell; //first: price, second's first: order ID, second's second: quantity
	//Resting Buy Queue
	map <real, list < pair<uinteger, ulinteger> >>* resting_buy; //first: price, second's first: order ID, second's second: quantity
	//Resting Order ID
	unordered_map <uinteger, tuple <uinteger, real, list < pair<uinteger, ulinteger> >::iterator > >* resting_order_ID; //first: order ID, second's first: side, second's second: price, second's third: iterator of the list
};