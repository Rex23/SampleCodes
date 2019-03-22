#include "Matching_Engine.h"

void Matching_Engine::cancel_order(const uinteger& orderid)
{
	auto iter = resting_order_ID->find(orderid);
	if (iter != resting_order_ID->end())
	{
		auto side = get<0>(iter->second);
		auto price = get<1>(iter->second);
		auto list_iter = get<2>(iter->second);

		if (side == 0) //Buy
		{
			(*resting_buy)[price].erase(list_iter);
		}
		else if (side == 1) //Sell
		{
			(*resting_sell)[price].erase(list_iter);
		}

		resting_order_ID->erase(iter);
	}
}

void Matching_Engine::book_keeping(const uinteger& side, const real& price, const uinteger& orderid, const ulinteger& quantity)
{
	if (side == 0) //Buy
	{
		(*resting_buy)[price].push_back(make_pair(orderid, quantity));
		(*resting_order_ID)[orderid] = make_tuple(side, price, next((*resting_buy)[price].end(),-1));
	}
	else if (side == 1) //Sell
	{
		(*resting_sell)[price].push_back(make_pair(orderid, quantity));
		(*resting_order_ID)[orderid] = make_tuple(side, price, next((*resting_sell)[price].end(),-1));
	}
}

void Matching_Engine::make_a_deal(const uinteger& msgtype, const uinteger& orderid, const uinteger& side, ulinteger& quantity, const real& price)
{
	if (side == 0) //Buy
	{
		if (resting_sell -> size() == 0)
			return;
		else
		{
			auto iter = resting_sell->begin();

			while (iter != resting_sell->end())
			{
				auto resting_lowest_price = iter->first;

				if (resting_lowest_price > price)
				{
					return;
				}
				else
				{
					auto& same_price_orders = iter->second;

					auto iter2 = same_price_orders.begin();

					while (iter2 != same_price_orders.end())
					{
						auto resting_order_id = iter2->first;
						auto& resting_quantity = iter2->second;

						if (resting_quantity > quantity)
						{
							resting_quantity -= quantity;

							print_trade_event(o_stream, quantity, resting_lowest_price);
							print_orderfullyfilled(o_stream, orderid);
							print_orderpartiallyfilled(o_stream, resting_order_id, resting_quantity);

							quantity = 0;
							return;
						}
						else if (resting_quantity == quantity)
						{
							same_price_orders.erase(iter2);
							resting_order_ID->erase(resting_order_id);

							print_trade_event(o_stream, quantity, resting_lowest_price);
							print_orderfullyfilled(o_stream, orderid);
							print_orderfullyfilled(o_stream, resting_order_id);

							quantity = 0;

							if (same_price_orders.size() == 0)
								resting_sell->erase(iter);
							return;
						}
						else
						{
							iter2 = same_price_orders.erase(iter2);
							resting_order_ID->erase(resting_order_id);
							quantity -= resting_quantity;

							print_trade_event(o_stream, resting_quantity, resting_lowest_price);
							print_orderpartiallyfilled(o_stream, orderid, quantity);
							print_orderfullyfilled(o_stream, resting_order_id);
						}
					}
				}

				if (iter->second.size() == 0)
					iter = resting_sell->erase(iter);
				else
					iter++;
			}
		}
	}
	else if (side == 1) //Sell
	{
		if (resting_buy->size() == 0)
			return;
		else
		{
			vector <real> to_be_erased;
			to_be_erased.reserve(100);

			auto iter = resting_buy->rbegin();

			while (iter != resting_buy->rend())
			{
				auto resting_highest_price = iter->first;

				if (resting_highest_price < price)
				{
					return;
				}
				else
				{
					auto& same_price_orders = iter->second;

					auto iter2 = same_price_orders.begin();

					while (iter2 != same_price_orders.end())
					{
						auto resting_order_id = iter2->first;
						auto& resting_quantity = iter2->second;

						if (resting_quantity > quantity)
						{
							resting_quantity -= quantity;

							print_trade_event(o_stream, quantity, resting_highest_price);
							print_orderfullyfilled(o_stream, orderid);
							print_orderpartiallyfilled(o_stream, resting_order_id, resting_quantity);

							quantity = 0;

							for (auto iterr = to_be_erased.begin(); iterr != to_be_erased.end(); iterr++)
							{
								resting_buy->erase(*iterr);
							}

							return;
						}
						else if (resting_quantity == quantity)
						{
							same_price_orders.erase(iter2);

							resting_order_ID->erase(resting_order_id);
							print_trade_event(o_stream, quantity, resting_highest_price);
							print_orderfullyfilled(o_stream, orderid);
							print_orderfullyfilled(o_stream, resting_order_id);

							quantity = 0;

							for (auto iterr = to_be_erased.begin(); iterr != to_be_erased.end(); iterr++)
							{
								resting_buy->erase(*iterr);
							}

							if (same_price_orders.size() == 0)
								resting_buy->erase(price);

							return;
						}
						else
						{
							iter2 = same_price_orders.erase(iter2);
							resting_order_ID->erase(resting_order_id);
							quantity -= resting_quantity;

							print_trade_event(o_stream, resting_quantity, resting_highest_price);
							print_orderpartiallyfilled(o_stream, orderid, quantity);
							print_orderfullyfilled(o_stream, resting_order_id);
						}
					}
				}

				if (iter->second.size() == 0)
					to_be_erased.push_back(resting_highest_price);

				iter++;
			}

			for (auto iterr = to_be_erased.begin(); iterr != to_be_erased.end(); iterr++)
			{
				resting_buy->erase(*iterr);
			}
		}
	}
}

void Matching_Engine::insert_order_into_book(tuple<uinteger, uinteger, uinteger, ulinteger, real>& new_order)
{
	//msgtype, orderid, side, quantity, price
	auto msgtype = get<0>(new_order);
	auto orderid = get<1>(new_order);

	if (msgtype == 0) //If msgtype == 0, add order request
	{
		if (resting_order_ID -> find(orderid) != resting_order_ID -> end())
		{
			print_bad_message(1, o_stream_err);
			return;
		}

		auto side = get<2>(new_order);
		auto quantity = get<3>(new_order);
		auto price = get<4>(new_order);

		if (side == 0) //Buy
		{
			make_a_deal(msgtype, orderid, side, quantity, price);

			if (quantity > 0)
				book_keeping(0, price, orderid, quantity);
		}
		else if (side == 1) //Sell
		{
			make_a_deal(msgtype, orderid, side, quantity, price);

			if (quantity > 0)
				book_keeping(1, price, orderid, quantity);
		}
	}
	else if (msgtype == 1) //If msgtype == 1, cancel order request
	{
		cancel_order(orderid);
	}
}

void Matching_Engine::process_order(string& tmp, istream& in_stream, ostream& o_stream, ostream& o_stream_err)
{
	tuple<uinteger, uinteger, uinteger, ulinteger, real> new_order;
	bool success = separate_order_items(tmp, new_order, o_stream_err);
	
	if (success) insert_order_into_book(new_order);
}

bool Matching_Engine::separate_order_items(string& tmp, tuple<uinteger, uinteger, uinteger, ulinteger, real>& new_order, ostream& o_stream)
{
	stringstream ss(tmp);

	string token;

	try
	{
		uinteger count = 0;
		uinteger msgtype = 0;
		uinteger orderid;
		uinteger side;

		//If the quantity of the order is larger than 2^31 - 1, define the macro "LARGE_QUANTITY" either in 'Headers.h' or on compile line.
		ulinteger quantity;

		real price;

		while (getline(ss, token, ','))
		{
			count++;

			switch (count)
			{
				case 1: 
				{
					msgtype = stoi(token);
					if (msgtype != 0 && msgtype != 1) 
					{
						print_bad_message(2, o_stream);
						return false;
					}
					break;
				}
				case 2: 
				{
					integer orderid_temp;

					orderid_temp = stoi(token);
					if (orderid_temp <= 0)
					{
						print_bad_message(3, o_stream);
						return false;
					}

					orderid = orderid_temp;

					break;
				}
				case 3:
				{
					side = stoi(token);

					if (side != 0 && side != 1)
					{
						print_bad_message(4, o_stream);
						return false;
					}
					break;
				}
				case 4:
				{
					linteger quantity_temp;

#ifdef LARGE_QUANTITY
					quantity_temp = sizeof(long) == 4 ? stoll(token) : stol(token);
#else
					quantity_temp = stoi(token);
#endif
					if (quantity_temp <= 0)
					{
						print_bad_message(5, o_stream);
						return false;
					}

					quantity = quantity_temp;

					break;
				}
				case 5:
				{
					price = stod(token);

					//if (price < 0.0)
					//{
					//	print_bad_message(9, o_stream);
					//	return false;
					//}

					break;
				}
				default: break;
			}
		}

		if (count != 2 && count != 5)
		{
			print_bad_message(6, o_stream);
			return false;
		}
		else if (msgtype == 0 && count != 5)
		{
			print_bad_message(7, o_stream);
			return false;
		}
		else if (msgtype == 1 && count != 2)
		{
			print_bad_message(8, o_stream);
			return false;
		}
		else
		{
			if (count == 2)
				new_order = make_tuple(msgtype, orderid, 2, 0, 0);
			else if (count == 5)
				new_order = make_tuple(msgtype, orderid, side, quantity, price);
		}
	}
	catch (exception& exp)
	{
		(void) exp;
		print_bad_message(0, o_stream);
		return false;
	}

	return true;
}

void Matching_Engine::print_bad_message(const unsigned type, ostream& o_stream) const
{ 
	if (type == 0)
		o_stream << "Unknown message type: BADMESSAGE\n";
	else if (type == 1)
		o_stream << "Error: the order ID has alreay been in the book.\n";
	else if (type == 2)
		o_stream << "Error: the input message type is neither 0 or 1.\n";
	else if (type == 3)
		o_stream << "Error: the order ID is smaller than or equal to 0.\n";
	else if (type == 4)
		o_stream << "Error: the order side is neither buy or sell.\n";
	else if (type == 5)
		o_stream << "Error: the order quantity is smaller than or equal to 0.\n";
	else if (type == 6)
		o_stream << "Error: the order information is neither of 2 or 5 items.\n";
	else if (type == 7)
		o_stream << "Error: the order for message type 0 is not of 5 items.\n";
	else if (type == 8)
		o_stream << "Error: the order for message type 1 is not of 2 items.\n";
	else if (type == 9)
		o_stream << "Error: the order price is less than 0.0.\n";
}

