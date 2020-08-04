#include <iostream>
#include <thread>
#include <string>
#include <algorithm>
#include<vector>
#include<mutex>
#include <condition_variable>
#include <future>
#include <chrono>

using namespace std;

//Example:
void a_fun(string& a_str)
{
	cout << "function is called.\n";
	cout << this_thread::get_id() << endl; //this_tread::get_id()
}

int example_hardware_concurrency()
{
	string a_str("dasf");

	thread t1(a_fun, ref(a_str));
	thread t2(a_fun, ref(a_str));

	t1.join();
	t2.join();

	cout << "finished.\n";

	cout << std::thread::hardware_concurrency() << endl; //Don’t forget the ()!

	return 0;
}

//Example:

void thread_function()
{
	std::cout << "Inside Thread :: ID  = " << std::this_thread::get_id() << std::endl;
}

int example_get_id()
{
	std::thread threadObj1(thread_function);
	std::thread threadObj2(thread_function);

	if (threadObj1.get_id() != threadObj2.get_id())
		std::cout << "Both Threads have different IDs" << std::endl;

	std::cout << "From Main Thread :: ID of Thread 1 = " << threadObj1.get_id() << std::endl;
	std::cout << "From Main Thread :: ID of Thread 2 = " << threadObj2.get_id() << std::endl;

	threadObj1.join();
	threadObj2.join();
	return 0;
}

//Example:

class WorkerThread
{
public:
	void operator()()
	{
		std::cout << "Worker Thread " << std::this_thread::get_id() << " is Executing" << std::endl;
	}
};
int example_for_each()
{
	std::vector<std::thread> threadList;
	for (int i = 0; i < 10; i++)
	{
		threadList.push_back(std::thread(WorkerThread()));
	}
	// Now wait for all the worker thread to finish i.e.
	// Call join() function on each of the std::thread object
	std::cout << "wait for all the worker thread to finish" << std::endl;
	std::for_each(threadList.begin(), threadList.end(), std::mem_fn(&std::thread::join));
	std::cout << "Exiting from Main Thread" << std::endl;
	return 0;
}

//Example:
class ThreadRAII
{
	std::thread & m_thread;
public:
	ThreadRAII(std::thread  & threadObj) : m_thread(threadObj)
	{

	}
	~ThreadRAII()
	{
		// Check if thread is joinable then detach the thread
		if (m_thread.joinable())
		{
			m_thread.detach();
		}
	}
};
void thread_function2()
{
	for (int i = 0; i < 10000; i++);
	std::cout << "thread_function Executing" << std::endl;
}

int example_RAII()
{
	std::thread threadObj(thread_function2);

	// If we comment this Line, then program will crash
	ThreadRAII wrapperObj(threadObj);
	return 0;
}

//Example:
void newThreadCallback(int * p)
{
	std::cout << "Inside Thread :  "" : p = " << p << std::endl;
	std::chrono::milliseconds dura(1000);
	std::this_thread::sleep_for(dura);
	*p = 19;
}
void startNewThread()
{
	int * p = new int();
	*p = 10;
	std::cout << "Inside Main Thread :  "" : *p = " << *p << std::endl;
	std::thread t(newThreadCallback, p);
	t.detach();
	delete p;
	p = NULL;
}
int example_pass_pointer_chrono_this_thread_sleep_for()
{
	startNewThread();
	std::chrono::milliseconds dura(2000);
	std::this_thread::sleep_for(dura);
	return 0;
}

//Example:
void threadCallback(int const & x)
{
	int & y = const_cast<int &>(x);
	y++;
	std::cout << "Inside Thread x = " << x << std::endl;
}
int example_pass_reference()
{
	int x = 9;
	std::cout << "In Main Thread : Before Thread Start x = " << x << std::endl;
	std::thread threadObj(threadCallback, std::ref(x));
	threadObj.join();
	std::cout << "In Main Thread : After Thread Joins x = " << x << std::endl;
	return 0;
}

//Example:
class DummyClass {
public:
	DummyClass()
	{}
	DummyClass(const DummyClass & obj)
	{}
	void sampleMemberFunction(int x)
	{
		std::cout << "Inside sampleMemberFunction " << x << std::endl;
	}
};

int example_call_class_method() {
	DummyClass dummyObj;
	int x = 10;
	std::thread threadObj(&DummyClass::sampleMemberFunction, &dummyObj, x);
	threadObj.join();
	return 0;
}

//Example:
class Wallet
{
	int mMoney;
	std::mutex mutex;
public:
	Wallet() :mMoney(0) {}
	int getMoney() { return mMoney; }
	void addMoney(int money)
	{
		mutex.lock();
		for (int i = 0; i < money; ++i)
		{
			mMoney++;
		}
		mutex.unlock();
	}
};
int testMultithreadedWallet()
{
	Wallet walletObject;
	std::vector<std::thread> threads;
	for (int i = 0; i < 5; ++i) {
		threads.push_back(std::thread(&Wallet::addMoney, &walletObject, 1000));
	}

	for (int i = 0; i < threads.size(); i++)
	{
		threads.at(i).join();
	}
	return walletObject.getMoney();
}

int example_push_back_mutex()
{

	int val = 0;
	for (int k = 0; k < 1000; k++)
	{
		if ((val = testMultithreadedWallet()) != 5000)
		{
			std::cout << "Error at count = " << k << "  Money in Wallet = " << val << std::endl;
			//break;
		}
	}
	return 0;
}

//Example:
class Wallet2
{
	int mMoney;
	std::mutex mutex;
public:
	Wallet2() :mMoney(0) {}
	int getMoney() { return mMoney; }
	void addMoney(int money)
	{
		std::lock_guard<std::mutex> lockGuard(mutex);
		// In constructor it locks the mutex

		for (int i = 0; i < money; ++i)
		{
			// If some exception occurs at this
			// poin then destructor of lockGuard
			// will be called due to stack unwinding.
			//
			mMoney++;
		}
		// Once function exits, then destructor
		// of lockGuard Object will be called.
		// In destructor it unlocks the mutex.
	}
};

int example_lock_guard()
{
	Wallet2 walletObject;

	std::vector<std::thread> threads;
	for (int i = 0; i < 5; ++i) {
		threads.push_back(std::thread(&Wallet2::addMoney, &walletObject, 1000));
	}

	for (int i = 0; i < threads.size(); i++)
	{
		threads.at(i).join();
	}

	return walletObject.getMoney();
}

//Example:
class Application
{
	std::mutex m_mutex;
	std::condition_variable m_condVar;
	bool m_bDataLoaded;

public:

	Application()
	{
		m_bDataLoaded = false;
	}
	void loadData()
	{
		// Make This Thread sleep for 1 Second
		std::this_thread::sleep_for(std::chrono::milliseconds(1000));
		std::cout << "Loading Data from XML" << std::endl;
		// Lock The Data structure
		std::lock_guard<std::mutex> guard(m_mutex);
		// Set the flag to true, means data is loaded
		m_bDataLoaded = true;
		// Notify the condition variable
		m_condVar.notify_one();
	}

	bool isDataLoaded()
	{
		return m_bDataLoaded;
	}

	void mainTask()
	{
		std::cout << "Do Some Handshaking" << std::endl;
		// Acquire the lock
		std::unique_lock<std::mutex> mlock(m_mutex);
		// Start waiting for the Condition Variable to get signaled
		// Wait() will internally release the lock and make the thread to block
		// As soon as condition variable get signaled, resume the thread and
		// again acquire the lock. Then check if condition is met or not
		// If condition is met then continue else again go in wait.
		m_condVar.wait(mlock, std::bind(&Application::isDataLoaded, this));
		std::cout << "Do Processing On loaded Data" << std::endl;
	}

};

int example_condition_variable()
{
	Application app;
	std::thread thread_1(&Application::mainTask, &app);
	std::thread thread_2(&Application::loadData, &app);
	thread_2.join();
	thread_1.join();
	std::cout << "Finished.\n";

	return 0;
}

//Example: Many times we encounter a situation where we want a thread to return a result
void initiazer(std::promise<int> * promObj)
{
	std::cout << "Inside Thread" << std::endl;     promObj->set_value(35);
}

int example_promise_future()
{
	std::promise<int> promiseObj;
	std::future<int> futureObj = promiseObj.get_future();
	std::thread th(initiazer, &promiseObj);
	std::cout << futureObj.get() << std::endl;
	th.join();
	return 0;
}

//Example:
using namespace std::chrono;

std::string fetchDataFromDB(std::string recvdData)
{
	// Make sure that function takes 5 seconds to complete
	std::this_thread::sleep_for(seconds(5));

	//Do stuff like creating DB Connection and fetching Data
	return "DB_" + recvdData;
}

std::string fetchDataFromFile(std::string recvdData)
{
	// Make sure that function takes 5 seconds to complete
	std::this_thread::sleep_for(seconds(5));

	//Do stuff like fetching Data File
	return "File_" + recvdData;
}

int example_async_system_clock()
{
	// Get Start Time
	system_clock::time_point start = system_clock::now();

	std::future<std::string> resultFromDB = std::async(std::launch::async, fetchDataFromDB, "Data");

	//Fetch Data from File
	std::string fileData = fetchDataFromFile("Data");

	//Fetch Data from DB
	// Will block till data is available in future<std::string> object.
	std::string dbData = resultFromDB.get();

	// Get End Time
	auto end = system_clock::now();

	auto diff = duration_cast < std::chrono::seconds > (end - start).count();
	std::cout << "Total Time Taken = " << diff << " Seconds" << std::endl;

	//Combine The Data
	std::string data = dbData + " :: " + fileData;

	//Printing the combined Data
	std::cout << "Data = " << data << std::endl;

	return 0;
}

//Example:
// Fetch some data from DB
std::string getDataFromDB(std::string token)
{
	// Do some stuff to fetch the data
	std::string data = "Data fetched from DB by Filter :: " + token;
	return data;
}

int example_packaged_task()
{
	// Create a packaged_task<> that encapsulated the callback i.e. a function
	std::packaged_task<std::string(std::string)> task(getDataFromDB);

	// Fetch the associated future<> from packaged_task<>
	std::future<std::string> result = task.get_future();

	// Pass the packaged_task to thread to run asynchronously
	std::thread th(std::move(task), "Arg");

	// Join the thread. Its blocking and returns when thread is finished.
	th.join();

	// Fetch the result of packaged_task<> i.e. value returned by getDataFromDB()
	std::string data = result.get();

	std::cout << data << std::endl;

	return 0;
}

//Example:
int example_packaged_task_lambda()
{
	// Create a packaged_task<> that encapsulated a lambda function
	std::packaged_task<std::string(std::string)> task([](std::string token) {
		// Do some stuff to fetch the data
		std::string data = "Data From " + token;
		return data;
	});

	// Fetch the associated future<> from packaged_task<>
	std::future<std::string> result = task.get_future();

	// Pass the packaged_task to thread to run asynchronously
	std::thread th(std::move(task), "Arg");

	// Join the thread. Its blocking and returns when thread is finished.
	th.join();

	// Fetch the result of packaged_task<> i.e. value returned by getDataFromDB()
	std::string data = result.get();

	std::cout << data << std::endl;

	return 0;
}

void test_thread()
{
	int case_index = 14;

	switch (case_index)
	{
	case 1:
	{
		cout << "Example 1:\n";
		example_hardware_concurrency();
		break;
	}
	case 2:
	{
		cout << "Example 2:\n";
		example_get_id();
		break;
	}
	case 3:
	{
		cout << "Example 3:\n";
		example_for_each();
		break;
	}
	case 4:
	{
		cout << "Example 4:\n";
		example_RAII();
		break;
	}
	case 5:
	{
		cout << "Example 5:\n";
		example_pass_pointer_chrono_this_thread_sleep_for();
		break;
	}
	case 6:
	{
		cout << "Example 6:\n";
		example_pass_reference();
		break;
	}
	case 7:
	{
		cout << "Example 7:\n";
		example_call_class_method();
		break;
	}
	case 8:
	{
		cout << "Example 8:\n";
		example_push_back_mutex();
		break;
	}
	case 9:
	{
		cout << "Example 9:\n";
		example_lock_guard();
	}
	case 10:
	{
		cout << "Example 10:\n";
		example_condition_variable();
	}
	case 11:
	{
		cout << "Example 11:\n";
		example_promise_future();
	}
	case 12:
	{
		cout << "Example 12:\n";
		example_async_system_clock();
	}
	case 13:
	{
		cout << "Example 13:\n";
		example_packaged_task();
	}
	case 14:
	{
		cout << "Example 14:\n";
		example_packaged_task_lambda();
	}
	default: break;
	}
}













