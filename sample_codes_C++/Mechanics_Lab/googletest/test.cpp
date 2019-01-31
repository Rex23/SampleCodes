#include "Utilities.h"
#include "gtest/gtest.h"
#include <vector>

using namespace std;

namespace
{
	TEST(ObtainPolygonArea, Area)
	{
		New_Class_Utilities::Class_Utilities a_class;

		vector <vector <double> > Points{{1,1,1},{1,1,1},{1,1,1}};

		EXPECT_EQ(0.0, a_class.ObtainPolygonArea(Points));
	}
}


