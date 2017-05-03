#include <algorithm>

#include "src/core/point.h"
#include "gtest/gtest.h"

using namespace std;

TEST(point_test, constructor_test)
{
    Point::setDimension(2);
    Point p("1 1 test");

    EXPECT_EQ(2,Point::dimension());
    EXPECT_EQ("1 1  test", p.toString());
}
