#include <algorithm>

#include "src/common/point.h"
#include "gtest/gtest.h"

using namespace std;

TEST(point_test, constructor_test)
{
    Point::setDimension(2);
    Point p1("1 1");
    Point p2("2 2 test");

    EXPECT_EQ((unsigned int)2, Point::dimension());
    EXPECT_EQ("1 1", p1.toString());
    EXPECT_EQ("2 2 test", p2.toString());
}
