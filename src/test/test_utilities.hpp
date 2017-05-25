
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>

#ifndef TEST_UTILITIES
#define TEST_UTILITIES

/*testing::AssertionResult assertJSONEqual(const char *expected_expr, 
										const char *actual_expr,
					 					Value& expected, Value& actual);

#define EXPECT_JSON_EQ(a, b) EXPECT_PRED_FORMAT2(assertJSONEqual, a, b)
#define ASSERT_JSON_EQ(a, b) ASSERT_PRED_FORMAT2(assertJSONEqual, a, b)*/


bool inline toleranceEquals (double a, double b, double tolerance) {	/**< Compare two doubles with a tolerance */
	double diff = a - b;
	if (fabs(diff) < tolerance) {
		return true;
	} else {
		std::cerr << a << " " << b << " " << diff << std::endl;
		return false;
	}
}

double inline drand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double inline test_max(double a, double b) {
	if (a > b) return a;
	return b;
}

double inline test_min(double a, double b) {
	if (a < b) return a;
	return b;
}

double inline nm2meters(double a) {return (a/1852);}
double inline meters2nm(double a) {return (a*1852);}

#define TESTLEN	(1000)

#endif /* TEST_UTILITIES */