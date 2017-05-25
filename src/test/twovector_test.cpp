#include <stdexcept>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <ofstream>
#include "geogeometry.hpp"
#include "rapidjson/rapidjson.h"
#include "test_utilities.hpp"

using namespace rapidjson;
using namespace std;

#define TOL 0.0000001	// Tolerance for floating point comparisons

class VectorTest : public ::testing::Test {
public:
	VectorTest () {
		srand(time(NULL));
		logfile.open("twovector_test.log", std::fstream::out | std::fstream::app);
	}
	ofstream logfile;
}

TEST(VectorTest, DefaultConstructor) {
	TwoVector v;
	EXPECT_FALSE(v.isValid());
	EXPECT_EQ(v.x(), NAN);
	EXPECT_EQ(v.y(), NAN);
}

TEST(VectorTest, NumericalConstructor) {
	for (int i = 0; i < TESTLEN; i++) {
		double x = drand(-100000.0, 100000.0);
		double y = drand(-100000.0, 100000.0);
		TwoVector v(x,y);
		TwoVector u(y,x);
		EXPECT_TRUE(v.isValid());
		EXPECT_TRUE(u.isValid());
		EXPECT_EQ(v.x(), x);
		EXPECT_EQ(v.y(), y);
		EXPECT_EQ(u.x(), y);
		EXPECT_EQ(u.y(), x);
	}
}

TEST(VectorTest, JSON_Constructor) {
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double x = drand(-100000.0, 100000.0);
		double y = drand(-100000.0, 100000.0);
		buf << "{\"x\":" << to_string(x) << ",\"y\":"<< to_string(y) << "}" << endl;
		Document d;
		d.Parse(buf.c_str());
		TwoVector v(d);
		EXPECT_TRUE(v.isValid());
		EXPECT_EQ(v.x(), x);
		EXPECT_EQ(v.y(), y);
	}
}

TEST(VectorTest, JSON_Object) {
	for (int i = 0; i < TESTLEN; i++) {
		Document d;
		double x = drand(-100000.0, 100000.0);
		double y = drand(-100000.0, 100000.0);
		TwoVector v(x,y);
		EXPECT_TRUE(v.isValid());
		Value j = v.pack(d);
		EXPECT_TRUE(j.hasMember("x"));
		EXPECT_TRUE(j.hasMember("y"));
		Value *xptr = Pointer("/x").Get(j);
		Value *yptr = Pointer("/y").Get(j);
		EXPECT_FALSE(xptr == nullptr);
		EXPECT_FALSE(yptr == nullptr);
		EXPECT_TRUE(xptr->IsDouble());
		EXPECT_TRUE(yptr->IsDouble());
		TwoVector u(j);
		EXPECT_TRUE(u.isValid());
		EXPECT_EQ(v.x(), u.x());
		EXPECT_EQ(v.y(), u.y());
	}
}

TEST (VectorTest, AngleMagOutput) {
	TwoVector v { 1 , 0 };
	TwoVector s { 3 , 4 };
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.mag(), 5.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.angleRad(), 0.9272952, TOL));
	EXPECT_TRUE(toleranceEquals(s.angleDeg(), 53.1301023, TOL));
}

TEST (VectorTest, PolarInterpretation) {
	TwoVector u;
	TwoVector v { 1 , 0 };

	EXPECT_TRUE(v.isValid());
	EXPECT_FALSE(u.isValid());

	v.x(3.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 3.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 3.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	
	v.y(4);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 5.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 3.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 4.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.9272952, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 53.1301023, TOL));
	
	u.x(1.0);
	EXPECT_FALSE(u.isValid());
	EXPECT_TRUE(toleranceEquals(u.x(), 1.0, TOL));
	EXPECT_TRUE(isnan(u.y()));
	EXPECT_TRUE(isnan(u.mag()));
	EXPECT_TRUE(isnan(u.angleRad()));
	EXPECT_TRUE(isnan(u.angleDeg()));
	
	u.y(1.0);
	EXPECT_TRUE(u.isValid());
	EXPECT_TRUE(toleranceEquals(u.x(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(u.y(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(u.mag(), sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(u.angleRad(), M_PI_4, TOL));
	EXPECT_TRUE(toleranceEquals(u.angleDeg(), 45.0, TOL));
}

TEST (VectorTest, AngleMagInput) {
	TwoVector v { 1 , 0 };
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	
	v.mag(5.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 5.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	
	v.angleRad(atan2(4, 3));
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 5.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 3.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 4.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), 0.9272952, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 53.1301023, TOL));
	
	v.angleRad(atan2(3, 4));
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 5.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 4.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 3.0, TOL));
}

TEST (VectorTest, GetVectorRad) {
	TwoVector v = TwoVector::getVectorRad(M_PI_4, 2.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleRad(), M_PI_4, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 45.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), sqrt(2.0), TOL));
}

TEST (VectorTest, VectorRotations) {
	TwoVector v {2, 0};
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 45.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), sqrt(2.0), TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 90.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 2.0, TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 135.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), sqrt(2.0), TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 180.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), -135.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), -sqrt(2.0), TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), -90.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), -2.0, TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), -45.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), sqrt(2.0), TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), -sqrt(2.0), TOL));
	v.rotateDeg(45.0);
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.angleDeg(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));
}

TEST (VectorTest, VectorAddition) {
	TwoVector u {2, 0};
	TwoVector v {0, 1};
	TwoVector s;

	EXPECT_TRUE(u.isValid());
	EXPECT_TRUE(v.isValid());
	EXPECT_FALSE(s.isValid());
	
	u += v;
	EXPECT_TRUE(u.isValid());
	EXPECT_TRUE(toleranceEquals(u.x(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(u.y(), 1.0, TOL));
	
	v -= u;
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.x(), -2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));	
	
	s = u + v;
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(toleranceEquals(u.x(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(u.y(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));	
	EXPECT_TRUE(toleranceEquals(s.x(), 0.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.y(), 1.0, TOL));	
	
	s = u - v;
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(toleranceEquals(u.x(), 2.0, TOL));
	EXPECT_TRUE(toleranceEquals(u.y(), 1.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -2.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 0.0, TOL));	
	EXPECT_TRUE(toleranceEquals(s.x(), 4.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.y(), 1.0, TOL));	
}

TEST (VectorTest, ScalarMultiplication) {
	TwoVector s;
	TwoVector v {3, 4};
	EXPECT_TRUE(v.isValid());
	EXPECT_FALSE(s.isValid());
	
	v *= 3;
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 15.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 9.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 12.0, TOL));	
	
	v *= -1;
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 15.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), -9.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), -12.0, TOL));

	v *= -(1.0/2.0); 
	EXPECT_TRUE(v.isValid());
	EXPECT_TRUE(toleranceEquals(v.x(), 4.5, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 6.0, TOL));
	EXPECT_TRUE(toleranceEquals(v.mag(), 7.5, TOL));
	
	s = v * 2;
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 7.5, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 4.5, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 6.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.mag(), 15.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.x(), 9.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.y(), 12.0, TOL));
	
	s = v / 2;
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(toleranceEquals(v.mag(), 7.5, TOL));
	EXPECT_TRUE(toleranceEquals(v.x(), 4.5, TOL));
	EXPECT_TRUE(toleranceEquals(v.y(), 6.0, TOL));
	EXPECT_TRUE(toleranceEquals(s.mag(), 3.75, TOL));
	EXPECT_TRUE(toleranceEquals(s.x(), 2.25, TOL));
	EXPECT_TRUE(toleranceEquals(s.y(), 3, TOL));
}

TEST (VectorTest, DotProduct) {
	TwoVector s {1, 0};
	TwoVector u {0, 1};
	TwoVector v {3, 4};
	double w;
	
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(u.isValid());
	EXPECT_TRUE(v.isValid());

	w  = v * s;
	EXPECT_TRUE(toleranceEquals(w, 3, TOL));
	
	w  = v * u;
	EXPECT_TRUE(toleranceEquals(w, 4, TOL));

	w = s * u;
	EXPECT_TRUE(toleranceEquals(w, 0, TOL));
}

TEST (VectorTest, PerProduct) {
	TwoVector s {1, 0};
	TwoVector u {0, 1};
	TwoVector v {3, 4};
	
	EXPECT_TRUE(s.isValid());
	EXPECT_TRUE(u.isValid());
	EXPECT_TRUE(v.isValid());

	EXPECT_TRUE(toleranceEquals(s.perp(u), 0, TOL));
	EXPECT_TRUE(toleranceEquals(s.perp(v), 3, TOL));
	EXPECT_TRUE(toleranceEquals(u.perp(s), 0, TOL));
	EXPECT_TRUE(toleranceEquals(u.perp(v), -4, TOL));
	EXPECT_TRUE(toleranceEquals(v.perp(s), 3, TOL));
	EXPECT_TRUE(toleranceEquals(v.perp(u), -4, TOL));
}
