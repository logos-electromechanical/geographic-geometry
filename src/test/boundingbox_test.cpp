#include <stdexcept>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <ofstream>
#include <utility>
#include "geogeometry.hpp"
#include "rapidjson/rapidjson.h"
#include "test_utilities.hpp"

using namespace rapidjson;
using namespace std;

#define TOL (0.00001)	// Tolerance for floating point comparisons

class BoundingBoxTest : public ::testing::Test {
public:
	BoundingBoxTest () {
		srand(time(NULL));
		logfile.open("position_test.log", std::fstream::out | std::fstream::app);
	}

	BoundingBox makeRandBox (double minlon, double maxlon, double minlat, double maxlat) {
		Position sw(drand(minlon, maxlon), drand(minlat, maxlat));
		Position ne(drand(getlon(sw), maxlon), drand(getlat(sw), maxlat));
		return makeBBox(sw,ne);
	}

	ofstream logfile;
}

TEST (BoundingBoxTest, PositionsValid) {
	BoundingBox b;
	EXPECT_FALSE(bboxValid(b));
	EXPECT_FALSE(PositionsValid(getne(b)));
	EXPECT_FALSE(PositionsValid(getsw(b)));
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-198.0,198.0,-99.0,99.0);
		if ((getlon(getne(bb)) < -180.0) || (getlon(getne(bb)) > 180.0) || 
			(getlat(getne(bb)) < -90.0) || (getlat(getne(bb)) > 90.0) || 
			getsw(getlon(getsw(bb)) < -180.0) || (getlon(getsw(bb)) > 180.0) || 
			(getlat(getsw(bb)) < -90.0) || (getlat(getsw(bb)) > 90.0)) {
			EXPECT_FALSE(bboxValid(bb));
		} else {
			EXPECT_TRUE(bboxValid(bb));
		}
	}
}

TEST (BoundingBoxTest, OrderInvalid) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-180.0,180.0,-90.0,90.0);
		swap(getne(bb), getsw(bb));
		EXPECT_FALSE(bboxValid(bb));
	}
}

TEST (BoundingBoxTest, JSON_Valid) {
	for (i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,90.0);
		Document d0, d1;
		Value v = packBBox(bb, d0);
		d0.Accept(writer);
		char* out = buffer.GetString();
		d1.Parse(out);
		BoundingBox bb1 = loadBBox(d1);
		EXPECT_TRUE(bboxValid(bb0));
		EXPECT_TRUE(bboxValid(bb1));
		EXPECT_TRUE(toleranceEquals(getlon(getne(bb0)), getlon(getne(bb1)), TOL));
		EXPECT_TRUE(toleranceEquals(getlat(getne(bb0)), getlat(getne(bb1)), TOL));
		EXPECT_TRUE(toleranceEquals(getlon(getsw(bb0)), getlon(getsw(bb1)), TOL));
		EXPECT_TRUE(toleranceEquals(getlat(getsw(bb0)), getlat(getsw(bb1)), TOL));
	}
}

TEST (BoundingBoxTest, Contains) {

}
