#include <stdexcept>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <utility>
#include "geogeometry.hpp"
#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/writer.h"
#include "test_utilities.hpp"

using namespace rapidjson;
using namespace std;
using namespace GeoGeometry;

#define TOL (0.00001)	// Tolerance for floating point comparisons

class BoundingBoxTest : public ::testing::Test {
public:
	BoundingBoxTest () {
		srand(time(NULL));
		logfile.open("position_test.log", std::fstream::out | std::fstream::app);
	}

	BoundingBox makeRandBox (double minlon, double maxlon, double minlat, double maxlat) {
		Position sw = makePosition(drand(minlon, maxlon), drand(minlat, maxlat));
		Position ne = makePosition(drand(getlon(sw), maxlon), drand(getlat(sw), maxlat));
		return makeBBox(sw,ne);
	}

	ofstream logfile;
};

TEST_F (BoundingBoxTest, PositionsValid) {
	BoundingBox b;
	EXPECT_FALSE(bboxValid(b));
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-198.0,198.0,-99.0,99.0);
		if ((getlon(getne(bb)) < -180.0) || (getlon(getne(bb)) > 180.0) || 
			(getlat(getne(bb)) < -90.0) || (getlat(getne(bb)) > 90.0) || 
			(getlon(getsw(bb)) < -180.0) || (getlon(getsw(bb)) > 180.0) || 
			(getlat(getsw(bb)) < -90.0) || (getlat(getsw(bb)) > 90.0)) {
			EXPECT_FALSE(bboxValid(bb));
		} else {
			EXPECT_TRUE(bboxValid(bb));
		}
	}
}

TEST_F (BoundingBoxTest, OrderInvalid) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-180.0,180.0,-90.0,90.0);
		BoundingBox bb0 = bb;
		swap(getne(bb), getsw(bb));						// out of order
		swap(getlon(getne(bb0)), getlon(getsw(bb0)));	// test the antimeridian crossing case
		EXPECT_FALSE(bboxValid(bb));
		EXPECT_TRUE(bboxValid(bb0));
	}
}

TEST_F (BoundingBoxTest, JSON_Valid) {
	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,90.0);
		Document d0, d1;
		Value v = packBBox(bb0, d0);
		v.Accept(writer);
		const char* out = buf.GetString();
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

TEST_F (BoundingBoxTest, Contains) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,90.0);

		// construct a point within the box
		Position p = makePosition(drand(getlon(getsw(bb0)), getlon(getne(bb0))),
			drand(getlat(getsw(bb0)), getlat(getne(bb0))));
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(bboxContains(bb0, p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainNorth) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,89.0);

		// construct a point north of the box
		Position p = makePosition(drand(-180.0,180.0), drand(getlat(getne(bb0)), 90.0));
		EXPECT_TRUE(positionValid(p));
		EXPECT_FALSE(bboxContains(bb0, p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainSouth) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-89.0,90.0);

		// construct a point south of the box
		Position p = makePosition(drand(-180.0,180.0), drand(-90.0, getlat(getsw(bb0))));
		EXPECT_TRUE(positionValid(p));
		EXPECT_FALSE(bboxContains(bb0, p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainEast) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,179.0,-90.0,90.0);

		// construct a point south of the box
		Position p = makePosition(drand(getlon(getne(bb0)),180.0), drand(-90.0, 90.0));
		EXPECT_TRUE(positionValid(p));
		EXPECT_FALSE(bboxContains(bb0, p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainWest) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-179.0,180.0,-90.0,90.0);

		// construct a point south of the box
		Position p = makePosition(drand(-180, getlon(getsw(bb0))), drand(-90.0, 90.0));
		EXPECT_TRUE(positionValid(p));
		EXPECT_FALSE(bboxContains(bb0, p));
	}
}

TEST_F (BoundingBoxTest, Antimeridian) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-179.0,179.0,-90.0,90.0);
		swap(getlon(getsw(bb)), getlon(getne(bb)));
		EXPECT_TRUE(bboxValid(bb));
		// test what would be in the box if it did not cross the meridian
		Position p = makePosition(drand(getlon(getsw(bb)), getlon(getne(bb))), 
			drand(getlat(getsw(bb)), getlat(getne(bb))));
		EXPECT_TRUE(positionValid(p));
		EXPECT_FALSE(bboxContains(bb, p));
		// test a point in east half of the box
		p = makePosition(drand(-180.0, getlon(getne(bb))), 
			drand(getlat(getsw(bb)), getlat(getne(bb))));
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(bboxContains(bb, p));
		// test a point in west half of the box
		p = makePosition(drand(getlon(getsw(bb)), 180.0), 
			drand(getlat(getsw(bb)), getlat(getne(bb))));
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(bboxContains(bb, p));
	}
}



	

