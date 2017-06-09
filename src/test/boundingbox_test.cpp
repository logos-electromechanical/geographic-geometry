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
		Position sw(drand(minlon, maxlon), drand(minlat, maxlat));
		Position ne(drand(sw.lon(), maxlon), drand(sw.lat(), maxlat));
		return BoundingBox(ne,sw);
	}

	ofstream logfile;
};

TEST_F (BoundingBoxTest, PositionsValid) {
	BoundingBox b;
	EXPECT_FALSE(b.isValid());
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-198.0,198.0,-99.0,99.0);
		if ((bb.ne().lon() < -180.0) || (bb.ne().lon() > 180.0) || 
			(bb.ne().lat() < -90.0) || (bb.ne().lat() > 90.0) || 
			(bb.sw().lon() < -180.0) || (bb.sw().lon() > 180.0) || 
			(bb.sw().lat() < -90.0) || (bb.sw().lat() > 90.0)) {
			EXPECT_FALSE(bb.isValid());
		} else {
			EXPECT_TRUE(bb.isValid());
		}
	}
}

TEST_F (BoundingBoxTest, OrderInvalid) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-180.0,180.0,-90.0,90.0);
		BoundingBox bb0 = bb;
		swap(bb.ne(), bb.sw());						// out of order
		swap(bb0.ne().lon(), bb0.sw().lon());	// test the antimeridian crossing case
		EXPECT_FALSE(bb.isValid());
		EXPECT_TRUE(bb0.isValid());
	}
}

TEST_F (BoundingBoxTest, JSON_Valid) {
	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,90.0);
		Document d0, d1;
		Value v = bb0.pack(d0);
		v.Accept(writer);
		const char* out = buf.GetString();
		d1.Parse(out);
		BoundingBox bb1(d1);
		EXPECT_TRUE(bb0.isValid());
		EXPECT_TRUE(bb1.isValid());
		EXPECT_TRUE(toleranceEquals(bb0.ne().lon(), bb1.ne().lon(), TOL));
		EXPECT_TRUE(toleranceEquals(bb0.ne().lat(), bb1.ne().lat(), TOL));
		EXPECT_TRUE(toleranceEquals(bb0.sw().lon(), bb1.sw().lon(), TOL));
		EXPECT_TRUE(toleranceEquals(bb0.sw().lat(), bb1.sw().lat(), TOL));
	}
}

TEST_F (BoundingBoxTest, Contains) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,90.0);

		// construct a point within the box
		Position p(drand(bb0.sw().lon(), bb0.ne().lon()),
			drand(bb0.sw().lat(), bb0.ne().lat()));
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(bb0.contains(p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainNorth) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-90.0,89.0);

		// construct a point north of the box
		Position p(drand(-180.0,180.0), drand(bb0.ne().lat(), 90.0));
		EXPECT_TRUE(p.isValid());
		EXPECT_FALSE(bb0.contains(p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainSouth) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,180.0,-89.0,90.0);

		// construct a point south of the box
		Position p(drand(-180.0,180.0), drand(-90.0, bb0.sw().lat()));
		EXPECT_TRUE(p.isValid());
		EXPECT_FALSE(bb0.contains(p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainEast) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-180.0,179.0,-90.0,90.0);

		// construct a point south of the box
		Position p(drand(bb0.ne().lon(),180.0), drand(-90.0, 90.0));
		EXPECT_TRUE(p.isValid());
		EXPECT_FALSE(bb0.contains(p));
	}
}

TEST_F (BoundingBoxTest, DoesNotContainWest) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb0 = makeRandBox(-179.0,180.0,-90.0,90.0);

		// construct a point south of the box
		Position p(drand(-180, bb0.sw().lon()), drand(-90.0, 90.0));
		EXPECT_TRUE(p.isValid());
		EXPECT_FALSE(bb0.contains(p));
	}
}

TEST_F (BoundingBoxTest, Antimeridian) {
	for (int i = 0; i < TESTLEN; i++) {
		BoundingBox bb = makeRandBox(-179.0,179.0,-90.0,90.0);
		swap(bb.sw().lon(), bb.ne().lon());
		EXPECT_TRUE(bb.isValid());
		// test what would be in the box if it did not cross the meridian
		Position p(drand(bb.sw().lon(), bb.ne().lon()), 
			drand(bb.sw().lat(), bb.ne().lat()));
		EXPECT_TRUE(p.isValid());
		EXPECT_FALSE(bb.contains(p));
		// test a point in east half of the box
		p = Position(drand(-180.0, bb.ne().lon()), 
			drand(bb.sw().lat(), bb.ne().lat()));
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(bb.contains(p));
		// test a point in west half of the box
		p = Position(drand(bb.sw().lon(), 180.0), 
			drand(bb.sw().lat(), bb.ne().lat()));
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(bb.contains(p));
	}
}



	

