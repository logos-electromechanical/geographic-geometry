#include <stdexcept>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
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

class PositionTest : public ::testing::Test {
	public:
		PositionTest () {
			srand(time(NULL));
			logfile.open("position_test.log", std::fstream::out | std::fstream::app);
		}

		ofstream logfile;
};

TEST_F (PositionTest, Validity) {
	logfile << "=================================================" << endl;
	logfile << "Testing validity for randomly generated positions" << endl;
	logfile << "Approximately 19 percent should be invalid" << endl;
	logfile << "=================================================" << endl;
	int i, invalid = 0;
	int invalid_detected = 0; 
	for (i = 0; i < TESTLEN; i++) {
		Position p(drand(-198.0, 198.0), drand(-99.0, 99.0), drand(-1000000,100000));
		logfile << "Test Position " << i << " : " << to_string(p.lon()) << ",";
		logfile << to_string(p.lat()) << "," << to_string(p.ele()); 
		logfile << "," << to_string(p.isValid()) << endl;
		if (!p.isValid()) invalid_detected++;
		if ((p.lon() < -180.0) || (p.lon() > 180.0) || 
			(p.lat() < -90.0) || (p.lat() > 90.0)) {
			EXPECT_FALSE(p.isValid());
			invalid++;
		} else {
			EXPECT_TRUE(p.isValid());
		}
	}
	for (i = 0; i < TESTLEN; i++) {
		Position p(drand(-198.0, 198.0), drand(-99.0, 99.0));
		logfile << "Test Position " << i << " : " << to_string(p.lon()) << ",";
		logfile << to_string(p.lat()) << "," << to_string(p.ele()); 
		logfile << "," << to_string(p.isValid()) << endl;
		if (!p.isValid()) invalid_detected++;
		if ((p.lon() < -180.0) || (p.lon() > 180.0) || 
			(p.lat() < -90.0) || (p.lat() > 90.0)) {
			EXPECT_FALSE(p.isValid());
			invalid++;
		} else {
			EXPECT_TRUE(p.isValid());
		}
	}
	logfile << "Total positions: " << to_string(i*2);
	logfile << " Invalid detected: " << to_string(invalid_detected);
	logfile << " Total invalid: " << to_string(invalid);
	EXPECT_FALSE(Position(NAN, 0, 0).isValid());
	EXPECT_FALSE(Position(0, NAN, 0).isValid());
	EXPECT_FALSE(Position(0, 0, NAN).isValid());
	EXPECT_FALSE(Position(NAN, 0).isValid());
	EXPECT_FALSE(Position(0, NAN).isValid());
}

TEST_F (PositionTest, JSON_Valid) {
	logfile << "=================================================" << endl;
	logfile << "Testing JSON pack/unpack" << endl;
	logfile << "=================================================" << endl;

	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		Position p(drand(-180.0, 180.0), drand(-90.0, 90.0), drand(-1000000,100000));
		logfile << "Test Position " << i << " : " << p.lon() << "," << p.lat() << "," << p.ele(); 
		logfile << "," << to_string(p.isValid()) << endl;
		Document d0, d1;
		Value v = p.pack(d0);
		v.Accept(writer);
		const char* out = buf.GetString();
		d1.Parse(out);
		Position p1;
		EXPECT_TRUE(p1.load(d1));
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(p1.isValid());
		EXPECT_TRUE(toleranceEquals(p.lon(), p1.lon(), TOL));
		EXPECT_TRUE(toleranceEquals(p.lat(), p1.lat(), TOL));
		EXPECT_TRUE(toleranceEquals(p.ele(), p1.ele(), TOL));
	}
	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		Position p(drand(-180.0, 180.0), drand(-90.0, 90.0));
		logfile << "Test Position " << i << " : " << p.lon() << "," << p.lat() << "," << p.ele(); 
		logfile << "," << to_string(p.isValid()) << endl;
		Document d0, d1;
		Value v = p.pack(d0);
		v.Accept(writer);
		const char* out = buf.GetString();
		d1.Parse(out);
		Position p1(d1);
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(p1.isValid());
		EXPECT_TRUE(toleranceEquals(p.lon(), p1.lon(), TOL));
		EXPECT_TRUE(toleranceEquals(p.lat(), p1.lat(), TOL));
		EXPECT_TRUE(toleranceEquals(p.ele(), p1.ele(), TOL));
	}
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double lon = drand(-180.0,180.0);
		double lat = drand(-90.0,90.0);
		buf << "[" << to_string(lon) << "," << to_string(lat) << "]";
		logfile << "Test Position " << i << " : " << buf.str() << endl;
		Document d;
		d.Parse(buf.str().c_str());
		Position p(d);
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(toleranceEquals(p.lon(), lon, TOL));
		EXPECT_TRUE(toleranceEquals(p.lat(), lat, TOL));
		EXPECT_TRUE(toleranceEquals(p.ele(), 0.0, TOL));
	}
}

TEST_F (PositionTest, IsNorthTest) {
	logfile << "=================================================" << endl;
	logfile << "Testing isNorth() function" << endl;
	logfile << "=================================================" << endl;
	for (int i = 0; i < TESTLEN; i++) {

		// construct a random point
		Position s0(drand(-180.0, 180.0), drand(-89.0, 89.0));
		EXPECT_TRUE(s0.isValid());

		// test for a point to the north of the first point
		Position p(drand(-180.0, 180.0), drand(s0.lat(), 90.0));
		EXPECT_TRUE(p.isValid());
		EXPECT_TRUE(s0.isNorth(p));

		// test for a point to the south of the first point
		Position p1(drand(-180.0, 180.0), drand(-90.0, s0.lat()));
		// cerr << to_string(p.getlon()) << "," << to_string(p.getlat()) << endl;
		// cerr << to_string(sb.angleDeg()) << "," << to_string(sb.mag()) << endl;
		EXPECT_TRUE(p1.isValid());
		EXPECT_FALSE(s0.isNorth(p1));
	}
}

TEST_F (PositionTest, IsEastTest) {
	logfile << "=================================================" << endl;
	logfile << "Testing isEast() function" << endl;
	logfile << "=================================================" << endl;
	for (int i = 0; i < TESTLEN; i++) {
		// construct a random point
		Position s0(drand(-180.0, 180.0), drand(-90.0, 90.0));
		double l0 = drand(s0.lon(), (s0.lon() + 350.0));
		if (l0 > 180.0) l0 -= 360.0;
		EXPECT_TRUE(s0.isValid());

		// test for a point to the east of the first point and west of the random meridian 
		// (expect true)
		Position s1(drand(s0.lon(), l0), drand(-90.0, 90.0)); 
		s1.norm();
		EXPECT_TRUE(s1.isValid());
		EXPECT_TRUE(s0.isEast(s1, l0));

		// test for a point to the east of the first point and east of the random meridian 
		// (expect false)
		Position s2(drand(l0, (l0 + 9.99)), drand(-90.0, 90.0)); 
		s2.norm();
		EXPECT_TRUE(s2.isValid());
		EXPECT_FALSE(s0.isEast(s2, l0));

		// test for a point to the west of the first point and east of the antimeridian and 
		// a random meridian (expect false)
		l0 = drand(s0.lon(), 180.0);
		Position s3(drand(-180.0, s0.lon()), drand(-90.0, 90.0));
		EXPECT_TRUE(s3.isValid());
		EXPECT_FALSE(s0.isEast(s3, l0));

		// test for a point to the west of the first point, west of the antimeridian, and 
		// the east of a random meridian (expect false)
	}
}
