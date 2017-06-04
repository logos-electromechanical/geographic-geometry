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
		Position p = makePosition(drand(-198.0, 198.0), drand(-99.0, 99.0), drand(-1000000,100000));
		logfile << "Test Position " << i << " : " << getlon(p) << "," << getlat(p) << "," << getele(p); 
		logfile << "," << to_string(positionValid(p)) << endl;
		if (!positionValid(p)) invalid_detected++;
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
			invalid++;
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
	for (i = 0; i < TESTLEN; i++) {
		Position p = makePosition(drand(-198.0, 198.0), drand(-99.0, 99.0));
		logfile << "Test Position " << i << " : " << getlon(p) << "," << getlat(p) << "," << getele(p); 
		logfile << "," << to_string(positionValid(p)) << endl;
		if (!positionValid(p)) invalid_detected++;
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
			invalid++;
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
	logfile << "Total positions: " << to_string(i*2);
	logfile << " Invalid detected: " << to_string(invalid_detected);
	logfile << " Total invalid: " << to_string(invalid);
	EXPECT_FALSE(positionValid(makePosition(NAN, 0, 0)));
	EXPECT_FALSE(positionValid(makePosition(0, NAN, 0)));
	EXPECT_FALSE(positionValid(makePosition(0, 0, NAN)));
	EXPECT_FALSE(positionValid(makePosition(NAN, 0)));
	EXPECT_FALSE(positionValid(makePosition(0, NAN)));
}

TEST_F (PositionTest, JSON_Valid) {
	logfile << "=================================================" << endl;
	logfile << "Testing JSON pack/unpack" << endl;
	logfile << "=================================================" << endl;

	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		Position p = makePosition(drand(-180.0, 180.0), drand(-90.0, 90.0), drand(-1000000,100000));
		logfile << "Test Position " << i << " : " << getlon(p) << "," << getlat(p) << "," << getele(p); 
		logfile << "," << to_string(positionValid(p)) << endl;
		Document d0, d1;
		Value v = packPosition(p, d0);
		v.Accept(writer);
		const char* out = buf.GetString();
		d1.Parse(out);
		Position p1 = loadPosition(d1);
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(positionValid(p1));
		EXPECT_TRUE(toleranceEquals(getlon(p), getlon(p1), TOL));
		EXPECT_TRUE(toleranceEquals(getlat(p), getlat(p1), TOL));
		EXPECT_TRUE(toleranceEquals(getele(p), getele(p1), TOL));
	}
	for (int i = 0; i < TESTLEN; i++) {
		StringBuffer buf;
		Writer<StringBuffer> writer(buf);
		Position p = makePosition(drand(-180.0, 180.0), drand(-90.0, 90.0));
		logfile << "Test Position " << i << " : " << getlon(p) << "," << getlat(p) << "," << getele(p); 
		logfile << "," << to_string(positionValid(p)) << endl;
		Document d0, d1;
		Value v = packPosition(p, d0);
		v.Accept(writer);
		const char* out = buf.GetString();
		d1.Parse(out);
		Position p1 = loadPosition(d1);
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(positionValid(p1));
		EXPECT_TRUE(toleranceEquals(getlon(p), getlon(p1), TOL));
		EXPECT_TRUE(toleranceEquals(getlat(p), getlat(p1), TOL));
		EXPECT_TRUE(toleranceEquals(getele(p), getele(p1), TOL));
	}
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double lon = drand(-180.0,180.0);
		double lat = drand(-90.0,90.0);
		buf << "[" << to_string(lon) << "," << to_string(lat) << "]";
		logfile << "Test Position " << i << " : " << buf.str() << endl;
		Document d;
		d.Parse(buf.str().c_str());
		Position p = loadPosition(d);
		EXPECT_TRUE(positionValid(p));
		EXPECT_TRUE(toleranceEquals(getlon(p), lon, TOL));
		EXPECT_TRUE(toleranceEquals(getlat(p), lat, TOL));
		EXPECT_TRUE(toleranceEquals(getele(p), 0.0, TOL));
	}
}

TEST_F (PositionTest, IsLeftTest) {
	logfile << "=================================================" << endl;
	logfile << "Testing isLeft() function" << endl;
	logfile << "=================================================" << endl;
	for (int i = 0; i < TESTLEN; i++) {

		// construct a random line segment 1000m long
		Point s0(drand(-179.0,179.0), drand(-80.0, 80.0));
		EXPECT_TRUE(s0.isValid());
		double bearing = drand(0.0,360.0);
		double linelen = 10000;
		TwoVector s(1,0);
		s.angleDeg(bearing);
		s.mag(linelen);
		Point s1 = s0.project(s,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s1.isValid());

		// cerr << to_string(i) << endl;
		// cerr << to_string(bearing) << endl;
		// cerr << to_string(s0.getlon()) << "," << to_string(s0.getlat()) << endl;
		// cerr << to_string(s1.getlon()) << "," << to_string(s1.getlat()) << endl;

		// test for a point on the line
		// Point p = s0.project((s/2),CourseTypeEnum::RhumbLine);
		// EXPECT_TRUE(p.isValid());
		// EXPECT_TRUE(toleranceEquals(isLeft(p.getPosition(), s0.getPosition(), s1.getPosition()), 0.0, TOL));

		// test for a point to the left of the line
		TwoVector sa = s/2;
		sa.rotateDeg(drand(-90.0,-10.0));
		Point p = s0.project(sa,CourseTypeEnum::RhumbLine);
		// cerr << to_string(p.getlon()) << "," << to_string(p.getlat()) << endl;
		// cerr << to_string(sa.angleDeg()) << "," << to_string(sa.mag()) << endl;
		EXPECT_TRUE(p.isValid());
		EXPECT_GT(isLeft(p.getPosition(), s0.getPosition(), s1.getPosition()), 0.0);

		// test for a point to the right of the line
		TwoVector sb = s/2;
		sb.rotateDeg(drand(10.0,90.0));
		p = s0.project(sb,CourseTypeEnum::RhumbLine);
		// cerr << to_string(p.getlon()) << "," << to_string(p.getlat()) << endl;
		// cerr << to_string(sb.angleDeg()) << "," << to_string(sb.mag()) << endl;
		EXPECT_TRUE(p.isValid());
		EXPECT_LT(isLeft(p.getPosition(), s0.getPosition(), s1.getPosition()), 0.0);
	}
}

TEST_F (PositionTest, InSegmentTest) {
	logfile << "=================================================" << endl;
	logfile << "Testing inSegment()" << endl;
	logfile << "=================================================" << endl;

	// test the general case
	for (int i = 0; i < TESTLEN; i++) {

		// construct a random line segment 1000m long
		Point s0(drand(-179.0,179.0), drand(-80.0, 80.0));
		EXPECT_TRUE(s0.isValid());
		double bearing = drand(0.0,360.0);
		double linelen = 1000;
		TwoVector s(1,0);
		s.angleDeg(bearing);
		s.mag(linelen);
		Point s1 = s0.project(s,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s1.isValid());

		// construct a random point on the segment
		Point p0 = s0.project((s/(drand(1,linelen))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p0.isValid());
		EXPECT_TRUE(inSegment(p0.getPosition(), s0.getPosition(), s1.getPosition()));

		// construct a random point past the end of the segment
		Point p1 = s0.project((s*(drand(1,3))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p1.isValid());
		if (!p1.isValid()) {
			cerr << to_string(p1.getlon()) << "," << to_string(p1.getlat()); 
			cerr << endl; 
		}
		EXPECT_FALSE(inSegment(p1.getPosition(), s0.getPosition(), s1.getPosition()));

		// construct a random point before the beginning of the segment
		s.rotateDeg(180);
		Point p2 = s0.project((s*(drand(1,3))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p2.isValid());
		if (!p2.isValid()) {
			cerr << to_string(p2.getlon()) << "," << to_string(p2.getlat()); 
			cerr << endl;
		}
		EXPECT_FALSE(inSegment(p2.getPosition(), s0.getPosition(), s1.getPosition()));

	}

	// test the E-W line case
	for (int i = 0; i < TESTLEN; i++) {

		// construct a random E-W line segment 1000m long
		Point s0(drand(-179.0,179.0), drand(-80.0, 80.0));
		EXPECT_TRUE(s0.isValid());
		double bearing = drand(0.0,360.0);
		if (bearing < 180) {
			bearing = 90;
		} else bearing = 270;
		double linelen = 1000;
		TwoVector s(1,0);
		s.angleDeg(bearing);
		s.mag(linelen);
		Point s1 = s0.project(s,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s1.isValid());

		// construct a random point on the segment
		Point p0 = s0.project((s/(drand(1,linelen))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p0.isValid());
		EXPECT_TRUE(inSegment(p0.getPosition(), s0.getPosition(), s1.getPosition()));

		// construct a random point past the end of the segment
		Point p1 = s0.project((s*(drand(1,linelen))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p1.isValid());
		EXPECT_FALSE(inSegment(p1.getPosition(), s0.getPosition(), s1.getPosition()));

		// construct a random point before the beginning of the segment
		s.rotateDeg(180);
		Point p2 = s0.project((s*(drand(1,linelen))), CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(p2.isValid());
		EXPECT_FALSE(inSegment(p2.getPosition(), s0.getPosition(), s1.getPosition()));
	}
}
