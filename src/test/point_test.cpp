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

#define TOL (0.00001)	// Tolerance for floating point comparisons
#define JFK {(40 + (38/60)), -(73 + (47/60))}
#define LAX {(33 + (57/60)), -(118 + (24/60))}

class PointTest : public ::testing::Test {
public:
	PointTest () {
		srand(time(NULL));
		logfile.open("point_test.log", std::fstream::out | std::fstream::app);
	}
	ofstream logfile;
}

TEST(PointTest, DefaultConstructor) {
	// test default constructor
	Point a;
	EXPECT_EQ(a.gettype(), GeometryType::Point);
	EXPECT_FALSE(a.isValid());
	EXPECT_EQ(a.getlon(), NAN);
	EXPECT_EQ(a.getlat(), NAN);
	EXPECT_EQ(a.getele(), 0);
}

TEST(PointTest, NumericalConstructor) {
	// test lon/lat/ele constructor
	// approx 19% of test points will be invalid 
	for (int i = 0; i < TESTLEN; i++) {
		double lon = drand(-198.0,198.0);
		double lat = drand(-99.0, 99.0);
		Point p(lon, lat);
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(p.isValid());
		} else {
			EXPECT_TRUE(p.isValid());
		}
	}
}

TEST(PointTest, PositionConstructor) {
	// test Position constructor
	// approx 19% of test points will be invalid 
	for (int i = 0; i < TESTLEN; i++) {
		double lon = drand(-198.0,198.0);
		double lat = drand(-99.0, 99.0);
		Point p(makePosition(lon, lat));
		EXPECT_EQ(lon, p.getlon());
		EXPECT_EQ(lat, p.getlat());
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
}

TEST(PointTest, JSON_Constructor) {
	// test Position constructor
	// approx 19% of test points will be invalid 
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double lon = drand(-198.0,198.0);
		double lat = drand(-99.0, 99.0);
		double ele = drand(-500.0, 25000.0);
		buf << "[" << to_string(lon) << "," << to_string(lat) << "," << to_string(ele) << "]";
		logfile << "Test Position " << i << " : " << buf << endl;
		Document d;
		d.Parse(buf.c_str());
		Point p(d);
		EXPECT_EQ(lon, p.getlon());
		EXPECT_EQ(lat, p.getlat());
		EXPECT_EQ(ele, p.getele());
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double lon = drand(-198.0,198.0);
		double lat = drand(-99.0, 99.0);
		buf << "[" << to_string(lon) << "," << to_string(lat) << "]";
		logfile << "Test Position " << i << " : " << buf << endl;
		Document d;
		d.Parse(buf.c_str());
		Point p(d);
		EXPECT_EQ(lon, p.getlon());
		EXPECT_EQ(lat, p.getlat());
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
	for (int i = 0; i < TESTLEN; i++) {
		ostringstream buf;
		double lon = drand(-198.0,198.0);
		double lat = drand(-99.0, 99.0);
		buf << "{\"type\":\"Point\",\"coordinates\":";
		buf << "[" << to_string(lon) << "," << to_string(lat) << "]}";
		logfile << "Test Position " << i << " : " << buf << endl;
		Document d;
		d.Parse(buf.c_str());
		Point p(d);
		EXPECT_EQ(lon, p.getlon());
		EXPECT_EQ(lat, p.getlat());
		if ((getlon(p) < -180.0) || (getlon(p) > 180.0) || 
			(getlat(p) < -90.0) || (getlat(p) > 90.0)) {
			EXPECT_FALSE(positionValid(p));
		} else {
			EXPECT_TRUE(positionValid(p));
		}
	}
}

TEST(PointTest, JSONobject) {
	for (int i = 0; i < TESTLEN; i++) {
		Document d;
		double lon = drand(-180.0,180.0);
		double lat = drand(-90.0, 90.0);
		Point p(lon,lat);
		EXPECT_TRUE(p.isValid());
		Value v = p.pack(d);
		EXPECT_TRUE(v.hasMember("type"));
		EXPECT_TRUE(v.hasMember("coordinates"));
		Value *typeptr = Pointer("/type").Get(v);
		Value *coordptr= Pointer("/coordinates").Get(v);
		EXPECT_FALSE(typeptr == nullptr);
		EXPECT_FALSE(coordptr == nullptr);
		EXPECT_TRUE(typeptr->IsString());
		EXPECT_TRUE(coordptr->IsArray());
		EXPECT_EQ(coordptr->Size(), 3);
		Point p1(v);
		EXPECT_EQ(p.getlon(), p1.getlon());
		EXPECT_EQ(p.getlat(), p1.getlat());
		EXPECT_EQ(p.getele(), p1.getele());
		EXPECT_TRUE(p1.isValid());
	}
}

TEST(PointTest, JSONarray) {
	for (int i = 0; i < TESTLEN; i++) {
		Document d;
		double lon = drand(-180.0,180.0);
		double lat = drand(-90.0, 90.0);
		Point p(lon,lat);
		EXPECT_TRUE(p.isValid());
		Value v = p.packArray(d);
		EXPECT_TRUE(v.IsArray());
		EXPECT_EQ(v.Size(), 3);
		Point p1(v);
		EXPECT_EQ(p.getlon(), p1.getlon());
		EXPECT_EQ(p.getlat(), p1.getlat());
		EXPECT_EQ(p.getele(), p1.getele());
		EXPECT_TRUE(p1.isValid());
	}
}

// Taken from the worked examples in the Aviation Formulary
TEST (PointTest, GreatCircleBearing) {
	Point lax LAX;
	Point jfk JFK;
	EXPECT_EQ((int)round(lax.bearing(jfk)), 66);
}

// Taken from the worked examples in the Aviation Formulary
TEST (PointTest, RhumbLineBearing) {
	Point lax LAX;
	Point jfk JFK;
	EXPECT_EQ((int)round(lax.bearing(jfk, CourseTypeEnum::RhumbLine)), 79);
	EXPECT_EQ((int)round(jfk.bearing(lax, CourseTypeEnum::RhumbLine)), (79-180));
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, GreatCircleDistance) {
	Point lax LAX;
	Point jfk JFK;
	EXPECT_TRUE(toleranceEquals(meters2nm(lax.distance(jfk)), 2193.0, 1));
	EXPECT_TRUE(toleranceEquals(meters2nm(jfk.distance(lax)), 2193.0, 1));
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, RhumbLineDistance) {
	Point lax LAX;
	Point jfk JFK;
	EXPECT_TRUE(toleranceEquals(meters2nm(lax.distance(jfk, CourseTypeEnum::RhumbLine)), 2214.0, 1));
	EXPECT_TRUE(toleranceEquals(meters2nm(jfk.distance(lax, CourseTypeEnum::RhumbLine)), 2214.0, 1));
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, GreatCircleTarget) {
	Position lax LAX;
	Position jfk JFK;
	TwoVector vec;
	vec = lax.target(jfk);
	EXPECT_TRUE(toleranceEquals(meters2nm(vec.mag()), 2193.0, 1));
	EXPECT_EQ((int)round(vec.angleDeg()), 66);
	vec = jfk.target(lax);
	EXPECT_TRUE(toleranceEquals(meters2nm(vec.mag()), 2193.0, 1));
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, RhumbLineTarget) {
	Point lax LAX;
	Point jfk JFK;
	TwoVector vec;
	vec = lax.target(jfk, CourseTypeEnum::RhumbLine);
	EXPECT_TRUE(toleranceEquals(meters2nm(vec.mag()), 2214.0, 1));
	EXPECT_EQ((int)round(vec.angleDeg()), 79);
	vec = jfk.target(lax, CourseTypeEnum::RhumbLine);
	EXPECT_TRUE(toleranceEquals(meters2nm(vec.mag()), 2214.0, 1));
	EXPECT_EQ((int)round(vec.angleDeg()), (79-180));
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, GreatCircleProject) {
	Location lax LAX;
	Location jfk JFK;
	Location loc;
	TwoVector vec = TwoVector::getVectorDeg(66, (nm2meters(100)));
	loc = lax.project(vec);
	EXPECT_TRUE(toleranceEquals(loc.lat, 33.666, 0.01));
	EXPECT_TRUE(toleranceEquals(loc.lon, -116.176, 0.01));	
}

// Taken from the worked examples in the Aviation Formulary, modified to get the right answers with GeographicLib
TEST (PointTest, RhumbLineProject) {
	Location lax LAX;
	Location jfk JFK;
	Location loc;
	TwoVector vec = TwoVector::getVectorDeg(79, nm2meters(2214));
	loc = lax.project(vec);
	EXPECT_TRUE(toleranceEquals(loc.lat, 32.1341, 0.1));
	EXPECT_TRUE(toleranceEquals(loc.lon, -74.0216, 0.1));
}

TEST(PointTest, Bearing_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		double lon = drand(-180.0,180.0);
		double lat = drand(-75.0, 75.0);
		double b = drand(0.0,360.0);
		double d = drand(100.0,10000.0);
		Point p(lon,lat);
		EXPECT_TRUE(p.isValid());
		Vector v {1,0};
		v.angleDeg(b);
		v.mag(d);
		EXPECT_TRUE(v.isValid());
		Point p1(p.project(v, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(p1.isValid());
		EXPECT_TRUE(toleranceEquals(p.bearing(p1, CourseTypeEnum::Approximate), b), 10.0);
	}
}

TEST(PointTest, Distance_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		double lon = drand(-180.0,180.0);
		double lat = drand(-75.0, 75.0);
		double b = drand(0.0,360.0);
		double d = drand(100.0,10000.0);
		Point p(lon,lat);
		EXPECT_TRUE(p.isValid());
		Vector v {1,0};
		v.angleDeg(b);
		v.mag(d);
		EXPECT_TRUE(v.isValid());
		Point p1(p.project(v, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(p1.isValid());
		EXPECT_TRUE(toleranceEquals(p.distance(p1, CourseTypeEnum::Approximate), d), 10.0);
	}
}

TEST(PointTest, Project_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-75.0, 75.0);
		double lon1 = lon0 + drand(-0.1,0.1);
		double lat1 = lat0 + drand(-0.1,0.1);
		Point p0(lon0, lat0);
		Point p1(lon1, lat1);
		EXPECT_TRUE(p0.isValid());
		EXPECT_TRUE(p1.isValid());
		Vector v = p0.target(p1.getPosition(),CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(v.isValid());
		Point p2(p0.project(v, CourseTypeEnum::Approximate));
		EXPECT_TRUE(p2.isValid());
		EXPECT_TRUE(toleranceEquals(p0.getlon(), p1.getlon(), 0.0001));
		EXPECT_TRUE(toleranceEquals(p0.getlat(), p1.getlat(), 0.0001));
	}
}

TEST(PointTest, Point2Segment_PointNearestLine_RhumbLine) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct closest position
		Vector v {1,0};
		v.angleDeg(b);
		v.mag(drand(TOL,d));
		EXPECT_TRUE(v.isValid());
		Point subp(s0.project(v, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(positionValid(subp));

		// construct point 
		double dist = drand(100.0,10000.0);
		double dir = (drand(-1.0,1.0) < 0) ? 270 : 90;
		Vector v1 {1,0};
		v1.angleDeg(dir);
		v1.mag(dist);
		EXPECT_TRUE(v1.isValid());
		Point p(subp.project(v1, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(subp.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(subp.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestEnd1_RhumbLine) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand((b+90), (b+270));
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s0.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s0.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s0.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestEnd2_RhumbLine) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand((b-90), (b+90));
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s1.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s1.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s1.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestDegenerate_RhumbLine) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		Point s0(lon0, lat0);
		EXPECT_TRUE(s0.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand(0, 360);
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s0.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s0.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s0.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestLine_Approximate) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct closest position
		Vector v {1,0};
		v.angleDeg(b);
		v.mag(drand(TOL,d));
		EXPECT_TRUE(v.isValid());
		Point subp(s0.project(v, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(positionValid(subp));

		// construct point 
		double dist = drand(100.0,10000.0);
		double dir = (drand(-1.0,1.0) < 0) ? 270 : 90;
		Vector v1 {1,0};
		v1.angleDeg(dir);
		v1.mag(dist);
		EXPECT_TRUE(v1.isValid());
		Point p(subp.project(v1, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::Approximate));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(subp.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(subp.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestEnd1_Approximate) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand((b+90), (b+270));
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s0.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::Approximate));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s0.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s0.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, Point2Segment_PointNearestEnd2_Approximate) {
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand((b-90), (b+90));
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s1.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s1, CourseTypeEnum::Approximate));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s1.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s1.getlat(),p1.getlat(),TOL));
	}	
}

TEST(PointTest, Point2Segment_PointNearestDegenerate_Approximate) {	
	for (int i=0; i < TESTLEN; i++) {
		// construct segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		Point s0(lon0, lat0);
		EXPECT_TRUE(s0.isValid());

		// construct point
		double dist = drand(100.0,10000.0);
		double dir = drand(0, 360);
		Vector v {1,0};
		v.angleDeg(dir);
		v.mag(dist);
		EXPECT_TRUE(v.isValid());
		Point p(s0.project(v, CourseTypeEnum::RhumbLine));

		// test point distance
		Point p1(point2segment(p, s0, s0, CourseTypeEnum::Approximate));
		EXPECT_TRUE(toleranceEquals(dist,p.distance(p1),1.0));
		EXPECT_TRUE(toleranceEquals(s0.getlon(),p1.getlon(),TOL));
		EXPECT_TRUE(toleranceEquals(s0.getlat(),p1.getlat(),TOL));
	}
}

TEST(PointTest, SegmentsIntersect_RhumbLine) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (180.0 - TOL));
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());
		v2.angleDeg(b + 180.0);
		v2.mag(drand(TOL, 50000));
		Point s3(intersect.project(v2, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s3.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,s2,s3,i0,i1,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 1);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}

TEST(PointTest, SegmentsDisjoint_RhumbLine) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (180.0 - TOL));
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());
		v2.angleDeg(b + 180.0);
		v2.mag(drand(TOL, 50000));
		Point s3(intersect.project(v2, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s3.isValid());

		// check non-intersection with disjoint lines
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s2,s1,s3,i0,i1,CourseTypeEnum::RhumbLine);
		EXPECT_FALSE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 0);
	}
}

TEST(PointTest, SegmentEndsOnLine_RhumbLine) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (360.0 - TOL));
		if (b == 180.0) b += TOL;
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,intersect,s2,i0,i1,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 1);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}

TEST(PointTest, SegmentsOverlap_RhumbLine) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		d = drand(d, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,intersect,s2,i0,i1,CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(i0.isValid());
		EXPECT_TRUE(i1.isValid());
		EXPECT_EQ(result, 2);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}

TEST(PointTest, SegmentsIntersect_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (180.0 - TOL));
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());
		v2.angleDeg(b + 180.0);
		v2.mag(drand(TOL, 50000));
		Point s3(intersect.project(v2, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s3.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,s2,s3,i0,i1,CourseTypeEnum::Approximate);
		EXPECT_TRUE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 1);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}

TEST(PointTest, SegmentsDisjoint_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (180.0 - TOL));
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());
		v2.angleDeg(b + 180.0);
		v2.mag(drand(TOL, 50000));
		Point s3(intersect.project(v2, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s3.isValid());

		// check non-intersection with disjoint lines
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s2,s1,s3,i0,i1,CourseTypeEnum::Approximate);
		EXPECT_FALSE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 0);
	}
}

TEST(PointTest, SegmentEndsOnLine_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		b = drand(TOL, (360.0 - TOL));
		if (b == 180.0) b += TOL;
		d = drand(TOL, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,intersect,s2,i0,i1,CourseTypeEnum::Approximate);
		EXPECT_TRUE(i0.isValid());
		EXPECT_FALSE(i1.isValid());
		EXPECT_EQ(result, 1);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}

TEST(PointTest, SegmentsOverlap_Approximate) {
	for (int i = 0; i < TESTLEN; i++) {
		// construct first segment
		double lon0 = drand(-180.0,180.0);
		double lat0 = drand(-85.0, 85.0);
		double lon1 = lon0 + drand(-1.0,1.0);
		double lat1 = lat0 + drand(-1.0,1.0);
		Point s0(lon0, lat0);
		Point s1(lon1, lat1);
		double b = s0.bearing(s1, CourseTypeEnum::RhumbLine);
		double d = s0.distance(s1, CourseTypeEnum::RhumbLine);
		EXPECT_TRUE(s0.isValid());
		EXPECT_TRUE(s1.isValid());

		// construct intersecting segment
		Vector v0 {1,0};
		Vector v1 {1,0};
		Vector v2 {1,0};
		v0.angleDeg(b);
		v0.mag(drand(TOL,d));
		EXPECT_TRUE(v0.isValid());	
		Point intersect(s0.project(v0, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(intersection.isValid());
		d = drand(d, 50000);
		v1.angleDeg(b);
		v1.mag(d);
		EXPECT_TRUE(v1.isValid());
		Point s2(intersect.project(v1, CourseTypeEnum::RhumbLine));
		EXPECT_TRUE(s2.isValid());

		// check intersection
		Point i0, i1;
		unsigned int result = segmentsIntersect(s0,s1,intersect,s2,i0,i1,CourseTypeEnum::Approximate);
		EXPECT_TRUE(i0.isValid());
		EXPECT_TRUE(i1.isValid());
		EXPECT_EQ(result, 2);
		EXPECT_TRUE(toleranceEquals(i0.getlon(), intersect.getlon(), TOL));
		EXPECT_TRUE(toleranceEquals(i0.getlat(), intersect.getlat(), TOL));
	}
}
