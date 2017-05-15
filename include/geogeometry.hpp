/*******************************************
 * geogeometry.hpp
 *
 * This library implements the geometry types of geojson. It provides
 * methods for performing such functions as nearest point, bearing,
 * inclusion, and the like. Works with loxodromes, geodesics, and 
 * approximate coordinates. 
 *
 * Coordinates are in decimal latitude and longitude on the WGS84
 * geoid. Distances and altitudes are in meters. Bearings are in 
 * degrees clockwise from north. 
 * 
 * Geometry algorithms sourced from http://geomalgorithms.com/
 * 
 * It requires the following libraries:
 * rapidjson
 * libGeographic
 *
 */

#ifndef GEOGEOMETRY_H
#define GEOGEOMETRY_H

#include <string>
#include <tuple>
#include <cmath>
#include <utility>
#include <rapidjson/rapidjson.h>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>

using namespace GeographicLib;
using namespace std;
using namespace rapidjson;

// Useful constants
const double smallnum = 0.0000000001;

/** @enum CourseType
 * 
 * @brief is a strongly type enum representing the type of a course
 * 
 * @var CourseTypeEnum::GreatCircle 
 * defines a course as a geodesic or great circle
 * @var CourseTypeEnum::RhumbLine
 * defines a course as a loxodrome or rhumb line
 * @var CourseTypeEnum::Approximate
 * defines a course as approximate, i.e. as a course on a cartesian grid.
 * It should only be used for distances less than 100 km
 *
 */
enum class CourseTypeEnum : uint8_t {
	GreatCircle,
	RhumbLine,
	Approximate
};

/** @enum GeometryType
 *
 * @brief is a strongly type enum representing the type of a geometry object
 * 
 * @var GeometryType::Point
 * is a single location with longitude, latitude, and elevation
 * @var GeometryType::MultiPoint
 * is an unordered collection of Points
 * @var GeometryType::LineString
 * is an ordered collection of Points representing a polyline
 * @var GeometryType::MultiLineString
 * is a collection of LineStrings
 * @var GeometryType::Polygon
 * is a LineString whose first and last points are identical. It represents the 
 * perimeter of an area. 
 * @var GeometryType::MultiPolygon
 * is an unordered collection of polygons
 * @var GeometryType::GeometryCollection
 * is a collection of geometry objects
 * @var GeometryType::Undefined 
 * is the default value
 *
 */
enum class GeometryType : uint8_t {
	Point, 
	MultiPoint,
	LineString,
	MultiLineString,
	Polygon,
	MultiPolygon,
	GeometryCollection,
	Undefined
};

/// @brief represents a position. Elements are lon in degrees, lat in degrees, and elevation in meters
typedef tuple<double, double, double> Position;
/// @brief represents a bounding box. Elements are the southwest and northeast corners
typedef tuple<Position, Position> BoundingBox;

/// @brief retrieve the the longitude from a Position tuple
inline double& getlon (const Position p) {return get<0>(p);};
/// @brief retrieve the the latitude from a Position tuple
inline double& getlat (const Position p) {return get<1>(p);};
/// @brief retrieve the the elevation from a Position tuple
inline double& getele (const Position p) {return get<2>(p);};

/// @brief returns a Position tuple with the given longitude, latitude, and elevation, in that order.
/// The ele parameter is optional and defaults to 0.0
inline Position makePosition (const double& lon, const double& lat, const double& ele = 0) {
	return make_tuple(lon, lat, ele);
};

/// @brief creates a Postion tuple from a GeoJSON array
inline Position loadPosition (const Value& v) {
	double lon = NAN;
	double lat = NAN;
	double ele = 0;
	Value *vp = nullptr;
	if (v.IsArray()) {
		vp = Pointer("/0").Get(v);
		if (vp && vp->IsDouble()) lon = vp->GetDouble();
		vp = Pointer("/1").Get(v);
		if (vp && vp->IsDouble()) lat = vp->GetDouble();
		vp = Pointer("/2").Get(v);
		if (vp && vp->IsDouble()) ele = vp->GetDouble();
	} 
	return makePosition(lon, lat, ele);
};

/// @brief creates a GeoJSON position array from the given Position in the given Document
inline Value packPosition (const Position& p, Document& d) {
	Value v(kArrayType);
	auto a = d.GetAllocator();
	v.PushBack(getlon(p), a);
	v.PushBack(getlat(p), a);
	v.PushBack(getele(p), a);
	return v;
}

/// @brief checks whether a Position refers to a real point
inline bool positionValid (const Position& p) {
	return (isfinite(getlon(p)) && 
			isfinite(getlat(p)) &&
			isfinite(getele(p)) &&
			(getlon(p) <= 180.0) &&
			(getlon(p) >= -180.0) &&
			(getlat(p) <= 90.0) &&
			(getlat(p) >= -90.0));
}

/// @brief retrieves the southwest corner of the given BoundingBox
inline Position getsw (const BoundingBox& bb) {return get<0>(bb);};
/// @brief retrieves the northeast corner of the given BoundingBox
inline Position getne (const BoundingBox& bb) {return get<1>(bb);};
/// @brief creates a bounding box from the given positions
inline BoundingBox makeBBox (const Position& sw, const Position& ne) {
	return make_tuple(sw, ne);
};

/// @brief creates a BoundBox from a GeoJSON bounding box array
inline BoundingBox loadBBox (const Value& v) {
	double swlon = NAN;
	double swlat = NAN;
	double nelon = NAN;
	double nelat = NAN;
	Value *vp = nullptr;
	if (v.IsArray()) {
		vp = Pointer("/0").Get();
		if (vp && vp->IsDouble()) swlon = vp->GetDouble();
		vp = Pointer("/1").Get();
		if (vp && vp->IsDouble()) swlat = vp->GetDouble();
		vp = Pointer("/2").Get();
		if (vp && vp->IsDouble()) nelon = vp->GetDouble();
		vp = Pointer("/3").Get();
		if (vp && vp->IsDouble()) nelat = vp->GetDouble();
	}
	return makeBBox(swlon, swlat, nelon, nelat);
};

/// @brief creates a GeoJSON bounding box array from the given BoundingBox in the given Document
inline Value packBBox (const BoundingBox& b, Document& d) {
	Value v(kArrayType);
	auto a = d.GetAllocator();
	v.PushBack(getlon(getsw(b)), a);
	v.PushBack(getlat(getsw(b)), a);
	v.PushBack(getlon(getne(b)), a);
	v.PushBack(getlat(getne(b)), a);
	return v;
};

/// @brief checks whether a BoundingBox refers to a real bounding box
inline bool bboxValid (const BoundingBox& b) {
	return(positionValid(getsw(b)) &&
			positionValid(getne(b)) &&
			(getlat(getne(b)) > getlat(getsw(b))));
};

// utility functions
inline double deg2rad (double deg) { return deg * ( M_PI / 180.0 ); }	/**< Convert degrees to radians */
inline double rad2deg (double rad) { return rad * ( 180.0 / M_PI ); }	/**< Convert radians to degrees */

/// @brief returns the number of meters per degree of latitude at given latitude
inline double metersPerDegreeLat (const double& lon) {
	phi = deg2rad(lat);
	return (111132.954 - 559.822 * cos(2 * phi) + 1.175 * cos(4 * phi));
};

/// @brief returns the number of meters per degree of longitude at the given latitude
inline double metersPerDegreeLon (const double& lon) {
	phi = deg2rad(lat);
	a = 6378137.0;	// WGS84 semi-major axis
	e = 0.0818192;	// WGS84 flattening
	return ((M_PI * a * cos(phi)) / (180 * sqrt(1 - (e * e * pow(sin(phi),2)))));
};

/// @brief return the point on a line segment closest to a point
/// Note that this function may behave strangely with great circle distances
/// This function was derived from http://geomalgorithms.com/a02-_lines.html
/// @param p is the point external to the line segment
/// @param s0 is one end of the line segment
/// @param s1 is the other end of the line segment
/// @param type the distance type used in this calculation
/// @returns the Position on the line segment closest to p
inline Position point2segment(const Point& p, const Point& s0, const Point& s1, const GeometryType type) {
	TwoVector v = s0.target(s1.getPosition(), type);
	TwoVector w = s0.target(p.getPosition(), type);
	double c1 = w * v;
	if (c1 <= 0) return p0;
	double c2 = v * v;
	if (c2 <= c1) return p1;
	double b = c1/c2;
	return s0.projection((b * v), type);
};

/// @brief determine whether a collinear point lies in a segment
/// This function was derived from http://geomalgorithms.com/a05-_intersect-1.html
/// @param p the point to test
/// @param s0 the start of the segment to test
/// @param s1 the end of the segment to test
/// @returns true if the point lies in the segment
inline bool inSegment(const Position& p, const Position& s0, const Position& s1) {
	if (fabs(getlat(s0) - getlat(s1)) > smallnum) {		// the segment is not vertical
		if ((getlat(s0) <= getlat(p)) && (getlat(p) <= getlat(s1))) return true;
		if ((getlat(s0) >= getlat(p)) && (getlat(p) >= getlat(s1))) return true;
	} else {											// segment is vertical
		if ((getlon(s0) <= getlon(p)) && (getlon(p) <= getlon(s1))) return true;
		if ((getlon(s0) >= getlon(p)) && (getlon(p) >= getlon(s1))) return true;
	}
	return false;
}

/// @brief determines whether two line segments intersect
/// Note that this function may behave strangely with great circle distances
/// This function was derived from http://geomalgorithms.com/a05-_intersect-1.html
/// @param s10 start of the first line segment
/// @param s11 end of the first line segment
/// @param s20 start of the second line segment
/// @param s21 end of the second line segment
/// @param intersect0 this Position is set to the intersection point, 
/// or untouched if there is no intersection
/// @param intersect1 if the segments overlap, intersect1 is set to the 
/// other end of the overlapping segment, otherwise untouched
/// @returns 0 if no intersection, 1 if they intersect at a point, and 2 if they overlap
inline unsigned int segmentsIntersect (const Point& s10, const Point& s11, 
										const Point& s20, const Point& s21,
										Point& intersect0, Point& intersect1,
										CourseTypeEnum type = CourseTypeEnum::RhumbLine) {
	TwoVector u = s10.target(s11.getPosition(), type);
	TwoVector v = s20.target(s21.getPosition(), type);
	TwoVector w = s20.target(s11.getPosition(), type);
	double D = u.perp(v);

	if (fabs(D) < smallnum)	{								// S1 & S2 are parallel
		if ((u.perp(v) != 0) || v.perp(w) != 0) return 0;	// They are not collinear
		// check if they are degenerate points
		double du = u * u;
		double dv = v * v;
		if ((du < smallnum) && (dv < smallnum)) {			// Both segments are points
			if (s10 != s20) return 0;						// They are distinct points
			intersect0 = s10;
			return 1;
		} else if (du < smallnum) {							// s1 is a single point
			if (inSegment(s10.getPosition(), s20.getPosition(), s21.getPosition())) {
				intersect0 = s10;
				return 1;
			} else return 0;
		} else if (dv < smallnum) {							// s2 is a single point
			if (inSegment(s20.getPosition(), s10.getPosition(), s11.getPosition())) {
				intersect0 = s20;
				return 1;
			} else return 0;
		} else {											// s1 and s2 are collinear and may overlap if it exists
			double t0, t1;									// endpoints of s1 in eqn for s2
			TwoVector w2 = s20.target(s11.getPosition(), type);
			if (fabs(v.x()) > smallnum) {
				t0 = w.x()/v.x();
				t1 = w2.x()/v.x();
			} else {
				t0 = w.y()/v.y();
				t1 = w2.y()/v.y();
			}
			if (t0 > t1) swap(t0,t1);						// t0 must be smaller than t1
			if ((t0 > 1) || (t1 < 0)) return 0;				// no overlap, return 
			t0 = t0 < 0 ? 0 : t0;							// clip to min 0
			t1 = t1 > 1 ? 1 : t1;							// clip to max 1
			if (fabs(t0-t1) < smallnum) {					// intersection is a point
				intersect0 = Point(s20.projection(v * t0, type));
				return 1;
			} else {
				intersect0 = Point(s20.projection(v * t0, type));
				intersect1 = Point(s20.projection(v * t1, type));
				return 2;
			}
		}
	} else {								// the segments are skew and map intersect
		double sI = v.perp(w)/D;			// get intersect parameter for s1
		if ((sI < 0) || (sI > 1)) return 0;	// no intersection
		double tI = u.perp(w)/D;			// get the intersect parameters for s2
		if ((tI < 0) || (tI > 1)) return 0;	// no intersection
		intersect0 = Point(s10.projection(u * sI, type));
		return 1;
	}
}

/** @class TwoVector
 *
 * @brief is a type for representing a two dimensional vector
 *
 */
class TwoVector {
	public:
		/// @brief creates a blank vector (both elements = NAN)
		TwoVector () {};
		/// @brief creates a vector
		/// @param x is the northing of the vector, in meters
		/// @param y is the easting of the vector, in meters
		TwoVector (double x, double y) : _x(x), _y(y) {};
		/// @brief creates a vector from the given json object
		/// @param input a reference to a json object with "x" & "y" members
		TwoVector (Value& input) {load(&input);};
		/// @brief creates a vector from the given json object
		/// @param input a pointer to a json object with "x" & "y" members
		TwoVector (Value* input) {load(input);};

		/// @brief populate this vector from the given json object
		/// @param d a pointer to a json object with "x" & "y" members
		/// @returns true if vector is successfully populated
		bool TwoVector::load(Value* d) {
			if (!d) return false;
			Value *xval = Pointer("/x").Get(d);
			Value *yval = Pointer("/y").Get(d);
			if (xval && yval && xval->IsDouble() && yval->IsDouble()) {
				_x = xval->GetDouble();
				_y = yval->GetDouble();
			} else return false;
			return isValid();
		}

		/// @brief creates a json object in Document d from this vector
		/// @param d a reference to a json Document to create the Value in
		/// @returns a json object in Document d
		Value TwoVector::pack (Document& d) const {
			Value v(kObject);
			v.AddMember("x", _x, d.GetAllocator());
			v.AddMember("y", _y, d.GetAllocator());
			return v;
		}

		/// @brief returns true if the vector is valid, i.e. has finite members
		bool isValid () {return ((isfinite(_x)) && (isfinite(_y)));};
		
		/// @brief create a vector with a given length and angle
		/// @param ang result vector angle, in radians
		/// @param mag result vector length, in meters
		/// @return resulting vector
		static TwoVector getVectorRad(double ang, double mag) {
			TwoVector result {1, 0};	// create a unit vector at angle 0
			result.rotateRad(ang);		// rotate that vector through the given ang
			result *= mag;				// multiply by the desired magnitude
			return result;
		}		
		
		/// @brief create a vector with a given length and angle
		/// @param ang result vector angle, in degrees
		/// @param mag result vector length, in meters
		/// @return resulting vector
		static TwoVector getVectorDeg(double ang, double mag) {return getVectorRad(deg2rad(ang),mag);};		/**< Get a vector with the given angle (in degrees) and magnitude */
		
		// get/set cartesian values, angle, or magnitude 
		double inline x () {return _x;};											/**< Get the x value */
		void inline x (double x) {_x=x;};											/**< Set the x value */
		double inline y () {return _y;};											/**< Get the y value */
		void inline y (double y) {_y=y;};											/**< Set the y value */
		double inline mag () {return sqrt((_x*_x)+(_y*_y));};						/**< Get the magnitude */
		void inline mag (double _mag) {*this *= _mag/mag();};						/**< Set the magnitude */
		double inline angleRad () {return atan2(_y,_x);}; 							/**< Get the angle in radians */
		void angleRad (double _ang) {this->rotateRad(_ang - this->angleRad());};	/**< Set the angle in radians */
		double inline angleDeg () {return rad2deg(atan2(_y,_x));}; 					/**< Get the angle in degrees */
		void inline angleDeg (double _ang) {angleRad(deg2rad(_ang));};				/**< Set the angle in radians */
		
		// rotations
		/// @brief rotates the vector through the given number of degrees
		/// @param _deg angle to rotate the vector through in degrees
		void inline rotateDeg (double _deg) {rotateRad(deg2rad(_deg));};	/**< Rotate vector through the given angle in degrees */
		
		/// @brief rotates the vector through the given number of radians
		/// @param _deg angle to rotate the vector through in radians
		void rotateRad (double _rad) {										/**< Rotate vector through the given angle in radians */
			double x, y;
			x = (_x * cos(_rad)) - (_y * sin(_rad));
			y = (_x * sin(_rad)) + (_y * cos(_rad));
			_x = x;
			_y = y;
		}											
		
		/// @brief calculates the perp product of two vectors 
		/// (useful for line intersections)
		/// @param v is the vector to perp with this vector
		double perp (TwoVector v) {return ((x() * v.x()) - (y() * v.y()));};

		// vector math
		inline TwoVector& operator+= (const TwoVector& r) {							/**< Vector addition */
			this->_x += r._x; 
			this->_y += r._y; 
			return *this; 
		};
		
		inline TwoVector& operator-= (const TwoVector& r) {							/**< Vector subtraction */ 
			this->_x -= r._x; 
			this->_y -= r._y; 
			return *this; 
		};
		
		inline TwoVector& operator*= (const double& r) {							/**< Scalar multiplication */ 
			this->_x *= r; 
			this->_y *= r; 
			return *this; 
		};
		
		friend double inline operator* (const TwoVector& l, const TwoVector& r) {	/**< Dot product */ 
			return ((l._x*r._x)+(l._y*r._y)); 
		};
		
		friend TwoVector inline operator+ (TwoVector l, const TwoVector& r) {		/**< Vector addition */
			l += r;
			return l;
		};
		
		friend TwoVector inline operator- (TwoVector l, const TwoVector& r) {		/**< Vector subtraction */ 
			l -= r;
			return l;
		}
		
		friend TwoVector inline operator* (TwoVector l, const double r) {			/**< Scalar multiplication */
			l *= r;
			return l;
		}
		
		friend TwoVector inline operator/ (TwoVector l, const double r) {			/**< Scalar division */
			l *= (1.0/r);
			return l;
		}
		
		// other vector functions
		TwoVector unit () {				/**< Get the corresponding unit vector */
			TwoVector result;
			result = *this * (1/this->mag());
			return result;
		}					
		
	protected:
		double _x = NAN;
		double _y = NAN;	 		
};

/** @class GeometryRoot
 *
 * @brief is a base class for representing GeoJSON geometry objects
 *
 */
class GeometryRoot {
	public:
		/// @brief creates a blank object
		virtual GeometryRoot () = 0;
		// @brief creates a geometry object from a GeoJSON object
		// @param val reference to a GeoJSON object
		GeometryRoot (const Value& val) {load(&val);};
		// @brief creates a geometry object from a GeoJSON object
		// @param val pointer to a GeoJSON object
		GeometryRoot (const Value* val) {load(val);};
		// @brief populates a geometry object from a GeoJSON object
		// @param val reference to a GeoJSON object
		// @returns true if object is successfully populated
		bool load (const Value& val) {return load(&val);};
		// @brief populated a geometry object from a GeoJSON object
		// @param val pointer to a GeoJSON object
		// @returns true if object is successfully populated
		virtual bool load (const Value* val) = 0;
		// @brief creates a GeoJSON object from this object
		// @param d the GeoJSON document to create the object in 
		// @returns the packed GeoJSON object
		virtual Value pack (Document& d) const = 0;
		// @brief creates a GeoJSON coordinates array from this object
		// @param d the GeoJSON document to create the array in 
		// @returns the packed GeoJSON coordinates array
		virtual Value packArray(Document& d) = 0;

		/// @brief generates a pointer to a new geometry object from the given GeoJSON object
		/// @param val pointer to a GeoJSON object
		/// @returns a pointer to a new GeoJSON object to a nullptr if the passed object is not valid GeoJSON
		static GeometryRoot *factory (Value* val) {
			Value *valtype = Pointer("/type").Get(val);
			if (valtype && valtype->IsString()) {
				string objType = valtype->GetString();
				if (objType == "Point") {
					return new Point(val);
				} else if (objType == "MultiPoint") {
					return new MultiPoint(val);
				} else if (objType == "LineString") {
					return new LineString(val);
				} else if (objType == "MultiLineString") {
					return new MultiLineString(val);
				} else if (objType == "Polygon") {
					return new Polygon(val);
				} else if (objType == "MultiPolygon") {
					return new MultiPolygon(val);
				} else if (objType == "GeometryCollection") {
					return new GeometryCollection(val);
				} else return nullptr;
			}
		}

		/// @brief fetches a constant reference to this object's type
		const GeometryType& gettype () const {return _mytype;};

		/// @brief fetches the bounding box of the current object, and generates it if it has not yet been generated
		const BoundingBox& getbox () {
			if (!_hasbox) _hasbox = makebox();
			return _bbox;
		}

		/// @brief checks whether the current object is valid
		virtual bool isValid () = 0;

	protected:
		virtual bool makebox () = 0;
		bool 			_hasbox = false;
		BoundingBox 	_bbox = ;
		GeometryType 	_mytype = GeometryType::Undefined;
};

/** @class Point
 *
 * @brief represents a GeoJSON Point type
 * It also has a variety of navigation related methods for points
 *
 */
class Point : public GeometryRoot {
	public:	
		/// @brief creates a blank Point
		Point () : _mytype(GeometryType::Point), 
			_posn(makePosition(NAN, NAN, 0)) {};

		/// @brief creates a Point at the given position
		/// @param posn the desired position
		Point (const Position& posn) : 
			_posn(posn), 
			_mytype(GeometryType::Point) {};

		/// @brief creates a Point at the given position
		/// @param lon the desired longitude (decimal degrees)
		/// @param lat the desired latitude (decimal degrees)
		/// @param ele the desired elevation (meters above datum)
		Point (const double lon, const double lat, const double ele = 0.0)  : 
			_mytype(GeometryType::Point),
			_posn(makePosition(lon, lat, ele)) {};

		/// @brief populated this point from a GeoJSON object
		/// @param val pointer to a GeoJSON object
		/// @returns true if object is successfully populated
		bool load (const Value* val) {
			_mytype = GeometryType::Point;
			if (!val) return false;
			if (val->IsObject()) {
				Value *typeptr = Pointer("/type").Get(val);
				Value *coordptr= Pointer("/coordinates").Get(val);
				if (typeptr && coordptr) {
					if (typeptr->IsString() && 
						(typeptr->GetString() == "Point") &&
						coordptr->IsArray()) {
						_posn = loadPosition(*coordptr);
					} else return false;
				} else return false;
			} else if (val->IsArray() && (val->Size() > 1)) {
				_posn = loadPosition(*val);
			} else return false;
			return isValid();
		};

		/// @brief creates a GeoJSON Point from this Point
		/// @param d the GeoJSON document to create the object in 
		/// @returns the packed GeoJSON Point
		Value pack (Document& d) {
			Value v(kObject);
			auto a = d.GetAllocator();
			v.AddMember("type", "Point", a);
			v.AddMember("coordinates", packPosition(_posn, d), a);
			return v;
		};

		/// @brief creates a GeoJSON Point coordinate array from this Point
		/// @param d the GeoJSON document to create the array in 
		/// @returns the packed GeoJSON Point coordinates array
		Value packArray (Document& d) {return packPosition(_posn, d);};

		/// @brief returns a reference to the embedded position
		Position& getPosition () {return _posn};

		/// @brief get a bearing from this Point to a distant position
		/// @param dest the Position of the target
		/// @param type the type of course to plot to the target. Note that the 
		/// a GreatCircle type will calculate the initial bearing. 
		/// @returns the bearing to target in degrees from north. Returns NAN if
		/// no valid bearing can be calculated
		double bearing (const Position& dest, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			if ((!this->isValid()) || (!positionValid(dest))) return NAN;
			double azi1, azi2, dist;
			switch (type) {
				case CourseTypeEnum::GreatCircle: {
					Geodesic::WGS84().Inverse(getlat(_posn), getlon(_posn), getlat(dest), getlon(dest), azi1, azi2);
					return azi1;
					break;
				}
				case CourseTypeEnum::RhumbLine:	{
					Rhumb::WGS84().Inverse(getlat(_posn), getlon(_posn), getlat(dest), getlon(dest), dist, azi1);
					return azi1;
					break;
				}
				case CourseTypeEnum::Approximate: {
					double deltalat = metersPerDegreeLat(getlat(_posn))*(getlat(_posn) - getlat(dest));
					double deltalon = metersPerDegreeLon(getlat(_posn))*(getlon(_posn) - getlon(dest));
					return rad2deg(atan2(deltalon, deltalat));
					break;
				}
				default:
					return NAN;
			}
			return NAN;
		};

		/// @brief get the distance from this Point to a distant position
		/// @param dest the Position of the target
		/// @param type the type of course to plot to the target. Note that the 
		/// a GreatCircle type will calculate the distance along the geodesic. 
		/// @returns the distance to target in meters
		double distance (const Position& dest, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			if ((!this->isValid()) || (!positionValid(dest))) return NAN;
			double azi1, azi2, dist;
			switch (type) {
				case CourseTypeEnum::GreatCircle: {
					Geodesic::WGS84().Inverse(getlat(_posn), getlon(_posn), getlat(dest), getlon(dest), azi1, azi2);
					return dist;
					break;
				}
				case CourseTypeEnum::RhumbLine:	{
					Rhumb::WGS84().Inverse(getlat(_posn), getlon(_posn), getlat(dest), getlon(dest), dist, azi1);
					return dist;
					break;
				}
				case CourseTypeEnum::Approximate: {
					double deltalat = metersPerDegreeLat(getlat(_posn))*(getlat(_posn) - getlat(dest));
					double deltalon = metersPerDegreeLon(getlat(_posn))*(getlon(_posn) - getlon(dest));
					return sqrt((deltalon*deltalon) + (deltalat*deltalat));
					break;
				}
				default:
					return NAN;
			}
			return NAN;
		};

		/// @brief get a vector from this Point to a distant position
		/// @param dest the Position of the target
		/// @param type the type of course to plot to the target. Note that the 
		/// a GreatCircle type will return the initial bearing of the geodesic and
		/// distance along the geodesic. 
		/// @returns a vector with a distance to the target in meters and an angle 
		/// in degress from north. 
		TwoVector target (const Position& dest, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			if ((!this->isValid()) || (!positionValid(dest))) return TwoVector {NAN,NAN};
			TwoVector t {1,0};
			t.rotateDeg(bearing(dest, type));
			t *= distance(dest, type);
			return t;
		};

		/// @brief get the Position defined by a vector from this Point 
		/// @param projection the vector to project along
		/// @param type the type of course to plot to the target. Note that the 
		/// a GreatCircle type will treat the projection vector as the initial bearing
		/// and distance of the geodesic
		/// @returns the projected position
		Position project (const TwoVector& projection, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			if ((!this->isValid()) || (!projection.isValid())) return makePosition(NAN,NAN);
			Position p = makePosition(NAN,NAN);
			switch (type) {
				case CourseTypeEnum::GreatCircle: {
					Geodesic::WGS84().Direct(getlat(_posn), getlon(_posn), projection.angleDeg(), projection.mag(), getlat(p), getlon(p));
					break;
				}
				case CourseTypeEnum::RhumbLine: {
					Rhumb::WGS84().Direct(getlat(_posn), getlon(_posn), projection.angleDeg(), projection.mag(), getlat(p), getlon(p));
					break;
				}
				case CourseTypeEnum::Approximate: {
					getlon(p) += projection.y()/metersPerDegreeLon(getlat(_posn));
					getlat(p) += projection.x()/metersPerDegreeLat(getlat(_posn));
					break;
				}
				default:
					break;
			}
			return p;
		};

		bool isValid () {return (positionValid(_posn));};

	private:
		Position _posn;

		bool makebox () {
			get<0>(_bbox) = _posn;
			get<1>(_bbox) = _posn;
		};
};

/** @class MultiPoint
 *
 * @brief represents a GeoJSON MultiPoint type
 * 
 */
class MultiPoint : public GeometryRoot {
	public:
		/// @brief creates a blank MultiPoint
		MultiPoint () : _mytype(GeometryType::MultiPoint) {};

		/// @brief creates a MultiPoint object from the given positions
		/// @param vector of Positions to create MultiPoint from 
		MultiPoint (const vector<Position>& v) : 
			_mytype(GeometryType::MultiPoint),
			_points(v) {};

		/// @brief creates a MultiPoint object from the given points
		/// @param vector of Points to create MultiPoint from
		MultiPoint (const vector<Point>& v) : _mytype(GeometryType::MultiPoint) {
			for (auto& a : v) {
				_points.emplace_back(Point(a));
			}
		};

		/// @brief is a const_iterator type for returning a target position
		typedef vector<Position>::const_iterator PositionConstItr;

		/// @brief populated this MultiPoint from a GeoJSON object
		/// @param val pointer to a GeoJSON object
		/// @returns true if object is successfully populated
		bool load (const Value* val) {
			_mytype = GeometryType::MultiPoint;
			if (!val) return false;
			Value *typeptr = nullptr;
			Value *coordptr = nullptr;
			Value *bboxptr = nullptr;
			if (val->IsObject()) {
				typeptr = Pointer("/type").Get(val);
				coordptr= Pointer("/coordinates").Get(val);
				bboxptr= Pointer("/bbox").Get(val);
				if (!typeptr || !coordptr || (!typeptr->IsString() || 
					!(typeptr->GetString() == "MultiPoint") ||
					!coordptr->IsArray())) return false;
				if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
			} else if (val->IsArray()) {
				coordptr = val;
			} else return false;
			_points.clear();
			for (auto& v : coordptr->GetArray()) {
				if (v.IsArray()) _points.emplace_back(loadPosition(v));
			}
			return isValid();
		};

		/// @brief creates a GeoJSON MultiPoint from this MultiPoint
		/// @param d the GeoJSON document to create the object in 
		/// @returns the packed GeoJSON MultiPoint
		Value pack (Document& d) {
			Value b(kArrayType);
			Value v(kObject);
			auto a = d.GetAllocator();
			v.AddMember("type", "MultiPoint", a);
			v.AddMember("coordinates", packArray(d), a);
			v.AddMember("bbox", packBBox(getbox(), d), a);
			return v;
		};

		/// @brief creates a GeoJSON MultiPoint coordinate array from this Point
		/// @param d the GeoJSON document to create the array in 
		/// @returns the packed GeoJSON MultiPoint coordinates array
		Value packArray (Document& d) {
			Value b(kArrayType);
			auto a = d.GetAllocator();
			for (auto& p : _points) {
				b.PushBack(packPosition(p, d), a);
			}
			return b;
		};

		/// @brief return a reference to the members of the MultiPoint
		vector<Position>& getPositions () const {return _points;};

		/// @brief returns a reference to a given member
		/// @params t the index to return
		/// @returns the position at index t
		Position& getPosition (const unsigned int& t) const {return _points[t];};

		/// @brief returns a Point created from a given member
		/// @params t the index to return
		/// @returns a Point created from the position at index t
		Points getPoint (const unsigned int& t) const {return Point(_points[t]);};

		/// @brief returns a vector of Points corresponding to the members of the MultiPoint
		vector<Point> getPoints () const {
			vector<Points> v;
			for (auto p&: _points) v.emplace_back(Point(p));
			return v;
		};

		/// @brief insert a new Point after the given index
		/// @param p the Point to insert
		/// @param posn the index to insert after
		/// @returns true if successful
		bool insert (const Point& p; const unsigned int& posn) {return insert(p.getPosition(), posn);};

		/// @brief insert a new Position after the given index
		/// @param p the Position to insert
		/// @param posn the index to insert after
		/// @returns true if successful
		bool insert (const Position& p; const unsigned int& posn) {
			if (posn >= _points.Size()) {
				_points.emplace_back(p);
			} else {
				_points.emplace((_points.begin()+posn), p);
			}
			return true;
		};

		/// @brief insert a new Point after the given const_interator
		/// @param i const_iterator of the insertion point
		/// @param p the Point to insert
		/// @returns true if successful
		bool insert (PositionConstItr& i, const Point& p) {
			return insert(i, p.getPosition());
		};

		/// @brief insert a new Position after the given const_interator
		/// @param i const_iterator of the insertion point
		/// @param p the Position to insert to insert
		/// @returns true if successful
		bool insert (PositionConstItr& i, const Position& p) {
			_points.emplace(i, p);
			return true;
		};

		/// @brief remove the point at the given index
		/// @param posn index to remove
		/// @returns true if successful
		bool deletePoint (const unsigned int& posn) {
			if (posn < _points.Size()) {
				_points.erase(_points.begin() + posn);
			} else return false;
			return true;
		};

		/// @brief remove the point at the given iterator
		/// @param p iterator to remove
		/// @returns true if successful
		bool deletePoint (PositionConstItr p) {
			_points.erase(p);
			return true;
		}

		/// @brief returns the size of the internal vector
		int size () {return _points.size();};
		/// @brief returns a const_iterator to the beginning of the internal vector
		PositionConstItr cbegin () {return _points.cbegin();};
		/// @brief returns a reverse const_iterator to the beginning of the internal vector
		PositionConstItr crbegin () {return _points.crbegin();};
		/// @brief returns a const_iterator to the end of the internal vector
		PositionConstItr cend () {return _points.cend();};
		/// @brief returns a reverse const_iterator to the end of the internal vector
		PositionConstItr crend () {return _points.crend();};

		/// @brief find the closest point to the given Point
		/// @param target point to find the closest point in the MultiPoint to 
		/// @param type type of the course to use for distance finding
		PositionConstItr closest (const Point& target, CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double dist = target.distance(_points.cbegin(), type);
			PositionConstItr p = _points.cbegin();
			for (auto& a : _points) {
				double newdist = target.distance(a.getPosition(), type);
				if (newdist < dist) {
					dist = newdist;
					p = a;
				} 
			}
			return p;
		};

		/// @brief check whether this MultiPoint is valid
		bool isValid () {
			bool valid = true;
			for (auto& p : _points) {
				valid &= positionValid(p);
			}
			return valid;
		}

	private:
		vector<Position> _points;
		bool makebox () {
			if (!isValid()) return false;
			double maxlat = getlat(_points.crbegin());
			double minlat = maxlat;
			double maxlon = getlon(_points.crbegin());
			double minlon = maxlon
			for (auto& a: _points) {

			}
			_bbox = makeBBox(makePosition(minlon, minlat), 
							makePosition(maxlon, maxlat));
			return true;
		};
};

/** @class LineString
 *
 * @brief represents a GeoJSON LineString type
 *
 */
class LineString : public GeometryRoot {
	public:
		/// @brief creates a blank LineString
		LineString () : _mytype(GeometryType::LineString) {};

		/// @brief creates a LineString object from the given positions
		/// @param vector of Positions to create LineString from 
		LineString (const vector<Position>& v) : 
			_mytype(GeometryType::LineString),
			_points(v) {};

		/// @brief creates a LineString object from the given points
		/// @param vector of Points to create LineString from
		LineString (const vector<Point>& v) : _mytype(GeometryType::LineString) {
			for (auto& a : v) {
				_points.emplace_back(Point(a));
			}
		};

		/// @brief creates a LineString object from the given MultiPoint
		/// @param v MultiPoint to create LineString from 
		LineString (const MultiPoint& v) : 
			_mytype(GeometryType::LineString),
			_points(v.getPositions()) {};

		/// @brief is a const_iterator type for returning a target vertex
		typedef vector<Position>::const_iterator VertexConstItr;

		/// @brief populated this LineString from a GeoJSON object
		/// @param val pointer to a GeoJSON object
		/// @returns true if object is successfully populated
		bool load (const Value* val) {
			_mytype = GeometryType::LineString;
			if (!val) return false;
			Value *typeptr = nullptr;
			Value *coordptr = nullptr;
			Value *bboxptr = nullptr;
			if (val->IsObject()) {
				typeptr = Pointer("/type").Get(val);
				coordptr= Pointer("/coordinates").Get(val);
				bboxptr= Pointer("/bbox").Get(val);
				if (!typeptr || !coordptr || (!typeptr->IsString() || 
					!(typeptr->GetString() == "LineString") ||
					!coordptr->IsArray())) return false;
				if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
			} else if (val->IsArray()) {
				coordptr = val;
			} else return false;
			_points.clear();
			for (auto& v : coordptr->GetArray()) {
				if (v.IsArray()) _points.emplace_back(loadPosition(v));
			}
			return isValid();
		};

		/// @brief creates a GeoJSON LineString from this LineString
		/// @param d the GeoJSON document to create the object in 
		/// @returns the packed GeoJSON LineString
		Value pack (Document& d) {
			Value b(kArrayType);
			Value v(kObject);
			auto a = d.GetAllocator();
			v.AddMember("type", "LineString", a);
			v.AddMember("coordinates", packArray(d), a);
			v.AddMember("bbox", packBBox(getbox(), d), a);
			return v;
		};

		/// @brief creates a GeoJSON LineString coordinate array from this LineString
		/// @param d the GeoJSON document to create the array in 
		/// @returns the packed GeoJSON LineString coordinates array
		Value packArray (Document& d) {
			Value b(kArrayType);
			auto a = d.GetAllocator();
			for (auto& p : _points) {
				b.PushBack(packPosition(p, d), a);
			}
			return b;
		};

		/// @brief return a reference to the members of the LineString
		vector<Position>& getVertexPositions () const {return _points;};

		/// @brief returns a Point created from a given vertex
		/// @params t the index to return
		/// @returns a Point created from the position at index t
		Point& getVertexPoint (const unsigned int& t) const {return Point(_points[t]);};

		/// @brief returns a reference to a given vertex
		/// @params t the index to return
		/// @returns the position at index t
		Position& getVertextPosition (const unsigned int& t) const {return _points[t];};

		/// @brief returns a vector of Points corresponding to the members of the LineString
		vector<Point> getVertexPoints () const {
			vector<Point> v;
			for (auto p&: _points) v.emplace_back(Point(p));
			return v;
		};

		/// @brief insert a new Point after the given index
		/// @param p the Point to insert
		/// @param posn the index to insert after
		/// @returns true if successful
		bool insert (const Point& p; const unsigned int& posn) {return insert(p.getPosition(), posn);};

		/// @brief insert a new Position after the given index
		/// @param p the Position to insert
		/// @param posn the index to insert after
		/// @returns true if successful
		bool insert (const Position& p; const unsigned int& posn) {
			if (posn >= _points.Size()) {
				_points.emplace_back(p);
			} else {
				_points.emplace((_points.begin()+posn), p);
			}
			return true;
		};

		/// @brief insert a new Point after the given const_interator
		/// @param i const_iterator of the insertion point
		/// @param p the Point to insert
		/// @returns true if successful
		bool insert (VertexConstItr& i, const Point& p) {
			return insert(i, p.getPosition());
		};

		/// @brief insert a new Position after the given const_interator
		/// @param i const_iterator of the insertion point
		/// @param p the Position to insert to insert
		/// @returns true if successful
		bool insert (VertexConstItr& i, const Position& p) {
			_points.emplace(i, p);
			return true;
		};

		/// @brief remove the point at the given index
		/// @param posn index to remove
		/// @returns true if successful
		bool deletePoint (const unsigned int& posn) {
			if (posn < _points.Size()) {
				_points.erase(_points.begin() + posn);
			} else return false;
			return true;
		};

		/// @brief remove the point at the given iterator
		/// @param p iterator to remove
		/// @returns true if successful
		bool deletePoint (VertexConstItr p) {
			_points.erase(p);
			return true;
		}

		/// @brief returns the size of the internal vector
		int size () {return _points.size();};
		/// @brief returns a const_iterator to the beginning of the internal vector
		VertexConstItr cbegin () {return _points.cbegin();};
		/// @brief returns a reverse const_iterator to the beginning of the internal vector
		VertexConstItr crbegin () {return _points.crbegin();};
		/// @brief returns a const_iterator to the end of the internal vector
		VertexConstItr cend () {return _points.cend();};
		/// @brief returns a reverse const_iterator to the end of the internal vector
		VertexConstItr crend () {return _points.crend();};

		/// @brief find the closest vertex to the given Point
		/// @param target point to find the closest vertex in the LineString to 
		/// @param type type of the course to use for distance finding
		/// @returns a const_iterator to the closest vertex to the target
		VertexConstItr closestVertex (const Point& target, CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double dist = target.distance(_points.cbegin(), type);
			VertexConstItr p = _points.cbegin();
			for (auto& a : _points) {
				double newdist = target.distance(a.getPosition(), type);
				if (newdist < dist) {
					dist = newdist;
					p = a;
				} 
			}
			return p;
		};

		/// @brief find the closest position on the LineString to the given Point
		/// @param target the given Point
		/// @param type the type of course to use for distanc measuring etc.
		/// @returns a tuple containing the desired position and an iterator to 
		/// to the start of the line segment that contains the closest position 
		tuple<Position, VertexConstItr> closestPoint (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double bestdist = target.distance(_points.cbegin(), type);
			VertexConstItr p = _points.cbegin();
			Position best, newposn;
			for (auto& a : _points) {
				if (a != _points.end()) {
					newposn = point2segment(target, Point(a), Point(a + 1), type);
					double newdist = target.distance(newposn, type);
					if (newdist < bestdist) {
						bestdist = newdist;
						best = newposn;
						p = a;
					}
				}
			}
			return make_tuple(best, p);
		};

		/// @brief find the line segments in the LineString that intersect with
		/// the given line segment.
		/// @param 
		tuple<MultiPoint, vector<VertexConstItr>> intersection (const Point& s0, const Point& s1, const CourseTypeEnum type = CourseTypeEnum::RhumbLine) {
			Point o0, o1;
			vector<VertexConstItr> outitr;
			MultiPoint out;
			for (auto& a : _points) {
				if (a != _points.end()) {
					if (segmentsIntersect(Point(a), Point(a+1), s0, s1, o0, o1, type)) {
						outitr.emplace_back(a);
						out.insert(out.cend(), a);
					}
				} 
			}
			return make_tuple(out, outitr);
		}
	
		/// @brief calculate the total length of the LineString
		/// @param type type of the course to use for distance finding
		/// @returns total line length in meters
		double length (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double result = 0;
			for (auto& a : _points) {
				if (a != _points.end()) {
					result += a.distance((a + 1), type);
				}
			}
			return result;
		};

		/// @brief check whether this LineString is valid
		bool isValid () {
			bool valid = true;
			for (auto& p : _points) {
				valid &= positionValid(p);
			}
			return valid;
		}

	private:
		vector<Position> _points; 
		bool makebox () {
			if (!isValid()) return false;
			double maxlat = getlat(_points.crbegin());
			double minlat = maxlat;
			double maxlon = getlon(_points.crbegin());
			double minlon = maxlon
			for (auto& a: _points) {
				if (maxlat < getlat(a)) maxlat = getlat(a);
				if (minlat > getlat(a)) minlat = getlat(a);
				if (maxlon < getlon(a)) maxlon = getlon(a);
				if (minlon > getlon(a)) minlon = getlon(a);
			}
			_bbox = makeBBox(makePosition(minlon, minlat), 
							makePosition(maxlon, maxlat));
			return true;
		};
};

/** @class MultiLineString
 *
 * @brief represents a GeoJSON MultiLineString type
 *
 */
class MultiLineString : public GeometryRoot {
	public:
		/// @brief creates a blank MultiLineString
		MultiLineString () : _mytype(GeometryType::MultiLineString) {};

		/// @brief creates a MultiLineString populated with the given LineStrings
		/// @param vector of LineString objects
		MultiLineString (const vector<LineString>& v) : 
			_mytype(GeometryType::MultiLineString),
			_lines(v) {};

		/// @brief creates a MultiLineString populated with the given MultiPoint objects
		/// Each MultiPoint is converted to a LineString
		/// @param vector of MultiPoint objects
		MultiLineString (const vector<MultiPoint>& v) : 
			_mytype(GeometryType::MultiLineString) {
			for (auto& a : v) _points.emplace_back(LineString(a));
		};

		/// @brief is a const_iterator type for returning a target LineString
		typedef vector<LineString>::const_iterator LinesConstItr;

		/// @brief populated this MultiLineString from a GeoJSON object
		/// @param val pointer to a GeoJSON object
		/// @returns true if object is successfully populated
		bool load (const Value* val) {
			_mytype = GeometryType::MultiLineString;
			if (!val) return false;
			Value *typeptr = nullptr;
			Value *coordptr = nullptr;
			if (val->IsObject()) {
				typeptr = Pointer("/type").Get(val);
				coordptr= Pointer("/coordinates").Get(val);
				bboxptr= Pointer("/bbox").Get(val);
				if (!typeptr || !coordptr || (!typeptr->IsString() || 
					!(typeptr->GetString() == "MultiLineString") ||
					!coordptr->IsArray())) return false;
				if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
			} else if (val->IsArray()) {
				coordptr = val;
			} else return false;
			_points.clear();
			for (auto& v : coordptr->GetArray()) {
				LineString l(v);
				if (v.IsArray()) _points.emplace_back(l);
			}
			return isValid();
		};

		/// @brief creates a GeoJSON MultiLineString from this MultiLineString
		/// @param d the GeoJSON document to create the object in 
		/// @returns the packed GeoJSON MultiLineString
		Value pack () const {
			Value b(kArrayType);
			Value v(kObject);
			auto a = d.GetAllocator();
			v.AddMember("type", "MultiLineString", a);
			v.AddMember("coordinates", packArray(d), a);
			v.AddMember("bbox", packBBox(getbox(), d), a);
			return v;
		};

		/// @brief creates a GeoJSON MultiLineString coordinate array from this Point
		/// @param d the GeoJSON document to create the array in 
		/// @returns the packed GeoJSON MultiLineString coordinates array
		Value packArray (Document& d) {
			Value b(kArrayType);
			auto a = d.GetAllocator();
			for (auto& p : _points) {
				b.PushBack(p.packArray(), a);
			}
			return b;
		};

		/// @brief returns a reference to the members of this MultiLineString
		vector<LineString>& getLines () const {return _lines;};

		/// @brief returns a reference to a given member
		/// @params t the index to return
		/// @returns the line at index t
		LineString& getLine (const unsigned int& l) const {return _lines[l];};

		/// @brief insert a new line after the given index. 
		/// Line can be either a LineString, a MultiPoint, or a vector of Positions
		/// @param t is the index of the position to add the new line
		/// @param l is the line to add
		/// @returns true if successful
		bool insert (const LineString& l, const unsigned int& t) {return insert()};
		bool insert (const MultiPoint& l, const unsigned int& t) {};
		bool insert (const vector<Position>& l, const unsigned int& t) {
			if (t >= _lines.Size()) {
				_lines.emplace_back(l);
			} else {
				_lines.emplace((_lines.begin()+t), l);
			}
			return true;
		};

		/// @brief insert a new line after the given iterator. 
		/// Line can be either a LineString, a MultiPoint, or a vector of Positions
		/// @param t is the iterator to the position to add the new line
		/// @param l is the line to add
		/// @returns true if successful
		bool insert (const LinesConstItr& t, const LineString& l) {return insert(i, l.getPositions());};
		bool insert (const LinesConstItr& t, const MultiPoint& l) {return insert(i, l.getPositions());};
		bool insert (const LinesConstItr& t, const vector<Position>& l) {
			_lines.emplace(t, l);
			return true;
		};

		/// @brief remove the line at the given index or iterator
		/// @param t the item to remove. This can be an unsigned integer or an iterator
		/// @returns true if successful
		bool deleteLine (const unsigned int& t) {
			if (t < _lines.Size()) {
				_lines.erase(_lines.begin() + t);
			} else return false;
			return true;
		};
		bool deleteLine (LinesConstItr& t) {
			_lines.erase(p);
			return true;
		};

		/// @brief returns the size of the internal vector
		int size () {return _lines.size();};
		/// @brief returns a const_iterator to the beginning of the internal vector
		LinesConstItr cbegin () {return _lines.cbegin();};
		/// @brief returns a reverse const_iterator to the beginning of the internal vector
		LinesConstItr crbegin () {return _lines.crbegin();};
		/// @brief returns a const_iterator to the end of the internal vector
		LinesConstItr cend () {return _lines.cend();};
		/// @brief returns a reverse const_iterator to the end of the internal vector
		LinesConstItr crend () {return _lines.crend();};

		/// @brief find the closest vertex to a given point and the line to which it belongs
		/// @param target the target point
		/// @param type the distance calculation to use
		/// @returns a tuple containing an iterator to the line containing the closest vertex
		/// and an iterator to the closest vertex
		tuple<LinesConstItr, LineString::VertexConstItr> closestVertex (const Point& target, 
			const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double dist = target.distance(_lines.cbegin().cbegin(), type);
			LinesConstItr l = _lines.cbegin();
			LineString::VertexConstItr p = _lines.cbegin().cbegin();
			for (auto& a : _lines) {
				for (auto& b : a.getPositions()) {
					double newdist = target.distance(b, type);
					if (newdist < dist) {
						dist = newdist;
						l = a;
						p = b;
					}
				}
			}
			return make_tuple<l, p>;
		};

		/// @brief find the closest point to a target point on a member of this MultiLineString
		/// @param target the desired target point
		/// @param type the distance calculation to use
		/// @returns a tuple containing the Position of the closest point, a const_iterator
		/// to the LineString containing the closest point, and a const_iterator to the 
		/// first vertex on the segment of that LineString that contains the closest point
		tuple<Position, LinesConstItr, LineString::VertexConstItr> closestPoint (const Point& target, 
			const CourseTypeEnum type = CourseTypeEnum::RhumbLine) const {
			double dist = target.distance(_lines.cbegin().cbegin(), type);
			LinesConstItr l = _lines.cbegin();
			LineString::VertexConstItr p = _lines.cbegin().cbegin();
			Position best, newposn;
			for (auto& a : _lines) {
				for (auto& b : a.getPositions()) {
					if (a != b.end()) {
						newposn = point2segment(target, Point(b), Point(b+1), type);
						double newdist = target.distance(newposn, type);
						if (newdist < dist) {
							dist = newdist;
							best = newposn;
							l = a;
							p = b;
						}
					}
				}
			}
			return make_tuple(best, l, p);
		};

		/// @brief find the lengths of all of the LineStrings in this MultiLineString
		/// @param type the distance measurement to use
		/// @returns a vector of doubles representing the lengths of each LineString in meters
		vector<double> lengths (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			vector<double> result;
			for (auto& a : _lines) result.emplace_back(a.length(type));
			return result;
		};
		
		/// @brief find the sum of the lengths of all of the LineStrings in this MultiLineString
		/// @param type the distance measurement to use
		/// @returns the total length of all the contained LineStrings in meters
		double totalLength (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			double result = 0;
			for (auto& a : _lines) result += a.length(type);
			return result;
		};

		/// @brief check whether all of the LineStrings are valid
		/// @returns true if all of the member LineStrings are valid
		bool isValid () {
			bool result = true;
			for (auto &a : _lines) {
				result &= a.isValid();
			}
			return result;
		};

	private:
		vector<LineString> _lines;
		bool makebox () {
			if (!isValid()) return false;
			double maxlat = getlat(_lines.crbegin().crbegin());
			double minlat = maxlat;
			double maxlon = getlon(_lines.crbegin().crbegin());
			double minlon = maxlon
			for (auto& a: _lines) {
				for (auto& b : a) {
					if (maxlat < getlat(b)) maxlat = getlat(b);
					if (minlat > getlat(b)) minlat = getlat(b);
					if (maxlon < getlon(b)) maxlon = getlon(b);
					if (minlon > getlon(b)) minlon = getlon(b);
				}
			}
			_bbox = makeBBox(makePosition(minlon, minlat), 
							makePosition(maxlon, maxlat));
			return true;
		};
};

/** @class Polygon
 *
 * @brief represents a GeoJSON Polygon type
 *
 */
class Polygon : public GeometryRoot {
	public:
		/// @brief creates a blank LineString
		Polygon () : _mytype(GeometryType::Polygon) {};

		/// @brief create a polygon from the given set of positions or points
		/// Note that this will not close the polygon if it is not already closed
		Polygon (const vector<Point>& v) : _mytype(GeometryType::Polygon) {
			for (auto& a : v) _vertices.emplace_back(a.getPosition());
		};

		/// @brief create a polygon from the given set of positions or points
		/// Note that this will not close the polygon if it is not already closed
		Polygon (const vector<Position>& v) : _mytype(GeometryType::Polygon),
			_vertices(v) {};

		/// @brief create a polygon from the given MultiPoint
		/// Note that this will not close the polygon if it is not already closed
		Polygon (const MultiPoint& v) : _mytype(GeometryType::Polygon),
			_vertices(v.getPositions()) {};
			
		/// @brief create a polygon from the given LineString
		/// Note that this will not close the polygon if it is not already closed
		Polygon (const LineString& v) : _mytype(GeometryType::Polygon),
			_vertices(v.getPositions()) {};

		/// @brief is a const_iterator type for returning a target vertex
		typedef vector<Position>::const_iterator VertexConstItr;

		/// @brief populated this Polygon from a GeoJSON object. It will
		/// close the polygon if it is not already closed. 
		/// @param val pointer to a GeoJSON object
		/// @returns true if object is successfully populated
		bool load (const Value* val) {
			_mytype = GeometryType::Polygon;
			if (!val) return false;
			Value *typeptr = nullptr;
			Value *coordptr = nullptr;
			Value *bboxptr = nullptr;
			if (val->IsObject()) {
				typeptr = Pointer("/type").Get(val);
				coordptr= Pointer("/coordinates").Get(val);
				bboxptr= Pointer("/bbox").Get(val);
				if (!typeptr || !coordptr || (!typeptr->IsString() || 
					!(typeptr->GetString() == "Polygon") ||
					!coordptr->IsArray())) return false;
				if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
			} else if (val->IsArray()) {
				coordptr = val;
			} else return false;
			_vertices.clear();
			for (auto& v : coordptr->GetArray()) {
				if (v.IsArray()) _vertices.emplace_back(loadPosition(v));
			}
			if (!isClosed()) close();
			return isValid();
		};

		/// @brief creates a GeoJSON Polygon from this Polygon
		/// @param d the GeoJSON document to create the object in 
		/// @returns the packed GeoJSON Polygon
		Value pack (Document& d) const {
			Value b(kArrayType);
			Value v(kObject);
			auto a = d.GetAllocator();
			v.AddMember("type", "Polygon", a);
			v.AddMember("coordinates", packArray(d), a);
			v.AddMember("bbox", packBBox(getbox(), d), a);
			return v;
		};

		/// @brief creates a GeoJSON Polygon coordinate array from this Polygon
		/// @param d the GeoJSON document to create the array in 
		/// @returns the packed GeoJSON Polygon coordinates array
		Value packArray (Document& d) const {
			Value b(kArrayType);
			auto a = d.GetAllocator();
			for (auto& p : _vertices) {
				b.PushBack(packPosition(p, d), a);
			}
			return b;
		};

		/// @brief closes this polygon if it is not already closed by 
		/// emplacing a copy of the first element to the back
		/// @returns true if the polygon is already closed or it is successfully 
		/// closed. Returns false if the first element is invalid
		bool close () {
			if (isClosed()) return true;
			if (positionValid(_vertices.cbegin())) {
				_vertices.emplace_back(_vertices.cbegin());
			} else return false;
		};

		vector<Position>& getVertices () const {return _vertices;};
		Position& getVertex (const unsigned int& t) {return _vertices[t];};
		Point& getVertexPoint (const unsigned int& t) {return Point(_vertices[t]);};

		/// @brief returns a vector of Points corresponding to the members of the LineString
		vector<Point> getVertexPoints () const {
			vector<Point> v;
			for (auto p&: _vertices) v.emplace_back(Point(p));
			return v;
		};

		Polygon makeApproximateHull (const unsigned int k) {};
		Polygon makeHull () {};
		bool includes (const Point& p) const {};
		bool includes (const Position& p) const {};
		bool deleteVertex (const unsigned int& posn) {};
		Point closestVertex (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {};
		Point closestPoint (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {};
		double perimeter () const {};
		double area () const {};

		/// @brief checks whether all vertices are valid positions and whether the polygon is closed
		bool isValid () {
			bool result = true;
			for (auto& a : _vertices) {
				result &= positionValid(a);
			}
			return (result & isClosed());
		};

		/// @brief checks whether the polygon is closed
		bool isClosed () {
			Position first = _vertices.cbegin();
			Position last = _vertices.cend();
			return (positionValid(first) && positionValid(last) && 
					(getlat(first) == getlat(last)) && 
					(getlon(first) == getlon(last)));
		};

	private:
		vector<Position> _vertices;
		bool makebox () {
			if (!isValid()) return false;
			double maxlat = getlat(_vertices.crbegin());
			double minlat = maxlat;
			double maxlon = getlon(_vertices.crbegin());
			double minlon = maxlon
			for (auto& a: _vertices) {
				if (maxlat < getlat(a)) maxlat = getlat(a);
				if (minlat > getlat(a)) minlat = getlat(a);
				if (maxlon < getlon(a)) maxlon = getlon(a);
				if (minlon > getlon(a)) minlon = getlon(a);
			}
			_bbox = makeBBox(makePosition(minlon, minlat), 
							makePosition(maxlon, maxlat));
			return true;
		};
};

class MultiPolygon : public GeometryRoot {
	public:
		MultiPolygon () {};
		MultiPolygon (const Value& val) {};
		MultiPolygon (const vector<LineString>& v) {};
		MultiPolygon (const vector<MultiPoint>& v) {};
		MultiPolygon (const vector<Polygon>& v) {};
		bool load (const Value& val) {};
		Value pack () const {};

		vector<Polygon>& getPolygons () const {return _polygons;};
		LineString& getPolygon (const unsigned int& l) const {return _polygons[l];};
		bool insert (const LineString& p, const unsigned int& t) const {};
		bool insert (const MultiPoint& p, const unsigned int& t) const {};
		bool insert (const Polygon& p, const unsigned int& t) const {};
		bool deleteLine (const unsigned int& t) {};
		bool includes (const Point& p) const {};
		bool includes (const Position& p) const {};
		Point closestVertex (const Point& target, Polygon& p, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {};
		Point closestPoint (const Point& target, Polygon& p, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {};
		vector<double> perimeters () const {};
		double totalPerimeters () const {};

	private:
		vector<Polygon> _polygons;
		bool makebox () {};
};

class GeometryCollection : public GeometryRoot {};

#endif /* GEOGEOMETRY_H */