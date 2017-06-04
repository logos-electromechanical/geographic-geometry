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
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>
#include "rapidjson/rapidjson.h"
#include "rapidjson/document.h"
#include "rapidjson/pointer.h"

using namespace GeographicLib;
using namespace std;
using namespace rapidjson;

namespace GeoGeometry {

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
	inline double& getlon (Position& p) {return get<0>(p);};
	inline const double& getlon (const Position& p) {return get<0>(p);};
	/// @brief retrieve the the latitude from a Position tuple
	inline double& getlat (Position& p) {return get<1>(p);};
	inline const double& getlat (const Position& p) {return get<1>(p);};
	/// @brief retrieve the the elevation from a Position tuple
	inline double& getele (Position& p) {return get<2>(p);};
	inline const double& getele (const Position& p) {return get<2>(p);};

	/// @brief returns a Position tuple with the given longitude, latitude, and elevation, in that order.
	/// The ele parameter is optional and defaults to 0.0
	inline Position makePosition (const double& lon, const double& lat, const double& ele = 0.0) {
		return make_tuple(lon, lat, ele);
	};

	/// @brief creates a Postion tuple from a GeoJSON array
	inline Position loadPosition (Value& v) {
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
		v.PushBack(Value().SetDouble(getlon(p)), d.GetAllocator());
		v.PushBack(Value().SetDouble(getlat(p)), d.GetAllocator());
		v.PushBack(Value().SetDouble(getele(p)), d.GetAllocator());
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

	/// @brief determine whether a collinear point lies in a segment
	/// This function was derived from http://geomalgorithms.com/a05-_intersect-1.html
	/// Note that this has only been tested to +/- 80 deg latitude. It will also fail
	/// for segments that cross the antimeridian due to ambiguity. 
	/// @param p the point to test
	/// @param s0 the start of the segment to test
	/// @param s1 the end of the segment to test
	/// @returns true if the point lies in the segment
	inline bool inSegment(const Position& p, const Position& s0, const Position& s1) {
		if (fabs(getlat(s0) - getlat(s1)) > smallnum) {		// the segment is not e-w
			if ((getlat(s0) <= getlat(p)) && (getlat(p) <= getlat(s1))) return true;
			if ((getlat(s0) >= getlat(p)) && (getlat(p) >= getlat(s1))) return true;
		} else {											// segment is e-w
			if ((getlon(s0) <= getlon(p)) && (getlon(p) <= getlon(s1))) return true;
			if ((getlon(s0) >= getlon(p)) && (getlon(p) >= getlon(s1))) return true;
		}
		return false;
	}

	/// @brief retrieves the southwest corner of the given BoundingBox
	inline Position& getsw (BoundingBox& bb) {return get<0>(bb);};
	inline const Position& getsw (const BoundingBox& bb) {return get<0>(bb);};
	/// @brief retrieves the northeast corner of the given BoundingBox
	inline Position& getne (BoundingBox& bb) {return get<1>(bb);};
	inline const Position& getne (const BoundingBox& bb) {return get<1>(bb);};
	/// @brief creates a bounding box from the given positions
	inline BoundingBox makeBBox (const Position& sw, const Position& ne) {
		return make_tuple(sw, ne);
	};

	/// @brief creates a BoundBox from a GeoJSON bounding box array
	inline BoundingBox loadBBox (Value& v) {
		double swlon = NAN;
		double swlat = NAN;
		double nelon = NAN;
		double nelat = NAN;
		Value *vp = nullptr;
		if (v.IsArray()) {
			vp = Pointer("/0").Get(v);
			if (vp && vp->IsDouble()) swlon = vp->GetDouble();
			vp = Pointer("/1").Get(v);
			if (vp && vp->IsDouble()) swlat = vp->GetDouble();
			vp = Pointer("/2").Get(v);
			if (vp && vp->IsDouble()) nelon = vp->GetDouble();
			vp = Pointer("/3").Get(v);
			if (vp && vp->IsDouble()) nelat = vp->GetDouble();
		}
		return makeBBox(makePosition(swlon, swlat), makePosition(nelon, nelat));
	};

	/// @brief creates a GeoJSON bounding box array from the given BoundingBox in the given Document
	inline Value packBBox (const BoundingBox& b, Document& d) {
		Value v(kArrayType);
		Value i(kNumberType);
		v.PushBack(i.SetDouble(getlon(getsw(b))), d.GetAllocator());
		v.PushBack(i.SetDouble(getlat(getsw(b))), d.GetAllocator());
		v.PushBack(i.SetDouble(getlon(getne(b))), d.GetAllocator());
		v.PushBack(i.SetDouble(getlat(getne(b))), d.GetAllocator());
		return v;
	};

	/// @brief checks whether a BoundingBox refers to a real bounding box
	inline bool bboxValid (const BoundingBox& b) {
		return(positionValid(getsw(b)) &&
				positionValid(getne(b)) &&
				(getlat(getne(b)) > getlat(getsw(b))));
	};

	///@brief checks whether a given Point or Location is within the giving bounding box
	inline bool bboxContains (const BoundingBox& b, const Position& p) {
		if (!bboxValid(b)) return false;
		if (getlon(getne(b)) > getlon(getsw(b))) {	// this checks that the bounding box does not cross the antimeridian
			return ((getlon(p) <= getlon(getne(b))) && (getlon(p) >= getlon(getsw(b))) &&
					(getlat(p) <= getlat(getne(b))) && (getlat(p) >= getlat(getsw(b))));
		} else { // this handles the case where the bounding box crosses the antimeridian
			return (((getlon(p) <= getlon(getne(b))) || (getlon(p) >= getlon(getsw(b)))) &&
					(getlat(p) <= getlat(getne(b))) && (getlat(p) >= getlat(getsw(b))));
		}
	}
 
	// utility functions
	inline double deg2rad (double deg) { return deg * ( M_PI / 180.0 ); }	/**< Convert degrees to radians */
	inline double rad2deg (double rad) { return rad * ( 180.0 / M_PI ); }	/**< Convert radians to degrees */

	/// @brief returns the number of meters per degree of latitude at given latitude
	inline double metersPerDegreeLat (const double& lat) {
		double phi = deg2rad(lat);
		return (111132.954 - 559.822 * cos(2 * phi) + 1.175 * cos(4 * phi));
	};

	/// @brief returns the number of meters per degree of longitude at the given latitude
	inline double metersPerDegreeLon (const double& lat) {
		double phi = deg2rad(lat);
		double a = 6378137.0;	// WGS84 semi-major axis
		double e = 0.0818192;	// WGS84 flattening
		return ((M_PI * a * cos(phi)) / (180 * sqrt(1 - (e * e * pow(sin(phi),2)))));
	};

	/// @brief compare whether two positions are equal
	inline bool operator== (Position& a, Position& b) {
		if ((getlat(a) == getlat(b)) && (getlon(a) == getlon(b))) return true;
		return false;
	}

	/// @brief compare whether two positions are unequal
	inline bool operator!= (Position& a, Position& b) {
		if (a == b) return false;
		return true;
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
			bool load(Value* d) {
				if (!d) return false;
				Value *xval = Pointer("/x").Get(*d);
				Value *yval = Pointer("/y").Get(*d);
				if (xval && yval && xval->IsDouble() && yval->IsDouble()) {
					_x = xval->GetDouble();
					_y = yval->GetDouble();
				} else return false;
				return isValid();
			}

			/// @brief creates a json object in Document d from this vector
			/// @param d a reference to a json Document to create the Value in
			/// @returns a json object in Document d
			Value pack (Document& d) const {
				Value v(kObjectType);
				v.AddMember("x", _x, d.GetAllocator());
				v.AddMember("y", _y, d.GetAllocator());
				return v;
			}

			/// @brief returns true if the vector is valid, i.e. has finite members
			bool isValid () const {return ((isfinite(_x)) && (isfinite(_y)));};
			
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
			double inline x () const {return _x;};											/**< Get the x value */
			void inline x (double x) {_x=x;};											/**< Set the x value */
			double inline y () const {return _y;};											/**< Get the y value */
			void inline y (double y) {_y=y;};											/**< Set the y value */
			double inline mag () const {return sqrt((_x*_x)+(_y*_y));};						/**< Get the magnitude */
			void inline mag (double _mag) {*this *= _mag/mag();};						/**< Set the magnitude */
			double inline angleRad () const {return atan2(_y,_x);}; 							/**< Get the angle in radians */
			void angleRad (double _ang) {this->rotateRad(_ang - this->angleRad());};	/**< Set the angle in radians */
			double inline angleDeg () const {return rad2deg(atan2(_y,_x));}; 					/**< Get the angle in degrees */
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
			
			friend TwoVector inline operator* (const double r, TwoVector l) {			/**< Scalar multiplication */
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

	// forward declarations of subclasses to make the factory work
	// class Point;
	// class MultiPoint;
	// class LineString;
	// class MultiLineString;
	// class Polygon;
	// class MultiPolygon;
	// class GeometryCollection;

	/** @class GeometryRoot
	 *
	 * @brief is a base class for representing GeoJSON geometry objects
	 *
	 */
	class GeometryRoot {
		public:
			/// @brief creates a blank object
			GeometryRoot () {};
			// @brief creates a geometry object from a GeoJSON object
			// @param val reference to a GeoJSON object
			GeometryRoot (Value& val) {load(&val);};
			// @brief creates a geometry object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			GeometryRoot (Value* val) {load(val);};
			// @brief populates a geometry object from a GeoJSON object
			// @param val reference to a GeoJSON object
			// @returns true if object is successfully populated
			bool load (Value& val) {return load(&val);};
			// @brief populated a geometry object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			// @returns true if object is successfully populated
			virtual bool load (Value* val) {return false;};
			// @brief creates a GeoJSON object from this object
			// @param d the GeoJSON document to create the object in 
			// @returns the packed GeoJSON object
			virtual Value pack (Document& d) {return Value();};
			// @brief creates a GeoJSON coordinates array from this object
			// @param d the GeoJSON document to create the array in 
			// @returns the packed GeoJSON coordinates array
			virtual Value packArray(Document& d) {return Value();};

			/// @brief fetches a constant reference to this object's type
			const GeometryType& gettype () const {return _mytype;};

			/// @brief fetches the bounding box of the current object, and generates it if it has not yet been generated
			const BoundingBox& getbox () {
				if (!_hasbox) _hasbox = makebox();
				return _bbox;
			}

			/// @brief checks whether the current object is valid
			virtual bool isValid () const = 0;

			static GeometryRoot* factory (Value* val);

		protected:
			GeometryRoot (GeometryType t) : _mytype(t) {};
			virtual bool makebox () = 0;
			bool 			_hasbox = false;
			BoundingBox 	_bbox;
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
			Point () : GeometryRoot(GeometryType::Point), 
				_posn(makePosition(NAN, NAN, 0)) {};

			/// @brief creates a Point at the given position
			/// @param posn the desired position
			Point (const Position& posn) : 
				_posn(posn), 
				GeometryRoot(GeometryType::Point) {};

			/// @brief creates a Point at the given position
			/// @param lon the desired longitude (decimal degrees)
			/// @param lat the desired latitude (decimal degrees)
			/// @param ele the desired elevation (meters above datum)
			Point (const double lon, const double lat, const double ele = 0.0)  : 
				GeometryRoot(GeometryType::Point),
				_posn(makePosition(lon, lat, ele)) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			Point (Value& val) {load(&val);};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			Point (Value* val) {load(val);};

			/// @brief populated this point from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			bool load (Value* val) {
				if (!val) return false;
				_mytype = GeometryType::Point;
				if (val->IsObject()) {
					Value *typeptr = Pointer("/type").Get(*val);
					Value *coordptr= Pointer("/coordinates").Get(*val);
					if (typeptr && coordptr) {
						if (typeptr->IsString() && 
							(string(typeptr->GetString()) == "Point") &&
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
				Value v(kObjectType);
				v.AddMember("type", "Point", d.GetAllocator());
				v.AddMember("coordinates", packPosition(_posn, d), d.GetAllocator());
				return v;
			};

			/// @brief creates a GeoJSON Point coordinate array from this Point
			/// @param d the GeoJSON document to create the array in 
			/// @returns the packed GeoJSON Point coordinates array
			Value packArray (Document& d) {return packPosition(_posn, d);};

			/// @brief returns a reference to the embedded position
			Position& getPosition () {return _posn;};

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
						Geodesic::WGS84().Inverse(getlat(), getlon(), GeoGeometry::getlat(dest), GeoGeometry::getlon(dest), azi1, azi2);
						return azi1;
						break;
					}
					case CourseTypeEnum::RhumbLine:	{
						Rhumb::WGS84().Inverse(getlat(), getlon(), GeoGeometry::getlat(dest), GeoGeometry::getlon(dest), dist, azi1);
						return azi1;
						break;
					}
					case CourseTypeEnum::Approximate: {
						double deltalat = metersPerDegreeLat(getlat())*(getlat() - GeoGeometry::getlat(dest));
						double deltalon = metersPerDegreeLon(getlat())*(getlon() - GeoGeometry::getlon(dest));
						double result = rad2deg(atan2(-deltalon, -deltalat));
						if (result > 360.0) result -= 360.0;
						if (result < 0.0) result += 360;
						return result;
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
						Geodesic::WGS84().Inverse(getlat(), getlon(), GeoGeometry::getlat(dest), GeoGeometry::getlon(dest), dist, azi1, azi2);
						return dist;
						break;
					}
					case CourseTypeEnum::RhumbLine:	{
						Rhumb::WGS84().Inverse(getlat(), getlon(), GeoGeometry::getlat(dest), GeoGeometry::getlon(dest), dist, azi1);
						return dist;
						break;
					}
					case CourseTypeEnum::Approximate: {
						double deltalat = metersPerDegreeLat(getlat())*(getlat() - GeoGeometry::getlat(dest));
						double deltalon = metersPerDegreeLon(getlat())*(getlon() - GeoGeometry::getlon(dest));
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
				if ((!this->isValid()) || (!projection.isValid())) return makePosition(NAN,NAN,0.0);
				Position p = makePosition(getlon(),getlat(),0.0);
				switch (type) {
					case CourseTypeEnum::GreatCircle: {
						Geodesic::WGS84().Direct(getlat(), getlon(), projection.angleDeg(), projection.mag(), GeoGeometry::getlat(p), GeoGeometry::getlon(p));
						break;
					}
					case CourseTypeEnum::RhumbLine: {
						Rhumb::WGS84().Direct(getlat(), getlon(), projection.angleDeg(), projection.mag(), GeoGeometry::getlat(p), GeoGeometry::getlon(p));
						break;
					}
					case CourseTypeEnum::Approximate: {
						GeoGeometry::getlon(p) += projection.y()/metersPerDegreeLon(getlat());
						GeoGeometry::getlat(p) += projection.x()/metersPerDegreeLat(getlat());
						break;
					}
					default:
						break;
				}
				// normalize longitude
				if (GeoGeometry::getlon(p) > 180.0) GeoGeometry::getlon(p) -= 360;
				if (GeoGeometry::getlon(p) < -180.0) GeoGeometry::getlon(p) += 360;
				return p;
			};

			bool isValid () const {
				if (_mytype == GeometryType::Point) return (positionValid(_posn));
				return false;
			};
			
			const double& getlon () const {return GeoGeometry::getlon(_posn);};
			const double& getlat () const {return GeoGeometry::getlat(_posn);};
			const double& getele () const {return GeoGeometry::getele(_posn);};

			friend bool inline operator== (Point& a, Point& b) {
				return (a.getPosition() == b.getPosition());
			}

			friend bool inline operator!= (Point& a, Point& b) {
				return (a.getPosition() != b.getPosition());
			}

		private:
			Position _posn;
			Point (GeometryType t) : GeometryRoot(t) {};
			bool makebox () {
				get<0>(_bbox) = _posn;
				get<1>(_bbox) = _posn;
				return true;
			};
	};

	/// @brief return the point on a line segment closest to a point
	/// Note that this function may behave strangely with great circle distances
	/// This function was derived from http://geomalgorithms.com/a02-_lines.html
	/// @param p is the point external to the line segment
	/// @param s0 is one end of the line segment
	/// @param s1 is the other end of the line segment
	/// @param type the distance type used in this calculation
	/// @returns the Position on the line segment closest to p
	Position point2segment(const Point& p, const Point& s0, const Point& s1, 
		const CourseTypeEnum type = CourseTypeEnum::RhumbLine);

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
	unsigned int segmentsIntersect (const Point& s10, const Point& s11, 
											const Point& s20, const Point& s21,
											Point& intersect0, Point& intersect1,
											CourseTypeEnum type = CourseTypeEnum::RhumbLine);

	/// @brief determine whether a point is to the left of a given line
	/// Note: this may not work at high latitudes (tested to +/- 80), and gives wrong 
	/// answers for segments that cross the antimeridian due to ambiguity. 
	/// @param p the Position of the point to test
	/// @param s0 a Position on the line to test
	/// @param s1 another Position on the line to test
	/// @returns > 0 if p is to the left of the line s
	/// == 0 if p is on the line s
	/// < 0 if p is to the right of line s
	inline double isLeft (const Position& p, const Position& s0, const Position& s1) {
		return (((getlon(s1) - getlon(s0)) * (getlat(p) - getlat(s0))) - 
				((getlon(p) - getlon(s0)) * (getlat(s1) - getlat(s0))));
	};

	/// @brief comparison function for sorting Positions by latitude
	inline bool latlt (const Position& a, const Position& b) {
		if (getlat(a) < getlat(b)) return true;
		return false;
	}

	/// @brief comparison function for sorting Positions by longitude
	inline bool lonlt (const Position& a, const Position& b) {
		if (getlon(a) < getlon(b)) return true;
		return false;
	}

	/** @class MultiPoint
	 *
	 * @brief represents a GeoJSON MultiPoint type
	 * 
	 */
	class MultiPoint : public GeometryRoot {
		public:
			/// @brief creates a blank MultiPoint
			MultiPoint () : GeometryRoot(GeometryType::MultiPoint) {};

			/// @brief creates a MultiPoint object from the given positions
			/// @param vector of Positions to create MultiPoint from 
			MultiPoint (const vector<Position>& v) : 
				GeometryRoot(GeometryType::MultiPoint),
				_points(v) {};

			/// @brief creates a MultiPoint object from the given points
			/// @param vector of Points to create MultiPoint from
			MultiPoint (vector<Point>& v) : GeometryRoot(GeometryType::MultiPoint) {
				for (auto& a : v) {
					_points.emplace_back(a.getPosition());
				}
			};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			MultiPoint (Value& val) : GeometryRoot(val) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			MultiPoint (Value* val) : GeometryRoot(val) {};

			/// @brief is a const_iterator type for returning a target position
			typedef vector<Position>::const_iterator PositionConstItr;

			/// @brief populated this MultiPoint from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			virtual bool load (Value* val) {
				if (!val) return false;
				_mytype = GeometryType::MultiPoint;
				Value *typeptr = nullptr;
				Value *coordptr = nullptr;
				Value *bboxptr = nullptr;
				if (val->IsObject()) {
					typeptr = Pointer("/type").Get(*val);
					coordptr = Pointer("/coordinates").Get(*val);
					bboxptr = Pointer("/bbox").Get(*val);
					if (!typeptr || !coordptr || (!typeptr->IsString() || 
						!(string(typeptr->GetString()) == "MultiPoint") ||
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
			virtual Value pack (Document& d) {
				Value b(kArrayType);
				Value v(kObjectType);
				v.AddMember("type", "MultiPoint", d.GetAllocator());
				v.AddMember("coordinates", packArray(d), d.GetAllocator());
				v.AddMember("bbox", packBBox(getbox(), d), d.GetAllocator());
				return v;
			};

			/// @brief creates a GeoJSON MultiPoint coordinate array from this Point
			/// @param d the GeoJSON document to create the array in 
			/// @returns the packed GeoJSON MultiPoint coordinates array
			Value packArray (Document& d) {
				Value b(kArrayType);
				for (auto& p : _points) {
					b.PushBack(packPosition(p, d), d.GetAllocator());
				}
				return b;
			};

			/// @brief return a reference to the members of the MultiPoint
			vector<Position>& getPositions () {return _points;};

			/// @brief returns a reference to a given member
			/// @params t the index to return
			/// @returns the position at index t
			Position& getPosition (const unsigned int& t) {return _points[t];};

			/// @brief returns a Point created from a given member
			/// @params t the index to return
			/// @returns a Point created from the position at index t
			Point getPoint (const unsigned int& t) {return Point(_points[t]);};

			/// @brief returns a vector of Points corresponding to the members of the MultiPoint
			vector<Point> getPoints () {
				vector<Point> v;
				for (auto& p: _points) v.emplace_back(Point(p));
				return v;
			};

			/// @brief insert a new Point after the given index
			/// @param p the Point to insert
			/// @param posn the index to insert after
			/// @returns true if successful
			bool insert (Point& p, const unsigned int& posn) {return insert(p.getPosition(), posn);};

			/// @brief insert a new Position after the given index
			/// @param p the Position to insert
			/// @param posn the index to insert after
			/// @returns true if successful
			bool insert (const Position& p, const unsigned int& posn) {
				if (posn >= _points.size()) {
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
			bool insert (PositionConstItr& i, Point& p) {
				return this->insert(i, p.getPosition());
			};

			/// @brief insert a new Position after the given const_interator
			/// @param i const_iterator of the insertion point
			/// @param p the Position to insert to insert
			/// @returns true if successful
			bool insert (const PositionConstItr& i, const Position& p) {
				_points.emplace(i, p);
				return true;
			};

			/// @brief remove the point at the given index
			/// @param posn index to remove
			/// @returns true if successful
			bool deletePoint (const unsigned int& posn) {
				if (posn < _points.size()) {
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

			/// @brief returns a const_iterator to the end of the internal vector
			PositionConstItr cend () {return _points.cend();};


			/// @brief find the closest point to the given Point
			/// @param target point to find the closest point in the MultiPoint to 
			/// @param type type of the course to use for distance finding
			PositionConstItr closestVertex (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				double dist = 1000000000; 	// One million km, larger than any distance on Earth
				PositionConstItr p = _points.cbegin();
				for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
					double newdist = target.distance(*a, type);
					if (newdist < dist) {
						dist = newdist;
						p = a;
					} 
				}
				return p;
			};

			/// @brief check whether this MultiPoint is valid
			virtual bool isValid () const {
				bool valid = true;
				for (auto& p : _points) {
					valid &= positionValid(p);
				}
				return valid;
			}

		protected:
			vector<Position> _points;
			MultiPoint (GeometryType t) : GeometryRoot(t) {};
			MultiPoint (const vector<Position>& v, GeometryType t) : _points(v), GeometryRoot(t) {};
			MultiPoint (vector<Point>& v, GeometryType t) : GeometryRoot(t) {
				for (auto& a : v) {
					_points.emplace_back(a.getPosition());
				}
			};
			bool makebox () {
				if (!isValid()) return false;
				// set the initial values such that all valid points are less than the mins
				// and greater than the maxes
				double maxlat = -100.0;
				double minlat = 100.0;
				double maxlon = -200.0;
				double minlon = 200.0;
				for (auto& a: _points) {
					if (getlat(a) > maxlat) maxlat = getlat(a);
					if (getlat(a) < minlat) minlat = getlat(a);
					if (getlon(a) > maxlon) maxlon = getlon(a);
					if (getlon(a) < minlon) minlon = getlon(a);
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
	class LineString : public MultiPoint {
		public:
			/// @brief creates a blank LineString
			LineString () : MultiPoint(GeometryType::LineString) {};

			/// @brief creates a LineString object from the given positions
			/// @param vector of Positions to create LineString from 
			LineString (vector<Position>& v) : 
				MultiPoint(v, GeometryType::LineString) {};

			/// @brief creates a LineString object from the given points
			/// @param vector of Points to create LineString from
			LineString (vector<Point>& v) : MultiPoint(GeometryType::LineString) {
				for (auto& a : v) {
					_points.push_back(a.getPosition());
				}
			};

			/// @brief creates a LineString object from the given MultiPoint
			/// @param v MultiPoint to create LineString from 
			LineString (MultiPoint& v) : 
				MultiPoint(v.getPositions(), GeometryType::LineString) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			LineString (Value& val) : MultiPoint(val) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			LineString (Value* val) : MultiPoint(val) {};

			/// @brief populated this LineString from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			bool load (Value* val) {
				MultiPoint::load(val);
				_mytype = GeometryType::LineString;
				return isValid();
			};

			/// @brief creates a GeoJSON LineString from this LineString
			/// @param d the GeoJSON document to create the object in 
			/// @returns the packed GeoJSON LineString
			Value pack (Document& d) {
				Value v = MultiPoint::pack(d);
				Pointer("/type").Set(d, "LineString");
				return v;
			};

			/// @brief find the closest position on the LineString to the given Point
			/// @param target the given Point
			/// @param type the type of course to use for distanc measuring etc.
			/// @returns a tuple containing the desired position and an iterator to 
			/// to the start of the line segment that contains the closest position 
			tuple<Position, PositionConstItr> closestPoint (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				double bestdist = 1000000000;
				PositionConstItr p = _points.cbegin();
				Position best, newposn;
				for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
					if (a != _points.end()) {
						newposn = point2segment(target, Point(*a), Point(*(a + 1)), type);
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
			/// @param s0 start point of the line segment to intersect
			/// @param s1 end point of the line segment to intersect
			tuple<MultiPoint, vector<PositionConstItr>> intersection (const Point& s0, const Point& s1, const CourseTypeEnum type = CourseTypeEnum::RhumbLine) {
				Point o0, o1;
				vector<PositionConstItr> outitr;
				MultiPoint out;
				for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
					if (a != _points.end()) {
						if (segmentsIntersect(Point(*a), Point(*(a+1)), s0, s1, o0, o1, type)) {
							outitr.emplace_back(a);
							out.insert(out.cend(), *a);
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
				for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
					if (a != _points.end()) {
						Point p(*a);
						result += p.distance(*(a + 1), type);
					}
				}
				return result;
			};

		protected:
			LineString (GeometryType t) : MultiPoint(t) {};
			LineString (const vector<Position> v, GeometryType t) : MultiPoint(v, t) {};
			LineString (vector<Point> v, GeometryType t) : MultiPoint(v, t) {};
			LineString (MultiPoint v, GeometryType t) : MultiPoint(v.getPositions(), t) {};

	};

	/** @class Polygon
	 *
	 * @brief represents a GeoJSON Polygon type
	 *
	 */
	class Polygon : public LineString {
		public:
			/// @brief creates a blank LineString
			Polygon () : LineString(GeometryType::Polygon) {};

			/// @brief create a polygon from the given set of positions or points
			/// Note that this will not close the polygon if it is not already closed
			Polygon (const vector<Point>& v) : LineString(v, GeometryType::Polygon) {};

			/// @brief create a polygon from the given set of positions or points
			/// Note that this will not close the polygon if it is not already closed
			Polygon (const vector<Position>& v) : LineString(v, GeometryType::Polygon) {};

			/// @brief create a polygon from the given MultiPoint
			/// Note that this will not close the polygon if it is not already closed
			Polygon (MultiPoint& v) : LineString(v, GeometryType::Polygon) {};
				
			/// @brief create a polygon from the given LineString
			/// Note that this will not close the polygon if it is not already closed
			Polygon (LineString& v) : LineString(v.getPositions(), GeometryType::Polygon) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			Polygon (Value& val) : LineString(val) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			Polygon (Value* val) : LineString(val) {};

			/// @brief populated this Polygon from a GeoJSON object. It will
			/// close the polygon if it is not already closed. 
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			bool load (Value* val) {
				MultiPoint::load(val);
				_mytype = GeometryType::Polygon;
				if (!isClosed()) close();
				return isValid();
			};

			/// @brief creates a GeoJSON Polygon from this Polygon
			/// @param d the GeoJSON document to create the object in 
			/// @returns the packed GeoJSON Polygon
			Value pack (Document& d) {
				Value v = MultiPoint::pack(d);
				Pointer("/type").Set(d, "Polygon");
				return v;
			};

			/// @brief closes this polygon if it is not already closed by 
			/// emplacing a copy of the first element to the back
			/// @returns true if the polygon is already closed or it is successfully 
			/// closed. Returns false if the first element is invalid
			bool close () {
				if (isClosed()) return true;
				if (positionValid(*(_points.cbegin()))) {
					_points.emplace_back(*(_points.cbegin()));
					return true;
				} 
				return false;
			};

			/// @brief Make an approximate convex hull with no more than k points using 
			/// the BFP algorithm described at http://geomalgorithms.com/a11-_hull-2.html.
			/// @param k the number of points in the resulting approximate convex hull.
			/// If k is larger than 4 times the number of points in this Polygon, it 
			/// it will be clipped to that number. 
			/// @returns the approximate convex hull
			Polygon makeApproximateHull (const unsigned int k) {
				return Polygon();
			};

			/// @brief Find a convex hull using Andrew's Monotone Chain algorithm 
			/// described at http://geomalgorithms.com/a10-_hull-1.html
			/// @returns the convex hull
			Polygon makeHull () {
				if (!isValid()) return Polygon();
				vector<Position> mypoly = _points;			// copy the points vector so we can sort it
				vector<Position> hull;
				sort(mypoly.begin(), mypoly.end(), latlt);
				sort(mypoly.begin(), mypoly.end(), lonlt);

				return Polygon();
			};

			bool includes (Point& p) {return includes(p.getPosition());};
			bool includes (Position& p) {
				if (bboxContains(_bbox, p)) {

				} 
				return false;
			};

			/// @brief checks whether all vertices are valid positions and whether the polygon is closed
			bool isValid () {
				return (MultiPoint::isValid() & isClosed());
			};

			/// @brief checks whether the polygon is closed
			bool isClosed () {
				Position first = *(_points.cbegin());
				Position last = *(_points.cend());
				return (positionValid(first) && positionValid(last) && 
						(getlat(first) == getlat(last)) && 
						(getlon(first) == getlon(last)));
			};
	};

	template<typename T> class Collection : public GeometryRoot {
		public:
			Collection () = delete;	

			/// @brief creates a GeoJSON MultiLineString coordinate array from this Point
			/// @param d the GeoJSON document to create the array in 
			/// @returns the packed GeoJSON MultiLineString coordinates array
			virtual Value packArray (Document& d) {
				Value b(kArrayType);
				for (auto p : _elements) {
					b.PushBack(p->packArray(d), d.GetAllocator());
				}
				return b;
			};

			/// @brief return a reference to the members of the Collection
			vector<T> getElements () {return _elements;};

			/// @brief returns a reference to a given member
			/// @params t the index to return
			/// @returns the element at index t
			T getElement (const unsigned int& t) {return _elements[t];};

			/// @brief insert a new element after the given index
			/// @param p the element to insert
			/// @param posn the index to insert after
			/// @returns true if successful
			bool insert (const T p, const unsigned int& posn) {
				if (posn >= _elements.size()) {
					_elements.emplace_back(p);
				} else {
					_elements.emplace((_elements.begin()+posn), p);
				}
				return true;
			};

			/// @brief insert a new element after the given const_interator
			/// @param i const_iterator of the insertion point
			/// @param p the element to insert
			/// @returns true if successful
			bool insert (const typename vector<T>::const_iterator i, T p) {
				_elements.emplace(i, p);
				return true;
			};

			/// @brief remove the element at the given index
			/// @param posn index to remove
			/// @returns true if successful
			bool deleteElement (const unsigned int& posn) {
				if (posn < _elements.size()) {
					_elements.erase(_elements.begin() + posn);
				} else return false;
				return true;
			};

			/// @brief remove the element at the given iterator
			/// @param p iterator to remove
			/// @returns true if successful
			bool deleteElement (typename vector<T>::const_iterator p) {
				_elements.erase(p);
				return true;
			}

			/// @brief returns the size of the internal vector
			int size () {return _elements.size();};
			/// @brief returns a const_iterator to the beginning of the internal vector
			typename vector<T>::const_iterator cbegin () {return _elements.cbegin();};

			/// @brief returns a const_iterator to the end of the internal vector
			typename vector<T>::const_iterator cend () {return _elements.cend();};

			/// @brief check whether this Collection is valid
			virtual bool isValid () const {
				bool valid = true;
				for (auto& p : _elements) {
					valid &= p->isValid();
				}
				return valid;
			}

			/// @brief find the closest vertex to a given point and the line to which it belongs
			/// @param target the target point
			/// @param type the distance calculation to use
			/// @returns a tuple containing an iterator to the line containing the closest vertex
			/// and an iterator to the closest vertex
			virtual tuple<MultiPoint::PositionConstItr, typename vector<T>::const_iterator> 
			closestVertex (const Point& target, 
				CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
					double dist = 1000000000;
					typename vector<T>::const_iterator best;
					MultiPoint::PositionConstItr bestPosn, newPosn;
					for (auto a = _elements.cbegin(); a != _elements.cend(); ++a) {
						newPosn = (*a)->closestVertex(target, type);
						double newdist = target.distance(*newPosn, type);
						if (newdist < dist) {
							bestPosn = newPosn;
							best = a;
						}
					}
					return make_tuple(bestPosn, best);
			};

			/// @brief find the closest point to a target point on a member of this MultiLineString
			/// @param target the desired target point
			/// @param type the distance calculation to use
			/// @returns a tuple containing the Position of the closest point, a const_iterator
			/// to the LineString containing the closest point, and a const_iterator to the 
			/// first vertex on the segment of that LineString that contains the closest point
			tuple<Position, MultiPoint::PositionConstItr, typename vector<T>::const_iterator> 
			closestPoint (const Point& target, 
				CourseTypeEnum type = CourseTypeEnum::RhumbLine) const {
					double dist = 1000000000;
					typename vector<T>::const_iterator best;
					tuple<Position, MultiPoint::PositionConstItr> bestPosn, newPosn;
					for (auto a = _elements.cbegin(); a != _elements.cend(); ++a) {
						newPosn = (*a)->closestPoint(target, type);
						double newdist = target.distance(get<0>(newPosn), type);
						if (newdist < dist) {
							bestPosn = newPosn;
							best = a;
						}
					}
					return make_tuple(get<0>(bestPosn), get<1>(bestPosn), best);
			}

		protected:
			Collection (GeometryType me) : GeometryRoot(me) {};
			Collection (Value& val, GeometryType me) : GeometryRoot(me) {load(&val);};
			Collection (Value* val, GeometryType me) : GeometryRoot(me) {load(val);};
			Collection (vector<T>& v, GeometryType me) : GeometryRoot(me) {
				for (auto& a : v) {
					_elements.emplace_back(a);
				}
			}
			vector<T> _elements;

			bool makebox () {
				if (!isValid()) return false;
				// set the initial values such that all valid points are less than the mins
				// and greater than the maxes
				double maxlat = -100.0;
				double minlat = 100.0;
				double maxlon = -200.0;
				double minlon = 200.0;
				for (auto& a: _elements) {
					BoundingBox b = a->getbox();
					if (getlat(getne(b)) > maxlat) maxlat = getlat(getne(b));
					if (getlat(getsw(b)) < minlat) minlat = getlat(getsw(b));
					if (getlon(getne(b)) > maxlon) maxlon = getlon(getne(b));
					if (getlon(getsw(b)) < minlon) minlon = getlon(getsw(b));
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
	class MultiLineString : public Collection<LineString*> {
		public:
			/// @brief creates a blank MultiLineString
			MultiLineString () : Collection<LineString*>(GeometryType::MultiLineString) {};

			/// @brief creates a MultiLineString populated with the given LineStrings
			/// @param vector of LineString objects
			MultiLineString (vector<LineString*>& v) : 
				Collection<LineString*>(v, GeometryType::MultiLineString) {};

			/// @brief creates a MultiLineString populated with the given LineStrings
			/// @param vector of LineString objects
			MultiLineString (vector<LineString>& v) : 
				Collection<LineString*>(GeometryType::MultiLineString) {
					for (auto& a : v) _elements.emplace_back(new LineString(a));
				};

			/// @brief creates a MultiLineString populated with the given MultiPoint objects
			/// Each MultiPoint is converted to a LineString
			/// @param vector of MultiPoint objects
			MultiLineString (vector<MultiPoint>& v) : 
				Collection<LineString*>(GeometryType::MultiLineString) {
					for (auto& a : v) _elements.emplace_back(new LineString(a));
				};

			MultiLineString (Value& val) : Collection<LineString*>(val, GeometryType::MultiLineString) {};
			MultiLineString (Value* val) : Collection<LineString*>(val, GeometryType::MultiLineString) {};

			/// @brief is a const_iterator type for returning a target LineString
			typedef vector<LineString*>::const_iterator LinesConstItr;

			/// @brief populated this MultiLineString from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			bool load (Value* val) {
				if (!val) return false;
				_mytype = GeometryType::MultiLineString;
				Value *typeptr = nullptr;
				Value *coordptr = nullptr;
				Value *bboxptr = nullptr;
				if (val->IsObject()) {
					typeptr = Pointer("/type").Get(*val);
					coordptr= Pointer("/coordinates").Get(*val);
					bboxptr= Pointer("/bbox").Get(*val);
					if (!typeptr || !coordptr || (!typeptr->IsString() || 
						!(string(typeptr->GetString()) == "MultiLineString") ||
						!coordptr->IsArray())) return false;
					if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
				} else if (val->IsArray()) {
					coordptr = val;
				} else return false;
				_elements.clear();
				for (auto& v : coordptr->GetArray()) {
					if (v.IsArray()) _elements.emplace_back(new LineString(v));
				}
				return isValid();
			};

			/// @brief creates a GeoJSON MultiLineString from this MultiLineString
			/// @param d the GeoJSON document to create the object in 
			/// @returns the packed GeoJSON MultiLineString
			Value pack (Document& d) {
				Value b(kArrayType);
				Value v(kObjectType);
				v.AddMember("type", "MultiLineString", d.GetAllocator());
				v.AddMember("coordinates", packArray(d), d.GetAllocator());
				v.AddMember("bbox", packBBox(getbox(), d), d.GetAllocator());
				return v;
			};

			/// @brief find the lengths of all of the LineStrings in this MultiLineString
			/// @param type the distance measurement to use
			/// @returns a vector of doubles representing the lengths of each LineString in meters
			vector<double> lengths (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				vector<double> result;
				for (auto& a : _elements) result.emplace_back(a->length(type));
				return result;
			};
			
			/// @brief find the sum of the lengths of all of the LineStrings in this MultiLineString
			/// @param type the distance measurement to use
			/// @returns the total length of all the contained LineStrings in meters
			double totalLength (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				double result = 0;
				for (auto& a : _elements) result += a->length(type);
				return result;
			};
	};

	class MultiPolygon : public Collection<Polygon*> {
		public:
			MultiPolygon () : Collection<Polygon*>(GeometryType::MultiPolygon) {};
			MultiPolygon (Value& val) : Collection<Polygon*>(val, GeometryType::MultiPolygon) {};
			MultiPolygon (Value* val) : Collection<Polygon*>(val, GeometryType::MultiPolygon) {};
			MultiPolygon (vector<LineString>& v) : Collection<Polygon*>(GeometryType::MultiPolygon) {
				for (auto& a : v) {
					_elements.emplace_back(new Polygon(a));
				}
			};

			MultiPolygon (vector<MultiPoint>& v) : Collection<Polygon*>(GeometryType::MultiPolygon) {
				for (auto& a : v) {
					_elements.emplace_back(new Polygon(a));
				}
			};

			MultiPolygon (vector<Polygon*>& v) : Collection<Polygon*>(v, GeometryType::MultiPolygon) {};

			/// @brief is a const_iterator type for returning a target LineString
			typedef vector<Polygon*>::const_iterator PolygonsConstItr;

			/// @brief populated this MultiPolygon from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			bool load (Value* val) {
				if (!val) return false;
				_mytype = GeometryType::MultiPolygon;
				Value *typeptr = nullptr;
				Value *coordptr = nullptr;
				Value *bboxptr = nullptr;
				if (val->IsObject()) {
					typeptr = Pointer("/type").Get(*val);
					coordptr= Pointer("/coordinates").Get(*val);
					bboxptr= Pointer("/bbox").Get(*val);
					if (!typeptr || !coordptr || (!typeptr->IsString() || 
						!(string(typeptr->GetString()) == "MultiPolygon") ||
						!coordptr->IsArray())) return false;
					if (bboxptr && bboxptr->IsArray()) _bbox = loadBBox(*bboxptr);
				} else if (val->IsArray()) {
					coordptr = val;
				} else return false;
				_elements.clear();
				for (auto& v : coordptr->GetArray()) {
					if (v.IsArray()) _elements.emplace_back(new Polygon(v));
				}
				return isValid();
			};

			Value pack (Document& d) {
				Value b(kArrayType);
				Value v(kObjectType);
				v.AddMember("type", "MultiPolygon", d.GetAllocator());
				v.AddMember("coordinates", packArray(d), d.GetAllocator());
				v.AddMember("bbox", packBBox(getbox(), d), d.GetAllocator());
				return v;
			};

			bool includes (Point& p) {return includes(p.getPosition());};
			bool includes (Position& p) {
				if (bboxContains(_bbox, p)) {
					for (auto a : _elements) {
						if (a->includes(p)) return true;
					}
				}
				return false;
			};

			vector<double> perimeters (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				vector<double> result;
				for (auto& a : _elements) result.emplace_back(a->length(type));
				return result;
			};

			double totalPerimeters (const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				double result = 0;
				for (auto& a : _elements) result += a->length(type);
				return result;
			};
	};

	class GeometryCollection : public Collection<GeometryRoot*> {};

	/// @brief return the point on a line segment closest to a point
	/// Note that this function may behave strangely with great circle distances
	/// This function was derived from http://geomalgorithms.com/a02-_lines.html
	/// @param p is the point external to the line segment
	/// @param s0 is one end of the line segment
	/// @param s1 is the other end of the line segment
	/// @param type the distance type used in this calculation
	/// @returns the Position on the line segment closest to p
	inline Position point2segment(Point& p, Point& s0, Point& s1, const CourseTypeEnum type) {
		TwoVector v = s0.target(s1.getPosition(), type);
		TwoVector w = s0.target(p.getPosition(), type);
		double c1 = w * v;
		if (c1 <= 0) return s0.getPosition();
		double c2 = v * v;
		if (c2 <= c1) return s1.getPosition();
		double b = c1/c2;
		return s0.project((b * v), type);
	};

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
	inline unsigned int segmentsIntersect (Point& s10, Point& s11, Point& s20, Point& s21,
											Point& intersect0, Point& intersect1,
											const CourseTypeEnum type) {
		TwoVector u = s10.target(s11.getPosition(), type);
		TwoVector v = s20.target(s21.getPosition(), type);
		TwoVector w = s20.target(s11.getPosition(), type);
		double D = u.perp(v);

		if (fabs(D) < smallnum)	{								// S1 & S2 are parallel
			if ((u.perp(w) != 0) || v.perp(w) != 0) return 0;	// They are not collinear
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
					intersect0 = Point(s20.project(v * t0, type));
					return 1;
				} else {
					intersect0 = Point(s20.project(v * t0, type));
					intersect1 = Point(s20.project(v * t1, type));
					return 2;
				}
			}
		} else {								// the segments are skew and map intersect
			double sI = v.perp(w)/D;			// get intersect parameter for s1
			if ((sI < 0) || (sI > 1)) return 0;	// no intersection
			double tI = u.perp(w)/D;			// get the intersect parameters for s2
			if ((tI < 0) || (tI > 1)) return 0;	// no intersection
			intersect0 = Point(s11.project((u * (-sI)), type));
			return 1;
		}
	}

	/// @brief generates a pointer to a new geometry object from the given GeoJSON object
	/// @param val pointer to a GeoJSON object
	/// @returns a pointer to a new GeoJSON object to a nullptr if the passed object is not valid GeoJSON
	inline GeometryRoot* GeometryRoot::factory (Value* val) {
		if (!val) return nullptr;
		Value *valtype = Pointer("/type").Get(*val);
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
			// } else if (objType == "GeometryCollection") {
			// 	return new GeometryCollection(val);
			} else return nullptr;
		}
		return nullptr;
	}
};

#endif /* GEOGEOMETRY_H */