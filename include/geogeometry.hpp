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
	const double LatLonTolerance = 0.00001;
	const double MeterTolerance = 2.0;

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

	class Position {
		public:
			Position () {};
			Position (const double& lon, const double& lat, const double& ele = 0) :
				_lat(lat), _lon(lon), _ele(ele) {};
			Position (Value& v) {load(v);};
			Position (Value* v) {load(v);};

			/// @brief populates this Position object from a GeoJSON array
			inline bool load (Value& v) {return load(&v);};
			inline bool load (Value* v) {
				if (v == nullptr) return false;
				Value *vp = nullptr;
				if (v->IsArray()) {
					vp = Pointer("/0").Get(*v);
					if (vp && vp->IsDouble()) _lon = vp->GetDouble();
					vp = Pointer("/1").Get(*v);
					if (vp && vp->IsDouble()) _lat = vp->GetDouble();
					vp = Pointer("/2").Get(*v);
					if (vp && vp->IsDouble()) _ele = vp->GetDouble();
				} 
				return this->isValid();
			};

			inline Value pack (Document& d) {
				Value v(kArrayType);
				v.PushBack(Value().SetDouble(this->lon()), d.GetAllocator());
				v.PushBack(Value().SetDouble(this->lat()), d.GetAllocator());
				v.PushBack(Value().SetDouble(this->ele()), d.GetAllocator());
				return v;
			};

			inline const double& lat () const {return _lat;};
			inline const double& lon () const {return _lon;};
			inline const double& ele () const {return _ele;};
			inline double& lat () {return _lat;};
			inline double& lon () {return _lon;};
			inline double& ele () {return _ele;};
			inline bool lat (const double& l) {
				if (!isfinite(l)) return false;
				if (fabs(l) > 90.0) return false;
				_lat = l;
				return isValid();
			};
			inline bool lon (const double& l) {
				if (!isfinite(l)) return false;
				if (fabs(l) > 180.0) return false;
				_lon = l;
				return isValid();
			};
			inline bool ele (const double& e) {
				if (!isfinite(e)) return false;
				_ele = e;
				return isValid();
			}

			inline void norm () {
				while (_lon > 180.0) _lon -= 360.0;
				while (_lon < -180.0) _lon += 360.0;
			}

			inline bool isValid () const {
				return (isfinite(_lat) && isfinite(_lon) && isfinite(_ele) &&
					(fabs(_lon) < 180.0) && (fabs(_lat) < 90.0));
			};

			/// @brief checks whether the Position p is north of this position
			inline bool isNorth (const Position& p) const {
				if (p.lat() > this->lat()) return true;
				return false;
			};

			/// @brief checks whether the Position p is to the east of this position and west of meridian l0
			inline bool isEast (const Position& p, const double& l0) const {
				// shift all the longitudes to 0-360 degrees
				double testlon = (p.lon() + 180.0);
				double thislon = (this->lon() + 180.0);
				double limitlon = (l0 + 180.0);
				// handle the cases where either testlon or limitlon are east over the antimeridian from thislon
				if (limitlon < thislon) {				// check whether limit is east over the antimeridian
					if (testlon < limitlon) {			// testlon is east of the antimeridian and west of limitlon
						return true;
					} 
				}
				if (thislon < testlon) return true;
				return false;

			};

			/// @brief compare whether two positions are equal
			inline bool operator== (const Position& other) const {
				return ((other.lat() == this->lat()) && 
					(other.lon() == this->lon()));
			}

			/// @brief compare whether two positions are unequal
			inline bool operator!= (const Position& other) const {
				return ((other.lat() != this->lat()) || 
					(other.lon() != this->lon()));
			}


		private:
			double _lat = NAN;
			double _lon = NAN;
			double _ele = NAN;
	};

	class BoundingBox {
		public:
			BoundingBox () {};
			BoundingBox (const Position& ne, const Position& sw) :
				_ne(ne), _sw(sw) {};
			BoundingBox (const double& ne_lon, const double& ne_lat, 
				const double& sw_lon, const double& sw_lat) :
				_ne(ne_lon, ne_lat), _sw(sw_lon, sw_lat) {};
			BoundingBox (Value& v) {load(v);};
			BoundingBox (Value* v) {load(v);};

			inline bool load (Value& v) {return load(&v);};
			inline bool load (Value* v) {
				if (v == nullptr) return false;
				Value *vp = nullptr;
				if (v->IsArray()) {
					vp = Pointer("/0").Get(*v);
					if (vp && vp->IsDouble()) _sw.lon(vp->GetDouble());
					vp = Pointer("/1").Get(*v);
					if (vp && vp->IsDouble()) _sw.lat(vp->GetDouble());
					vp = Pointer("/2").Get(*v);
					if (vp && vp->IsDouble()) _ne.lon(vp->GetDouble());
					vp = Pointer("/3").Get(*v);
					if (vp && vp->IsDouble()) _ne.lat(vp->GetDouble());
				}
				return isValid();
			}
	
			/// @brief creates a GeoJSON bounding box array from the given BoundingBox in the given Document
			inline Value pack (Document& d) {
				Value v(kArrayType);
				Value i(kNumberType);
				v.PushBack(i.SetDouble(_sw.lon()), d.GetAllocator());
				v.PushBack(i.SetDouble(_sw.lat()), d.GetAllocator());
				v.PushBack(i.SetDouble(_ne.lon()), d.GetAllocator());
				v.PushBack(i.SetDouble(_ne.lat()), d.GetAllocator());
				return v;
			};

			/// @brief checks whether a BoundingBox refers to a real bounding box
			inline bool isValid () {
				return(_sw.isValid() && _ne.isValid() &&
						(_ne.lat() > _sw.lat()));
			};

			/// @brief checks whether a given location is within this bounding box
			inline bool contains (Position& p) {
				if (!(this->isValid())) return false;
				if (_ne.lon() > _sw.lon()) {	// this checks that the bounding box does not cross the antimeridian
					return ((p.lon() <= _ne.lon()) && (p.lon() >= _sw.lon()) &&
							(p.lat() <= _ne.lat()) && (p.lat() >= _sw.lat()));
				} else { // this handles the case where the bounding box crosses the antimeridian
					return (((p.lon() <= _ne.lon()) || (p.lon() >= _sw.lon())) &&
							(p.lat() <= _ne.lat()) && (p.lat() >= _sw.lat()));
				}
			}

			inline Position& ne () {return _ne;};
			inline Position& sw () {return _sw;};
			inline const Position& ne () const {return _ne;};
			inline const Position& sw () const {return _sw;};

		private:
			Position _ne;
			Position _sw;
	};
 
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
			virtual bool load (Value& val) {return load(&val);};
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
			const BoundingBox& bbox () {
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
			Point () : GeometryRoot(GeometryType::Point) {};

			/// @brief creates a Point at the given position
			/// @param posn the desired position
			Point (const Position& posn) : _posn(posn), 
				GeometryRoot(GeometryType::Point) {};

			/// @brief creates a Point at the given position
			/// @param lon the desired longitude (decimal degrees)
			/// @param lat the desired latitude (decimal degrees)
			/// @param ele the desired elevation (meters above datum)
			Point (const double lon, const double lat, const double ele = 0.0)  : 
				GeometryRoot(GeometryType::Point),
				_posn(lon, lat, ele) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			Point (Value& val) {load(&val);};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			Point (Value* val) {load(val);};

			/// @brief populated this point from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			inline bool load (Value& val) {return load(&val);};
			inline bool load (Value* val) {
				if (!val) return false;
				_mytype = GeometryType::Point;
				if (val->IsObject()) {
					Value *typeptr = Pointer("/type").Get(*val);
					Value *coordptr= Pointer("/coordinates").Get(*val);
					if (typeptr && coordptr) {
						if (typeptr->IsString() && 
							(string(typeptr->GetString()) == "Point") &&
							coordptr->IsArray()) {
							_posn.load(*coordptr);
						} else return false;
					} else return false;
				} else if (val->IsArray() && (val->Size() > 1)) {
					_posn.load(*val);
				} else return false;
				return isValid();
			};

			/// @brief creates a GeoJSON Point from this Point
			/// @param d the GeoJSON document to create the object in 
			/// @returns the packed GeoJSON Point
			Value pack (Document& d) {
				Value v(kObjectType);
				v.AddMember("type", "Point", d.GetAllocator());
				v.AddMember("coordinates", _posn.pack(d), d.GetAllocator());
				return v;
			};

			/// @brief creates a GeoJSON Point coordinate array from this Point
			/// @param d the GeoJSON document to create the array in 
			/// @returns the packed GeoJSON Point coordinates array
			Value packArray (Document& d) {return _posn.pack(d);};

			/// @brief returns a reference to the embedded position
			Position& position () {return _posn;};
			const Position& position () const {return _posn;};

			/// @brief get a bearing from this Point to a distant position
			/// @param dest the Position of the target
			/// @param type the type of course to plot to the target. Note that the 
			/// a GreatCircle type will calculate the initial bearing. 
			/// @returns the bearing to target in degrees from north. Returns NAN if
			/// no valid bearing can be calculated
			double bearing (const Position& dest, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
				if ((!this->isValid()) || (!dest.isValid())) return NAN;
				double azi1, azi2, dist;
				switch (type) {
					case CourseTypeEnum::GreatCircle: {
						Geodesic::WGS84().Inverse(lat(), lon(), dest.lat(), dest.lon(), azi1, azi2);
						return azi1;
						break;
					}
					case CourseTypeEnum::RhumbLine:	{
						Rhumb::WGS84().Inverse(lat(), lon(), dest.lat(), dest.lon(), dist, azi1);
						return azi1;
						break;
					}
					case CourseTypeEnum::Approximate: {
						double deltalat = metersPerDegreeLat(lat())*(lat() - dest.lat());
						double deltalon = metersPerDegreeLon(lat())*(lon() - dest.lon());
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
				if ((!this->isValid()) || (!dest.isValid())) return NAN;
				double azi1, azi2, dist;
				switch (type) {
					case CourseTypeEnum::GreatCircle: {
						Geodesic::WGS84().Inverse(lat(), lon(), dest.lat(), dest.lon(), dist, azi1, azi2);
						return dist;
						break;
					}
					case CourseTypeEnum::RhumbLine:	{
						Rhumb::WGS84().Inverse(lat(), lon(), dest.lat(), dest.lon(), dist, azi1);
						return dist;
						break;
					}
					case CourseTypeEnum::Approximate: {
						double deltalat = metersPerDegreeLat(lat())*(lat() - dest.lat());
						double deltalon = metersPerDegreeLon(lat())*(lon() - dest.lon());
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
				if ((!this->isValid()) || (!dest.isValid())) return TwoVector {NAN,NAN};
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
				if ((!this->isValid()) || (!projection.isValid())) return Position(NAN,NAN,0.0);
				Position p = _posn;
				switch (type) {
					case CourseTypeEnum::GreatCircle: {
						Geodesic::WGS84().Direct(lat(), lon(), projection.angleDeg(), projection.mag(), p.lat(), p.lon());
						break;
					}
					case CourseTypeEnum::RhumbLine: {
						Rhumb::WGS84().Direct(lat(), lon(), projection.angleDeg(), projection.mag(), p.lat(), p.lon());
						break;
					}
					case CourseTypeEnum::Approximate: {
						p.lon() += projection.y()/metersPerDegreeLon(lat());
						p.lat() += projection.x()/metersPerDegreeLat(lat());
						break;
					}
					default:
						break;
				}
				// normalize longitude
				if (p.lon() > 180.0) p.lon() -= 360;
				if (p.lon() < -180.0) p.lon() += 360;
				return p;
			};

			bool isValid () const {
				if (_mytype == GeometryType::Point) return (_posn.isValid());
				return false;
			};
			
			const double& lon () const {return _posn.lon();};
			const double& lat () const {return _posn.lat();};
			const double& ele () const {return _posn.ele();};

			bool inline operator== (const Point& a) const {
				return (a.position() == this->position());
			}

			bool inline operator!= (const Point& a) const {
				return (a.position() != this->position());
			}

		private:
			Position _posn;
			Point (GeometryType t) : GeometryRoot(t) {};
			bool makebox () {
				_bbox.ne() = _posn;
				_bbox.sw() = _posn;
				return true;
			};
	};

	/// @brief comparison function for sorting Positions by latitude
	inline bool latlt (const Position& a, const Position& b) {
		return (a.lat() < b.lat());
	}

	/// @brief comparison function for sorting Positions by longitude
	inline bool lonlt (const Position& a, const Position& b) {
		return (a.lon() < b.lon());
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
					_points.emplace_back(a.position());
				}
			};

			// @brief creates a Point object from a GeoJSON object
			// @param val reference to a GeoJSON object
			MultiPoint (Value& val) : GeometryRoot(GeometryType::MultiPoint) {};

			// @brief creates a Point object from a GeoJSON object
			// @param val pointer to a GeoJSON object
			MultiPoint (Value* val) : GeometryRoot(GeometryType::MultiPoint) {};

			/// @brief is a const_iterator type for returning a target position
			typedef vector<Position>::const_iterator PositionConstItr;

			/// @brief populated this MultiPoint from a GeoJSON object
			/// @param val pointer to a GeoJSON object
			/// @returns true if object is successfully populated
			virtual bool load (Value& val) {return load(&val);};
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
					if (bboxptr && bboxptr->IsArray()) _bbox.load(*bboxptr);
				} else if (val->IsArray()) {
					coordptr = val;
				} else return false;
				_points.clear();
				for (auto& v : coordptr->GetArray()) {
					if (v.IsArray()) _points.emplace_back(Position(v));
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
				bbox();
				v.AddMember("bbox", _bbox.pack(d), d.GetAllocator());
				return v;
			};

			/// @brief creates a GeoJSON MultiPoint coordinate array from this Point
			/// @param d the GeoJSON document to create the array in 
			/// @returns the packed GeoJSON MultiPoint coordinates array
			Value packArray (Document& d) {
				Value b(kArrayType);
				for (auto& p : _points) {
					b.PushBack(p.pack(d), d.GetAllocator());
				}
				return b;
			};

			/// @brief return a reference to the members of the MultiPoint
			vector<Position>& positions () {return _points;};

			/// @brief returns a reference to a given member
			/// @params t the index to return
			/// @returns the position at index t
			Position& position (const unsigned int& t) {return _points[t];};

			/// @brief returns a Point created from a given member
			/// @params t the index to return
			/// @returns a Point created from the position at index t
			Point point (const unsigned int& t) {return Point(_points[t]);};

			/// @brief returns a vector of Points corresponding to the members of the MultiPoint
			vector<Point> points () {
				vector<Point> v;
				for (auto& p: _points) v.emplace_back(Point(p));
				return v;
			};

			/// @brief insert a new Point after the given index
			/// @param p the Point to insert
			/// @param posn the index to insert after
			/// @returns true if successful
			bool insert (Point& p, const unsigned int& posn) {return insert(p.position(), posn);};

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
				return this->insert(i, p.position());
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
					valid &= p.isValid();
				}
				return valid;
			}

		protected:
			vector<Position> _points;
			MultiPoint (GeometryType t) : GeometryRoot(t) {};
			MultiPoint (const vector<Position>& v, GeometryType t) : _points(v), GeometryRoot(t) {};
			MultiPoint (vector<Point>& v, GeometryType t) : GeometryRoot(t) {
				for (auto& a : v) {
					_points.emplace_back(a.position());
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
					if (a.lat() > maxlat) maxlat = a.lat();
					if (a.lat() < minlat) minlat = a.lat();
					if (a.lon() > maxlon) maxlon = a.lon();
					if (a.lon() < minlon) minlon = a.lon();
				}
				_bbox = BoundingBox(maxlon, maxlat, minlon, minlat);
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
					_points.push_back(a.position());
				}
			};

			/// @brief creates a LineString object from the given MultiPoint
			/// @param v MultiPoint to create LineString from 
			LineString (MultiPoint& v) : 
				MultiPoint(v.positions(), GeometryType::LineString) {};

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
			// tuple<Position, PositionConstItr> closestPoint (const Point& target, const CourseTypeEnum type = CourseTypeEnum::GreatCircle) const {
			// 	double bestdist = 1000000000;
			// 	PositionConstItr p = _points.cbegin();
			// 	Position best, newposn;
			// 	for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
			// 		if (a != _points.end()) {
			// 			newposn = point2segment(target, Point(*a), Point(*(a + 1)), type);
			// 			double newdist = target.distance(newposn, type);
			// 			if (newdist < bestdist) {
			// 				bestdist = newdist;
			// 				best = newposn;
			// 				p = a;
			// 			}
			// 		}
			// 	}
			// 	return make_tuple(best, p);
			// };

			/// @brief find the line segments in the LineString that intersect with
			/// the given line segment.
			/// @param s0 start point of the line segment to intersect
			/// @param s1 end point of the line segment to intersect
			// tuple<MultiPoint, vector<PositionConstItr>> intersection (const Point& s0, const Point& s1, const CourseTypeEnum type = CourseTypeEnum::RhumbLine) {
			// 	Point o0, o1;
			// 	vector<PositionConstItr> outitr;
			// 	MultiPoint out;
			// 	for (auto a = _points.cbegin(); a != _points.cend(); ++a) {
			// 		if (a != _points.end()) {
			// 			if (segmentsIntersect(Point(*a), Point(*(a+1)), s0, s1, o0, o1, type)) {
			// 				outitr.emplace_back(a);
			// 				out.insert(out.cend(), *a);
			// 			}
			// 		} 
			// 	}
			// 	return make_tuple(out, outitr);
			// }
		
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
			LineString (MultiPoint v, GeometryType t) : MultiPoint(v.positions(), t) {};

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
			Polygon (LineString& v) : LineString(v.positions(), GeometryType::Polygon) {};

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
				if ((_points.cbegin()->isValid())) {
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

			bool includes (Point& p) {return includes(p.position());};
			bool includes (Position& p) {
				if (_bbox.contains(p)) {

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
				return (first.isValid() && last.isValid() && (first == last));
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
					BoundingBox b = a->bbox();
					if (b.ne().lat() > maxlat) maxlat = b.ne().lat();
					if (b.sw().lat() < minlat) minlat = b.sw().lat();
					if (b.ne().lon() > maxlon) maxlon = b.ne().lon();
					if (b.sw().lon() < minlon) minlon = b.sw().lon();
				}
				_bbox = BoundingBox(maxlon, maxlat, minlon, minlat);
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
					if (bboxptr && bboxptr->IsArray()) _bbox.load(*bboxptr);
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
				v.AddMember("bbox", _bbox.pack(d), d.GetAllocator());
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
					if (bboxptr && bboxptr->IsArray()) _bbox.load(*bboxptr);
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
				v.AddMember("bbox", _bbox.pack(d), d.GetAllocator());
				return v;
			};

			bool includes (Point& p) {return includes(p.position());};
			bool includes (Position& p) {
				if (_bbox.contains(p)) {
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
	// inline Position point2segment(Point& p, Point& s0, Point& s1, const CourseTypeEnum type) {
	// 	TwoVector v = s0.target(s1.position(), type);
	// 	TwoVector w = s0.target(p.position(), type);
	// 	double c1 = w * v;
	// 	if (c1 <= 0) return s0.position();
	// 	double c2 = v * v;
	// 	if (c2 <= c1) return s1.position();
	// 	double b = c1/c2;
	// 	return s0.project((b * v), type);
	// };

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