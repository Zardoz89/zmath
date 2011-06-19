/**
* Module that defines a N-Dimensional Vector where N = 2, 3 or 4
*/
module zmath.vector;

import zmath.aux;

import std.math;
import std.conv;

version(unittest) {
	import std.stdio;
}

unittest {
	writeln("Unit test of Vector :");
	Vec2r v1 = Vec2r(0,2);
	auto v2 = Vec2r(2,0);
	auto v3 = Vec2r(1,1);
	
	assert (v2.x == 2);
	assert (v2.coor[0] == 2);
	assert (v2[0] == 2);
	
	// Check equality
	assert (v2 != v3);
	assert (v3 == v3);
	assert (v3.equal(v3));
	writeln("Equality Operator : OK");
	
	// Check change of sign
	auto nv3 = -v3;
	assert (nv3.x == -1);
	assert (nv3.y == -1);
	writeln("Change of sign : OK");
	
	// Check addition
	auto v12 = v1 + v2;
	auto v21 = v2 + v1;

	// Symetry
	assert(v12 == v21);
	// Value
	assert(v12.x == 2.0L);
	assert(v12.y == 2.0L);
	writeln("Addition : OK");
	
	// Check subtraction
	auto v1m3 = v1 - v3;
	assert(v1m3.x == -1.0L);
	assert(v1m3.y == 1.0L);
	// Subtraction asimetry
	auto v3m1 = v3 - v1;
	assert(v3m1.x == 1.0L);
	assert(v3m1.y == -1.0L);
	writeln("Subtraction : OK");
	
	// Product by Scalar
	auto v2ten = v2 *10.0L;
	assert(v2ten.x == (10.0L * v2.x));
	assert(v2ten.y == (10.0L * v2.y));
	v2ten = v2ten * (1/2.0L);
	assert(v2ten.x == (5.0L * v2.x));
	assert(v2ten.y == (5.0L * v2.y));
	writeln("Product by scalar : OK");
	
	
	// Dot Product and length
	auto v2dotv1 = v2 * v1;
	assert(v2dotv1 == 0.0L);
	assert(v2ten.length == 10.0L);
	assert(v2ten.sq_length == 100.0L);
	writeln("Dot Product and length : OK");
	
	// Test Distance
	auto dis = v1.distance(v2);
	assert (dis == 2* SQRT2 );
	writeln("Distance : OK");
	
	// Unitary Vector
	auto uv2ten = v2ten; 
	uv2ten.normalize();
	auto uv2 =  v2;
	uv2.normalize();
	assert (uv2 == uv2ten);
	assert (uv2.x == 1.0L);
	assert (uv2.y == 0.0L);
	assert (v2ten.x != 1.0L);
	assert (v2ten.y == 0.0L);
	assert (!v2ten.isUnit());
	assert (uv2ten.isUnit());
	writeln("Unitary Vector : OK");
	
	// Cross product
	auto r3v1 = Vec3r(0,0,10);
	auto r3v2 = Vec3r(2,0,2);
	auto cross0 = r3v1 & r3v1;
	auto cross12 = r3v1 & r3v2;
	auto cross21 = r3v2 & r3v1;
	
	assert (cross0.length == 0);
	assert (cross12.length == 20.0);
	assert (cross21.length == 20.0);
	
	assert (cross12 == Vec3r(0,20.0L, 0));
	assert (cross21 == Vec3r(0,-20.0L, 0));
	writeln("Cross Product : OK");
	
	// Check Rotation in R2
	auto rot90 = v2.rotate(PI_2);
	auto rotn90 = v2.rotate(-PI_2);
	assert (rot90.length == v2.length);
	assert (rotn90.length == v2.length);
	assert (rot90.equal(Vec2r(0,2)));
	assert (rotn90.equal(Vec2r(0,-2)));
	writeln("Rotation in R2 : OK");
	
	// Check proyection
	auto proy_v3_in_v2 = v3.projectOnTo(v2);
	assert (proy_v3_in_v2 == Vec2r.X_AXIS );
	writeln("Proyection of a Vector : OK");
	
	// Check isOk() and isFinite()
	assert(r3v1.isOk());
	assert(r3v1.isFinite());
	Vec4r nain;
	assert(! nain.isOk());
	assert(! nain.isFinite());
	writeln("isOk() and isFinite() : OK");
	
	writeln();
}

alias Vector!(real,2) Vec2r; /// Alias of a 2d Vector with reals
alias Vector!(real,3) Vec3r; /// Alias of a 3d Vector with reals
alias Vector!(real,4) Vec4r; /// Alias of a 4d Vector with reals

alias Vector!(double,2) Vec2d; /// Alias of a 2d Vector with doubles
alias Vector!(double,3) Vec3d; /// Alias of a 3d Vector with doubles
alias Vector!(double,4) Vec4d; /// Alias of a 4d Vector with doubles

alias Vector!(float,2) Vec2f; /// Alias of a 2d Vector with floats
alias Vector!(float,3) Vec3f; /// Alias of a 3d Vector with floats
alias Vector!(float,4) Vec4f; /// Alias of a 4d Vector with floats

/**
* N-Dimensional Vector over a FloatPoint type, where N must be 2,3 or 4
*/
public struct Vector(T, size_t dim_)
if (is(T == real) || is(T == double) || is(T == float) ) 
{
	
	enum size_t dim = dim_; /// Vector Dimension
	
	static assert (dim >= 2 && dim <= 4);
	static assert (is(T : real));
	
	union {
		T[dim] coor; /// Vector coords like Array
		
		struct {
			static if( dim >= 1) T x;
			static if( dim >= 2) T y;
			static if( dim >= 3) T z;
			static if( dim >= 4) T w;
		}
	}
	
	// Consts
	static if (dim == 2) { // R2
		public static enum Vector!(T,2) ZERO = Vector!(T,2)(0, 0);	 /// Origin
		public static enum Vector!(T,2) X_AXIS = Vector!(T,2)(1, 0); /// X Axis in R2
		public static enum Vector!(T,2) Y_AXIS = Vector!(T,2)(0, 1); /// Y Axis in R2
	}
	
	static if (dim == 3) { // R3
		public static enum Vector!(T,3) ZERO = Vector!(T,3)(0, 0, 0);	 	/// Origin
		public static enum Vector!(T,3) X_AXIS = Vector!(T,3)(1, 0, 0); /// X Axis in R3	
		public static enum Vector!(T,3) Y_AXIS = Vector!(T,3)(0, 1, 0); /// Y Axis in R3
		public static enum Vector!(T,3) Z_AXIS = Vector!(T,3)(0, 0, 1); /// Z Axis in R3
	}
	
	static if (dim == 4) { // R4
		public static enum Vector!(T,4) ZERO = Vector!(T,4)(0, 0, 0, 0);	 /// Origin
		public static enum Vector!(T,4) X_AXIS = Vector!(T,4)(1, 0, 0, 0); /// X Axis in R4	
		public static enum Vector!(T,4) Y_AXIS = Vector!(T,4)(0, 1, 0, 0); /// Y Axis in R4
		public static enum Vector!(T,4) Z_AXIS = Vector!(T,4)(0, 0, 1, 0); /// Z Axis in R4
		public static enum Vector!(T,4) W_AXIS = Vector!(T,4)(0, 0, 0, 1); /// W Axis in R4
	}
	
	
	/**
	* Build a new Vector from a set of initial values
	* If no there values for z and w, will be set to 0
	* Params:
	*	x = X coord
	*	y = Y coord
	*	z = Z coord
	*	w = W coord
	*/
	this(in T x, in T y, in T z = 0, in T w = 0) {
		this.x = x;
		this.y = y;
		static if (dim >= 3)
			this.z = z;
		static if (dim >= 4) 
			this.w = w;
	}
	
	/**
	* Build a new Vector from a array
	* If no there values for y, z and w, will be set to 0
	* Params:
	*	xs = Array with coords
	*/
	this(in T[] xs) {
		size_t i;
		for (; i< dim && i< xs.length ; i++) {
			coor[i] = xs[i];
		}
		for (;i< dim; i++) {
			coor[i] = 0;
		}
	}
	
	/+
	/**
	* Post-Blit
	*/
	this(this) {
		coor = coor.dup;
	}+/
	
	// Basic Properties ***************************************	
	
	/**
	* Returns i coord of this vector
	*/
	T opIndex(size_t i) const { return coor[i];}
	
	/**
	* Assigns a value to a coord
	*/
	void opIndexAssign(T c, size_t i) {
		coor[i] = c;
	}
	
	/**
	* Returns the actual length of this Vector
	*/
	@property T length() const {
		return sqrt(sq_length);
	}
	
	/**
	* Returns the actual squared length of this Vector
	*/
	@property T sq_length() const {
		return this * this;
	}
	
	// Operations ***************************************	
	
	/**
	* Define Equality 
	*/
	bool opEquals(ref const Vector rhs) const {
		if (x != rhs.x) return false;
		static if (dim >= 2) 
			if (y != rhs.y) return false;
		static if (dim >= 3) 
			if (z != rhs.z) return false;
		static if (dim == 4) 
			if (w != rhs.w) return false;
		return true;	
	}
	
	
	/**
	* Approximated equality
	*/
	bool equal(ref const Vector rhs) const {
		static if (dim == 2) 
			return approxEqual(x, rhs.x) && approxEqual(y, rhs.y);
		static if (dim == 3) 
			return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) && approxEqual(z, rhs.z);
		static if (dim == 4) 
			return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) && approxEqual(z, rhs.z) && approxEqual(w, rhs.w);
	}

		
	/**
	* Approximated equality with controlable precision
	*/
	bool equal(ref const Vector rhs, T maxDiff) const {
		static if (dim == 2) 
			return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff);
		static if (dim == 3) 
			return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff) && approxEqual(z, rhs.z, maxDiff);
		static if (dim == 4) 
			return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff) && approxEqual(z, rhs.z, maxDiff) && approxEqual(w, rhs.w, maxDiff);
	}
	
	/**
	* Approximated equality with controlable precision
	*/
	const bool equal(ref const Vector rhs, T maxRelDiff, T maxAbsDiff = 1e-05) const {
		static if (dim == 2) 
			return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff);
		static if (dim == 3) 
			return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff);
		static if (dim == 4) 
			return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff) && approxEqual(w, rhs.w, maxRelDiff, maxAbsDiff);
	}
	
	/**
	* Define unary operators + and -
	*/
	Vector opUnary(string op) () const
		if (op == "+" || op == "-")
	{
		static if (dim == 2) 
			return Vector( mixin(op~"x"), mixin(op~"y"));
		static if (dim == 3) 
			return Vector( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"));	
		static if (dim == 4) 
			return Vector( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"), mixin(op~"w"));
	}
	
	
	/**
	* Define binary operator + and -
	*/
	Vector opBinary(string op) (ref const Vector rhs) const
		if (op == "+" || op == "-")
	{
		static if (dim == 2) 
			return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"));
		static if (dim == 3) 
			return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"), mixin("z"~op~"rhs.z"));	
		static if (dim == 4) 
			return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"), mixin("z"~op~"rhs.z"), mixin("w"~op~"rhs.w"));
	}
	
	/**
	* Define Scalar multiplication
	*/
	Vector opBinary(string op) (in T rhs) const
		if (op == "*" )
	{
		static if (dim == 2) 
			return Vector(x *rhs, y *rhs);
		static if (dim == 3) 
			return Vector(x *rhs, y *rhs, z *rhs);	
		static if (dim == 4) 
			return Vector(x *rhs, y *rhs, z *rhs, w *rhs);
	}
	
	
	/**
	* Define Dot Product
	*/
	T opBinary(string op) (ref const Vector rhs) const
		if (op == "*" )
	{
		T tmp = 0;
		foreach (i, x; coor) {
			tmp += x * rhs.coor[i];
		}
		return tmp;
	}
	
	
	/**
	* Define Cross Product for R3 (operation c = a & b )
	*/
	Vector opBinary(string op) (ref const Vector rhs) const 
		if (op == "&" && dim == 3)
	{
		return Vector( coor[1]*rhs.coor[2] - coor[2]*rhs.coor[1],
			 						 coor[2]*rhs.coor[0] - coor[0]*rhs.coor[2],
			 						 coor[0]*rhs.coor[1] - coor[1]*rhs.coor[0]);
	}
	
	/**
	* It's a unitary vector (length == 1)
	*/
	@property bool isUnit() const {
		return approxEqual (abs( this.sq_length - 1.0), 0 );
	}
	
	/**
	* Normalize this vector 
	*/
	void normalize() {
		if (!isUnit()) {
			real l = 1 / this.length;
			static if (dim >= 2) {
				x = x *l; y = y *l;
			}	
			static if (dim >= 3) 
				z = z *l;
			static if (dim == 4) 
				w = w *l;
		}
	}
	
	/**
	* Returns the unit vector of this vector
	*/
	@property Vector unit() {
		Vector ret = this;
		ret.normalize();
		return ret;
	}
	
	/**
	* Obtains the proyection two vectors
	* Params:
	*	a = Vector to proyect
	* b = Vector where proyect vector a
	*/
	static Vector projectOnTo (Vector a,Vector b) {
		b.normalize();
		Vector ret = void; 
		auto s = a*b;
		ret = b * (a*b);
		return ret ;
	}
	
	/**
	* Obtains the proyection of this vector over other vector
	* Params:
	* b = Vector where proyect vector a
	*/
	Vector projectOnTo(Vector b) {
		return this.projectOnTo(this,b);
	}
	
	/**
	* Calculate the distance between two points pointed by this vector and other
	* Params :
	*	b = Vector B
	*/
	T distance( in Vector b) {
		auto d = this - b;
		return d.length;
	}
	
	/**
	* Calculate the squared distance between two points pointed by this vector and other
	* Params :
	*	b = Vector B
	*/
	T sq_distance( in Vector b) {
		auto d = this - b;
		return d.sq_length;
	}
	
	/**
	* Rotation in R2
	* Params:
	* 	axis = Rotation axis
	*/
	static if (dim == 2)
	Vector rotate (real angle) const {
		return Vector( x * cos(angle) - y * sin(angle), x * sin(angle) + y *cos(angle) );
	}
	
	// Misc ***************************************************************************
	
	/**
	* Checks that the vector not have a weird NaN value
	*/
	bool isOk() {
		if (isNaN(x) ==0 || isNaN(y) ==0) return true;
		static if (dim >= 3)
			if (isNaN(z) ==0) return true;
		static if (dim >= 4)
			if (isNaN(w) ==0) return true;	
		return false;	
	}
	
	/**
	* Checks that the vector have finite values
	*/
	bool isFinite() {
		if (std.math.isFinite(x) && std.math.isFinite(y)) return true;
		static if (dim >= 3)
			if (std.math.isFinite(z)) return true;
		static if (dim >= 4)
			if (std.math.isFinite(w)) return true;	
		return false;	
	}
	
	/**
	* Returns a string representation of this vector
	*/
	string toString() {
		string ret = "["~ to!string(x) ~", "~ to!string(y);
		static if (dim >=3)
			ret ~= ", "~ to!string(z);
		static if (dim >=4)
			ret ~= ", "~ to!string(w);
		ret ~= "]";
		return ret;		
	}
	
}

