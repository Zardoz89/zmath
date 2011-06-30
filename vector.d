/**
Defines a Vector of float point typefrom 2 to 4 dimension

License: $(LINK2 http://www.gnu.org/licenses/lgpl.txt, LGPL 3).

Authors: Luis Panadero Guardeño $(LINK http://zardoz.es)
*/
module zmath.vector;

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
	
	assert (isVector!(typeof(v1))); // Yep it's a Vector
	
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
	
	// Check casting
	auto v2f = cast(Vec2f) v2;
	auto v2d = cast(Vec2d) v2;
	auto vecf = Vec2f(1,1);
	//auto vec4d = toImpl!(Vec4d, Vec2f) (vecf);
	//auto vec4d = to!Vec4d (vecf);
	auto vec4d = cast(Vec4d) vecf;
	assert(is(typeof(v2f) == Vec2f));
	assert(is(typeof(v2d) == Vec2d));
	assert(is(typeof(vec4d) == Vec4d));
	
	writeln("Vector Casting : OK");
	
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
if (__traits(isFloating, T) ) {
	static enum size_t dim = dim_; 		/// Vector Dimension
	
	static assert (dim >= 2 && dim <= 4, "Not valid dimension size.");
	static assert (__traits(isFloating, T), "Type not is a Float Point type.");
	static assert (is(T : real), "Type not is like a Float Point type.");
	
	union {
		private T[dim] coor; /// Vector coords like Array
		
		struct {
			static if( dim >= 1) T x;
			static if( dim >= 2) T y;
			static if( dim >= 3) T z;
			static if( dim >= 4) T w;
		}
	}
	
	// Consts
	static if (dim == 2) { // for R2
		public static enum Vector!(T,2) ZERO = Vector!(T,2)(0, 0);	 /// Origin
		public static enum Vector!(T,2) X_AXIS = Vector!(T,2)(1, 0); /// X Axis in R2
		public static enum Vector!(T,2) Y_AXIS = Vector!(T,2)(0, 1); /// Y Axis in R2
	}
	
	static if (dim == 3) { // for R3
		public static enum Vector!(T,3) ZERO = Vector!(T,3)(0, 0, 0);	 	/// Origin
		public static enum Vector!(T,3) X_AXIS = Vector!(T,3)(1, 0, 0); /// X Axis in R3	
		public static enum Vector!(T,3) Y_AXIS = Vector!(T,3)(0, 1, 0); /// Y Axis in R3
		public static enum Vector!(T,3) Z_AXIS = Vector!(T,3)(0, 0, 1); /// Z Axis in R3
	}
	
	static if (dim == 4) { // for R4
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
	* Assigns a value to a i coord
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
	* Params:
	* rhs = Vector at rigth of '=='
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
	* Params:
	* rhs = Vector to compare with this vector
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
	* Params:
	* rhs = Vector to compare with this vector
	* maxDiff = 
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
	* Params:
	* rhs = Vector to compare with this vector
	* maxRelDiff = Maximun relative difference
	* maxAbsDiff = Maximun absolute difference
	*/
	bool equal(ref const Vector rhs, T maxRelDiff, T maxAbsDiff = 1e-05) const {
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
	* Returns : True if length approxEqual to 1.0
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
	* Returns : A Vector that it's proyection of Vector a over Vector b
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
	* b = Vector where proyect this vector
	* Returns : A Vector that it's proyection of this Vector over Vector b
	*/
	Vector projectOnTo(Vector b) {
		return this.projectOnTo(this,b);
	}
	
	/**
	* Calculate the distance between two points pointed by this vector and other
	* Params :
	*	b = Vector B
	* Returns : Distance between the point pointed by this vector and other point
	*/
	T distance( in Vector b) {
		auto d = this - b;
		return d.length;
	}
	
	/**
	* Calculate the squared distance between two points pointed by this vector and other
	* Params :
	*	b = Vector B
	* Returns : Squared distance between the point pointed by this vector and other point
	*/
	T sq_distance( in Vector b) {
		auto d = this - b;
		return d.sq_length;
	}
	
	/**
	* Rotation in R2
	* Params:
	* 	axis = Rotation axis
	* Returns : A vector that is the rotation of this vector
	*/
	static if (dim == 2)
	Vector rotate (real angle) const {
		return Vector( x * cos(angle) - y * sin(angle), x * sin(angle) + y *cos(angle) );
	}
	
	// Misc ***************************************************************************
	
	/**
	* Checks that the vector not have a weird NaN value
	* Returns : True if this vector not have a NaN value
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
	* Returns : True if this vector have finite values (not infinite value or NaNs)
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
	* Casting method to convert to other vector types
	*/
	Tout opCast( Tout ) () 
	if (isVector!(Tout) && Tout.dim >= dim) {
		static assert (isVector!(Tout), "This type not is a Vector");
		static assert (Tout.dim >= dim, "Original Vector bigger that destiny Vector");
		
		Tout newVector; auto i = 0;
		static if (is (typeof(newVector.x) == typeof(this.x)))  {
			for (; i < dim; i++)
				newVector.coor[i] =  coor[i];
		} else {
			for (; i < dim; i++)	
				newVector.coor[i] =  to!(typeof(newVector.x))(coor[i]); 
		}
			
		for (; i< Tout.dim; i++) // Expands a small vector to a bigger dimension with 0 value in extra dimension
			newVector.coor[i] = 0;
		
		return newVector;
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

/**
* Say if a thing it's a Vector
*/
template isVector(T)
{
	//immutable bool isVector = is(T == Vector);
  immutable bool isVector = __traits(compiles,
        (){  
            T t;
            static assert(T.dim >= 2 && T.dim <= 4);
            auto coor = t.coor;
            auto x = t.x;
            auto y = t.y;
            static if(t.dim >= 3)
                auto z = t.z;
            static if(t.dim >= 4)
                auto w = t.w;
            // TODO : Should test for methods ?        
        }
    );
}

/+
/**
* to converto a vector to other vector
*/
T toImpl(T, S)(S s) 
if (!implicitlyConverts!(S, T) && isVector!T && isVector!S )
{
    static assert (isVector!T, "This type not is a Vector"); // Ok, redundant
		static assert (T.dim >= S.dim, "Original Vector bigger that destiny Vector");
		
		T newVector; auto i = 0;
		static if (is (typeof(newVector.x) == typeof(s.x)))  {
			for (; i < S.dim; i++)
				newVector.coor[i] =  s.coor[i];
		} else {
			for (; i < S.dim; i++)	
				newVector.coor[i] =  to!(typeof(newVector.x))(s.coor[i]); 
		}
			
		for (; i< T.dim; i++) // Expands a small vector to a bigger dimension with 0 value in extra dimension
			newVector.coor[i] = 0;
		
		return newVector;
}+/