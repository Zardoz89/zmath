/**
* Module that defines a Quaternion
*/
module zmath.quaternion;

import zmath.vector;
import zmath.matrix;

import std.math;
import std.conv;

version(unittest) {
	import std.stdio;
}

unittest{
	writeln("Unit Test of Quaternion:");
	
	auto q1 = Qua_d(1,1,1,1);
	auto qk = Qua_d(0,0,1,0);
	auto qj = Qua_d(0,1,0,0);
	auto qi = Qua_d(1,0,0,0);
	
	assert (isQuaternion!(typeof(qj)));
	
	// Check equality
	assert( q1 == q1);
	assert( q1 != qj);
	assert( qj != q1);
	assert( q1.equal(q1) && ! q1.equal(qk) );
	writeln("Equality : OK");
	
	// Check addition / subtraction
	auto qki = qk + qi;
	auto q1j = q1 + qj;
	auto q1k = q1 - qk;
	assert (qki == Qua_d(1,0,1,0));
	assert (q1j == Qua_d(1,2,1,1));
	assert (q1k == Qua_d(1,1,0,1));
	writeln("Addition /subtraction : OK");
	
	// Check Hamilton multiplication
	auto q1i = q1 * qi;
	auto qi1 = qi * q1;
	assert(q1i == Qua_d(1,1,-1,-1)); // TODO Check this values
	assert (qi1 == Qua_d(1,-1,1,-1));
	writeln("Hamilton multiplication : OK");
	
	// Check Conjugate
	auto conj_q1 = q1.conj();
	assert(conj_q1 == Qua_d(-1,-1,-1,1));
	auto qr = q1 * q1.conj();
	assert(qr == Qua_d(0,0,0,4));
	writeln("Conjugate : OK");
	
	// Check conversion from Axis/angle 
	auto q_AxisX = Qua_d(Vec3r.X_AXIS,0);
	assert (q_AxisX == Qua_d(0,0,0,1)); // TODO Do more checks
	assert (q_AxisX.isUnit());
	writeln("Conversion from Axis & angle : OK");
	
	// Check conversion to Axis/angle
	double angle; 
	auto axis = q_AxisX.toAxisAngle(angle);
	assert(axis == Vec3d(1,0,0));
	assert(angle == 0);
	
	auto q_Axis = Qua_d(Vec3r(1,2,3), PI_4);
	axis = q_Axis.toAxisAngle(angle);
	auto cmp = Vec3d(1,2,3);
	cmp.normalize();
	assert(axis.equal(cmp) );
	assert(approxEqual(angle, PI_4));
	writeln("Conversion to Axis/angle : OK");
	
	// Check conversion from/to Euler angles
	auto q_bank = Qua_d(0,0, PI_2);
	assert (q_bank.equal( Qua_d(0.707107, 0, 0, 0.707107) ) );
	assert (q_bank.isUnit() );
	assert (q_bank.toEuler().equal(Vec3d(0,0, PI_2)));
	
	auto q_head = Qua_d(PI_4,0, 0);
	assert (q_head.isUnit() );
	assert (q_head.toEuler().equal(Vec3d(PI_4,0, 0)));
	
	auto q_elev = Qua_d(0, -PI_4, 0);
	assert (q_elev.isUnit() );
	assert (q_elev.toEuler().equal(Vec3d(0, -PI_4, 0)));
	
	auto q_r = Qua_d(1, -PI_4, PI_4);
	assert (q_r.isUnit() );
	assert (q_r.toEuler().equal(Vec3d(1, -PI_4, PI_4)));
	writeln("Conversion from/to Euler angles : OK");
	
	// Check rotation with quaternion
	Vec3d v1 = Vec3d(1,1,1);
	Vec3d v2 = q_bank.rotate(v1);
	assert (v2.equal( Vec3d(1, -1, 1) ));
	Vec3d v3 = q_head.rotate(v1);
	assert (v3.equal( Vec3d(1.41421, 1, 0) ));
	Vec3d v4 = q_elev.rotate(v1);
	assert (v4.equal( Vec3d(1.41421, 0, 1) ));
	
	writeln("Rotating a Vector with a Quaternion : OK");
	
	auto m = q_bank.toMatrix!(Mat4d)();
	// TODO Check this matrix and others
	
	writeln("To rotaton Matrix : OK");
	
	// check conversion between Quaternions
	auto qf = toImpl!(Qua_f, Qua_d) (q1);
	auto qreal = toImpl!(Qua_r, Qua_d) (q1);
	auto qd = toImpl!(Qua_d, Qua_f) (qf);
	
	assert(is(typeof(qf) == Qua_f));
	assert(is(typeof(qreal) == Qua_r));
	assert(is(typeof(qd) == Qua_d));
	for (auto i=0; i< 4; i++ ) {
		assert( qf[i] == 1);
		assert( qreal[i] == 1);
		assert( qd[i] == 1);
	}	
	writeln("Conversion between Quaternions : OK");
	
	writeln();
}


alias Quaternion!float	Qua_f;	/// Alias of a Quaternion with floats
alias Quaternion!double	Qua_d;	/// Alias of a Quaternion with doubles
alias Quaternion!real		Qua_r;	/// Alias of a Quaternion with reals

/**
* Quaternion over a FloatPoint type,
*/
public struct Quaternion(T)
if (is(T == real) || is(T == double) || is(T == float) ) 
{
	static assert (is(T : real));
	
	static enum size_t dim = 4; /// Quaternion Dimension
	
	union {
		private T[4] coor;				/// Quaternion coords in a rray
		struct {
			T i;						/// i complex component
			T j;						/// j complex component
			T k;						/// k complex component
			T w;						/// w real component
		}
		
		struct {
			T x;						/// 1 complex component
			T y;						/// j complex component
			T z;						/// k complex component
		}
	}
	
	/**
	* Build a new Quaternion from a set of initial values
	* If no there values for z and w, will be set to 0
	* Params:
	*	i = i imaginary component
	*	j = j imaginary component
	*	k = k imaginary component
	*	w = i real component
	*/
	this(in T i, in T j, in T k = 0, in T w = 0) {
		this.i = i;
		this.j = j;
		this.k = k;
		this.w = w;
	}
	
	/**
	* Build a new Quaternion from a array
	* If no there values for j, k and w, will be set to 0
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
	
	/**
	* Build a new Quaternion from a 3D 
	* w will be set to 0
	* Params:
	*	v = Vector 3d 
	*/
	this(in Vector!(T, 3) v) {
		size_t i;
		for (; i< dim && i< v.dim ; i++) {
			coor[i] = v[i];
		}
		for (;i< dim; i++) {
			coor[i] = 0;
		}
	}
	
	/**
	* Build a new Quaternion from a 4D 
	* Params:
	*	v = Vector 4d
	*/
	this(in Vector!(T, 4) v) {
		size_t i;
		for (; i< dim && i< v.dim ; i++) {
			coor[i] = v[i];
		}
	}
	/**
	* Creates a Quaternion from a 3d Vector and angle
	*	Params:
	*	v = Rotation axis
	* angle = Rotation angle in radians
	*	Returns:
	*	Quaternion that represents a rotation
	*/
	this( Vec3r v, in real angle) {
		v.normalize();
		auto sin_2 = sin(angle / 2.0L);
		i = v.x * sin_2;
		j = v.y * sin_2;
		k = v.z * sin_2;
		w = cos(angle /2.0L);
	}
	
	/**
	* Creates a new Quaternion from  Euler angles
	* Params :
	* heading = Rotation around Z axis in radians
	*	elevation = Rotation around Y axis in radians
	*	bank = Rotation around X axis in radians
	*	Returns:
	*	Quaternion that represents a rotation
	*/
	this( in real heading, in real elevation, in real bank)	{
		auto c1 = cos(heading/2);     // x =s1s2*c3  + c1c2 *s3
	  auto s1 = sin(heading/2);			// y =s1*c2*c3 + c1*s2*s3
	  auto c2 = cos(elevation/2);		// z =c1*s2*c3 - s1*c2*s3
	  auto s2 = sin(elevation/2);		// w =c1c2*c3 - s1s2*s3
	  auto c3 = cos(bank/2);
	  auto s3 = sin(bank/2);
	  auto c1c2 = c1*c2;
	  auto s1s2 = s1*s2;
	  i = s1s2*c3  + c1c2 *s3;
	  j = s1*c2*c3 + c1*s2*s3;
	  k = c1*s2*c3 - s1*c2*s3;
	  w = c1c2*c3 - s1s2*s3;
	}
	
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
	* Returns the actual length of this quaternion
	*/
	@property T length() const {
		return sqrt(sq_length);
	}
	
	/**
	* Returns the actual squared length of this quaternion
	*/
	@property T sq_length() const {
		return x*x + y*y + z*z + w*w;
	}
	
	// Operations ***************************************	
	
	/**
	* Define Equality 
	*/
	bool opEquals(ref const Quaternion rhs) const {
		if (x != rhs.x) return false;
		if (y != rhs.y) return false;
		if (z != rhs.z) return false;
		if (w != rhs.w) return false;
		return true;	
	}
	
	/**
	* Approximated equality
	*/
	bool equal(ref const Quaternion rhs) const {
			return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) && approxEqual(z, rhs.z) && approxEqual(w, rhs.w);
	}

		
	/**
	* Approximated equality with controlable precision
	*/
	bool equal(ref const Quaternion rhs, T maxDiff) const {
			return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff) && approxEqual(z, rhs.z, maxDiff) && approxEqual(w, rhs.w, maxDiff);
	}
	
	/**
	* Approximated equality with controlable precision
	*/
	const bool equal(ref const Quaternion rhs, T maxRelDiff, T maxAbsDiff = 1e-05) const {
			return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff) && approxEqual(w, rhs.w, maxRelDiff, maxAbsDiff);
	}
	
	/**
	* Define unary operators + and -
	*/
	Quaternion opUnary(string op) () const
		if (op == "+" || op == "-")
	{
			return Quaternion( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"), mixin(op~"w"));
	}
	
	
	/**
	* Define binary operator + and -
	*/
	Quaternion opBinary(string op) (ref const Quaternion rhs) const
		if (op == "+" || op == "-")
	{
			return Quaternion( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"), mixin("z"~op~"rhs.z"), mixin("w"~op~"rhs.w"));
	}
	
	/**
	* Define Scalar multiplication
	*/
	Quaternion opBinary(string op) (in T rhs) const
		if (op == "*" )
	{
			return Quaternion(x *rhs, y *rhs, z *rhs, w *rhs);
	}
	
	/**
	* Defines Hamilton product for Quaternions 
	* Params:
	* rhs = Quaternion at rigth of operation
	*	Return:
	*	Quaternion result of hamilton product of <strong>this</strong> and rhs (this * rhs)
	*/
	Quaternion opBinary(string op) (in Quaternion rhs) const
		if (op == "*" )
	{
			auto _w = w* rhs.w - x* rhs.x - y* rhs.y - z* rhs.z; // (a1a2 − b1b2 − c1c2 − d1d2)w
			auto _i = w* rhs.x + x* rhs.w + y* rhs.z - z* rhs.y; // (a1b2 + b1a2 + c1d2 − d1c2)i
			auto _j = w* rhs.y - x* rhs.z + y* rhs.w + z* rhs.x; // (a1c2 − b1d2 + c1a2 + d1b2)j
			auto _k = w* rhs.z + x* rhs.y - y* rhs.x + z* rhs.w; // (a1d2 + b1c2 − c1b2 + d1a2)k
	
	return Quaternion(_i, _j ,_k, _w);
	}
	
	/**
	* Obtain conjugate of a Quaternion
	* Params:
	*	q = Quaternion from were calc
	* Return:
	* Conjugate Quaternion of q
	*/
	Quaternion conj () const {
		return Quaternion(- x, -y, -z, w);
	}
	
	/**
	* It's a unitary Quaternion (length == 1)
	*/
	@property bool isUnit() const {
		return approxEqual (abs( this.sq_length - 1.0), 0 );
	}
	
	/**
	* Normalize this Quaternion
	*/
	void normalize() {
		if (!isUnit()) {
			real l = 1 / this.length;
			x = x *l;
			y = y *l;
			z = z *l;
			w = w *l;
		}
	}
	
	/**
	* Apply a rotation Quaternion over a 3d Vector
	*/
	Vector!(T, 3) rotate (in Vector!(T, 3) v) const
	in {
		assert (this.isUnit());
	}
	body {
		Quaternion vecQua = Quaternion (v);
		
		Quaternion resQua = vecQua * this.conj();
		resQua = this * resQua;
		
		return Vector!(T, 3)( resQua.x, resQua.y, resQua.z ); 
	}
	
	// Misc ***************************************	
	/**
	* Calc the Axis and angle of rotation of this Quaternion
	* Params:
	*	q = Quaternion that represents a rotation
	*	angle = Set to angle of rotation
	*	Returns:
	*	Axis of rotation vector
	*/
	Vector!(T, 3) toAxisAngle (out T angle)	 const
	in {
		assert (this.isUnit());
	}
	body {
		angle = 2* acos(w);
		if (approxEqual(angle, 0))
			return Vector!(T, 3)(1, 0, 0);
		auto inv_qw = 1/ sqrt(1 - w*w);
		return Vector!(T, 3)(i*inv_qw, j*inv_qw, k*inv_qw );
	}

	/**
	* Returns a Vector with x = Heading, y= Elevation and z = Bank
	*	Params:
	*	q = Quaternion that represents a rotation
	* Returns:
	*	Vector with Euler angles of rotation in radians
	*/
	Vector!(T, 3) toEuler() const	
	in {
		assert (this.isUnit());
	}
	body {
		T head = void, elev = void , bank = void;
		auto att = 2* x* y + 2* z* w;
		if (att >= 1) { // North Pole
			head = 2* atan2(x, w);
			elev = PI_2;
			bank = 0;
			
		} else if (att <= -1)  { // South Pole
			head = -2* atan2(x, w);
			elev = -PI_2;
			bank = 0;
		} else {
			head = atan2(2* y* w - 2* x* z , 1 - 2* y* y - 2* z* z);
			elev = asin(att);
			bank = atan2(2* x* w - 2* y* z , 1 - 2* x* x - 2* z* z);
		}
		
		return Vector!(T, 3)(head,elev,bank);
	}

	/**
	* Returns a Matrix that represents the rotation stored in a quaternion
	* Returns :
	* Rotation matrix in a 3 or 4 dimension matrix
	*/
	U toMatrix(U) () const
	if (is (U == Mat4r) || is (U == Mat3r) || is (U == Mat4d) || is (U == Mat3d) || is (U == Mat4f) || is (U == Mat3f) )
	in {
		assert (this.isUnit());
	}
	body {
		static assert (U.dim == 3 || U.dim == 4 );
		
		
		auto x2 = x*x;	auto y2 = y*y;	auto z2 = z*z;
		auto xy = x*y;	auto xz = x*z;	auto yz = y*z;
		auto wx = w*x;	auto wy = w*y;	auto wz = w*z;
		
		static if (U.dim == 4) {
			return U([ 1.0L - 2.0L * (y2 +z2), 2.0L * (xy - wz), 2.0L * (xz + wy), 0.0L,
				 				2.0L * (xy + wz), 1.0L - 2.0L * (x2 + z2), 2.0L * (yz - wx), 0.0L,
				 				2.0L * (xz - wy), 2.0L * (yz + wx), 1.0L - 2.0L * (x2 + y2), 0.0L,
				 				0.0L						, 0.0L						, 0.0L									 , 1.0L]);
		} else {
			return U([ 1.0L - 2.0L * (y2 +z2), 2.0L * (xy - wz)	, 2.0L * (xz + wy),
				 				 2.0L * (xy + wz), 1.0L - 2.0L * (x2 + z2), 2.0L * (yz - wx),
				 				 2.0L * (xz - wy), 2.0L * (yz + wx), 1.0L - 2.0L * (x2 + y2)]);
		}
	}

	/**
	* Returns a string representation of this Quaternion
	*/
	string toString() {
		string ret = "["~ to!string(x) ~", "~ to!string(y) ~", "~ to!string(z) ~", "~ to!string(w) ~"]";
		return ret;		
	}
	
}

/**
* Say if a thing it's a Quaternion
*/
template isQuaternion(T)
{
	//immutable bool isQuaternion = is(T == Vector);
  immutable bool isQuaternion = __traits(compiles,
        (){  
            T t;
            static assert(T.dim == 4);
            auto coor = t.coor;
            auto x = t.x;
            auto y = t.y;
            auto z = t.z;
            auto w = t.w;
            T conj = t.conj();
            // TODO : Should test for methods ?        
        }
    );
}

/**
* to converto a Quaternion to other Quaternion
*/
T toImpl(T, S)(S s) 
if (!implicitlyConverts!(S, T) && isQuaternion!T && isQuaternion!S )
{
    static assert (isQuaternion!T, "This type not is a Quaternion"); // Ok, redundant
		
		T newQ;
		static if (is (typeof(newQ.x) == typeof(s.x)))  {
			foreach (i ,se ;s.coor)
				newQ.coor[i] =  se;
		} else {
			foreach (i ,se ;s.coor)
				newQ.coor[i] =  to!(typeof(newQ.x))(se); 
		}
		
		return newQ;
}
