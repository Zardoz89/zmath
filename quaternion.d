/**
A Quaternion it's used to store a space orientation / rotation, and makes easy 
to interpolate rotations, at same time that avoid grimbal locks that have using
Euler angles

License: $(LINK2 http://www.gnu.org/licenses/lgpl.txt, LGPL 3).

Authors: Luis Panadero Guardeño $(LINK http://zardoz.es)
*/
module zmath.quaternion;

import zmath.vector;
import zmath.matrix;

import std.math;
import std.conv;

alias Quaternion!float  Qua_f;  /// Alias of a Quaternion with floats
alias Quaternion!double Qua_d;  /// Alias of a Quaternion with doubles
alias Quaternion!real   Qua_r;  /// Alias of a Quaternion with reals

/**
* Quaternion over a FloatPoint type,
*/
public struct Quaternion(T=float)
if (__traits(isFloating, T))  {
  static assert (__traits(isFloating, T));
  
  static enum size_t dim = 4; /// Quaternion Dimension
  
  union {
    private T[4] coor;  /// Quaternion coords in a rray
    struct {
      T i;  /// i complex component
      T j;  /// j complex component
      T k;  /// k complex component
      T w;  /// w real component
    }
    
    struct {
      T x;  /// alias of i component
      T y;  /// alias of j component
      T z;  /// alias of k component
    }
  }
  
  /**
  * Build a new Quaternion from a set of initial values
  * If no there values will be set to 0,0,0,1
  * Params:
  *	i = i imaginary component
  *	j = j imaginary component
  *	k = k imaginary component
  *	w = i real component
  */
  this(in T i=0, in T j=0, in T k = 0, in T w = 1) {
    this.i = i;
    this.j = j;
    this.k = k;
    this.w = w;
  }
  
  /**
  * Build a new Quaternion from a array
  * If no there values for j, and k, will be set to 0. w is set to 1
  * Params:
  *	xs = Array with coords
  */
  this(in T[] xs) {
    size_t l = xs.length > dim? dim : xs.length;
    coor[0..l] = xs[0..l].dup;
    if (l <dim) {
      coor[l..$-1] = 0; // imaginary part
      w = 1; // real part
    }
  }
  
  /**
  * Build a new Quaternion from a 3D Vector
  * w will be set to 1
  * Params:
  *	v = Vector 3d 
  */
  this(in Vector!(T, 3) v) {
    coor[0..3] = v.coor.dup;
    w = 1;
  }
  
  /**
  * Build a new Quaternion from a 4D Vector
  * Params:
  *	v = Vector 4d
  */
  this(in Vector!(T, 4) v) {
    coor = v.coor.dup;
  }
  
  unittest {
    auto q = Qua_r(1,2,3,4);
    assert(q.coor == [1,2,3,4]);
    real[] arr = [4.0,3.0,2.0,1.0, 20.0, 30.0];
    q = Qua_r(arr);
    arr[0] = 7.7;
    assert(q.coor == [4,3,2,1]);
    q = Qua_r(Vec4r(1,2,3,4));
    assert(q.coor == [1,2,3,4]);
  }
  
  /**
  * Creates a Quaternion from a 3d Vector and angle
  *	Params:
  *	v = Rotation axis
  * angle = Rotation angle in radians
  *	Returns:
  *	Quaternion that represents a rotation
  */
  this( Vector!(T,3) v, in T angle) {
    v.normalize();
    auto sin_2 = sin(angle / 2.0L);
    i = v.x * sin_2;
    j = v.y * sin_2;
    k = v.z * sin_2;
    w = cos(angle /2.0L);
  }
  
  unittest {
    auto q_AxisX = Qua_d(Vec3d.X_AXIS, 0);
    assert (q_AxisX == Qua_d(0,0,0,1)); // TODO Do more checks
    assert (q_AxisX.isUnit());
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
  this( in T heading, in T elevation, in T bank)	{
    auto c1 = cos(heading/2);     // x =s1s2*c3  + c1c2 *s3
    auto s1 = sin(heading/2);     // y =s1*c2*c3 + c1*s2*s3
    auto c2 = cos(elevation/2);   // z =c1*s2*c3 - s1*c2*s3
    auto s2 = sin(elevation/2);   // w =c1c2*c3 - s1s2*s3
    auto c3 = cos(bank/2);
    auto s3 = sin(bank/2);
    auto c1c2 = c1*c2;
    auto s1s2 = s1*s2;
    i = s1s2*c3  + c1c2 *s3;
    j = s1*c2*c3 + c1*s2*s3;
    k = c1*s2*c3 - s1*c2*s3;
    w = c1c2*c3 - s1s2*s3;
  }
  
  unittest {
    auto q_bank = Qua_d(0,0, PI_2);
    assert (q_bank.equal( Qua_d(0.707107, 0, 0, 0.707107) ) );
    assert (q_bank.isUnit() );
    auto q_head = Qua_d(PI_4,0, 0);
    assert (q_head.isUnit() );
    assert (q_head.toEuler().equal(Vec3d(PI_4,0, 0)));
    auto q_elev = Qua_d(0, -PI_4, 0);
    assert (q_elev.isUnit() );
    assert (q_elev.toEuler().equal(Vec3d(0, -PI_4, 0)));
    auto q_r = Qua_d(1, -PI_4, PI_4);
    assert (q_r.isUnit() );
    assert (q_r.toEuler().equal(Vec3d(1, -PI_4, PI_4)));
  }
  
  // Basic Properties **********************************************************
  
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
  
  unittest {
    auto q = Qua_f(10,0,0,0);
    assert(q.sq_length == 100);
    q = Qua_f(1,1,1,1);
    assert(q.length == 2);
  }
  
  // Operations ****************************************************************
  
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
    return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) 
        && approxEqual(z, rhs.z) && approxEqual(w, rhs.w);
  }
  
  /**
  * Approximated equality with controlable precision
  */
  bool equal(ref const Quaternion rhs, T maxDiff) const {
    return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff)
        && approxEqual(z, rhs.z, maxDiff) && approxEqual(w, rhs.w, maxDiff);
  }
  
  /**
  * Approximated equality with controlable precision
  */
  bool equal(ref const Quaternion rhs, T maxRelDiff,
                    T maxAbsDiff = 1e-05) const {
    return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) 
        && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) 
        && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff) 
        && approxEqual(w, rhs.w, maxRelDiff, maxAbsDiff);
  }
  
  unittest {
    auto q1 = Qua_d(1,2,3,4);
    auto qj = Qua_d(0,1,0,0);
    assert( q1 == q1);
    assert( q1 != qj);
    assert( qj != q1);
    assert( q1.equal(q1) && ! q1.equal(qj) );
  }
  
  /**
  * Define unary operators + and -
  */
  Quaternion opUnary(string op) () const
  if (op == "+" || op == "-") {
    return Quaternion( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"),
              mixin(op~"w"));
  }
  
  /**
  * Define binary operator + and -
  */
  Quaternion opBinary(string op) (ref const Quaternion rhs) const
  if (op == "+" || op == "-") {
    return Quaternion( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"),
              mixin("z"~op~"rhs.z"), mixin("w"~op~"rhs.w"));
  }
  
  unittest {
    // Check addition / subtraction
    auto qk = Qua_d(0,0,1,0);
    auto qj = Qua_d(0,1,0,0);
    auto qi = Qua_d(1,0,0,0);
    auto q1 = Qua_d(1,1,1,1);
    auto qki = qk + qi;
    auto q1j = q1 + qj;
    auto q1k = q1 - qk;
    assert (qki == Qua_d(1,0,1,0));
    assert (q1j == Qua_d(1,2,1,1));
    assert (q1k == Qua_d(1,1,0,1));
  }
  
  /**
  * Define Scalar multiplication
  */
  Quaternion opBinary(string op) (in T rhs) const
  if (op == "*" ) {
    return Quaternion(x *rhs, y *rhs, z *rhs, w *rhs);
  }
  
  /**
  * Defines Hamilton product for Quaternions 
  * Params:
  * rhs = Quaternion at rigth of operation
  *	Return:
  *	Quaternion result of hamilton product of <strong>this</strong> and rhs (this
* rhs)
  */
  Quaternion opBinary(string op) (in Quaternion rhs) const
    if (op == "*" ) {
    // (a1a2 − b1b2 − c1c2 − d1d2)w
    auto _w = w* rhs.w - x* rhs.x - y* rhs.y - z* rhs.z;
    // (a1b2 + b1a2 + c1d2 − d1c2)i
    auto _i = w* rhs.x + x* rhs.w + y* rhs.z - z* rhs.y;
    // (a1c2 − b1d2 + c1a2 + d1b2)j
    auto _j = w* rhs.y - x* rhs.z + y* rhs.w + z* rhs.x;
    // (a1d2 + b1c2 − c1b2 + d1a2)k
    auto _k = w* rhs.z + x* rhs.y - y* rhs.x + z* rhs.w;  
    return Quaternion(_i, _j ,_k, _w);
  }
  
  unittest {
    // Check Hamilton multiplication
    auto qi = Qua_d(1,0,0,0);
    auto q1 = Qua_d(1,1,1,1);
    auto q1i = q1 * qi;
    auto qi1 = qi * q1;
    assert(q1i == Qua_d(1,1,-1,-1)); // TODO Check this values
    assert (qi1 == Qua_d(1,-1,1,-1));
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
  
  unittest {
    auto q1 = Qua_d(1,1,1,1);
    auto conj_q1 = q1.conj();
    assert(conj_q1 == Qua_d(-1,-1,-1,1));
    auto qr = q1 * q1.conj();
    assert(qr == Qua_d(0,0,0,4));
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
      T l = 1 / this.length;
      x = x *l;
      y = y *l;
      z = z *l;
      w = w *l;
    }
  }
  
  unittest {
    auto q = Qua_r(10,10,1,5);
    assert(! q.isUnit);
    q.normalize();
    assert(q.isUnit);
  }
  
  /**
  * Apply a 3d rotation Quaternion over a 3d/4d Vector
  */
  K rotate (K) (K v) const 
  if (isVector!K && K.dim >= 3){    
    Quaternion vecQua = Quaternion (v);    
    Quaternion resQua = vecQua * this.conj();    
    resQua = this * resQua;
    return K( resQua.coor );
  }
  
  unittest {
    // Check rotation with quaternion
    auto q_bank = Qua_d(0,0, PI_2);
    auto q_head = Qua_d(PI_4,0, 0);
    auto q_elev = Qua_d(0, -PI_4, 0);
    auto v1 = Vec3d(1,1,1);
    auto v4d = cast(Vec4d) v1;
    auto v2 = q_bank.rotate(v1);    
    assert (v2.equal( Vec3d(1, -1, 1) ));
    auto v3 = q_head.rotate(v1);
    assert (v3.equal( Vec3d(1.41421, 1, 0) ));
    auto v4 = q_elev.rotate(v1);
    assert (v4.equal( Vec3d(1.41421, 0, 1) ));
    auto v5 = q_elev.rotate(v4d);
    assert (v5.equal( Vec4d(1.41421, 0, 1) ));
  }
  
  // Misc **********************************************************************
  
  /**
  * Calc the Axis and angle of rotation of this Quaternion
  * Params:
  *	q = Quaternion that represents a rotation
  *	angle = Set to angle of rotation
  *	Returns:
  *	Axis of rotation vector
  */
  Vector!(T, 3) toAxisAngle (out T angle) const
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
  
  unittest {
    // Check conversion to Axis/angle
    auto q_head = Qua_d(PI_4,0, 0);
    double angle; 
    auto axis = q_head.toAxisAngle(angle);
    assert(axis.equal(Vec3d(0,1,0)) );
    assert(approxEqual(angle, PI_4));
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
    
    return Vector!(T, 3)(head,elev,bank); // TODO: Return a tuple with 3
                                          // values ?
  }
  
  unittest {
    auto q = Qua_d(PI_4,PI_4, PI_4);
    auto rot = q.toEuler();
    assert(rot.equal(Vec3d(PI_4,PI_4, PI_4)) );
  }
    
  /**
  * Casting for convert between Quaternions
  */
  Tout opCast( Tout ) () 
  if (isQuaternion!(Tout)) {    
    Tout newQ;
    static if (is (typeof(newQ.x) == typeof(this.x)))  {
      foreach (i ,s ;coor)
        newQ.coor[i] =  s;
    } else {
      foreach (i ,s ;coor)
        newQ.coor[i] =  to!(typeof(newQ.x))(s); 
    }    
    return newQ;
  }
  
  unittest {
    auto q = Qua_d(1,1,1,1);
    auto qf = cast(Qua_f) q;
    auto qreal = cast(Qua_r) q;
    auto qd = cast(Qua_d) qf;
    
    assert(is(typeof(qf) == Qua_f));
    assert(is(typeof(qreal) == Qua_r));
    assert(is(typeof(qd) == Qua_d));
    for (auto i=0; i< 4; i++ ) {
      assert( qf[i] == 1);
      assert( qreal[i] == 1);
      assert( qd[i] == 1);
    }	
  }
  
  /**
  * Casting for make rotation Matrix
  */
  Tout opCast( Tout ) () 
  if (isMatrix!(Tout) && Tout.dim >= 3)
  in {
    assert (this.isUnit());
  } body {    
    auto x2 = x*x;	auto y2 = y*y;	auto z2 = z*z;
    auto xy = x*y;	auto xz = x*z;	auto yz = y*z;
    auto wx = w*x;	auto wy = w*y;	auto wz = w*z;    
    static if (Tout.dim == 4) {
      return Tout([ 1.0L - 2.0L * (y2 +z2), 2.0L * (xy - wz), 2.0L * (xz + wy),
        0.0L, 2.0L * (xy + wz), 1.0L - 2.0L * (x2 + z2), 2.0L * (yz - wx), 
        0.0L, 2.0L * (xz - wy), 2.0L * (yz + wx), 1.0L - 2.0L * (x2 + y2),
        0.0L, 0.0L, 0.0L, 0.0L, 1.0L]);
    } else {
      return Tout([ 1.0L - 2.0L * (y2 +z2), 2.0L * (xy - wz), 2.0L * (xz + wy),
        2.0L * (xy + wz), 1.0L - 2.0L * (x2 + z2), 2.0L * (yz - wx),
        2.0L * (xz - wy), 2.0L * (yz + wx), 1.0L - 2.0L * (x2 + y2)]);
    }
  }
  
  unittest {
    // Generation of Rotation matrix
    auto q_bank = Qua_d(0,0, PI_2);
    auto q_head = Qua_d(PI_4,0, 0);
    auto q_elev = Qua_d(0, -PI_4, 0);
    auto m = cast(Mat4d) q_bank;
    assert (approxEqual(m.determinant , 1) );
    assert (m.equal( Mat4d([1, 0, 0, 0,
                            0, 0,-1, 0,
                            0, 1, 0, 0,
                            0, 0, 0, 1]) ));
    
    m = cast(Mat4d) q_head;
    assert (approxEqual(m.determinant , 1) );
    assert (m.equal( Mat4d([0.707107, 0, 0.707107, 0,
                            0,        1, 0,        0,
                           -0.707107, 0, 0.707107, 0,
                            0,        0, 0,        1]) ));
    
    m = cast(Mat4d) q_elev;
    assert (approxEqual(m.determinant , 1) );
    assert (m.equal( Mat4d([0.707107, 0.707107, 0, 0,
                           -0.707107, 0.707107, 0, 0,
                            0,        0,        1, 0,
                            0,        0,        0, 1]) ));
  }
  
  /**
  * Returns a string representation of this Quaternion
  */
  string toString() {
    string ret = "["~ to!string(x) ~", "~ to!string(y) ~", "~ to!string(z)
                  ~", " ~ to!string(w) ~"]";
    return ret;
  }
  
}

/**
* Say if a thing it's a Quaternion
*/
template isQuaternion(T)  {
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

unittest {
  auto q = Qua_f(1,2,3,4);
  assert (isQuaternion!(typeof(q)));
  assert (! isQuaternion!(int));
}
