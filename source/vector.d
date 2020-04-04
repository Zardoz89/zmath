/**
Defines a Vector of float point typefrom 2 to 4 dimension
*/
module zmath.vector;

import std.math;
import std.conv;

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
  static enum size_t dim = dim_;    /// Vector Dimension  
  
  static assert (dim >= 2 && dim <= 4, "Not valid dimension size.");
  static assert (__traits(isFloating, T), "Type not is a Float Point type.");
  static assert (is(T : real), "Type not is like a Float Point type.");
  
  union {
    package T[dim] coor; /// Vector coords like Array
    
    struct {
      static if( dim >= 1) T x;
      static if( dim >= 2) T y;
      static if( dim >= 3) T z;
      static if( dim >= 4) T w;
    }
  }
  
  // Consts
  static if (dim == 2) { // for R2
    public static enum Vector!(T,2) ZERO = Vector!(T,2)(0, 0);  ///Origin
    public static enum Vector!(T,2) X_AXIS = Vector!(T,2)(1, 0); ///X Axis in R2
    public static enum Vector!(T,2) Y_AXIS = Vector!(T,2)(0, 1); ///Y Axis in R2
  }
  
  static if (dim == 3) { // for R3
    public static enum Vector!(T,3) ZERO = Vector!(T,3)(0, 0, 0); ///Origin
    /**X Axis in R3 */
    public static enum Vector!(T,3) X_AXIS = Vector!(T,3)(1, 0, 0);
    /**Y Axis in R3 */
    public static enum Vector!(T,3) Y_AXIS = Vector!(T,3)(0, 1, 0);
    /**Z Axis in R3 */
    public static enum Vector!(T,3) Z_AXIS = Vector!(T,3)(0, 0, 1);
  }
  
  static if (dim == 4) { // for R4
    public static enum Vector!(T,4) ZERO = Vector!(T,4)(0, 0, 0, 0); /// Origin
    /**X Axis in R4 */
    public static enum Vector!(T,4) X_AXIS = Vector!(T,4)(1, 0, 0, 0);
    /**Y Axis in R4 */
    public static enum Vector!(T,4) Y_AXIS = Vector!(T,4)(0, 1, 0, 0);
    /**Z Axis in R4 */
    public static enum Vector!(T,4) Z_AXIS = Vector!(T,4)(0, 0, 1, 0);
    /**W Axis in R4 (used like a scale factor with 3d Maths with 4d Matrixes)*/
    public static enum Vector!(T,4) W_AXIS = Vector!(T,4)(0, 0, 0, 1);
  }
  
  
  /**
  * Build a new Vector from a set of initial values
  * If no there values for z and w, will be set to 0
  * Params:
  *	x = X coord
  *	y = Y coord
  *	z = Z coord
  *	w = W coord (scale factor in 3d math)
  */
  this(in T x, in T y, in T z = 0, in T w = 1) {
    this.x = x;
    this.y = y;
    static if (dim >= 3)
      this.z = z;
    static if (dim >= 4) 
      this.w = w;
  }
  
  /**
  * Build a new Vector from a array
  * If no there values for y and z, will be set to 0. w sill beset to 1
  * Params:
  *	xs = Array with coords
  */
  this(in T[] xs) {
    size_t l = xs.length > dim? dim : xs.length;
    coor[0..l] = xs[0..l].dup;
    if (l <dim) {
      static if (dim <4) {
        coor[l..dim] = 0;
      } else {
        coor[l..dim-1] = 0;
        w = 1;
      }
    }
  }
  
  /+
  /**
  * Post-Blit
  */
  this(this) {
    coor = coor.dup;
  }+/
  
  unittest {
    Vec2r v1 = Vec2r(0,2);
    auto v2 = Vec2r(2,0);
    auto v3 = Vec2r([1,1]);
    float[] arr = [5,20.1,0,10,20];
    auto v4 = Vec4f(arr);
    assert (v3 == Vec2r(1,1));
    assert (v4 == Vec4f(5,20.1,0,10));
    
    assert (v2.x == 2);   // By Coordinate name
    assert (v2.coor[0] == 2); // By direct access to array (only internal use)
    assert (v2[0] == 2);   // By opIndex
  }
  
  // Basic Properties **********************************************************
  
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
  
  unittest {
    auto v = Vec3f(10,0,0);
    auto v2 = Vec4d(1,1,1,1);
    assert (v.length == 10.0);
    assert (v2.sq_length == 4.0);
  }
  
  // Operations ****************************************************************
  
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
  
  unittest {
    // Check change of sign
    auto v = Vec2r([1,1]);
    auto n = -v;
    assert (n.x == -1);
    assert (n.y == -1);
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
      return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) 
          && approxEqual(z, rhs.z);
    static if (dim == 4) 
      return approxEqual(x, rhs.x) && approxEqual(y, rhs.y) 
          && approxEqual(z, rhs.z) && approxEqual(w, rhs.w);
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
      return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff) 
              && approxEqual(z, rhs.z, maxDiff);
    static if (dim == 4) 
      return approxEqual(x, rhs.x, maxDiff) && approxEqual(y, rhs.y, maxDiff) 
              && approxEqual(z, rhs.z, maxDiff) 
              && approxEqual(w, rhs.w, maxDiff);
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
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) 
              && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff);
    static if (dim == 3) 
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) 
              && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) 
              && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff);
    static if (dim == 4) 
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) 
              && approxEqual(y, rhs.y, maxRelDiff, maxAbsDiff) 
              && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff) 
              && approxEqual(w, rhs.w, maxRelDiff, maxAbsDiff);
  }
  
  /**
  * Define unary operators + and -
  */
  Vector opUnary(string op) () const
    if (op == "+" || op == "-") {
    static if (dim == 2) 
      return Vector( mixin(op~"x"), mixin(op~"y"));
    static if (dim == 3) 
      return Vector( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"));
    static if (dim == 4) 
      return Vector( mixin(op~"x"), mixin(op~"y"), mixin(op~"z"),
                      mixin(op~"w"));
  }
  
  unittest {
    // Check equality
    auto v =  Vec4f(4, 2, 1, 0);
    auto v2 = Vec4f(4, 2, 3, 1);
    assert (v2 != v);
    assert (v == v);
    assert (v.equal(v));
    assert (! v.equal(v2));
  }
  
  /**
  * Define binary operator + and -
  */
  Vector opBinary(string op) (ref const Vector rhs) const
    if (op == "+" || op == "-") {
    static if (dim == 2) 
      return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"));
    static if (dim == 3) 
      return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"),
                      mixin("z"~op~"rhs.z"));	
    static if (dim == 4) 
      return Vector( mixin("x"~op~"rhs.x"), mixin("y"~op~"rhs.y"),
                      mixin("z"~op~"rhs.z"), mixin("w"~op~"rhs.w"));
  }
  
  unittest {
    // Check addition
    auto v1 = Vec4d(1,2,3,4);
    auto v2 = Vec4d(1,1,1,1);
    auto v12 = v1 + v2;
    auto v21 = v2 + v1;
    
    // Symetry
    assert(v12 == v21);
    // Value
    assert(v12 == Vec4d(2,3,4,5));
    
    // Check subtraction
    auto v1m2 = v1 - v2;
    assert(v1m2 == Vec4d(0,1,2,3));
    // Subtraction asimetry
    auto v2m1 = v2 - v1;
    assert(v2m1 == Vec4d(0,-1,-2,-3));
  }
  
  /**
  * Define Scalar multiplication
  */
  Vector opBinary(string op) (in T rhs) const
    if (op == "*" ) {
    static if (dim == 2) 
      return Vector(x *rhs, y *rhs);
    static if (dim == 3) 
      return Vector(x *rhs, y *rhs, z *rhs);	
    static if (dim == 4) 
      return Vector(x *rhs, y *rhs, z *rhs, w *rhs);
  }
  
  unittest {
    // Product by Scalar
    auto vten = Vec4r(0,1,2,3) *10.0L;
    assert(vten == Vec4r(0,10,20,30));
    vten = vten * (1/2.0L);
    assert(vten == Vec4r(0,5,10,15));		
  }
  
  /**
  * Define Dot Product
  */
  T opBinary(string op) (ref const Vector rhs) const
    if (op == "*" ) {
    T tmp = 0;
    foreach (i, x; coor) {
      tmp += x * rhs.coor[i];
    }
    return tmp;
  }
  
  unittest {
    // Dot Product and length
    auto v1 = Vec4f(10,0,0,0);
    auto v2 = Vec4f(1,1,1,0);
    auto v2dotv1 = v2 * v1;
    assert(v2dotv1 == 10.0);
    v1 = Vec4f(10,0,0,0);
    v2 = Vec4f(0,1,1,0);
    v2dotv1 = v2 * v1;
    assert(v2dotv1 == 0.0);
  }
  
  /**
  * Define Cross Product for R3 (operation c = a & b )
  */
  Vector opBinary(string op) (ref const Vector rhs) const 
    if (op == "&" && dim == 3)  {
    return Vector( coor[1]*rhs.coor[2] - coor[2]*rhs.coor[1],
                  coor[2]*rhs.coor[0] - coor[0]*rhs.coor[2],
                  coor[0]*rhs.coor[1] - coor[1]*rhs.coor[0]);
  }
  
  unittest {
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
  
  unittest {
    // Unitary Vector
    auto v = Vec4f(10,10,10,10);
    Vec4f uv = v; 
    uv.normalize();
    auto uv2 =  Vec4f(3,3,3,3);
    uv2.normalize();
    assert (uv == uv2);
    assert (uv == Vec4f(.5,.5,.5,.5));
    assert (!v.isUnit());
    assert (uv.isUnit());
  }
  
  /**
  * Obtains the projection two vectors
  * Params:
  *	a = Vector to project
  * b = Vector where project vector a
  * Returns : A Vector that it's projection of Vector a over Vector b
  */
  static Vector projectOnTo (Vector a,Vector b) {
    b.normalize();
    Vector ret = void; 
    auto s = a*b;
    ret = b * (a*b);
    return ret ;
  }
  
  /**
  * Obtains the projection of this vector over other vector
  * Params:
  * b = Vector where project this vector
  * Returns : A Vector that it's projection of this Vector over Vector b
  */
  Vector projectOnTo(Vector b) {
    return this.projectOnTo(this,b);
  }
  
  unittest {
    // Check projection
    auto b = Vec3f.X_AXIS;
    auto a = Vec3f(1,20,-30);
    auto proy = a.projectOnTo(b);
    assert (proy == Vec3f.X_AXIS );
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
  * Calculate the squared distance between two points pointed by this vector 
  * and other
  * Params :
  *	b = Vector B
  * Returns : Squared distance between the point pointed by this vector and
  *   other point
  */
  T sq_distance( in Vector b) {
    auto d = this - b;
    return d.sq_length;
  }
  
  unittest {
    // Test Distance
    auto v1 = Vec4d(1,1,1,1);
    auto v2 = Vec4d(-1,-1,-1,-1);
    auto dis = v1.distance(v2);
    auto dis_sq = v1.sq_distance(v2);
    assert (dis == 4 );
    assert (dis_sq == 16 );
  }
  
  /**
  * Rotation in R2
  * Params:
  * axis = Rotation axis
  * Returns : A vector that is the rotation of this vector
  */
  static if (dim == 2)
  Vector rotate (real angle) const {
    return Vector( x * cos(angle) - y * sin(angle),
                   x * sin(angle) + y * cos(angle) );
  }
  
  unittest {
    // Check Rotation in R2
    auto v = Vec2f(10,0);
    auto rot90 = v.rotate(PI_2);
    auto rotn90 = v.rotate(-PI_2);
    assert (rot90.length == v.length);
    assert (rotn90.length == v.length);
    assert (rot90.equal(Vec2f(0,10)));
    assert (rotn90.equal(Vec2f(0,-10)));
  }
  
  // Misc **********************************************************************
  
  /**
  * Checks that the vector not have a weird NaN value
  * Returns : True if this vector not have a NaN value
  */
  @property bool isOk() {
    if (isNaN(x) || isNaN(y)) return false;
    static if (dim >= 3)
      if (isNaN(z)) return false;
    static if (dim >= 4)
      if (isNaN(w)) return false;	
    return true;	
  }
  
  /**
  * Checks that the vector have finite values
  * Returns : True if this vector have finite values (not infinite value or
  *   NaNs)
  */
  @property bool isFinite() {
    if (std.math.isFinite(x) && std.math.isFinite(y)) {
      static if (dim >= 3)
        if (! std.math.isFinite(z)) return false;
      static if (dim >= 4)
        if (! std.math.isFinite(w)) return false;	
      return true;	
    } else {
      return false;
    }
  }
  
  unittest {
    // Check isOk() and isFinite()
    Vec4d vNaN;
    assert(!vNaN.isOk);
    assert(!vNaN.isFinite);
    Vec4d vOK = Vec4d(1,2,3,4);
    assert(vOK.isOk);
    assert(vOK.isFinite);
    Vec4d vInf = Vec4d(1.0 /0.0, 10, 0, 1);
    assert (!vInf.isFinite);
    assert (vInf.isOk);
  }
  
  /**
  * Casting method to convert to other vector types
  */
  Tout opCast( Tout ) () 
  if (isVector!(Tout)) {
    static assert (isVector!(Tout), "This type not is a Vector");
     
    Tout newVector = void; auto i = 0;
    static if (is (typeof(newVector.x) == typeof(this.x)))  {
      for (; i < dim && i< Tout.dim; i++)
        newVector.coor[i] =  coor[i];
    } else {
      for (; i < dim && i< Tout.dim; i++)
        newVector.coor[i] =  to!(typeof(newVector.x))(coor[i]); 
    }
    
    static if (Tout.dim >=3 && dim < 3) // Expands a small vector to a bigger
      newVector.z = 0; // Z by default to 0
    
    static if (Tout.dim == 4 && dim < 4) {
      newVector.w = 1; 
      // w is usually a scale factor used in 3d math with 4d matrixes
    }
    
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
  
  unittest {
    // Check casting
    auto v2 = Vec2r(10,-10);
    auto v2f = cast(Vec2f) v2;
    auto v2d = cast(Vec2d) v2;
    auto vecf = Vec2f(10,10);
    //auto vec4d = toImpl!(Vec4d, Vec2f) (vecf);
    //auto vec4d = to!Vec4d (vecf); // Not yet
    auto vec4d = cast(Vec4d) vecf;
    assert(is(typeof(v2f) == Vec2f));
    assert(is(typeof(v2d) == Vec2d));
    assert(is(typeof(vec4d) == Vec4d));
    auto r2 = Vec2d(5,4);
    auto r4 = cast(Vec4d) r2;
    assert (is(typeof(r4) == Vec4d));
    assert (r4 == Vec4d(5,4,0,1));
    r2 = cast(Vec2d) r4;
    assert (r2 == Vec2d(5,4));
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

unittest {
  auto v1 = Vec2f(5,7);
  assert (isVector!(typeof(v1))); // Yep it's a Vector
  assert (! isVector!(int));
}
