/**
Defines a Vector of float point typefrom 2 to 4 dimension
*/
module zmath.vector;

alias Vec2r = Vector!(real, 2); /// Alias of a 2d Vector with reals
alias Vec3r = Vector!(real, 3); /// Alias of a 3d Vector with reals
alias Vec4r = Vector!(real, 4); /// Alias of a 4d Vector with reals

alias Vec2d = Vector!(double, 2); /// Alias of a 2d Vector with doubles
alias Vec3d = Vector!(double, 3); /// Alias of a 3d Vector with doubles
alias Vec4d = Vector!(double, 4); /// Alias of a 4d Vector with doubles

alias Vec2f = Vector!(float, 2); /// Alias of a 2d Vector with floats
alias Vec3f = Vector!(float, 3); /// Alias of a 3d Vector with floats
alias Vec4f = Vector!(float, 4); /// Alias of a 4d Vector with floats

/**
 * N-Dimensional Vector over a FloatPoint type, where N must be 2,3 or 4
 */
public struct Vector(T, size_t dim_) if (__traits(isFloating, T)) {
nothrow:
  static enum size_t dim = dim_; /// Vector Dimension

  static assert(dim >= 2 && dim <= 4, "Not valid dimension size.");
  static assert(__traits(isFloating, T), "Type not is a Float Point type.");
  static assert(is(T : real), "Type not is like a Float Point type.");

  union {
    package T[dim] coor; /// Vector coords like Array

    struct {
      static if (dim >= 1) {
        T x; /// X coord
        alias r = x; /// R component
        alias roll = x; /// Euler roll angle
        alias bank = x; /// Euler roll angle
      }
      static if (dim >= 2) {
        T y; /// Y coord
        alias g = y; /// G component
        alias pitch = y; /// Euler pith angle
        alias attidue = y; /// Euler pith angle
      }
      static if (dim >= 3) {
        T z; /// Z coord
        alias b = z; /// B component
        alias yaw = z; /// Euler yaw angle
        alias heading = z; /// Euler yaw angle
      }
      static if (dim >= 4) {
        T w; /// W coord
        alias a = w; /// Alpha component
      }
    }
  }
  // Consts
  static if (dim == 2) { // for R2
    public static enum Vector!(T, 2) ZERO = Vector!(T, 2)(0, 0); /// Origin
    public static enum Vector!(T, 2) X_AXIS = Vector!(T, 2)(1, 0); /// X Axis in R2
    public static enum Vector!(T, 2) Y_AXIS = Vector!(T, 2)(0, 1); /// Y Axis in R2

  } else static if (dim == 3) { // for R3
    public static enum Vector!(T, 3) ZERO = Vector!(T, 3)(0, 0, 0); /// Origin
    public static enum Vector!(T, 3) X_AXIS = Vector!(T, 3)(1, 0, 0); /// X Axis in R3
    public static enum Vector!(T, 3) Y_AXIS = Vector!(T, 3)(0, 1, 0); /// Y Axis in R3
    public static enum Vector!(T, 3) Z_AXIS = Vector!(T, 3)(0, 0, 1); /// Z Axis in R3

  } else static if (dim == 4) { // for R4
    public static enum Vector!(T, 4) ZERO = Vector!(T, 4)(0, 0, 0, 0); /// Origin
    public static enum Vector!(T, 4) X_AXIS = Vector!(T, 4)(1, 0, 0, 0); /// X Axis in R4
    public static enum Vector!(T, 4) Y_AXIS = Vector!(T, 4)(0, 1, 0, 0); /// Y Axis in R4
    public static enum Vector!(T, 4) Z_AXIS = Vector!(T, 4)(0, 0, 1, 0); /// Z Axis in R4
    /** W Axis in R4 (used like a scale factor with 3d Maths with 4d Matrixes) */
    public static enum Vector!(T, 4) W_AXIS = Vector!(T, 4)(0, 0, 0, 1);
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
  @nogc this(in T x, in T y = 0, in T z = 0, in T w = 1) pure {
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
  this(in T[] xs) pure {
    size_t l = xs.length > dim ? dim : xs.length;
    this.coor[0 .. l] = xs[0 .. l].dup;
    if (l < dim) {
      static if (dim < 4) {
        this.coor[l .. dim] = 0;
      } else {
        this.coor[l .. dim - 1] = 0;
        w = 1;
      }
    }
  }

  // Basic Properties **********************************************************

  /**
   * Returns i coord of this vector
   */
  @nogc T opIndex(size_t i) pure const {
    return coor[i];
  }

  /**
   * Assigns a value to a i coord
   */
  void opIndexAssign(T c, size_t i) {
    coor[i] = c;
  }

  /**
   * Returns the actual length of this Vector
   */
  @property @nogc T length() pure const {
    import std.math : sqrt;

    return sqrt(sq_length);
  }

  /**
   * Returns the actual squared length of this Vector
   */
  @property @nogc T sq_length() pure const {
    return this * this;
  }

  // Operations ****************************************************************

  /**
   * Define Equality
   * Params:
   * rhs = Vector at rigth of '=='
   */
  @nogc bool opEquals()(auto ref const Vector rhs) pure const {
    if (x != rhs.x)
      return false;
    static if (dim >= 2) {
      if (y != rhs.y) {
        return false;
      }
    }
    static if (dim >= 3) {
      if (z != rhs.z) {
        return false;
      }
    }
    static if (dim == 4) {
      if (w != rhs.w) {
        return false;
      }
    }
    return true;
  }

  /**
   * Approximated equality with controlable precision
   * Params:
   * rhs = Vector to compare with this vector
   * maxRelDiff = Maximun relative difference
   * maxAbsDiff = Maximun absolute difference
   *
   * See: std.math : approxEqual
   */
  @nogc bool approxEqual()(auto ref const Vector rhs, T maxRelDiff = 1e-2, T maxAbsDiff = 1e-05) pure const {
    import std.math : approxEqual;

    static if (dim == 2) {
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y,
          rhs.y, maxRelDiff, maxAbsDiff);
    } else static if (dim == 3) {
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y,
          maxRelDiff, maxAbsDiff) && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff);
    } else static if (dim == 4) {
      return approxEqual(x, rhs.x, maxRelDiff, maxAbsDiff) && approxEqual(y, rhs.y, maxRelDiff,
          maxAbsDiff) && approxEqual(z, rhs.z, maxRelDiff, maxAbsDiff)
        && approxEqual(w, rhs.w, maxRelDiff, maxAbsDiff);
    }
  }

  /**
   * Define unary operators + and -
   */
  @nogc Vector opUnary(string op)() pure const if (op == "+" || op == "-") {
    static if (dim == 2) {
      return Vector(mixin(op ~ "x"), mixin(op ~ "y"));
    } else static if (dim == 3) {
      return Vector(mixin(op ~ "x"), mixin(op ~ "y"), mixin(op ~ "z"));
    } else static if (dim == 4) {
      return Vector(mixin(op ~ "x"), mixin(op ~ "y"), mixin(op ~ "z"), mixin(op ~ "w"));
    }
  }

  /**
   * Define binary operator + and -
   */
  @nogc Vector opBinary(string op)(auto ref const Vector rhs) pure const if (op == "+" || op == "-") {
    static if (dim == 2) {
      return Vector(mixin("x" ~ op ~ "rhs.x"), mixin("y" ~ op ~ "rhs.y"));
    } else static if (dim == 3) {
      return Vector(mixin("x" ~ op ~ "rhs.x"), mixin("y" ~ op ~ "rhs.y"), mixin("z" ~ op ~ "rhs.z"));
    } else static if (dim == 4) {
      return Vector(mixin("x" ~ op ~ "rhs.x"), mixin("y" ~ op ~ "rhs.y"),
          mixin("z" ~ op ~ "rhs.z"), mixin("w" ~ op ~ "rhs.w"));
    }
  }

  /**
   * Define Scalar multiplication
   */
  @nogc Vector opBinary(string op)(in T rhs) pure const if (op == "*" || op == "/") {
    static if (dim == 2) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"));
    } else static if (dim == 3) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"), mixin("z" ~ op ~ "rhs"));
    } else static if (dim == 4) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"),
          mixin("z" ~ op ~ "rhs"), mixin("w" ~ op ~ "rhs"));
    }
  }

  /**
   * Define Scalar multiplication
   */
  @nogc Vector opBinaryRight(string op)(in T rhs) pure const if(op == "*" || op == "/") {
    static if (dim == 2) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"));
    } else static if (dim == 3) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"), mixin("z" ~ op ~ "rhs"));
    } else static if (dim == 4) {
      return Vector(mixin("x" ~ op ~ "rhs"), mixin("y" ~ op ~ "rhs"),
          mixin("z" ~ op ~ "rhs"), mixin("w" ~ op ~ "rhs"));
    }
  }

  /**
   * Define Dot Product
   */
  @nogc T opBinary(string op)(auto ref const Vector rhs) pure const if (op == "*") {
    T tmp = 0;
    foreach (i, x; this.coor) {
      tmp += x * rhs.coor[i];
    }
    return tmp;
  }

  /**
   * Define Cross Product for R3 (operation c = a & b )
   */
  @nogc Vector opBinary(string op)(auto ref const Vector rhs) pure const
  if (op == "&" && dim == 3) {
    // dfmt off
    return Vector(
        coor[1] * rhs.coor[2] - coor[2] * rhs.coor[1],
        coor[2] * rhs.coor[0] - coor[0] * rhs.coor[2],
        coor[0] * rhs.coor[1] - coor[1] * rhs.coor[0]
        );
    // dfmt on
  }

  /**
  * It's a unitary vector (length == 1)
  * Returns : True if length approxEqual to 1.0
  */
  @property @nogc bool isUnit() pure const {
    import std.math : approxEqual, abs;

    return approxEqual(abs(this.sq_length - 1.0), 0);
  }

  /**
  * Normalize this vector
  */
  @nogc void normalize() {
    if (!isUnit()) {
      const T invL = 1 / this.length;
      static if (dim >= 2) {
        x = x * invL;
        y = y * invL;
      }
      static if (dim >= 3)
        z = z * invL;
      static if (dim == 4)
        w = w * invL;
    }
  }

  /**
  * Returns the unit vector of this vector
  */
  @property @nogc Vector unit() pure const {
    Vector ret = this;
    ret.normalize;
    return ret;
  }

  /**
  * Obtains the projection of two vectors
  * Params:
  *	a = Vector to project
  * b = Vector where project vector a
  * Returns : A Vector that it's projection of Vector a over Vector b
  */
  @nogc static Vector projectOnTo()(auto ref const Vector a, auto ref const Vector b) pure {
    const bNormalized = b.unit();
    Vector ret = b * (a * bNormalized);
    return ret;
  }

  /**
  * Obtains the projection of this vector over other vector
  * Params:
  * b = Vector where project this vector
  * Returns : A Vector that it's projection of this Vector over Vector b
  */
  @nogc Vector projectOnTo()(auto ref const Vector b) pure {
    return this.projectOnTo(this, b);
  }

  /**
  * Calculate the distance between two points pointed by this vector and other
  * Params :
  *	b = Vector B
  * Returns : Distance between the point pointed by this vector and other point
  */
  @nogc T distance()(auto ref const Vector b) pure {
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
  @nogc T sq_distance()(auto ref const Vector b) pure {
    auto d = this - b;
    return d.sq_length;
  }

  /**
  * Rotation in R2
  * Params:
  * angle = Rotation angle
  * Returns : A new vector that is the rotation of this vector
  */
  static if (dim == 2) {
    @nogc Vector rotate(real angle) pure const {
      import std.math : sin, cos;

      return Vector(x * cos(angle) - y * sin(angle), x * sin(angle) + y * cos(angle));
    }
  }

  // Misc **********************************************************************

  /**
   * Return a pointer of the internal array
   */
  @property @nogc T* ptr() pure {
    return this.coor.ptr;
  }

  /**
  * Checks that the vector not have a weird NaN value
  * Returns : True if this vector not have a NaN value
  */
  @property @nogc bool isOk() pure const {
    import std.math : isNaN;

    if (isNaN(x) || isNaN(y))
      return false;
    static if (dim >= 3)
      if (isNaN(z))
        return false;
    static if (dim >= 4)
      if (isNaN(w))
        return false;
    return true;
  }

  /**
  * Checks that the vector have finite values
  * Returns : True if this vector have finite values (not infinite value or
  *   NaNs)
  */
  @property @nogc bool isFinite() pure const {
    import std.math : isFinite;

    if (isFinite(x) && isFinite(y)) {
      static if (dim >= 3)
        if (!isFinite(z))
          return false;
      static if (dim >= 4)
        if (!isFinite(w))
          return false;
      return true;
    } else {
      return false;
    }
  }

  /**
  * Casting method to convert to other vector types
  */
  @nogc Tout opCast(Tout)() pure const if (isVector!(Tout)) {
    import std.conv : to;

    Tout newVector = void;
    auto i = 0;
    static if (is(typeof(newVector.x) == typeof(this.x))) {
      for (; i < dim && i < Tout.dim; i++)
        newVector.coor[i] = coor[i];
    } else {
      for (; i < dim && i < Tout.dim; i++)
        newVector.coor[i] = to!(typeof(newVector.x))(coor[i]);
    }

    static if (Tout.dim >= 3 && dim < 3) // Expands a small vector to a bigger
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
    import std.conv : to;
    try {
      string ret = "[" ~ to!string(x) ~ ", " ~ to!string(y);
      static if (dim >= 3)
        ret ~= ", " ~ to!string(z);
      static if (dim >= 4)
        ret ~= ", " ~ to!string(w);
      ret ~= "]";
      return ret;
    } catch (Exception ex) {
      assert(false); // This never should happen
    }
  }
}

/**
* Say if a thing it's a Vector
*/
template isVector(T) {
  immutable bool isVector = __traits(compiles, () {
    T t;
    static assert(T.dim >= 2 && T.dim <= 4);
    auto coor = t.coor;
    auto x = t.x;
    auto y = t.y;
    static if (t.dim >= 3)
      auto z = t.z;
    static if (t.dim >= 4)
      auto w = t.w;
    // TODO : Should test for methods ?
  });
}

static assert(Vec2f.sizeof == 8);
static assert(Vec3d.sizeof == 24);
static assert(Vec4f.sizeof == 16);

unittest {
  Vec2r v1 = Vec2r(0, 2);
  auto v2 = Vec2r(2, 0);
  auto v3 = Vec2r([1, 1]);
  float[] arr = [5, 20.1, 0, 10, 20];
  auto v4 = Vec4f(arr);

  assert(v2.x == 2); // By Coordinate name
  assert(v2.coor[0] == 2); // By direct access to array (only internal use)
}

unittest {
  auto v = Vec3f(10, 0, 0);
  auto v2 = Vec4d(1, 1, 1, 1);
  assert(v2[0] == 1); // By opIndex
  v2[0] = 123;
  assert(v2[0] == 123); // Assign by opIndex
  assert(v.length == 10.0);
  assert(Vec4d(1, 1, 1, 1).sq_length == 4);
}

unittest {
  // Check equality
  auto v1 = Vec4f(4, 2, 1, 0);
  auto v2 = Vec4f(4, 2, 3, 1);
  auto v3 = Vec4f(4, 2, 1, 0);
  auto v4 = Vec4f(4, 2, 1.00001, 0);
  assert(v2 != v1);
  assert(v1 == v1);
  assert(v1 == Vec4f(4, 2, 1, 0));
  assert(v3 == v1);
  assert(v1 != v4);
  assert(v1.approxEqual(v1));
  assert(!v1.approxEqual(v2));
  assert(v1.approxEqual(v4));
  assert(v1.approxEqual(Vec4f(4, 2, 1, 0)));
}

unittest {
  // Check change of sign
  auto v = Vec2r([1, 1]);
  auto n = -v;
  assert(n.x == -1);
  assert(n.y == -1);
}

unittest {
  // Check addition
  auto v1 = Vec4d(1, 2, 3, 4);
  auto v2 = Vec4d(1, 1, 1, 1);
  auto v12 = v1 + v2;
  auto v21 = v2 + v1;

  // Symetry
  assert(v12 == v21);
  // Value
  assert(v12 == Vec4d(2, 3, 4, 5));

  // Check subtraction
  auto v1m2 = v1 - v2;
  assert(v1m2 == Vec4d(0, 1, 2, 3));
  // Subtraction asimetry
  auto v2m1 = v2 - v1;
  assert(v2m1 == Vec4d(0, -1, -2, -3));
}

unittest {
  // Product by Scalar
  auto vten = Vec4r(0, 1, 2, 3) * 10.0L;
  assert(vten == Vec4r(0, 10, 20, 30));
  vten = vten / 2.0L;
  assert(vten == Vec4r(0, 5, 10, 15));
  vten = 10 * vten ;
  assert(vten == Vec4r(0, 50, 100, 150));
}

unittest {
  // Dot Product and length
  auto v1 = Vec4f(10, 1, 0, 0);
  auto v2 = Vec4f(1, 1, 1, 0);
  auto v2dotv1 = v2 * v1;
  assert(v2dotv1 == 11.0);
  v2 = Vec4f(0, 1, 1, 0);
  v2dotv1 = v2 * Vec4f(10, 0, 1, 0);
  assert(v2dotv1 == 1.0);
}

unittest {
  // Cross product
  auto r3v1 = Vec3r(0, 0, 10);
  auto r3v2 = Vec3r(2, 0, 2);
  auto cross0 = r3v1 & r3v1;
  auto cross12 = r3v1 & r3v2;
  auto cross21 = r3v2 & r3v1;

  assert(cross0.length == 0);
  assert(cross12.length == 20.0);
  assert(cross21.length == 20.0);

  assert(cross12 == Vec3r(0, 20.0L, 0));
  assert(cross21 == Vec3r(0, -20.0L, 0));
}

unittest {
  // Unitary Vector
  auto v = Vec4f(10, 10, 10, 10);
  Vec4f uv = v;
  uv.normalize();
  auto uv2 = Vec4f(3, 3, 3, 3);
  uv2.normalize();
  assert(uv == uv2);
  assert(uv == Vec4f(.5, .5, .5, .5));
  assert(!v.isUnit());
  assert(uv.isUnit());
}

unittest {
  // Check projection
  auto b = Vec3f.X_AXIS;
  auto a = Vec3f(1, 20, -30);
  auto proy = a.projectOnTo(b);
  assert(proy == Vec3f.X_AXIS);
}

unittest {
  // Test Distance
  auto v1 = Vec4d(1, 1, 1, 1);
  auto v2 = Vec4d(-1, -1, -1, -1);
  auto dis = v1.distance(v2);
  auto dis_sq = v1.sq_distance(v2);
  assert(dis == 4);
  assert(dis_sq == 16);

  dis = v1.distance(Vec4d(-1, -1, -1, -1));
  assert(dis == 4);
}

unittest {
  import std.math : PI_2;

  // Check Rotation in R2
  auto v = Vec2f(10, 0);
  auto rot90 = v.rotate(PI_2);
  auto rotn90 = v.rotate(-PI_2);
  assert(rot90.length == v.length);
  assert(rotn90.length == v.length);
  assert(rot90.approxEqual(Vec2f(0, 10)));
  assert(rotn90.approxEqual(Vec2f(0, -10)));
}

unittest {
  // Check isOk() and isFinite()
  Vec4d vNaN;
  assert(!vNaN.isOk);
  assert(!vNaN.isFinite);
  Vec4d vOK = Vec4d(1, 2, 3, 4);
  assert(vOK.isOk);
  assert(vOK.isFinite);
  Vec4d vInf = Vec4d(1.0 / 0.0, 10, 0, 1);
  assert(!vInf.isFinite);
  assert(vInf.isOk);
}

unittest {
  // Check casting
  auto v2 = Vec2r(10, -10);
  auto v2f = cast(Vec2f) v2;
  auto v2d = cast(Vec2d) v2;
  auto vecf = Vec2f(10, 10);
  //auto vec4d = toImpl!(Vec4d, Vec2f) (vecf);
  //auto vec4d = to!Vec4d (vecf); // Not yet
  auto vec4d = cast(Vec4d) vecf;
  assert(is(typeof(v2f) == Vec2f));
  assert(is(typeof(v2d) == Vec2d));
  assert(is(typeof(vec4d) == Vec4d));
  auto r2 = Vec2d(5, 4);
  auto r4 = cast(Vec4d) r2;
  assert(is(typeof(r4) == Vec4d));
  assert(r4 == Vec4d(5, 4, 0, 1));
  r2 = cast(Vec2d) r4;
  assert(r2 == Vec2d(5, 4));
}

unittest {
  const v1 = Vec2f(5, 7);
  assert(isVector!(typeof(v1))); // Yep it's a Vector
  assert(!isVector!(int));
}

