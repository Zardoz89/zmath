/**
Defines a squared Matrix with column-major ordering from 2x2 to 4x4 size
*/
module zmath.matrix;

import zmath.vector;

alias Mat2f = Matrix!(float, 2); /// 2D squared matrix of floats
alias Mat3f = Matrix!(float, 3); /// 3D squared matrix of floats
alias Mat4f = Matrix!(float, 4); ///	4D squared matrix of floats

alias Mat2d = Matrix!(double, 2); /// 2D squared matrix of doubles
alias Mat3d = Matrix!(double, 3); /// 3D squared matrix of doubles
alias Mat4d = Matrix!(double, 4); ///	4D squared matrix of doubles

alias Mat2r = Matrix!(real, 2); /// 2D squared matrix of reals
alias Mat3r = Matrix!(real, 3); /// 3D squared matrix of reals
alias Mat4r = Matrix!(real, 4); ///	4D squared matrix of reals

/**
 * Defines a squared Matrix of n = 2, 3 or 4 size, like a linear array of numbers in a row-major order
 */
public struct Matrix(T, size_t dim_)
      if (__traits(isFloating, T) && dim_ >= 2 && dim_ <= 4)
{
nothrow:
  static enum size_t dim = dim_; /// Matrix Dimension
  static enum size_t cells = dim * dim; /// Matrix number of cells

  private alias VRow = Vector!(T, dim);
  private alias VCol = Vector!(T, dim);

  union {
    private T[cells] cell; /// Matrix like of a array of cells in row-major
    private T[dim][dim] c; /// Matrix like of a 2d array of cells in row-major
    private VRow[dim_] row; /// Matrix like of a array of rows vectors
  }

  // Consts
  static if (dim == 2) { // 2x2
    // dfmt off
    /** Zero R2 matrix */
    public static enum Matrix!(T,2) ZERO = Matrix!(T,2)([
        0, 0,
        0, 0
    ]);
    /** Identity R2 matrix */
    public static enum Matrix!(T,2) IDENTITY = Matrix!(T,2)([
        1, 0,
        0, 1
    ]);
    // dfmt on
  } else static if (dim == 3) { // 3x3
    // dfmt off
    /** Zero R3 matrix */
    public static enum Matrix!(T,3) ZERO = Matrix!(T,3)([
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    ]);
    /** Identity R3 matrix */
    public static enum Matrix!(T,3) IDENTITY = Matrix!(T,3)([
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    ]);
    // dfmt on
  } else static if (dim == 4) { // 4x4
    // dfmt off
    /** Zero R4 matrix */
    public static enum Matrix!(T,4) ZERO = Matrix!(T,4)([
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
    ]);
    /** Identity R4 matrix */
    public static enum Matrix!(T,4) IDENTITY = Matrix!(T,4)([
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]);
    // dfmt on
  }

  // Basic properties *****************************************************

  /**
   * Returns i, j cell
   */
  @nogc T opIndex(size_t row, size_t col) pure const {
    return this.row[row][col];
  }

  /**
   * Assigns a new cell value
   */
  @nogc void opIndexAssign(K)(K c, size_t row, size_t col) if (is(K : real)) {
    this.row[row][col] = c;
  }

  /**
   * Returns i row vector
   */
  @nogc VRow opIndex(size_t i) pure const {
    return this.row[i];
  }

  /**
   * Assigns a new row vector
   */
  @nogc void opIndexAssign(K)(K v, size_t i) if (isVector!(K)) {
    static assert(v.dim <= dim, "Vector has a bigger dimension that this matrix.");
    static if (!is(K == VRow)) {
      this.row[i] = cast(VRow) v;
    } else {
      this.row[i] = v;
    }
  }

  // Operations ***************************************

  /**
   * Define Equality
   */
  @nogc bool opEquals()(auto ref const Matrix rhs) pure const {
    if (this.row[0] != rhs.row[0] || this.row[1] != rhs.row[1]) {
      return false;
    }
    static if (dim >= 3) {
      if (this.row[2] != rhs.row[2]) {
        return false;
      }
    } else static if (dim >= 4) {
      if (this.row[3] != rhs.row[3]) {
        return false;
      }
    }
    return true;
  }

  /**
  * Approximated equality
  */
  @nogc bool approxEqual()(auto ref const Matrix rhs, T maxRelDiff = 1e-2, T maxAbsDiff = 1e-05) pure const {
    if (!this.row[0].approxEqual(rhs.row[0], maxRelDiff, maxAbsDiff)
        || !this.row[1].approxEqual(rhs.row[1], maxRelDiff, maxAbsDiff)) {
      return false;
    }
    static if (dim >= 3) {
      if (!this.row[2].approxEqual(rhs.row[2], maxRelDiff, maxAbsDiff)) {
        return false;
      }
    } else static if (dim >= 4) {
      if (!this.row[3].approxEqual(rhs.row[3], maxRelDiff, maxAbsDiff)) {
        return false;
      }
    }
    return true;
  }

  /**
  * Define unary operators + and -
  */
  @nogc Matrix opUnary(string op)() pure const if (op == "+" || op == "-") {
    import std.conv : to;

    string makeMix(size_t d, string op) {
      string ret = "return Matrix( [";
      foreach (i; 0 .. (d * d)) {
        ret ~= op ~ "cell[" ~ to!string(i) ~ "] ";
        if (i != (d * d - 1))
          ret ~= ",";
      }
      ret ~= "] );";
      return ret;
    }

    static if (dim == 2) {
      mixin(makeMix(2, op));
    } else static if (dim == 3) {
      mixin(makeMix(3, op));
    } else static if (dim == 4) {
      mixin(makeMix(4, op));
    }
  }

  /**
  * Define binary operator + and -
  */
  @nogc Matrix opBinary(string op)(auto ref const Matrix rhs) pure const
      if ((op == "+" || op == "-") && dim >= 2 && dim <= 4) {
    import std.conv : to;

    string makeMixOp(size_t d, string op) {
      string ret = "return Matrix( [";
      foreach (i; 0 .. (d * d)) {
        ret ~= "cell[" ~ to!string(i) ~ "] " ~ op ~ "rhs.cell[" ~ to!string(i) ~ "] ";
        if (i != (d * d - 1))
          ret ~= ",";
      }
      ret ~= "] );";
      return ret;
    }

    mixin(makeMixOp(dim, op));
  }

  /**
  * Define Scalar multiplication
  */
  @nogc Matrix opBinary(string op)(in real rhs) pure const if (op == "*" || op == "/") {
    T[this.cells] ret = mixin("this.cell[]" ~ op ~ "cast(T) rhs");
    return Matrix(ret);
  }

  /**
  * Define Scalar multiplication
  */
  @nogc Matrix opBinaryRight(string op)(in real rhs) pure const if (op == "*" || op == "/") {
    T[this.cells] ret = mixin("this.cell[]" ~ op ~ "cast(T) rhs");
    return Matrix(ret);
  }


  /**
  * Define Matrix Product
  */
  @nogc Matrix opBinary(string op)(auto ref const Matrix rhs) pure const if (op == "*") {
    Matrix mat;
    foreach (size_t i; 0 .. dim) { // Runs along result rows
      foreach (size_t j; 0 .. dim) { // Runs along result columns
        static if (dim == 2) {
          mat[i, j] = this[i, 0] * rhs[0, j] + this[i, 1] * rhs[1, j];
        } else static if (dim == 3) {
          mat[i, j] = this[i, 0] * rhs[0, j] + this[i, 1] * rhs[1, j] + this[i, 2] * rhs[2, j];
        } else static if (dim == 4) {
          mat[i, j] = this[i, 0] * rhs[0, j] + this[i, 1] * rhs[1, j] + this[i, 2] * rhs[2, j] + this[i, 3] * rhs[3, j];
        }
      }
    }

    return mat;
  }

  /**
  * Define Matrix x Vector
  */
  @nogc VCol opBinary(string op)(auto ref const VRow rhs) pure const if (op == "*") {
    static assert(rhs.dim == dim, "The vector dimension must be equal to Matrix dimmension");

    VCol ret;
    foreach (size_t i; 0 .. dim) { // Runs along result coords
      static if (dim == 2) {
        ret[i] = this[i, 0] * rhs[0] + this[i, 1] * rhs[1];
      } else static if (dim == 3) {
        ret[i] = this[i, 0] * rhs[0] + this[i, 1] * rhs[1] + this[i, 2] * rhs[2];
      } else static if (dim == 4) {
        ret[i] = this[i, 0] * rhs[0] + this[i, 1] * rhs[1] + this[i, 2] * rhs[2] + this[i, 3] * rhs[3];
      }
    }
    return ret;
  }

  /**
  * Define Vector x Matrix
  */
  @nogc VCol opBinaryRight(string op)(auto ref const VCol lhs) pure const if (op == "*") {
    static assert(lhs.dim == dim, "The vector dimension must be equal to Matrix dimmension");
    VCol ret;
    foreach (size_t j; 0 .. dim) { // Runs along result coords
      static if (dim == 2) {
        ret[j] = this[0, j] * lhs[0] + this[1, j] * lhs[1];
      } else static if (dim == 3) {
        ret[j] = this[0, j] * lhs[0] + this[1, j] * lhs[1] + this[2, j] * lhs[2];
      } else static if (dim == 4) {
        ret[j] = this[0, j] * lhs[0] + this[1, j] * lhs[1] + this[2, j] * lhs[2] + this[3, j] * lhs[3];
      }
    }
    return ret;
  }

  /**
  * Returns transposed matrix
  *
  * It can be used to "convert" Matrix to a row ordered matrix
  */
  @nogc Matrix transpose() pure const {
    Matrix mat;
    foreach (i; 0 .. dim) { // Runs along result rows
      foreach (j; 0 .. dim) { // Runs along result columns
        mat[j, i] = this[i, j];
      }
    }
    return mat;
  }

  /**
  * Returns Determinant of this matrix
  */
  @nogc T determinant() pure const {
    static if (dim == 2) {
      return this.c[0][0] * this.c[1][1] - (this.c[1][0] * this.c[0][1]);
    } else static if (dim == 3) {
      T aei = this.c[0][0] * this.c[1][1] * this.c[2][2];
      T bfg = this.c[0][1] * this.c[1][2] * this.c[2][0];
      T cdh = this.c[0][2] * this.c[1][0] * this.c[2][1];
      T afh = this.c[0][0] * this.c[1][2] * this.c[2][1];
      T bdi = this.c[0][1] * this.c[1][0] * this.c[2][2];
      T ceg = this.c[0][2] * this.c[1][1] * this.c[2][0];
      return aei + bfg + cdh - afh - bdi - ceg;
    } else {
      // dfmt off
      return (
           this.c[0][0] * this.c[1][1] - this.c[0][1] * this.c[1][0])
        * (this.c[2][2] * this.c[3][3] - this.c[2][3] * this.c[3][2])
        - (this.c[0][0] * this.c[1][2] - this.c[0][2] * this.c[1][0])
        * (this.c[2][1] * this.c[3][3] - this.c[2][3] * this.c[3][1])
        + (this.c[0][0] * this.c[1][3] - this.c[0][3] * this.c[1][0])
        * (this.c[2][1] * this.c[3][2] - this.c[2][2] * this.c[3][1])
        + (this.c[0][1] * this.c[1][2] - this.c[0][2] * this.c[1][1])
        * (this.c[2][0] * this.c[3][3] - this.c[2][3] * this.c[3][0])
        - (this.c[0][1] * this.c[1][3] - this.c[0][3] * this.c[1][1])
        * (this.c[2][0] * this.c[3][2] - this.c[2][2] * this.c[3][0])
        + (this.c[0][2] * this.c[1][3] - this.c[0][3] * this.c[1][2])
        * (this.c[2][0] * this.c[3][1] - this.c[2][1] * this.c[3][0]);
      // dfmt on
    }
  }

  /**
  * Calcs this matrix inverse
  * Returns: A matrix full of NaNs if this matrix isn't inverible (!isOk =1)
  */
  @nogc Matrix inverse() pure const {
    import std.math : approxEqual;

    Matrix mat;
    auto det = this.determinant();
    if (approxEqual(det, 0)) // Not invertible matrix
      return mat; // At this point Mat it's full of NaNs and a not valid Matrix

    static if (dim == 2) {
      mat[0, 0] = this[1, 1];
      mat[1, 1] = this[0, 0];
      mat[0, 1] = -this[0, 1];
      mat[1, 0] = -this[1, 0];
    } else if (dim == 3) {
      mat[0, 0] = this[1, 1] * this[2, 2] - this[1, 2] * this[2, 1];
      mat[0, 1] = this[0, 2] * this[2, 1] - this[0, 1] * this[2, 2];
      mat[0, 2] = this[0, 1] * this[1, 2] - this[0, 2] * this[1, 1];

      mat[1, 0] = this[1, 2] * this[2, 0] - this[1, 0] * this[2, 2];
      mat[1, 1] = this[0, 0] * this[2, 2] - this[0, 2] * this[2, 0];
      mat[1, 2] = this[0, 2] * this[1, 0] - this[0, 0] * this[1, 2];

      mat[2, 0] = this[1, 0] * this[2, 1] - this[1, 1] * this[2, 2];
      mat[2, 1] = this[0, 1] * this[2, 0] - this[0, 0] * this[2, 1];
      mat[2, 2] = this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0];
    } else { // dim = 4
      // dftm off
      mat[0, 0] = this[1, 1] * (this[2, 2] * this[3, 3] - this[2, 3] * this[3, 2])
        + this[1, 2] * (this[2, 3] * this[3, 1] - this[2, 1] * this[3, 3])
        + this[1, 3] * (this[2, 1] * this[3, 2] - this[2, 2] * this[3, 1]);
      mat[0, 1] = this[2, 1] * (this[0, 2] * this[3, 3] - this[0, 3] * this[3, 2])
        + this[2, 2] * (this[0, 3] * this[3, 1] - this[0, 1] * this[3, 3])
        + this[2, 3] * (this[0, 1] * this[3, 2] - this[0, 2] * this[3, 1]);
      mat[0, 2] = this[3, 1] * (this[0, 2] * this[1, 3] - this[0, 3] * this[1, 2])
        + this[3, 2] * (this[0, 3] * this[1, 1] - this[0, 1] * this[1, 3])
        + this[3, 3] * (this[0, 1] * this[1, 2] - this[0, 2] * this[1, 1]);
      mat[0, 3] = this[0, 1] * (this[1, 3] * this[2, 2] - this[1, 2] * this[2, 3])
        + this[0, 2] * (this[1, 1] * this[2, 3] - this[1, 3] * this[2, 1])
        + this[0, 3] * (this[1, 2] * this[2, 1] - this[1, 1] * this[2, 2]);

      mat[1, 0] = this[1, 2] * (this[2, 0] * this[3, 3] - this[2, 3] * this[3, 0])
        + this[1, 3] * (this[2, 2] * this[3, 0] - this[2, 0] * this[3, 2])
        + this[1, 0] * (this[2, 3] * this[3, 2] - this[2, 2] * this[3, 3]);
      mat[1, 1] = this[2, 2] * (this[0, 0] * this[3, 3] - this[0, 3] * this[3, 0])
        + this[2, 3] * (this[0, 2] * this[3, 0] - this[0, 0] * this[3, 2])
        + this[2, 0] * (this[0, 3] * this[3, 2] - this[0, 2] * this[3, 3]);
      mat[1, 2] = this[3, 2] * (this[0, 0] * this[1, 3] - this[0, 3] * this[1, 0])
        + this[3, 3] * (this[0, 2] * this[1, 0] - this[0, 0] * this[1, 2])
        + this[3, 0] * (this[0, 3] * this[1, 2] - this[0, 2] * this[1, 3]);
      mat[1, 3] = this[0, 2] * (this[1, 3] * this[2, 0] - this[1, 0] * this[2, 3])
        + this[0, 3] * (this[1, 0] * this[2, 2] - this[1, 2] * this[2, 0])
        + this[0, 0] * (this[1, 2] * this[2, 3] - this[1, 3] * this[2, 2]);

      mat[2, 0] = this[1, 3] * (this[2, 0] * this[3, 1] - this[2, 1] * this[3, 0])
        + this[1, 0] * (this[2, 1] * this[3, 3] - this[2, 3] * this[3, 1])
        + this[1, 1] * (this[2, 3] * this[3, 0] - this[2, 0] * this[3, 3]);
      mat[2, 1] = this[2, 3] * (this[0, 0] * this[3, 1] - this[0, 1] * this[3, 0])
        + this[2, 0] * (this[0, 1] * this[3, 3] - this[0, 3] * this[3, 1])
        + this[2, 1] * (this[0, 3] * this[3, 0] - this[0, 0] * this[3, 3]);
      mat[2, 2] = this[3, 3] * (this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0])
        + this[3, 0] * (this[0, 1] * this[1, 3] - this[0, 3] * this[1, 1])
        + this[3, 1] * (this[0, 3] * this[1, 0] - this[0, 0] * this[1, 3]);
      mat[2, 3] = this[0, 3] * (this[1, 1] * this[2, 0] - this[1, 0] * this[2, 1])
        + this[0, 0] * (this[1, 3] * this[2, 1] - this[1, 1] * this[2, 3])
        + this[0, 1] * (this[1, 0] * this[2, 3] - this[1, 3] * this[2, 0]);

      mat[3, 0] = this[1, 0] * (this[2, 2] * this[3, 1] - this[2, 1] * this[3, 2])
        + this[1, 1] * (this[2, 0] * this[3, 2] - this[2, 2] * this[3, 0])
        + this[1, 2] * (this[2, 1] * this[3, 0] - this[2, 0] * this[3, 1]);
      mat[3, 1] = this[2, 0] * (this[0, 2] * this[3, 1] - this[0, 1] * this[3, 2])
        + this[2, 1] * (this[0, 0] * this[3, 2] - this[0, 2] * this[3, 0])
        + this[2, 2] * (this[0, 1] * this[3, 0] - this[0, 0] * this[3, 1]);
      mat[3, 2] = this[3, 0] * (this[0, 2] * this[1, 1] - this[0, 1] * this[1, 2])
        + this[3, 1] * (this[0, 0] * this[1, 2] - this[0, 2] * this[1, 0])
        + this[3, 2] * (this[0, 1] * this[1, 0] - this[0, 0] * this[1, 1]);
      mat[3, 3] = this[0, 0] * (this[1, 1] * this[2, 2] - this[1, 2] * this[2, 1])
        + this[0, 1] * (this[1, 2] * this[2, 0] - this[1, 0] * this[2, 2])
        + this[0, 2] * (this[1, 0] * this[2, 1] - this[1, 1] * this[2, 0]);
      // dfmt on
    }

    return mat / det;
  }

  // Misc***********************************************************************

  /**
  * Checks that the matrix not have a weird NaN value
  */
  @nogc bool isOk() pure const {
    foreach (c; this.row) {
      if (!c.isOk())
        return false;
    }
    return true;
  }

  /**
  * Checks that the matrix have finite values
  */
  @nogc bool isFinite() pure const {
    foreach (c; this.row) {
      if (!c.isFinite())
        return false;
    }
    return true;
  }

  /**
   * Return a pointer of the internal array
   */
  @property @nogc T* ptr() pure nothrow {
    return this.cell.ptr;
  }

  unittest {
    auto m = Mat3f.IDENTITY;
    auto pt = m.ptr;
    pt[0] = 5;
    assert(m[0, 0] == 5);
    assert(pt[0] == 5);
  }

  /**
  * Casting operation that allow convert between Matrix types
  */
  @nogc Tout opCast(Tout)() pure const if (isMatrix!(Tout) && Tout.dim >= dim) {
    static assert(isMatrix!(Tout), "This type not is a Matrix");
    static assert(Tout.dim >= dim, "Original Matrix bigger that destiny" ~ " Vector");
    Tout newMat;
    auto i = 0;
    static if (is(typeof(newMat[0, 0]) == typeof(this[0, 0]))) {
      for (; i < dim; i++)
        newMat[i] = this[i];
    } else {
      for (; i < dim; i++) //Vector Casting auto expands rows
        newMat[i] = cast(newMat.VCol) this[i];
    }

    // Expands a small matrix to a bigger dimension with a full 0's columns
    for (; i < Tout.dim; i++)
      newMat[i] = cast(newMat.VCol) newMat.VCol.ZERO;

    return newMat;
  }

  /**
  * Returns a visual representation of this matrix
  */
  string toString() const {
    import std.conv : to;

    string ret; // I -> row, j->col
    try {
      foreach (i; 0 .. dim) {
        ret ~= "|";
        foreach (j; 0 .. dim) {
          ret ~= to!string(this.c[i][j]);
          if (j < (dim - 1))
            ret ~= ", ";
        }
        ret ~= "|";
        if (i < (dim - 1))
          ret ~= "\n";
      }
    } catch (Exception ex) {
      assert(false); // This never should happen
    }
    return ret;
  }
}

/**
* Say if a thing it's a Matrix
*/
template isMatrix(T) {
  //immutable bool isMatrix = is(T == Matrix);
  immutable bool isMatrix = __traits(compiles, () {
    T t;
    static assert(T.dim >= 2 && T.dim <= 4);
    auto cell = t.cell;
    auto cell00 = t[0, 0];
    auto col0 = t[0];
    auto coln = t[T.dim];
    t.isOk();

    // TODO : Should test for methods ?
  });
}

static assert(Mat2f.sizeof == 4 * 2*2);
static assert(Mat3d.sizeof == 8 * 3*3);
static assert(Mat4f.sizeof == 4 * 4*4);

unittest {
  // Fundamental Matrixes
  auto z = Mat4d.ZERO;
  auto ide = Mat4d.IDENTITY;

  foreach (ind; 0 .. z.cells) {
    assert(z.cell[ind] == 0.0L);
  }

  foreach (i; 0 .. ide.dim) {
    foreach (j; 0 .. ide.dim) {
      if (i == j) {
        assert(ide[i, j] == 1.0L);
      } else {
        assert(ide[i, j] == 0.0L);
      }
    }
  }
}

unittest {
  auto m = Mat4f.IDENTITY;
  assert(m[0, 0] == 1);
  assert(m[1, 1] == 1);
  assert(m[2, 2] == 1);
  assert(m[3, 3] == 1);
  assert(m[3, 0] == 0);
  m[3, 0] = 5;
  assert(m[3, 0] == 5);
}

unittest {
  // Testing opIndexAssign of rows
  auto m = Mat3r.IDENTITY;
  const v = Vec2r(-1, -2);
  const v2 = Vec3f(1, 2, 3);
  m[1] = v; // Ha! internal vector casting
  m[2] = v2;
  assert(m[1, 0] == -1);
  assert(m[1, 1] == -2);
  assert(m[1, 2] == 0);
  assert(m[2, 0] == 1);
  assert(m[2, 1] == 2);
  assert(m[2, 2] == 3);
}

unittest {
  // Check == and !=
  auto z = Mat4r.ZERO;
  auto ide = Mat4r.IDENTITY;
  assert(z != ide);
  assert(ide == ide);
  assert(ide.approxEqual(ide));
  assert(!ide.approxEqual(z));
}

unittest {
  // Check + - unary operators
  auto n_ide = -Mat4r.IDENTITY;
  foreach (i; 0 .. Mat4r.dim) {
    foreach (j; 0 .. Mat4r.dim) {
      if (i == j) {
        assert(n_ide[i, j] == -1.0L);
      } else {
        assert(n_ide[i, j] == 0.0L);
      }
    }
  }
}

unittest {
  // Check Addition and subtraction
  // dfmt off
  auto a = Mat4r([1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 1, 1]);
  // dfmt on
  auto b = a + Mat4r.IDENTITY;
  auto c = Mat4r.ZERO - a;

  foreach (i; 0 .. Mat4r.dim) {
    foreach (j; 0 .. Mat4r.dim) {
      assert(c[i, j] == -1.0);
      if (i == j) {
        assert(b[i, j] == 2.0L);
      } else {
        assert(b[i, j] == 1.0L);
      }
    }
  }
}

unittest {
  // Scalar multiplication
  // dfmt off
  auto a = Mat4f([2, 1, 1, 1,
                1, 2, 1, 1,
                1, 1, 2, 1,
                1, 1, 1, 2]);
  // dfmt on
  a = a * 2;
  foreach (i; 0 .. Mat4f.dim) {
    foreach (j; 0 .. Mat4f.dim) {
      if (i == j) {
        assert(a[i, j] == 4.0);
      } else {
        assert(a[i, j] == 2.0);
      }
    }
  }

  // dfmt off
  auto b = Mat4f([
      2, 1, 1, 1,
      1, 2, 1, 1,
      1, 1, 2, 1,
      1, 1, 1, 2]);
  // dfmt on

  b = b * 2;
  assert(a == b);
}

unittest {
  // Check Matrix multiplication - Neutral element
  // dfmt off
  const a = Mat4f([
      2, 1, 1, 1,
      1, 2, 1, 1,
      1, 1, 2, 1,
      1, 1, 1, 2]);
  // dfmt on
  const axi = a * Mat4f.IDENTITY;
  const ixa = Mat4f.IDENTITY * a;
  assert(a == axi, "Multiplication by identity isn't correct");
  assert(axi == ixa);

  // dfmt off
  const a3 = Mat3r([
      1,4,7,
      2,5,8,
      3,6,9]);
  const b3 = Mat3r([
      9,6,3,
      8,5,2,
      7,4,1]);
  const shouldBeC3 = Mat3r([
       90, 54, 18,
      114, 69, 24,
      138, 84, 30]);
  const shouldBeD3 = Mat3r([
      30, 84, 138,
      24, 69, 114,
      18, 54,  90]);
  // dfmt on
  const c3 = a3 * b3;
  const d3 = b3 * a3;
  assert(c3 == shouldBeC3, "A * B should be \n" ~ shouldBeC3.toString() ~ " but was\n" ~ c3.toString());
  assert(c3 != d3);
  assert(d3 == shouldBeD3, "B * A should be \n" ~ shouldBeD3.toString() ~ " but was\n" ~ d3.toString());
}

unittest {
  // dfmt off
  const mat = Mat2f([
      0, -1,
      1,  0]);
  // dfmt on
  const vec = Vec2f(1, 2);
  const vecl = vec * mat;
  const vecr = mat * vec;
  assert(vec != vecl);
  assert(vec != vecr);
  assert(vecl != vecr);
  assert(vecr == Vec2f(-2, 1));
  assert(vecl == Vec2f(2, -1));

  // TODO : Do it with other type of matrix
}

unittest { // Check Matrix transposition
  auto ide = Mat4f.IDENTITY;
  auto ide_t = ide.transpose();
  assert(ide == ide_t);
  auto a = Mat3r([1, 4, 7, 2, 5, 8, 3, 6, 9]);
  auto a_t = a.transpose();
  auto should_a_t = Mat3r([1, 2, 3, 4, 5, 6, 7, 8, 9]);
  assert(a_t == should_a_t);
}

unittest { // Check Matrix Determinant
  auto ide = Mat4f.IDENTITY;
  auto ide_t = ide.transpose;
  assert(ide_t.determinant == 1);
  assert(ide.determinant == 1);
  // dfmt off
  auto n = Mat4r([1, 20, 1, 0,
                1,  5, 1, 1,
                1, -1, 3, 1,
                1, 10, 1, 1]);
  // dfmt on
  auto d_n = n.determinant;
  assert(d_n == -10);
  auto b = Mat4r.IDENTITY;
  b = b * 2;
  b = n + b;
  auto d_b = b.determinant;
  auto prod = n * b;
  auto d_prod = prod.determinant;
  assert(d_n * d_b == d_prod);

  auto ide3 = Mat3f.IDENTITY;
  assert(ide3.determinant == 1);
  auto ide2 = Mat2f.IDENTITY;
  assert(ide2.determinant == 1);
}

unittest { // Check matrix inversion
  // dfmt off
  auto n = Mat4r([
      1, 20, 1, 0,
      1,  5, 1, 1,
      1, -1, 3, 1,
      1, 10, 1, 1
  ]);
  // dfmt on
  auto n_inv = n.inverse;
  // dfmt off
  assert(n_inv == Mat4r([
      1,  5.1, -0.5, -4.6,
      0, -0.2,    0,  0.2,
      0, -1.1,  0.5,  0.6,
    -1,   -2,    0,    3
  ]), "n_inv\n" ~ n_inv.toString());
  // dfmt on
  auto a = Mat3r([1, 4, 7, 2, 5, 8, 3, 6, 9]);
  auto a_inv = a.inverse;
  assert(!a_inv.isOk()); //not have inverse, so returns a matrix full of NaN
  // dfmt off
  auto m = Mat2r([
      10, 5,
      5, 2
  ]);
  // dfmt on
  auto m_inv = m.inverse;
  // dfmt off
  assert (m_inv == Mat2r([
        -.4, 1,
        1, -2
  ]));
  // dfmt on
}

unittest {
  auto ide = Mat4d.IDENTITY;
  assert(ide.isOk());
  assert(ide.isFinite());
  Mat3r full_nan;
  assert(!full_nan.isOk());
  assert(!full_nan.isFinite());
}

unittest {
  // Check casting
  auto mat2f = Mat2f.IDENTITY;
  mat2f[1, 0] = 5;
  auto mat2r = cast(Mat2r) mat2f;
  auto mat2d = cast(Mat2d) mat2f;
  auto mat4d = cast(Mat4d) mat2r;
  assert(is(typeof(mat2r) == Mat2r));
  assert(is(typeof(mat2d) == Mat2d));
  assert(is(typeof(mat4d) == Mat4d));

  assert(mat2f[0, 0] == mat4d[0, 0]);
  assert(mat2f[0, 1] == mat4d[0, 1]);
  assert(mat2f[1, 0] == mat4d[1, 0]);
  assert(mat2f[1, 1] == mat4d[1, 1]);
}

unittest {
  auto ide = Mat4d.IDENTITY;
  auto v = Vec3f.X_AXIS;
  assert(isMatrix!(typeof(ide)));
  assert(!isMatrix!(typeof(v)));
}
