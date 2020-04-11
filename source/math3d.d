/**
Usefull functions to work in 3d math
*/
module zmath.math3d;

import zmath.matrix;
import zmath.vector;
import zmath.quaternion;

import std.math, std.traits;

/**
 * Checks if a type is a valid Matrix for 3d math over 4d
 */
template is4dMat(M) {
  immutable bool is4dMat = (isMatrix!M && __traits(compiles, () {
      M t;
      static assert(M.dim == 4);
    }));
}

unittest {
  assert(is4dMat!Mat4f);
  assert(!is4dMat!float);
  assert(is4dMat!Mat4d);
  assert(!is4dMat!Mat2f);
}

/**
 * Creates a projection matrix (similar to gluPerspective)
 * Params:
 * fov = Field of View in radians (PI_2 -> 90º)
 * aspect = Aspect ratio (x/y)
 * zMin = Near clipping plane
 * zMax = Far clipping plane
 * Returns a perspective projection matrix
 */
M perspectiveMat(M = Mat4f, T = float, U = float, V = float, W = float)(T fov,
    U aspect, V zMin, W zMax)
    if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V && isNumeric!W)
in {
  import std.math : PI;

  assert(zMin > 0, "Zmin equal or less that zero");
  assert(zMax > 0, "Zmax equal or less that zero");
  assert(fov > 0 && fov <= 2 * PI, "Not valid FOV");
  assert(aspect > 0, "Not valid aspect ratio");
}
body {
  import std.math : tan;

  M mat = M.ZERO; // all set to 0
  real yMax = zMin * tan(fov); // Max precision
  real yMin = -yMax;
  real xMax = yMax * aspect;
  real xMin = -xMax;

  mat[0, 0] = (2.0 * zMin) / (xMax - xMin);
  mat[1, 1] = (2.0 * zMin) / (yMax - yMin);
  mat[0, 2] = (xMax + xMin) / (xMax - xMin);
  mat[1, 2] = (yMax + yMin) / (yMax - yMin);
  mat[2, 2] = -((zMax + zMin) / (zMax - zMin));
  mat[3, 2] = -1.0;
  mat[2, 3] = -((2.0 * (zMax * zMin)) / (zMax - zMin));
  mat[3, 3] = 0;

  return mat;
}

unittest {
  const auto proy = perspectiveMat!Mat4d(PI_4, 800.0 / 600.0, 1.0, 100.0);
  // dfmt off
  const auto expected = Mat4d([
        0.75,  0,       0,  0,
           0,  1,       0,  0,
           0,  0, -1.0202, -1,
           0,  0, -2.0202,  0
  ]);
  // dfmt on
  assert(proy.approxEqual(expected),
      "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ proy.toString());
  // Check value take from glOrtho
}

/**
 * Creates a orthographic projection matrix
 * Params:
 * xMin = Left clipping plane
 * xMax = Right clipping plane
 * yMin = Bottom clipping plane
 * yMax = Top clipping plane
 * zMin = Near clipping plane
 * zMax = Far clipping plane
 * Returns a orthographic projection matrix
 */
M orthoMat(M = Mat4f, T = float, U = float, V = float, W = float, X = float, Y = float)(
    T xMin, U xMax, V yMin, W yMax, X zMin = 0, Y zMax = 1)
    if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V
      && isNumeric!W && isNumeric!X && isNumeric!Y)
in {
  assert(xMin != xMax);
  assert(yMin != yMax);
  assert(zMin != zMax);
}
body {
  M mat = M.ZERO;
  mat[0, 0] = 2.0 / (xMax - xMin);
  mat[1, 1] = 2.0 / (yMax - yMin);
  mat[2, 2] = -2.0 / (zMax - zMin);
  mat[0, 3] = -((xMax + xMin) / (xMax - xMin));
  mat[1, 3] = -((yMax + yMin) / (yMax - yMin));
  mat[2, 3] = -((zMax + zMin) / (zMax - zMin));
  mat[3, 3] = 1.0;
  return mat;
}

/**
 * Creates a orthographic projection matrix
 * Params:
 * width = Width of visible zone
 * heigth = Heigth of visible zone
 * deep = Deep of visible zone
 * Returns a orthographic projection matrix
 */
M orthoMat(M = Mat4f, T = float, U = float, V = float)(T width = 1, U height = 1, V deep = 1)
    if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V)
in {
  assert(width > 0);
  assert(height > 0);
  assert(deep > 0);
}
body {
  real x = -width / 2.0L; // Max precision
  real y = -height / 2.0L;
  return orthoMat!M(x, -x, y, -y, 0f, deep);
}

unittest {
  const ortho = orthoMat(100, 75, 100);
  // dfmt off
  const expected = Mat4f([
      0.02,         0,     0, 0,
         0, 0.0266667,     0, 0,
         0,         0, -0.02, 0,
        -0,        -0,    -1, 1
  ]);
  // dfmt on
  assert(ortho.approxEqual(expected),
      "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ ortho.toString());
  // Values from glOrtho
}

/**
 * Creates a trasnlation matrix from a 2d/3d Vector
 * Params:
 * v = 2d/3d Vector
 * Returns a translation matrix
 */
M translateMat(M = Mat4f, V = Vec3f)(V v) if (isVector!V && V.dim <= 3 && is4dMat!M) {
  M m = M.IDENTITY;
  m[3] = v; // internal cast in opIndexAssign set w=1 by default
  return m;
}

/**
 * Creates a translation matrix from xyz coords
 * Params:
 * x = X coord
 * y = Y coord
 * z = Z coord
 * Returns a translation matrix
 */
M translateMat(M = Mat4f, T = float)(T x, T y, T z) if (is4dMat!M && is(T : real)) {
  return translateMat!(M, Vector!(T, 3))(Vector!(T, 3)(x, y, z));
}

unittest {
  auto m1 = translateMat(1f, 2f, 3f);
  auto m2 = translateMat(Vec3f(1, 2, 3));
  assert(m1 == m2);
  assert(m1[0, 3] == 1);
  assert(m1[1, 3] == 2);
  assert(m1[2, 3] == 3);
  assert(m1[3, 3] == 1);
}

/**
 * Creates a uniform 3d scale matrix
 * Params:
 * s = scale factor
 * Returns a scale matrix
 */
M scaleMat(M = Mat4f, T = float)(T s) if (is4dMat!M && is(T : real)) {
  return M.IDENTITY * s;
}

/**
 * Creates a not uniform 3d scale matrix
 * Params:
 * x = X axis scale
 * y = Y axis scale
 * z = Z axis scale
 * Returns a scale matrix
 */
M scaleMat(M = Mat4f, T = float)(T x, T y, T z) if (is4dMat!M && is(T : real)) {
  M m = M.IDENTITY;
  m[0, 0] = x;
  m[1, 1] = y;
  m[2, 2] = z;
  return m;
}

/**
 * Creates a not uniform 3d scale matrix
 * Params:
 * v = Vector with scale factor. If is a 2d Vector, z axis scale factor is 1
 * Returns a scale matrix
 */
M scaleMat(M = Mat4f, V = Vec3f)(V v) if (isVector!V && V.dim <= 3 && is4dMat!M) {
  M m = M.IDENTITY;
  m[0, 0] = v.x;
  m[1, 1] = v.y;
  static if (V.dim == 2) {
    m[2, 2] = 1;
  } else {
    m[2, 2] = v.z;
  }
  return m;
}

unittest {
  import std.conv : to;

  auto m = scaleMat(10);
  assert(m[0, 0] == 10);
  assert(m[1, 1] == 10);
  assert(m[2, 2] == 10);
  assert(m[3, 3] == 10);

  m = scaleMat(10, 10, .5);
  assert(m[0, 0] == 10);
  assert(m[1, 1] == 10);
  assert(m[2, 2] == .5);
  assert(m[3, 3] == 1);

  m = scaleMat(Vec2f(4, 4));
  assert(m[0, 0] == 4);
  assert(m[1, 1] == 4);
  assert(m[2, 2] == 1);
  assert(m[3, 3] == 1);
}

/+
/**
 * Creates a rotation matrix from a axis and angle in radians
 * Params:
 * v = Rotation axis
 * angle = angle in radians
 * Returns a 3d rotation matrix
 */
M rotMat(M=Mat4f, V=Vec3f, T=float) (V v, T angle)
if (isVector!V && V.dim <= 3 && is4dMat!M && is(T : real)) {
  import std.math : approxEqual, sin, cos;
  auto mag = v.sq_length;
  if (mag == 0) {
    return M.IDENTITY;
  } else if (! approxEqual(mag, 1)) {
    v.normalize;
  }

  auto s = sin(angle);
  auto c = cos(angle);

  auto xx = v.x * v.x;
  auto yy = v.y * v.y;
  auto zz = v.z * v.z;
  auto xy = v.x * v.y;
  auto yz = v.y * v.z;
  auto zx = v.z * v.x;
  auto xs = v.x * v.s;
  auto ys = v.y * v.s;
  auto zs = v.z * v.s;
  auto one_c = 1.0 - c;

  M m = M.IDENTITY;
  m[0,0] = (one_c * xx) + c;
  m[0,1] = (one_c * xy) - zs;
  m[0,2] = (one_c * zx) + ys;

  m[1,0] = (one_c * xy) + zs;
  m[1,1] = (one_c * yy) + c;
  m[1,2] = (one_c * yz) - xs;

  m[2,0] = (one_c * zx) - ys;
  m[2,1] = (one_c * yz) + xs;
  m[2,2] = (one_c * zz) + c;

  return m;
}

/**
 * Creates a rotation matrix from a axis and angle in radians
 * Params:
 * x = X coord of rotation axis
 * y = Y coord of rotation axis
 * z = Z coord of rotation axis
 * angle = angle in radians
 * Returns a 3d rotation matrix
 */
M rotMat(M=Mat4f,T=float, U=float) (T x, T y, T z, U angle)
if (is4dMat!M && is(T : real) && is(U :real)) {
  return rotMat!(M,Vector!(3,T),T) (Vector!(3,T)(x,y,z), angle);
}

/**
 * Creates a rotation matrix from a Quaternion
 * Params:
 * q = Quaternion that represents a rotation
 * Returns a 3d rotation matrix
 */
M rotMat(M=Mat4f,Q=Qua_f) (Q q)
if (is4dMat!M && isQuaternion!Q) {
  return cast(M) q;
}

unittest {
  auto m = rotMat(Vec3f.X_AXIS, 0);
  assert (approxEqual(m.determinant , 1) );
  assert (m.equal(Mat4f.IDENTITY));

  auto q = Qua_f(0,0, PI_2);
  m = rotMat(q);
  assert (approxEqual(m.determinant , 1) );
  // dfmt off
  assert (m.equal( Mat4d([1, 0, 0, 0,
                         0, 0,-1, 0,
                         0, 1, 0, 0,
                         0, 0, 0, 1]) ));
  // dfmt on
}
+/
// TODO More usefull functions, etc...
