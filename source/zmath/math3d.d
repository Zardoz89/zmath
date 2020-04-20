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
  immutable bool is4dMat = (isMatrix!M && __traits(compiles,
        (){
            M t;
            static assert(M.dim == 4);
        }
    ));
}

unittest {
  assert(is4dMat!Mat4f);
  assert(! is4dMat!float);
  assert(is4dMat!Mat4d);
  assert(! is4dMat!Mat2f);
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
@nogc M perspectiveMat(M=Mat4f, T=float, U=float, V=float, W=float)
                          (T fov, U aspect, V zMin, W zMax) pure nothrow
if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V && isNumeric!W)
in {
  import std.math : PI;

  assert (zMin > 0, "Zmin equal or less that zero");
  assert (zMax > 0, "Zmax equal or less that zero");
  assert (fov > 0 && fov <= 2 * PI, "Not valid FOV");
  assert (aspect > 0, "Not valid aspect ratio");
} body  {
  import std.math : tan;

  M mat = M.ZERO; // all set to 0

  auto f = 1.0f / tan(fov / 2.0f);
  auto d = 1.0f / (zMin - zMax);

  mat[0,0] = f / aspect;
  mat[1,1] = f;
  mat[2,2] = (zMax + zMin) * d;
  mat[2,3] = 2.0f * d * zMax * zMin;
  mat[3,2] = -1;

  return mat;
}

unittest {
  const auto proy = perspectiveMat(PI_4, 800.0 / 600.0, 1.0, 100.0);
  // dfmt off
  const auto expected = Mat4f([
     1.81066,       0,       0,       0,
           0, 2.41421,       0,       0,
           0,       0, -1.0202, -2.0202,
           0,       0,      -1,       0
  ]);
  // dfmt on
  assert (proy.approxEqual(expected), "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ proy.toString() );
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
@nogc M orthoMat(M=Mat4f,T=float, U=float, V=float, W=float, X=float, Y=float)
                    (T xMin, U xMax, V yMin, W yMax, X zMin = 0, Y zMax = 1) pure nothrow
if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V && isNumeric!W
    && isNumeric!X && isNumeric!Y)
in {
  assert (xMin < xMax, "Invalid values for xMin and xMax");
  assert (yMin < yMax, "Invalid values for yMin and yMax");
  assert (zMin != zMax, "Invalid values for zMin and zMax");
} body {
  M mat = M.ZERO;

  const widthRatio = 1.0L / (xMax - xMin);
  const heightRatio = 1.0L / (yMax- yMin);
  const depthRatio = 1.0L / (zMax - zMin);

  const sx = 2.0 * widthRatio;
  const sy = 2.0 * heightRatio;
  const sz = 2.0 * depthRatio;

  const tx = -(xMax + xMin) * widthRatio;
  const ty = -(yMax + yMin) * heightRatio;
  const tz = -(zMax + zMin) * depthRatio;

  mat[0,0] = sx;
  mat[1,1] = sy;
  mat[2,2] = sz;
  mat[3,3] = 1;

  mat[0,3] = tx;
  mat[1,3] = ty;
  mat[2,3] = tz;
  return mat;
}

/**
 * Creates a orthographic projection matrix
 * Params:
 * width = Width of visible zone
 * height = Height of visible zone
 * deep = Deep of visible zone
 * Returns a orthographic projection matrix
 */
@nogc M orthoMat(M=Mat4f,T=float, U=float, V=float)
                  (T width =1, U height = 1, V deep = 1) pure nothrow
if (is4dMat!M && isNumeric!T && isNumeric!U && isNumeric!V)
in {
  assert (width > 0);
  assert (height > 0);
  assert (deep > 0);
} body  {
  auto x = width / 2.0L; // Max precision
  auto y = height / 2.0L;
  auto z = deep / 2.0L;
  return orthoMat!M(-x, x, -y, y, z, -z);
}


unittest {
  auto ortho = orthoMat(-2, 2, -2, 2, 0, 10);
  // dfmt off
  auto expected = Mat4f([
      0.5,   0,    0,  0,
        0, 0.5,    0, -0,
        0,   0,  0.2, -1,
        0,   0,    0,  1
  ]);
  // dfmt on
  assert (ortho.approxEqual(expected), "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ ortho.toString() );

  ortho = orthoMat(-10, 10, -10, 10, -1, 100f);
  // dfmt off
  expected = Mat4f([
      0.1,   0,        0,         0,
        0, 0.1,        0,         0,
        0,   0, 0.019802, -0.980198,
        0,   0,        0,         1
  ]);
  // dfmt on
  assert (ortho.approxEqual(expected), "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ ortho.toString() );

  ortho = orthoMat(2, 2, 10f);
  // dfmt off
  expected = Mat4f([
      1, 0,    0,  0,
      0, 1,    0, -0,
      0, 0,  0.2, -0,
      0, 0,    0,  1
  ]);
  // dfmt on
  assert (ortho.approxEqual(expected), "Expected\n" ~ expected.toString() ~ " but obtained\n" ~ ortho.toString() );
}

/**
 * Creates a translation matrix from a 2d/3d Vector
 * Params:
 * v = 2d/3d Vector
 * Returns a translation matrix
 */
@nogc M translateMat(M=Mat4f, V=Vec3f) (V v) pure nothrow
if (isVector!V && V.dim <= 3 && is4dMat!M) {
  M m = M.IDENTITY;
  m[0, 3] = v.x;
  m[1, 3] = v.y;
  m[2, 3] = v.z;
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
@nogc M translateMat(M=Mat4f, T=float) (T x, T y, T z) pure nothrow
if (is4dMat!M && is(T : real)) {
  return translateMat!(M, Vector!(T, 3))(Vector!(T, 3) (x,y,z));
}

unittest {
  auto m1 = translateMat(1f, 2f, 3f);
  auto m2 = translateMat(Vec3f(1, 2, 3));
  assert(m1 == m2);
  assert(m2 == Mat4f([
        1, 0, 0, 1,
        0, 1, 0, 2,
        0, 0, 1, 3,
        0, 0, 0, 1,
  ]), "Failed with\n" ~ m2.toString);
}

/**
 * Creates a uniform 3d scale matrix
 * Params:
 * s = scale factor
 * Returns a scale matrix
 */
@nogc M scaleMat(M=Mat4f, T=float) (T s) pure nothrow
if (is4dMat!M && is(T : real)) {
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
@nogc M scaleMat(M=Mat4f, T=float) (T x, T y, T z) pure nothrow
if (is4dMat!M && is(T : real)) {
  M m = M.IDENTITY;
  m[0,0] = x;
  m[1,1] = y;
  m[2,2] = z;
  return m;
}

/**
 * Creates a not uniform 3d scale matrix
 * Params:
 * v = Vector with scale factor. If is a 2d Vector, z axis scale factor is 1
 * Returns a scale matrix
 */
@nogc M scaleMat(M=Mat4f, V=Vec3f) (V v) pure nothrow
if (isVector!V && V.dim <= 3 && is4dMat!M) {
  M m = M.IDENTITY;
  m[0,0] = v.x;
  m[1,1] = v.y;
  static if (V.dim == 2) {
    m[2,2] = 1;
  } else {
    m[2,2] = v.z;
  }
  return m;
}

unittest {
  import std.conv : to;
  auto m = scaleMat(10);
  assert (m[0,0] == 10);
  assert (m[1,1] == 10);
  assert (m[2,2] == 10);
  assert (m[3,3] == 10);

  m = scaleMat(10, 10, .5);
  assert (m[0,0] == 10);
  assert (m[1,1] == 10);
  assert (m[2,2] == .5);
  assert (m[3,3] == 1);

  m = scaleMat(Vec2f(4,4));
  assert (m[0,0] == 4);
  assert (m[1,1] == 4);
  assert (m[2,2] == 1);
  assert (m[3,3] == 1);
}

/**
 * Creates a rotation matrix from a axis and angle in radians
 * Params:
 * v = Rotation axis
 * angle = angle in radians
 * Returns a 3d rotation matrix
 */
@nogc M rotMat(M=Mat4f, V=Vec3f, T=float) (V v, T angle) pure nothrow
    if (isVector!V && V.dim <= 3 && is4dMat!M && __traits(isFloating, T)) {
  import std.math : approxEqual, sin, cos;

  auto mag = v.sq_length;
  if (mag == 0) {
    return M.IDENTITY;
  } else if (! approxEqual(mag, 1)) {
    v.normalize;
  }

  const s = sin(angle);
  const c = cos(angle);
  const oneMinusC = 1.0 - c;

  const xx = v.x * v.x;
  const yy = v.y * v.y;
  const zz = v.z * v.z;

  const xy = v.x * v.y;
  const yz = v.y * v.z;
  const zx = v.z * v.x;

  const xs = v.x * s;
  const ys = v.y * s;
  const zs = v.z * s;

  M m = M.IDENTITY;
  m[0,0] = (oneMinusC * xx) + c;
  m[0,1] = (oneMinusC * xy) - zs;
  m[0,2] = (oneMinusC * zx) + ys;

  m[1,0] = (oneMinusC * xy) + zs;
  m[1,1] = (oneMinusC * yy) + c;
  m[1,2] = (oneMinusC * yz) - xs;

  m[2,0] = (oneMinusC * zx) - ys;
  m[2,1] = (oneMinusC * yz) + xs;
  m[2,2] = (oneMinusC * zz) + c;

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
@nogc M rotMat(M=Mat4f,T=float, U=float) (T x, T y, T z, U angle) pure nothrow
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
  import std.conv : to;

  auto m = rotMat(Vec3f.X_AXIS, 0f);
  assert (approxEqual(m.determinant , 1), "Expected rotation Matrix determinant to be 1, but obtained " ~
      to!string(m.determinant) ~ "\n" ~ m.toString());
  assert (m.approxEqual(Mat4f.IDENTITY));

  auto q = Qua_f(0, 0, PI_2);
  m = rotMat(q);
  assert (approxEqual(m.determinant , 1) );
  // dfmt off
  assert (m.approxEqual( Mat4f([
          0, -1, 0, 0,
          1,  0, 0, 0,
          0,  0, 1, 0,
          0,  0, 0, 1
  ])));
  // dfmt on
}

/**
 * Look At projection
 * Based on gfm code
 */
@nogc M lookAt(M=Mat4f, V=Vec3f, T=float) (V eye, V target, V up) pure nothrow
if (is4dMat!M && isNumeric!T && isVector!V && V.dim <= 3)
{
    auto z = (eye - target);
    z.normalize();
    auto x = (-up & z);
    x.normalize();
    auto y = z & -x; // Cross product

    return M([-x.x,        -x.y,        -x.z,      x * eye,
               y.x,         y.y,         y.z,   -(y * eye),
               z.x,         z.y,         z.z,   -(z * eye),
                0f,          0f,          0f,           1f]);
}

// TODO More usefull functions, etc...
