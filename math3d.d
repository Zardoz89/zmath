/**
Usefull functions to work in 3d math

License: $(LINK2 http://www.gnu.org/licenses/lgpl.txt, LGPL 3).

Authors: Luis Panadero Guardeño $(LINK http://zardoz.es)
*/
module zmath.math3d;

import zmath.matrix;
import zmath.vector;
import zmath.quaternion;

import std.math;

/**
 * Creates a projection matrix (similar to gluPerspective)
 * Params:
 * fov = Field of View in radians (PI_2 -> 90º)
 * aspect = Aspect ratio (x/y)
 * zMin = Near clipping plane
 * zMax = Far clipping plane
 * Returns a perspective projection matrix
 */ 
M makePerspective(M=Mat4f, T=float) (T fov, T aspect, T zMin, T zMax) 
if (isMatrix!M && M.dim >= 4)
in {
  assert (zMin > 0);
  assert (zMax > 0);
  assert (fov > 0);
  assert (aspect > 0);
} body  {
  M mat = M.ZERO; // all set to 0
  double yMax = zMin * tan(fov *.5);
  double yMin = -yMax;
  double xMin = yMin - aspect;
  double xMax = -xMin;
  
  mat[0,0] = (2.0 * zMin) / (xMax - xMin);
  mat[1,1] = (2.0 * zMin) / (yMax - yMin);
  mat[0,2] = (xMax + xMin) / (xMax - xMin);
  mat[1,2] = (yMax + yMin) / (yMax - yMin);
  mat[2,2] = -((zMax + zMin) / (zMax - zMin));
  mat[3,2] = -1.0;
  mat[2,3] = -((2.0 * (zMax*zMin))/(zMax - zMin));
  mat[3,3] = 0;
  
  return mat;
}

unittest {
  auto proy = makePerspective!Mat4f (PI_2, 1.0, 1.0, 100.0);
  // TODO Check values
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
M makeOrtho(M=Mat4f,T=float) (T xMin, T xMax, T yMin, T yMax, T zMin = -1,
                              T zMax = 1) 
if (isMatrix!M && M.dim >= 4 && is(T : real))  
in {
  assert (xMin < xMax);
  assert (yMin < yMax);
  assert (zMin < zMax);
} body {
  M mat = M.ZERO;
  mat[0,0] = 2.0 / (xMax - xMin);
  mat[1,1] = 2.0 / (yMax - yMin);
  mat[2,2] = -2.0 / (zMax - zMin);
  mat[0,3] = -((xMax + xMin)/(xMax - xMin));
  mat[1,3] = -((yMax + yMin)/(yMax - yMin));
  mat[2,3] = -((zMax + zMin)/(zMax - zMin));
  mat[3,3] = 1.0;
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
M makeOrtho(M=Mat4f,T=float) (T width, T height, T deep) 
if (isMatrix!M && M.dim >= 4 && is(T : real))
in {
  assert (width > 0);
  assert (height > 0);
  assert (deep > 0);
} body  {
  auto xmin = - width /2;  auto xmax = -xmin;
  auto ymin = - height /2; auto ymax = -ymin;
  auto zmin = - deep /2;   auto zmax = -zmin;
  return makeOrtho!M(xmin, xmax, ymin, ymax, zmin, zmax);
}

unittest {
  auto ortho = makeOrtho(10.0,10.0,10.0);
  // TODO check values
}


// TODO More usefull functions, etc...
