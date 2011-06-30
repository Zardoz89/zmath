/**
Some auxiliar math funtions

License: $(LINK2 http://www.gnu.org/licenses/lgpl.txt, LGPL 3).

Authors: Luis Panadero Guardeño $(LINK http://zardoz.es)
*/
module zmath.aux;

import std.math;

version (unittest) {
	import std.stdio;
}

enum M_1_180 = 1 / 180.0L; /// Inverse of 180

/**
* Clamps a float point between -1.0 and +1.0
* Params:
*	x = Float point number to clamp 
*/
@safe pure nothrow T clamp (T) (in T x) 
if (__traits(isFloating, T)) {
	if ( x > 1.0L) 
		return 1.0L;
	if ( x < 1.0L)
		return -1.0L;
	return x;	
}

unittest {
	writeln("Unit test of Aux :");
	double d = 25;
	assert(clamp(d) == 1.0);
	real r = -25;
	assert(clamp(r) == -1.0);
	
	writeln("Clamp float : OK");
}

/**
* Converts angles in degrees to radians
* Params:
* x = Angle in grades
*/
@safe pure nothrow T toRadians (T) (in T x) 
if (__traits(isScalar, T)) {
	return x * PI * M_1_180;
}

unittest {
	assert(approxEqual(toRadians(180.0), PI));
	assert(approxEqual(toRadians(360.0), PI * 2.0));
	
	writeln("toRadians : OK");
}

/**
* Converts angles in radians to degrees
* Params:
* x = Angle in radians
*/
@safe pure nothrow T toDegrees (T) (in T x) 
if (__traits(isScalar, T)) {
	return x * 180.0L * M_1_PI; // x * 180 / PI
}

unittest {
	assert(approxEqual(toDegrees(PI),180.0));
	
	writeln("toDegrees : OK");
}

/**
* Compare two float point numbers with assuming that are equal if are in range of maxAbsDiff
* Params:
* a = A float point number
* b = Other float point
* maxRelDiff = Max relative difference 
* maxAbsDiff = Max absoulte difference
* Returns:
* If _a are aproximated equal that _b, returns 0. Otherwise, if _a > _b, returns 1 and if _a < _b , returns -1;
*/
int cmpFloat (T, U) ( in T a, in U b, T maxRelDiff = 1e-2, T maxAbsDiff = 1e-5 ) {
	static assert(__traits(isFloating, T) && is(U : real) );
	
	if (approxEqual(a, b, maxRelDiff, maxAbsDiff))
		return 0;
  if ( a < b ) return -1 ;
  return  1 ;
}

unittest {
	assert(cmpFloat(3.0, 0.0) == 1);
	assert(cmpFloat(0.0, 3) == -1);
	assert(cmpFloat(0.0, 0.00001) == 0);
	
	writeln("cmpFloat : OK");
	writeln();
}
