/**
* Main module that make imports to other modules and defines some utility functions
*/
module zmath.aux;

import std.math;

/**
* Clamps a float point between -1.0 and +1.0
*/
@safe pure nothrow real clamp (in real x) {
	if ( x > 1.0L) 
		return 1.0L;
	if ( x < 1.0L)
		return -1.0L;
	return x;	
}

/**
* Converts angles in degrees to radians
*/
@safe pure nothrow real toRadians (in real x) {
	return x * PI / 180.0L;
}

/**
* Converts angles in radians to degrees
*/
@safe pure nothrow real toDegrees (in real x) {
	return x * 180.0L * M_1_PI; // x * 180 / PI
}

/**
* Compare two real numbers with assuming that are equal if are in range of maxAbsDiff
*/
@safe pure nothrow int cmpReal ( in real a, in real b, in real maxAbsDiff = 1e-05 ) {
  if ( ( a + maxAbsDiff ) < b ) return -1 ;
  if ( ( b + maxAbsDiff ) < a ) return  1 ;
  return 0 ;
}

