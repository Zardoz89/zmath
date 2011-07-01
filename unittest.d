public import zmath.aux;
public import zmath.vector;
public import zmath.matrix;
public import zmath.quaternion;

import std.stdio;

int main (string[] args) {
	
	version (all)	{
		// Bring in unit test for module by referencing function in it
		clamp(3.0f);							// aux
		auto v = Vec3f(1,2,3);		// Vector
		auto m = Mat2f.IDENTITY;	// Matrix
		auto q = Qua_f();					// Quaternion
		
	}
	
	writeln("Success!\n");
  return 0;
}