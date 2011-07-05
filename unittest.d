public import zmath.aux;
public import zmath.vector;
public import zmath.matrix;
public import zmath.quaternion;
public import zmath.math3d;

import std.stdio;

int main (string[] args) {
	
	version (all)	{
		// Bring in unit test for module by referencing function in it
		clamp(3.0f);							// aux
		auto v = Vec3f(1,2,3);		// Vector
		auto m = Mat2f.IDENTITY;	// Matrix
		auto q = Qua_f();					// Quaternion
		auto ortho = orthoMat(100,100,100); //Math3d
	}
	
	writeln("Success!\n");
  return 0;
}