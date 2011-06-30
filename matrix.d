/**
* Module that defines a squared Matrix compatible with OpenGL
*/
module zmath.matrix;

import zmath.aux;
import zmath.vector;

import std.math;
import std.parallelism;
import std.conv;

version (unittest) {
	import std.stdio;
}

unittest {
	writeln("Unit test of Matrix :");
	
	// Fundamental Matrixes
	auto z = Mat4r.ZERO;
	auto ide = Mat4r.IDENTITY;
	//writeln("Zero:");
	//writeln(z); writeln();
	
	foreach (ind; 0..z.cells) {
		assert(z.cell[ind] == 0.0L);
	}
	writeln("Zero Matrix : OK");
	
	//writeln("Identity:");
	//writeln(ide); writeln();
	foreach (i; 0..Mat4r.dim) {
		foreach (j; 0..Mat4r.dim) {
			if (i == j) {
				assert (ide[i,j] == 1.0L);
			} else {
				assert (ide[i,j] == 0.0L);
			}	
		}
	}
	writeln("Identity Matrix : OK");
	
	// Check isMatrix
	assert(isMatrix!(typeof(ide)));
	writeln("isMatrix : OK");
	
	// Check == and !=
	assert(z != ide);
	assert(ide == ide);
	assert(ide.equal(ide));
	assert(! ide.equal(z));
	writeln("Equality Operator : OK");
	
	// Check isOk() and is Finite()
	assert(ide.isOk());
	assert(ide.isFinite());
	Mat3r full_nan;
	assert(!full_nan.isOk());
	assert(!full_nan.isFinite());
	writeln("isOk and isFinite : OK");
	
	// Check Column assignation
	auto tcol = ide;
	auto v = cast(tcol.VCol) Vec3f(1, 2 ,3 );
	tcol[2] = v;
	assert (tcol[2] == v);
	writeln("Column vector assignation : OK");
	
	// Check + - unary operators
	auto n_ide = -ide;
	foreach (i; 0..Mat4r.dim) {
		foreach (j; 0..Mat4r.dim) {
			if (i == j) {
				assert (n_ide[i,j] == -1.0L);
			} else {
				assert (n_ide[i,j] == 0.0L);
			}	
		}
	}
	writeln("Unary Operator : OK");
	
	// Check Addition and subtraction
	auto a = Mat4r([1, 1, 1, 1,
										1, 1, 1, 1,
										1, 1, 1, 1,
										1, 1, 1, 1]);
	auto b = a + ide;
	auto c = z - a;
	
	foreach (i; 0..Mat4r.dim) {
		foreach (j; 0..Mat4r.dim) {
			if (i == j) {
				assert (b[i,j] == 2.0L);
			} else {
				assert (b[i,j] == 1.0L);
			}	
		}
	}
	writeln("Adittion Operator : OK");
	
	foreach (i; 0..Mat4r.dim) {
		foreach (j; 0..Mat4r.dim) {
			assert (c[i,j] == -1.0L);
		}
	}
	writeln("Subtraction Operator : OK");
	
	// Scalar multiplication
	b = b *2;
	foreach (i; 0..Mat4r.dim) {
		foreach (j; 0..Mat4r.dim) {
			if (i == j) {
				assert (b[i,j] == 4.0L);
			} else {
				assert (b[i,j] == 2.0L);
			}	
		}
	}
	writeln("Scalar Multiplication : OK");
	
	// Check Matrix multiplication
	auto axi = a * ide;
	auto ixa = ide * a;
	assert (axi == ixa);
	assert (a == axi);
	writeln("Multiplication - Neutral Element : OK");
	
	auto a3 = Mat3r([1,4,7 , 2,5,8, 3,6,9]);
	auto b3 = Mat3r([9,6,3 , 8,5,2, 7,4,1]);
	auto shouldBeC3 = Mat3r([30,84,138 , 24,69,114, 18,54,90]);
	auto c3 = a3*b3;
	auto d3 = b3*a3;
	assert(c3 == shouldBeC3);
	assert(c3 != d3);
	writeln("Multiplication : OK");
	
	auto vec = ide.VCol(1,2 ,3 ,4);
	auto vecl = vec * ide;
	auto vecr = ide * vec;
	assert (vec == vecl);
	assert (vec == vecr);
	
	writeln("Matrix x Vector : OK");
	
	// Check Matrix transposition
	auto ide_t = ide.transposed();
	assert (ide == ide_t);
	
	auto a3_t = a3.transposed();
	auto should_a3_t = Mat3r([1,2,3 ,4,5,6, 7,8,9]);
	assert (a3_t == should_a3_t);
	writeln("Transposition : OK");
	
	// Check Matrix Determinant
	assert (ide_t.determinant == 1);
	assert (ide.determinant == 1);
	auto n = Mat4r([1, 20, 1, 0,
									1,  5, 1, 1,
									1, -1, 3, 1,
									1, 10, 1, 1]);
	auto d_n = n.determinant;
	assert (d_n == -10 );
	auto d_b = b.determinant;
	auto prod = n *b;
	auto d_prod = prod.determinant;
	assert(d_n* d_b == d_prod);
	assert( a3_t.determinant == 0);
	writeln("Determinant : OK");
	
	// Check matrix inversion
	auto n_inv = n.inverse;
	assert(n_inv == Mat4r([   1,  5.1, -0.5, -4.6,
												 		0, -0.2,    0,  0.2,
												 		0, -1.1,  0.5,  0.6,
												 	 -1,   -2,    0,    3]) );
	auto a3_inv = a3.inverse;
	assert(!a3_inv.isOk() ); // not have inverse, so returns a matrix full os NaNs
	auto m2 = Mat2r([10, 5,  5, 2]);
	auto m2_inv = m2.inverse;
	assert (m2_inv == Mat2r([-.4, 1,  1, -2]) );
	writeln("Inversion : OK");
	
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
	writeln("Matrix Casting : OK");
	
	writeln();
}

alias Matrix!(float,2) 	Mat2f;		/// 2D squared matrix of floats
alias Matrix!(float,3) 	Mat3f;		/// 3D squared matrix of floats
alias Matrix!(float,4) 	Mat4f;		///	4D squared matrix of floats

alias Matrix!(double,2) Mat2d;		/// 2D squared matrix of doubles
alias Matrix!(double,3) Mat3d;		/// 3D squared matrix of doubles
alias Matrix!(double,4) Mat4d;		///	4D squared matrix of doubles

alias Matrix!(real,2) 	Mat2r;		/// 2D squared matrix of reals
alias Matrix!(real,3) 	Mat3r;		/// 3D squared matrix of reals
alias Matrix!(real,4) 	Mat4r;		///	4D squared matrix of reals

/**
* Defines a squared Matrix of n = 2, 3 or 4 size, like a linear array of numbers
*/
struct Matrix(T, size_t dim_)
if (is(T == real) || is(T == double) || is(T == float)  
		&& (( dim_ >= 2 ) && ( dim_ <= 4 ) ))   
{ 
	static enum size_t dim = dim_; 			/// Matrix Dimension
	static enum size_t cells = dim*dim; 	/// Matrix number of cells
	
	alias Vector!(T,dim_) VCol;
	
	union {
		private T cell[cells]; 					/// Matrix like of a array of cells of major-column
		private VCol col[dim_];					/// Matrix like of a array of column vectors
	}
	
	// Consts
	static if (dim == 2) { // 2x2
		public static enum Matrix!(T,2) ZERO = {[0, 0, 0 ,0]};	 
		public static enum Matrix!(T,2) IDENTITY = {[1, 0, 0 ,1]};	 
	}
	static if (dim == 3) { // 3x3
		public static enum Matrix!(T,3) ZERO = {[0, 0, 0,  0 ,0, 0,  0 ,0, 0]};	 
		public static enum Matrix!(T,3) IDENTITY = {[1, 0, 0,  0, 1, 0,  0, 0, 1]};	 
	}
	static if (dim == 4) { // 4x4
		public static enum Matrix!(T,4) ZERO = {[0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0]};	 
		public static enum Matrix!(T,4) IDENTITY = {[1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]};	 
	}

	// Basic properties *****************************************************
	
	/**
	* Returns i, j cell
	*/
	T opIndex(size_t row, size_t cl) const { return col[cl][row];}
	
	/**
	* Assigns a new cell value
	*/
	void opIndexAssign(T c, size_t row, size_t cl) {
		col[cl][row] = c;
	}
	
	/**
	* Returns j column vector
	*/
	ref VCol opIndex(size_t j) const { return col[j];}
	
	/**
	* Assigns a new column vector
	*/
	void opIndexAssign(VCol v, size_t j) {
		//static assert (v.dim <= dim, "Vector have a bigger dimension that this maxtrix.");
		//static if (v.dim != dim) {
		//	col[j] = cast(VCol) v;
		//} else {
			col[j] = v;
		//}
	}
	
	/**
	* Calcs offset to locate a particular a(i,j), where i -> row ; j -> col
	*/
	@safe private pure nothrow size_t offset (size_t j) const {
		return dim * j;
	}
	
	/**
	* Calcs position in cell array to locate a particular a(i,j), where i -> row ; j -> col
	*/
	@safe private pure nothrow size_t pos(size_t i, size_t j) const {return i + offset(j);}
	
	// Operations ***************************************	
	
	/**
	* Define Equality 
	*/
	bool opEquals(ref const Matrix rhs) const {
		if (col[0] != rhs.col[0] || col[1] != rhs.col[1]) return false;
		static if (dim >= 3)
			if (col[2] != rhs.col[2]) return false;
		static if (dim >= 4)
			if (col[3] != rhs.col[3]) return false;
		return true;	
	}
	
	
	/**
	* Approximated equality
	*/
	bool equal(ref const Matrix rhs) const {
		if (! col[0].equal(rhs.col[0]) || ! col[1].equal(rhs.col[1]) ) return false;
		static if (dim >= 3)
			if (! col[2].equal(rhs.col[2])) return false;
		static if (dim >= 4)
			if (! col[3].equal(rhs.col[3])) return false;
		return true;
	}

	/**
	* Approximated equality with controlable precision
	*/
	bool equal(ref const Matrix rhs, T maxDiff) const {
		if (! col[0].equal(rhs.col[0], maxDiff) || ! col[1].equal(rhs.col[1], maxDiff) ) return false;
		static if (dim >= 3)
			if (! col[2].equal(rhs.col[2], maxDiff)) return false;
		static if (dim >= 4)
			if (! col[3].equal(rhs.col[3], maxDiff)) return false;
		return true;
	}
	
	/**
	* Approximated equality
	*/
	bool equal(ref const Matrix rhs, T maxRelDiff, T maxAbsDiff = 1e-05) const {
		if (! col[0].equal(rhs.col[0], maxRelDiff, maxAbsDiff) || ! col[1].equal(rhs.col[1], maxRelDiff, maxAbsDiff) ) return false;
		static if (dim >= 3)
			if (! col[2].equal(rhs.col[2], maxRelDiff, maxAbsDiff)) return false;
		static if (dim >= 4)
			if (! col[3].equal(rhs.col[3], maxRelDiff, maxAbsDiff)) return false;
		return true;
	}
	
	/**
	* Define unary operators + and -
	*/
	Matrix opUnary(string op) () const
		if (op == "+" || op == "-")
	{
		string makeMix(size_t d, string op) {
			string ret = "return Matrix( [";
			foreach (i; 0..(d*d)) {
				ret ~= op ~ "cell["~to!string(i)~"] ";
				if (i != (d*d -1))
					ret ~= ",";
			}
			ret ~= "] );";
			return ret;
		}

		static if (dim == 2) 
			mixin(makeMix(2, op));
		static if (dim == 3) 
			mixin(makeMix(3, op));
		static if (dim == 4) 
			mixin(makeMix(4, op));
	}
	
	/**
	* Define binary operator + and -
	*/
	Matrix opBinary(string op) (ref const Matrix rhs) const
		if (op == "+" || op == "-")
	{
		string makeMixOp(size_t d, string op) {
			string ret = "return Matrix( [";
			foreach (i; 0..(d*d)) {
				ret ~= "cell["~to!string(i)~"] " ~op ~ "rhs.cell["~to!string(i)~"] ";
				if (i != (d*d -1))
					ret ~= ",";
			}
			ret ~= "] );";
			return ret;
		}
		
		static if (dim == 2) 
			mixin(makeMixOp(2, op));
		static if (dim == 3) 
			mixin(makeMixOp(3, op));
		static if (dim == 4) 
			mixin(makeMixOp(4, op));
	}
	
	/**
	* Define Scalar multiplication
	*/
	Matrix opBinary(string op) (in real rhs) const
		if (op == "*" )
	{
		T[this.cells] ret = cell[] * rhs;
		return Matrix(ret);
	}
	
	
	/**
	* Define Matrix Product
	*/
	Matrix opBinary(string op) (ref const Matrix rhs) const
		if (op == "*" )
	{
		Matrix mat;/+ = Matrix.ZERO;
		
		// naive method
		foreach (size_t j; 0..dim) { // Runs along result columns
			foreach (size_t i; 0..dim) { // Runs along result rows
				foreach (size_t u; 0..dim) { // Calcs result cell matrix value
					// value += Aiu * Buj
					mat[i,j] = mat[i,j] + this[i,u] * rhs[u,j] ; 
				}
			}
		}
		+/
		
		// direct method
		foreach (size_t j; 0..dim) { // Runs along result columns
			foreach (size_t i; 0..dim) { // Runs along result rows
				static if (dim == 2) {
					mat[i, j] = this[i,0] * rhs[0,j] + this[i,1] * rhs[1,j] ;
				} else if (dim == 3) {
					mat[i, j] = this[i,0] * rhs[0,j] + this[i,1] * rhs[1,j] + this[i,2] * rhs[2,j] ;
				} else if (dim == 4) {
					mat[i, j] = this[i,0] * rhs[0,j] + this[i,1] * rhs[1,j] + this[i,2] * rhs[2,j] + this[i,3] * rhs[3,j] ;
				}
			}
		}
		
		return mat;
	}
	
	/**
	* Define Matrix x Vector
	*/
	VCol opBinary(string op) (ref const VCol rhs) const
		if (op == "*" )
	{
		static assert (rhs.dim == dim, "The vector dimension must be smae that Matrix dimmension");
		VCol ret;
		foreach (size_t j; 0..dim) { // Runs along result coords	
			static if (dim == 2) {
				ret[j] = this[j,0] * rhs[0] + this[j,1] * rhs[1] ;
			} else if (dim == 3) {
				ret[j] = this[j,0] * rhs[0] + this[j,1] * rhs[1] + this[j,2] * rhs[2] ;
			} else if (dim == 4) {
				ret[j] = this[j,0] * rhs[0] + this[j,1] * rhs[1] + this[j,2] * rhs[2] + this[j,3] * rhs[3] ;
			}
		}
		return ret;
	}
	
	/**
	* Define Vector x Matrix
	*/
	VCol opBinaryRight(string op) (ref const VCol lhs) const
		if (op == "*" )
	{
		static assert (lhs.dim == dim, "The vector dimension must be smae that Matrix dimmension");
		VCol ret;
		foreach (size_t j; 0..dim) { // Runs along result coords	
			static if (dim == 2) {
				ret[j] = this[0,j] * lhs[0] + this[1,j] * lhs[1] ;
			} else if (dim == 3) {
				ret[j] = this[0,j] * lhs[0] + this[1,j] * lhs[1] + this[2,j] * lhs[2] ;
			} else if (dim == 4) {
				ret[j] = this[0,j] * lhs[0] + this[1,j] * lhs[1] + this[2,j] * lhs[2] + this[3,j] * lhs[3] ;
			}
		}
		return ret;
	}
	
	/**
	* Returns transposed matrix
	*/
	Matrix transposed () const{
		Matrix mat;
		foreach (i; 0..dim) { // Runs along result rows 
			foreach (j; 0..dim) { // Runs along result columns	
				mat[j,i] = this.cell[pos(i,j)];
			}
		}
		return mat;
	}
	
	/**
	* Returns Determinant of this matrix
	*/
	pure T determinant () const {
		
		static if (dim == 2)
			return cell[pos(0,0)] * cell[pos(1,1)] - (cell[pos(1,0)] * cell[pos(0,1)]) ;
		else if (dim == 3) {
			T aei = cell[pos(0,0)] * cell[pos(1,1)] * cell[pos(2,2)];
			T bfg = cell[pos(0,1)] * cell[pos(1,2)] * cell[pos(2,0)];
			T cdh = cell[pos(0,2)] * cell[pos(1,0)] * cell[pos(2,1)];
			T afh = cell[pos(0,0)] * cell[pos(1,2)] * cell[pos(2,1)];
			T bdi = cell[pos(0,1)] * cell[pos(1,0)] * cell[pos(2,2)];
			T ceg = cell[pos(0,2)] * cell[pos(1,1)] * cell[pos(2,0)];
			return aei + bfg + cdh -afh - bdi -ceg;
		}else {
			return 	 (cell[pos(0,0)] * cell[pos(1,1)] - cell[pos(0,1)] * cell[pos(1,0)]) * (cell[pos(2,2)] * cell[pos(3,3)] - cell[pos(2,3)] * cell[pos(3,2)])
             - (cell[pos(0,0)] * cell[pos(1,2)] - cell[pos(0,2)] * cell[pos(1,0)]) * (cell[pos(2,1)] * cell[pos(3,3)] - cell[pos(2,3)] * cell[pos(3,1)])
             + (cell[pos(0,0)] * cell[pos(1,3)] - cell[pos(0,3)] * cell[pos(1,0)]) * (cell[pos(2,1)] * cell[pos(3,2)] - cell[pos(2,2)] * cell[pos(3,1)])
             + (cell[pos(0,1)] * cell[pos(1,2)] - cell[pos(0,2)] * cell[pos(1,1)]) * (cell[pos(2,0)] * cell[pos(3,3)] - cell[pos(2,3)] * cell[pos(3,0)])
             - (cell[pos(0,1)] * cell[pos(1,3)] - cell[pos(0,3)] * cell[pos(1,1)]) * (cell[pos(2,0)] * cell[pos(3,2)] - cell[pos(2,2)] * cell[pos(3,0)])
             + (cell[pos(0,2)] * cell[pos(1,3)] - cell[pos(0,3)] * cell[pos(1,2)]) * (cell[pos(2,0)] * cell[pos(3,1)] - cell[pos(2,1)] * cell[pos(3,0)]);
		}
	}
	
	/**
	* Calcs this maxtrix inverse
	* Returns: A !isOk() matrix full of NaNs if this matrix isn't inverible
	*/
	Matrix inverse() const {
		Matrix mat;
		auto det = this.determinant();
		if (approxEqual(det, 0)) // Not invertible matrix
			return mat; // At this point Mat it's full of NaNs and a not valid Matrix
			
		static if (dim == 2) {
			mat[0,0] = this[1,1];
			mat[1,1] = this[0,0];
			mat[0,1] = - this[0,1];
			mat[1,0] = - this[1,0];
		} else if (dim == 3) {
			mat[0,0] = this[1,1] * this[2,2] - this[1,2] * this[2,1];
			mat[0,1] = this[0,2] * this[2,1] - this[0,1] * this[2,2];
			mat[0,2] = this[0,1] * this[1,2] - this[0,2] * this[1,1];
			
			mat[1,0] = this[1,2] * this[2,0] - this[1,0] * this[2,2];
			mat[1,1] = this[0,0] * this[2,2] - this[0,2] * this[2,0];
			mat[1,2] = this[0,2] * this[1,0] - this[0,0] * this[1,2];
			
			mat[2,0] = this[1,0] * this[2,1] - this[1,1] * this[2,2];
			mat[2,1] = this[0,1] * this[2,0] - this[0,0] * this[2,1];
			mat[2,2] = this[0,0] * this[1,1] - this[0,1] * this[1,0];
		} else { // dimm = 4
			mat[0,0] = this[1, 1] * (this[2, 2] * this[3, 3] - this[2, 3] * this[3, 2] ) + this[1, 2] *( this[2, 3] * this[3, 1] - this[2, 1] * this[3, 3] ) + this[1, 3] *( this[2, 1] * this[3, 2] - this[2, 2] * this[3, 1]);
      mat[0,1] = this[2, 1] * (this[0, 2] * this[3, 3] - this[0, 3] * this[3, 2] ) + this[2, 2] *( this[0, 3] * this[3, 1] - this[0, 1] * this[3, 3] ) + this[2, 3] *( this[0, 1] * this[3, 2] - this[0, 2] * this[3, 1]);          
			mat[0,2] = this[3, 1] * (this[0, 2] * this[1, 3] - this[0, 3] * this[1, 2] ) + this[3, 2] *( this[0, 3] * this[1, 1] - this[0, 1] * this[1, 3] ) + this[3, 3] *( this[0, 1] * this[1, 2] - this[0, 2] * this[1, 1]);
      mat[0,3] = this[0, 1] * (this[1, 3] * this[2, 2] - this[1, 2] * this[2, 3] ) + this[0, 2] *( this[1, 1] * this[2, 3] - this[1, 3] * this[2, 1] ) + this[0, 3] *( this[1, 2] * this[2, 1] - this[1, 1] * this[2, 2]);
      
      mat[1,0] = this[1, 2] * (this[2, 0] * this[3, 3] - this[2, 3] * this[3, 0] ) + this[1, 3] * (this[2, 2] * this[3, 0] - this[2, 0] * this[3, 2] ) + this[1, 0] * (this[2, 3] * this[3, 2] - this[2, 2] * this[3, 3]);
      mat[1,1] = this[2, 2] * (this[0, 0] * this[3, 3] - this[0, 3] * this[3, 0] ) + this[2, 3] * (this[0, 2] * this[3, 0] - this[0, 0] * this[3, 2] ) + this[2, 0] * (this[0, 3] * this[3, 2] - this[0, 2] * this[3, 3]);
      mat[1,2] = this[3, 2] * (this[0, 0] * this[1, 3] - this[0, 3] * this[1, 0] ) + this[3, 3] * (this[0, 2] * this[1, 0] - this[0, 0] * this[1, 2] ) + this[3, 0] * (this[0, 3] * this[1, 2] - this[0, 2] * this[1, 3]);
      mat[1,3] = this[0, 2] * (this[1, 3] * this[2, 0] - this[1, 0] * this[2, 3] ) + this[0, 3] * (this[1, 0] * this[2, 2] - this[1, 2] * this[2, 0] ) + this[0, 0] * (this[1, 2] * this[2, 3] - this[1, 3] * this[2, 2]);
      
      mat[2,0] = this[1, 3] * (this[2, 0] * this[3, 1] - this[2, 1] * this[3, 0] ) + this[1, 0] * (this[2, 1] * this[3, 3] - this[2, 3] * this[3, 1] ) + this[1, 1] * (this[2, 3] * this[3, 0] - this[2, 0] * this[3, 3]);
      mat[2,1] = this[2, 3] * (this[0, 0] * this[3, 1] - this[0, 1] * this[3, 0] ) + this[2, 0] * (this[0, 1] * this[3, 3] - this[0, 3] * this[3, 1] ) + this[2, 1] * (this[0, 3] * this[3, 0] - this[0, 0] * this[3, 3]);
      mat[2,2] = this[3, 3] * (this[0, 0] * this[1, 1] - this[0, 1] * this[1, 0] ) + this[3, 0] * (this[0, 1] * this[1, 3] - this[0, 3] * this[1, 1] ) + this[3, 1] * (this[0, 3] * this[1, 0] - this[0, 0] * this[1, 3]);
      mat[2,3] = this[0, 3] * (this[1, 1] * this[2, 0] - this[1, 0] * this[2, 1] ) + this[0, 0] * (this[1, 3] * this[2, 1] - this[1, 1] * this[2, 3] ) + this[0, 1] * (this[1, 0] * this[2, 3] - this[1, 3] * this[2, 0]);
      
      mat[3,0] = this[1, 0] * (this[2, 2] * this[3, 1] - this[2, 1] * this[3, 2] ) + this[1, 1] * (this[2, 0] * this[3, 2] - this[2, 2] * this[3, 0] ) + this[1, 2] * (this[2, 1] * this[3, 0] - this[2, 0] * this[3, 1]);
			mat[3,1] = this[2, 0] * (this[0, 2] * this[3, 1] - this[0, 1] * this[3, 2] ) + this[2, 1] * (this[0, 0] * this[3, 2] - this[0, 2] * this[3, 0] ) + this[2, 2] * (this[0, 1] * this[3, 0] - this[0, 0] * this[3, 1]);
			mat[3,2] = this[3, 0] * (this[0, 2] * this[1, 1] - this[0, 1] * this[1, 2] ) + this[3, 1] * (this[0, 0] * this[1, 2] - this[0, 2] * this[1, 0] ) + this[3, 2] * (this[0, 1] * this[1, 0] - this[0, 0] * this[1, 1]);
			mat[3,3] = this[0, 0] * (this[1, 1] * this[2, 2] - this[1, 2] * this[2, 1] ) + this[0, 1] * (this[1, 2] * this[2, 0] - this[1, 0] * this[2, 2] ) + this[0, 2] * (this[1, 0] * this[2, 1] - this[1, 1] * this[2, 0]);
 		}
		
		det = 1 / det;
		return mat * det;
	}
	
	// Misc *********************************************************************************
	
	/**
	* Checks that the matrix not have a weird NaN value
	*/
	bool isOk() {
		foreach (c ; col) {
			if (!c.isOk()) return false;
		}
		return true;
	}
	
	/**
	* Checks that the matrix have finite values
	*/
	bool isFinite() {
		foreach (c ; col) {
			if (! c.isFinite()) return false;
		}
		return true;
	}
	
	/**
	* to!type method to convert to other vector types
	*/
	Tout opCast( Tout ) () 
	if (isMatrix!(Tout) && Tout.dim >= dim) {
		static assert (isMatrix!(Tout), "This type not is a Matrix");
		static assert (Tout.dim >= dim, "Original Matrix bigger that destiny Vector");
		
		Tout newMat; auto i = 0;
		static if (is (typeof(newMat[0, 0]) == typeof(this[0, 0])))  {
			for (; i < dim; i++)
				newMat[i] =  this[i];
		} else {
			for (; i < dim; i++)	
				newMat[i] =  cast(newMat.VCol) this[i]; //Vector Casting auto expands rows
		}
			
		for (; i< Tout.dim; i++) // Expands a small matrix to a bigger dimension with a full 0's columns 
			newMat[i] = cast(newMat.VCol) newMat.VCol.ZERO;
		
		return newMat;
	}
	
	/**
	* Returns a visual representation of this matrix
	*/
	string toString() const {
		string ret; // I -> row, j->col
		foreach (i; 0..dim) {
			ret ~= "|"; 
			foreach (j; 0..dim) { 
				ret ~= to!string( cell[pos(i,j)] );
				if (j < (dim-1))
					ret ~= ", ";
			}
			ret ~= "|"; 
			if (i < (dim-1))
				ret ~= "\n"; 
		}
		return ret;
	}
}

/**
* Say if a thing it's a Matrix
*/
template isMatrix(T)	{
	//immutable bool isMatrix = is(T == Matrix);
  immutable bool isMatrix = __traits(compiles,
        (){  
            T t;
            static assert(T.dim >= 2 && T.dim <= 4);
            auto cell = t.cell;
            auto cell00 = t[0,0];
            auto col0 = t[0];
            auto coln = t[T.dim];
            
            // TODO : Should test for methods ?
        }
    );
}