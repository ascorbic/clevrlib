/* AS3
	Copyright 2007 Sphex LLP
	Ported from Jama: http://math.nist.gov/javanumerics/jama/
	
	This work is licensed under the Creative Commons Attribution 2.0 UK: England & Wales License. 
	To view a copy of this license, visit http://creativecommons.org/licenses/by/2.0/uk/ or send a 
	letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

*/

package com.clevr.matrixalgebra {

	/*
		*	Matrix class ported from Jama.
		*
		*	@langversion ActionScript 3.0
		*	@playerversion Flash 9.0
		*
		*	@author Matt Kane
		*	@since  04.06.2007
		*/
	public class RealMatrix {

			/** Array for internal storage of elements.
			*/
			private var A:Array;

		/** Row and column dimensions.
			*/
			private var m:int, n:int;

		/* ------------------------
			Constructors
			* ------------------------ */

			/** Construct an m-by-n matrix of zeros. 
			@param m    Number of rows.
			@param n    Number of colums.
			*/

		public function RealMatrix (mm:*, nn:int = 0, s:Number = NaN) {
			if (mm is Array) {
				m = mm.length;
				n = mm[0].length;
				for (i = 0; i < m; i++) {
				   if (mm[i].length != n) {
				      throw new Error("All rows must have the same length.");
				   }
				}
				this.A = mm;
		      	return;
			}
			m = mm;
			n = nn;
			A = new Array(m);
			var i:int;
			var j:int;
			for (i = 0; i < m; i++) {
				A[i] = new Array(n);
			}
			if (!isNaN(s)) {
				for (i = 0; i < m; i++) {
		         for (j = 0; j < n; j++) {
		            A[i][j] = s;
		         }
		      }
			}
		}


			/** Construct a matrix from a copy of a 2-D array.
			@param arr    Two-dimensional array of doubles.
			*/

		public static function constructWithCopy(arr:Array):RealMatrix {
			var mm:int = arr.length;
			var nn:int = arr[0].length;
			var X:RealMatrix = new RealMatrix(mm,nn);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < mm; i++) {
				if (arr[i].length != nn) {
					throw new Error	("All rows must have the same length.");
				}
				for (j = 0; j < nn; j++) {
					C[i][j] = arr[i][j];
					/*trace("C[" + i + "][" + j + "]: " + C[i][j]);*/
				}
			}
			return X;
		}


		/** Access the internal two-dimensional array.
			@return     Pointer to the two-dimensional array of matrix elements.
			*/

		public function getArray ():Array {
			return A;
		}

		/** Copy the internal two-dimensional array.
			@return     Two-dimensional array copy of matrix elements.
			*/

		public function getArrayCopy ():Array {
			var C:Array = new Array(m);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				C[i] = new Array(n);
				for (j = 0; j < n; j++) {
					C[i][j] = A[i][j];
				}
			}
			return C;
		}

		/** Make a one-dimensional column packed copy of the internal array.
			@return     Matrix elements packed in a one-dimensional array by columns.
			*/

		public function getColumnPackedCopy ():Array {
			var vals:Array = new Array(m*n);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					vals[i+j*m] = A[i][j];
				}
			}
			return vals;
		}

		/** Make a one-dimensional row packed copy of the internal array.
			@return     Matrix elements packed in a one-dimensional array by rows.
			*/

		public function getRowPackedCopy ():Array {
			var vals:Array = new Array(m*n);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					vals[i*n+j] = A[i][j];
				}
			}
			return vals;
		}

		/** Get row dimension.
			@return     m, the number of rows.
			*/

		public function getRowDimension ():int {
			return m;
		}

		/** Get column dimension.
			@return     n, the number of columns.
			*/

		public function getColumnDimension ():int {
			return n;
		}

		/** Get a single element.
			@param i    Row index.
			@param j    Column index.
			@return     A(i,j)
			@exception  ArrayIndexOutOfBoundsException
			*/

		public function get (i:int, j:int):Number {
			return A[i][j];
		}

		/** Get a submatrix.
			@param i0   Initial row index
			@param i1   Final row index
			@param j0   Initial column index
			@param j1   Final column index
			@return     A(i0:i1,j0:j1)
			*/

		public function getMatrix (i0:int, i1:int, j0:int, j1:int):RealMatrix {
			var X:RealMatrix = new RealMatrix(i1-i0+1,j1-j0+1);
			var B:Array = X.getArray();
			for (var i:int = i0; i <= i1; i++) {
				for (var j:int = j0; j <= j1; j++) {
					B[i-i0][j-j0] = A[i][j];
				}
			}
			return X;
		}


		/** Set a single element.
			@param i    Row index.
			@param j    Column index.
			@param s    A(i,j).
		*/

		public function set (i:int, j:int, s:Number):void {
			A[i][j] = s;
		}


		/** Matrix transpose.
			@return    A
		*/

		public function transpose ():RealMatrix {
			var X:RealMatrix = new RealMatrix(n,m);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[j][i] = A[i][j];
				}
			}
			return X;
		}

		/** One norm
			@return    maximum column sum.
			*/

		public function norm1 ():Number {
			var f:Number = 0;
			var i:int;
			var s:Number;
			for (var j:int = 0; j < n; j++) {
				s = 0;
				for (i = 0; i < m; i++) {
					s += Math.abs(A[i][j]);
				}
				f = Math.max(f,s);
			}
			return f;
		}

		/** Two norm
			@return    maximum singular value.
		*/

		public function norm2 ():Number {
			return (new SingularValueDecomposition(this).norm2());
		}

		/** Infinity norm
			@return    maximum row sum.
		*/

		public function normInf ():Number {
			var f:Number = 0;
			var j:int;
			var s:Number;
			for (var i:int = 0; i < m; i++) {
				s = 0;
				for (j = 0; j < n; j++) {
					s += Math.abs(A[i][j]);
				}
				f = Math.max(f,s);
			}
			return f;
		}

		/** Frobenius norm
			@return    sqrt of sum of squares of all elements.
			*/

		public function normF ():Number {
			var f:Number = 0;
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					f = Maths.hypot(f,A[i][j]);
				}
			}
			return f;
		}

		/**  Unary minus
			@return    -A
		*/

		public function uminus ():RealMatrix {
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = -A[i][j];
				}
			}
			return X;
		}

		/** C = A + B
			@param B    another matrix
			@return     A + B
		*/

		public function plus (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = A[i][j] + B.A[i][j];
				}
			}
			return X;
		}

		/** A = A + B
			@param B    another matrix
			@return     A + B
		*/

		public function plusEquals (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] + B.A[i][j];
				}
			}
			return this;
		}

		/** C = A - B
			@param B    another matrix
			@return     A - B
			*/

		public function minus (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = A[i][j] - B.A[i][j];
				}
			}
			return X;
		}

		/** A = A - B
			@param B    another matrix
			@return     A - B
			*/

		public function minusEquals (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] - B.A[i][j];
				}
			}
			return this;
		}

		/** Element-by-element multiplication, C = A.*B
			@param B    another matrix
			@return     A.*B
			*/

		public function arrayTimes (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = A[i][j] * B.A[i][j];
				}
			}
			return X;
		}

		/** Element-by-element multiplication in place, A = A.*B
			@param B    another matrix
			@return     A.*B
			*/

		public function arrayTimesEquals (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] * B.A[i][j];
				}
			}
			return this;
		}

		/** Element-by-element right division, C = A./B
			@param B    another matrix
			@return     A./B
			*/

		public function arrayRightDivide (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var X:RealMatrix = new RealMatrix(m,n);
			var j:int;
			var C:Array = X.getArray();
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = A[i][j] / B.A[i][j];
				}
			}
			return X;
		}

		/** Element-by-element right division in place, A = A./B
			@param B    another matrix
			@return     A./B
			*/

		public function arrayRightDivideEquals (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = A[i][j] / B.A[i][j];
				}
			}
			return this;
		}

		/** Element-by-element left division, C = A.\B
			@param B    another matrix
			@return     A.\B
			*/

		public function arrayLeftDivide (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = B.A[i][j] / A[i][j];
				}
			}
			return X;
		}

		/** Element-by-element left division in place, A = A.\B
			@param B    another matrix
			@return     A.\B
			*/

		public function arrayLeftDivideEquals (B:RealMatrix):RealMatrix {
			checkMatrixDimensions(B);
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = B.A[i][j] / A[i][j];
				}
			}
			return this;
		}

		/** Multiply a matrix by a scalar, C = s*A
			@param s    scalar
			@return     s*A
			*/

		public function times (s:*):RealMatrix {
			if (s is RealMatrix) {
				return timesMatrix(s);
			}
			var X:RealMatrix = new RealMatrix(m,n);
			var C:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					C[i][j] = s*A[i][j];
				}
			}
			return X;
		}

		/** Multiply a matrix by a scalar in place, A = s*A
			@param s    scalar
			@return     replace A by s*A
			*/

		public function timesEquals (s:Number):RealMatrix {
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					A[i][j] = s*A[i][j];
				}
			}
			return this;
		}

		/** Linear algebraic matrix multiplication, A * B
			@param B    another matrix
			@return     Matrix product, A * B
			*/

		public function timesMatrix (B:RealMatrix):RealMatrix {
			if (B.m != n) {
				throw new Error("Matrix inner dimensions must agree.");
			}
			var X:RealMatrix = new RealMatrix(m,B.n);
			var C:Array = X.getArray();
			var Bcolj:Array = new Array(n);
			var Arowi:Array;
			var s:Number;
			var k:int;
			var i:int;
			for (var j:int = 0; j < B.n; j++) {
				for (k = 0; k < n; k++) {
					Bcolj[k] = B.A[k][j];
				}
				for (i = 0; i < m; i++) {
					Arowi = A[i];
					s = 0;
					for (k = 0; k < n; k++) {
						s += Arowi[k]*Bcolj[k];
					}
					C[i][j] = s;
				}
			}
			return X;
		}

		/** LU Decomposition
			@return     LUDecomposition
			@see LUDecomposition
			*/

		/*public function lu ():LUDecomposition {
			return new LUDecomposition(this);
		}*/

		/** QR Decomposition
			@return     QRDecomposition
			@see QRDecomposition
			*/

		public function qr ():QRDecomposition {
			return new QRDecomposition(this);
		}

		/** Cholesky Decomposition
			@return     CholeskyDecomposition
			@see CholeskyDecomposition
			*/

		/*public CholeskyDecomposition chol () {
			return new CholeskyDecomposition(this);
		}*/

		/** Singular Value Decomposition
			@return     SingularValueDecomposition
			@see SingularValueDecomposition
			*/

		public function svd ():SingularValueDecomposition {
			return new SingularValueDecomposition(this);
		}

		/** Eigenvalue Decomposition
			@return     EigenvalueDecomposition
			@see EigenvalueDecomposition
			*/

		/*public EigenvalueDecomposition eig () {
			return new EigenvalueDecomposition(this);
		}*/

		/** Solve A*X = B
			@param B    right hand side
			@return     solution if A is square, least squares solution otherwise
		*/

		public function solve (B:RealMatrix):RealMatrix {
			/*return (m == n ? (new LUDecomposition(this)).solve(B) :*/
			return (new QRDecomposition(this)).solve(B);
		}

		/** Solve X*A = B, which is also A'*X' = B'
			@param B    right hand side
			@return     solution if A is square, least squares solution otherwise.
		*/

		public function solveTranspose (B:RealMatrix):RealMatrix {
			return transpose().solve(B.transpose());
		}

		/** Matrix inverse or pseudoinverse
			@return     inverse(A) if A is square, pseudoinverse otherwise.
			*/

		public function inverse ():RealMatrix {
			return solve(identity(m,m));
		}

		/** Matrix determinant
			@return     determinant
			*/

		/*public function det ():Number {
			return new LUDecomposition(this).det();
		}*/

		/** Matrix rank
			@return     effective numerical rank, obtained from SVD.
			*/

		public function rank ():int {
			return new SingularValueDecomposition(this).rank();
		}

		/** Matrix condition (2 norm)
			@return     ratio of largest to smallest singular value.
			*/

		public function cond ():Number {
			return new SingularValueDecomposition(this).cond();
		}

		/** Matrix trace.
			@return     sum of the diagonal elements.
			*/

		public function trace ():Number {
			var t:Number = 0;
			for (var i:int = 0; i < Math.min(m,n); i++) {
				t += A[i][i];
			}
			return t;
		}

		/** Generate matrix with random elements
			@param m    Number of rows.
			@param n    Number of colums.
			@return     An m-by-n matrix with uniformly distributed random elements.
			*/

		public static function random (m:int, n:int):RealMatrix {
			var A:RealMatrix = new RealMatrix(m,n);
			var X:Array = A.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					X[i][j] = Math.random();
				}
			}
			return A;
		}

		/** Generate identity matrix
			@param m    Number of rows.
			@param n    Number of colums.
			@return     An m-by-n matrix with ones on the diagonal and zeros elsewhere.
			*/

		public static function identity (m:int, n:int):RealMatrix {
			var A:RealMatrix = new RealMatrix(m,n);
			var X:Array = A.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					X[i][j] = (i == j ? 1.0 : 0.0);
				}
			}
			return A;
		}



			/** Check if size(A) == size(B) **/

		private function checkMatrixDimensions (B:RealMatrix):void {
			if (B.m != m || B.n != n) {
				throw new Error("Matrix dimensions must agree.");
			}
		}

	}

}
