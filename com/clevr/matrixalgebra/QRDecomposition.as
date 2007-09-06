/* AS3
	Copyright 2007 Sphex LLP
	Ported from Jama: http://math.nist.gov/javanumerics/jama/
	
	This work is licensed under the Creative Commons Attribution 2.0 UK: England & Wales License. 
	To view a copy of this license, visit http://creativecommons.org/licenses/by/2.0/uk/ or send a 
	letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

*/

package com.clevr.matrixalgebra {

	/** QR Decomposition.

		For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
		orthogonal matrix Q and an n-by-n upper triangular matrix R so that
		A = Q*R.
		The QR decompostion always exists, even if the matrix does not have
		full rank, so the constructor will never fail.  The primary use of the
		QR decomposition is in the least squares solution of nonsquare systems
		of simultaneous linear equations.  This will fail if isFullRank()
		*/

	public class QRDecomposition {

		/* ------------------------
			Class variables
			* ------------------------ */

			/** Array for internal storage of decomposition.
			*/
			private var QR:Array;

		/** Row and column dimensions.
			*/
			private var m:int, n:int;

		/** Array for internal storage of diagonal of R.
			*/
			private var Rdiag:Array;

		/* ------------------------
			Constructor
			* ------------------------ */

			/** QR Decomposition, computed by Householder reflections.
			@param A    Rectangular matrix
			*/

		public function QRDecomposition (A:RealMatrix) {
			// Initialize.
			QR = A.getArrayCopy();
			m = A.getRowDimension();
			n = A.getColumnDimension();
			Rdiag = new Array(n);

			// Main loop.
			var i:int;
			var j:int;
			var nrm:Number;
			var s:Number;
			for (var k:int = 0; k < n; k++) {
				// Compute 2-norm of k-th column without under/overflow.
				nrm = 0;
				for (i = k; i < m; i++) {
					nrm = Maths.hypot(nrm,QR[i][k]);
				}

				if (nrm != 0.0) {
					// Form k-th Householder vector.
					if (QR[k][k] < 0) {
						nrm = -nrm;
					}
					for (i = k; i < m; i++) {
						QR[i][k] /= nrm;
					}
					QR[k][k] += 1.0;

					// Apply transformation to remaining columns.
					for (j = k+1; j < n; j++) {
						s = 0.0; 
						for (i = k; i < m; i++) {
							s += QR[i][k]*QR[i][j];
						}
						s = -s/QR[k][k];
						for (i = k; i < m; i++) {
							QR[i][j] += s*QR[i][k];
						}
					}
				}
				Rdiag[k] = -nrm;
			}
		}

		/* ------------------------
			Public Methods
		* ------------------------ */

		/** Is the matrix full rank?
			@return     true if R, and hence A, has full rank.
		*/

		public function isFullRank ():Boolean {
			for (var j:int = 0; j < n; j++) {
				if (Rdiag[j] == 0)
					return false;
			}
			return true;
		}

		/** Return the Householder vectors
			@return     Lower trapezoidal matrix whose columns define the reflections
		*/

		public function getH ():RealMatrix {
			var X:RealMatrix = new RealMatrix(m,n);
			var H:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < m; i++) {
				for (j = 0; j < n; j++) {
					if (i >= j) {
						H[i][j] = QR[i][j];
					} else {
						H[i][j] = 0.0;
					}
				}
			}
			return X;
		}

		/** Return the upper triangular factor
			@return     R
		*/

		public function getR ():RealMatrix {
			var X:RealMatrix = new RealMatrix(n,n);
			var R:Array = X.getArray();
			var j:int
			for (var i:int = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (i < j) {
						R[i][j] = QR[i][j];
					} else if (i == j) {
						R[i][j] = Rdiag[i];
					} else {
						R[i][j] = 0.0;
					}
				}
			}
			return X;
		}

		/** Generate and return the (economy-sized) orthogonal factor
			@return     Q
		*/

		public function getQ ():RealMatrix {
			var X:RealMatrix = new RealMatrix(m,n);
			var Q:Array = X.getArray();
			var i:int;
			var k:int;
			var j:int;
			var s:Number;
			for (k = n-1; k >= 0; k--) {
				for (i = 0; i < m; i++) {
					Q[i][k] = 0.0;
				}
				Q[k][k] = 1.0;
				for (j = k; j < n; j++) {
					if (QR[k][k] != 0) {
						s = 0.0;
						for (i = k; i < m; i++) {
							s += QR[i][k]*Q[i][j];
						}
						s = -s/QR[k][k];
						for (i = k; i < m; i++) {
							Q[i][j] += s*QR[i][k];
						}
					}
				}
			}
			return X;
		}

		/** Least squares solution of A*X = B
			@param B    A Matrix with as many rows as A and any number of columns.
			@return     X that minimizes the two norm of Q*R*X-B.
			*/

		public function solve (B:RealMatrix):RealMatrix {
			if (B.getRowDimension() != m) {
				throw new Error("Matrix row dimensions must agree.");
			}
			if (!this.isFullRank()) {
				throw new Error("Matrix is rank deficient.");
			}

			// Copy right hand side
			var nx:int = B.getColumnDimension();
			var X:Array = B.getArrayCopy();
			
			var k:int;
			var j:int;
			var i:int;
			var s:Number;
			
			// Compute Y = transpose(Q)*B
			for (k = 0; k < n; k++) {
				for (j = 0; j < nx; j++) {
					s = 0; 
					for (i = k; i < m; i++) {
						s += QR[i][k]*X[i][j];
					}
					s = -s/QR[k][k];
					for (i = k; i < m; i++) {
						X[i][j] += s*QR[i][k];
					}
				}
			}
			// Solve R*X = Y;
			for (k = n-1; k >= 0; k--) {
				for (j = 0; j < nx; j++) {
					X[k][j] /= Rdiag[k];
				}
				for (i = 0; i < k; i++) {
					for (j = 0; j < nx; j++) {
						X[i][j] -= X[k][j]*QR[i][k];
					}
				}
			}
			return (new RealMatrix(X,n,nx).getMatrix(0,n-1,0,nx-1));
		}
	}

}
