/* AS3
	Copyright 2007 Sphex LLP.
	Ported from Jama: http://math.nist.gov/javanumerics/jama/
	
	This work is licensed under the Creative Commons Attribution 2.0 UK: England & Wales License. 
	To view a copy of this license, visit http://creativecommons.org/licenses/by/2.0/uk/ or send a 
	letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

*/

package com.clevr.matrixalgebra {

	/** Singular Value Decomposition.
		For an m-by-n matrix A with m >= n, the singular value decomposition is
		an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
		an n-by-n orthogonal matrix V so that A = U*S*V'.
		The singular values, sigma[k] = S[k][k], are ordered so that
		sigma[0] >= sigma[1] >= ... >= sigma[n-1].
		The singular value decompostion always exists, so the constructor will
		never fail.  The matrix condition number and the effective numerical
		rank can be computed from this decomposition.
		*/

	public class SingularValueDecomposition {

			/** Arrays for internal storage of U and V.
			*/
			private var U:Array, V:Array;

		/** Array for internal storage of singular values.
			*/
			private var s:Array;

		/** Row and column dimensions.
			*/
			private var m:int, n:int;

		/* ------------------------
			Constructor
			* ------------------------ */

			/** Construct the singular value decomposition
			@param Arg    Rectangular matrix
			*/

		public function SingularValueDecomposition (Arg:RealMatrix) {

			// Derived from LINPACK code.
			// Initialize.
			var yyyy:int;
			
			var A:Array = Arg.getArrayCopy();
			m = Arg.getRowDimension();
			n = Arg.getColumnDimension();

			var nu:int = Math.min(m,n);
			s = new Array (Math.min(m+1,n));
			U = new Array (m);
			var i:int;
			for (i = 0; i < U.length; i++) 
				U[i] = new Array(nu)
				
			V = new Array(n);
			for (i = 0; i < V.length; i++) 
				V[i] = new Array(m);
				
			var e:Array = new Array(n);
			var work:Array = new Array(m);
			var wantu:Boolean = true;
			var wantv:Boolean = true;

			// Reduce A to bidiagonal form, storing the diagonal elements
			// in s and the super-diagonal elements in e.

			var nct:int = Math.min(m-1,n);
			var nrt:int = Math.max(0,Math.min(n-2,m));
			var j:int;
			var t:Number;
			
			for (var k:int = 0; k < Math.max(nct,nrt); k++) {
				if (k < nct) {

					// Compute the transformation for the k-th column and
					// place the k-th diagonal in s[k].
					// Compute 2-norm of k-th column without under/overflow.
					s[k] = 0;
					for (i = k; i < m; i++) {
						s[k] = Maths.hypot(s[k],A[i][k]);
					}
					if (s[k] != 0.0) {
						if (A[k][k] < 0) {
							s[k] = -s[k];
						}
						for (i = k; i < m; i++) {
							A[i][k] /= s[k];
						}
						A[k][k] += 1;
					}
					s[k] = -s[k];
				}
				for (j = k+1; j < n; j++) {
					if ((k < nct) && (s[k] != 0.0))  {

						// Apply the transformation.

						t = 0;
						for (i = k; i < m; i++) {
							t += A[i][k]*A[i][j];
						}
						t = -t/A[k][k];
						for (i = k; i < m; i++) {
							A[i][j] += t*A[i][k];
						}
					}

					// Place the k-th row of A into e for the
					// subsequent calculation of the row transformation.

					e[j] = A[k][j];
					/*trace ('e[' + j + ']: ' + e[j]);*/
				}
				if (wantu && (k < nct)) {

					// Place the transformation in U for subsequent back
					// multiplication.

					for (i = k; i < m; i++) {
						U[i][k] = A[i][k];
					}
				}
				if (k < nrt) {

					// Compute the k-th row transformation and place the
					// k-th super-diagonal in e[k].
					// Compute 2-norm without under/overflow.
					e[k] = 0;
					for (i = k+1; i < n; i++) {
						e[k] = Maths.hypot(e[k],e[i]);
					}
					if (e[k] != 0.0) {
						if (e[k+1] < 0.0) {
							e[k] = -e[k];
						}
						for (i = k+1; i < n; i++) {
							e[i] /= e[k];
						}
						e[k+1] += 1.0;
					}
					e[k] = -e[k];
					if ((k+1 < m) && (e[k] != 0.0)) {

						// Apply the transformation.

						for (i = k+1; i < m; i++) {
							work[i] = 0.0;
						}
						for (j = k+1; j < n; j++) {
							for (i = k+1; i < m; i++) {
								work[i] += e[j]*A[i][j];
							}
						}
						for (j = k+1; j < n; j++) {
							t = -e[j]/e[k+1];
							for (i = k+1; i < m; i++) {
								A[i][j] += t*work[i];
							}
						}
					}
					if (wantv) {

						// Place the transformation in V for subsequent
						// back multiplication.

						for (i = k+1; i < n; i++) {
							V[i][k] = e[i];
						}
					}
				}
			}

			// Set up the final bidiagonal matrix or order p.

			var p:int = Math.min(n,m+1);
			if (nct < n) {
				s[nct] = A[nct][nct];
			}
			if (m < p) {
				s[p-1] = 0.0;
			}
			if (nrt+1 < p) {
				e[nrt] = A[nrt][p-1];
			}
			e[p-1] = 0.0;

			// If required, generate U.

			if (wantu) {
				for (j = nct; j < nu; j++) {
					for (i = 0; i < m; i++) {
						U[i][j] = 0.0;
					}
					U[j][j] = 1.0;
				}
				for (k = nct-1; k >= 0; k--) {
					if (s[k] != 0.0) {
						for (j = k+1; j < nu; j++) {
							t = 0;
							for (i = k; i < m; i++) {
								t += U[i][k]*U[i][j];
							}
							t = -t/U[k][k];
							for (i = k; i < m; i++) {
								U[i][j] += t*U[i][k];
							}
						}
						for (i = k; i < m; i++ ) {
							U[i][k] = -U[i][k];
						}
						U[k][k] = 1.0 + U[k][k];
						for (i = 0; i < k-1; i++) {
							U[i][k] = 0.0;
						}
					} else {
						for (i = 0; i < m; i++) {
							U[i][k] = 0.0;
						}
						U[k][k] = 1.0;
					}
				}
			}

			// If required, generate V.

			if (wantv) {
				for (k = n-1; k >= 0; k--) {
					if ((k < nrt) && (e[k] != 0.0)) {
						for (j = k+1; j < nu; j++) {
							t = 0;
							for (i = k+1; i < n; i++) {
								t += V[i][k]*V[i][j];
							}
							t = -t/V[k+1][k];
							for (i = k+1; i < n; i++) {
								V[i][j] += t*V[i][k];
							}
						}
					}
					for (i = 0; i < n; i++) {
						V[i][k] = 0.0;
					}
					V[k][k] = 1.0;
				}
			}

			// Main iteration loop for the singular values.

			var pp:int = p-1;
			var iter:int = 0;
			var eps:Number = Math.pow(2.0,-52.0);
			var tiny:Number = Math.pow(2.0,-966.0);
			var ks:int;
			var ttt:Number;
			var f:Number;
			var kase:int;
			var iteration:int = 0;
			var debug:Boolean = false;
			while (p > 0) {
				/*if (iteration++ % 100 == 0) {
					trace('iteration: ' + iteration + ' p: ' + p);
					debug = true;
				} else debug = false;*/
				// Here is where a test for too many iterations would go.

				// This section of the program inspects for
				// negligible elements in the s and e arrays.  On
				// completion the variables kase and k are set as follows.

				// kase = 1     if s(p) and e[k-1] are negligible and k<p
				// kase = 2     if s(k) is negligible and k<p
				// kase = 3     if e[k-1] is negligible, k<p, and
				//              s(k), ..., s(p) are not negligible (qr step).
				// kase = 4     if e(p-1) is negligible (convergence).

				for (k = p-2; k >= -1; k--) {
					if (debug) {
						trace('k: ' + k);
					}
					if (k == -1) {
						break;
					}
					if (Math.abs(e[k]) <= tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))) {
						e[k] = 0;
						break;
					} else if (debug) {
						trace('e[k]: ' + Math.abs(e[k]) + ' sth: ' + (tiny + eps*(Math.abs(s[k]) + Math.abs(s[k+1]))));
					}
				}
				if (k == p-2) {
					kase = 4;
				} else {
					for (ks = p-1; ks >= k; ks--) {
						if (ks == k) {
							break;
						}
						ttt = (ks != p ? Math.abs(e[ks]) : 0.) + 
							(ks != k+1 ? Math.abs(e[ks-1]) : 0.);
						if (Math.abs(s[ks]) <= tiny + eps*ttt)  {
							s[ks] = 0.0;
							break;
						}
					}
					if (ks == k) {
						kase = 3;
					} else if (ks == p-1) {
						kase = 1;
					} else {
						kase = 2;
						k = ks;
					}
				}
				k++;

				// Perform the task indicated by kase.

				switch (kase) {

					// Deflate negligible s(p).

					case 1: {
						f = e[p-2];
						e[p-2] = 0.0;
						for (j = p-2; j >= k; j--) {
							t = Maths.hypot(s[j],f);
							cs = s[j]/t;
							sn = f/t;
							s[j] = t;
							if (j != k) {
								f = -sn*e[j-1];
								e[j-1] = cs*e[j-1];
							}
							if (wantv) {
								for (i = 0; i < n; i++) {
									t = cs*V[i][j] + sn*V[i][p-1];
									V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
									V[i][j] = t;
								}
							}
						}
					}
					break;

					// Split at negligible s(k).

					case 2: {
						f = e[k-1];
						e[k-1] = 0.0;
						for (j = k; j < p; j++) {
							ttt = Maths.hypot(s[j],f);
							var cs:Number = s[j]/ttt;
							var sn:Number = f/ttt;
							s[j] = ttt;
							f = -sn*e[j];
							e[j] = cs*e[j];
							if (wantu) {
								for (i = 0; i < m; i++) {
									ttt = cs*U[i][j] + sn*U[i][k-1];
									U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
									U[i][j] = ttt;
								}
							}
						}
					}
					break;

					// Perform one qr step.

					case 3: {

						// Calculate the shift.

						var scale:Number = Math.max(Math.max(Math.max(Math.max(
							Math.abs(s[p-1]),Math.abs(s[p-2])),Math.abs(e[p-2])), 
							Math.abs(s[k])),Math.abs(e[k]));
						var sp:Number = s[p-1]/scale;
						var spm1:Number = s[p-2]/scale;
						var epm1:Number = e[p-2]/scale;
						var sk:Number = s[k]/scale;
						var ek:Number = e[k]/scale;
						var b:Number = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
						var c:Number = (sp*epm1)*(sp*epm1);
						var shift:Number = 0.0;
						if ((b != 0.0) || (c != 0.0)) {
							shift = Math.sqrt(b*b + c);
							if (b < 0.0) {
								shift = -shift;
							}
							shift = c/(b + shift);
						}
						f = (sk + sp)*(sk - sp) + shift;
						var g:Number = sk*ek;
						// Chase zeros.
						
						for (j = k; j < p-1; j++) {
							ttt = Maths.hypot(f,g);
							cs = f/ttt;
							sn = g/ttt;
							if (j != k) {
								e[j-1] = ttt;
							}
							f = cs*s[j] + sn*e[j];
							e[j] = cs*e[j] - sn*s[j];
							g = sn*s[j+1];
							s[j+1] = cs*s[j+1];
							if (wantv) {
								for (i = 0; i < n; i++) {
									ttt = cs*V[i][j] + sn*V[i][j+1];
									V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
									V[i][j] = ttt;
								}
							}
							ttt = Maths.hypot(f,g);
							cs = f/ttt;
							sn = g/ttt;
							s[j] = ttt;
							f = cs*e[j] + sn*s[j+1];
							s[j+1] = -sn*e[j] + cs*s[j+1];
							g = sn*e[j+1];
							e[j+1] = cs*e[j+1];
							if (wantu && (j < m-1)) {
								for (i = 0; i < m; i++) {
									ttt = cs*U[i][j] + sn*U[i][j+1];
									U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
									U[i][j] = ttt;
								}
							}
						}
						e[p-2] = f;
						iter = iter + 1;
					}
					break;

					// Convergence.

					case 4: {

						// Make the singular values positive.

						if (s[k] <= 0.0) {
							s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
							if (wantv) {
								for (i = 0; i <= pp; i++) {
									V[i][k] = -V[i][k];
								}
							}
						}

						// Order the singular values.

						while (k < pp) {
							if (s[k] >= s[k+1]) {
								break;
							}
							ttt = s[k];
							s[k] = s[k+1];
							s[k+1] = ttt;
							if (wantv && (k < n-1)) {
								for (i = 0; i < n; i++) {
									ttt = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = ttt;
								}
							}
							if (wantu && (k < m-1)) {
								for (i = 0; i < m; i++) {
									ttt = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = ttt;
								}
							}
							k++;
						}
						iter = 0;
						p--;
					}
					break;
				}
			}
		}

		/* ------------------------
			Public Methods
			* ------------------------ */

			/** Return the left singular vectors
			@return     U
			*/

		public function getU ():RealMatrix {
			return new RealMatrix(U,m,Math.min(m+1,n));
		}

		/** Return the right singular vectors
			@return     V
			*/

		public function getV ():RealMatrix {
			return new RealMatrix(V,n,n);
		}

		/** Return the one-dimensional array of singular values
			@return     diagonal of S.
			*/

		public function getSingularValues ():Array {
			return s;
		}

		/** Return the diagonal matrix of singular values
			@return     S
			*/

		public function getS ():RealMatrix {
			var X:RealMatrix = new RealMatrix(n,n);
			var S:Array = X.getArray();
			var j:int;
			for (var i:int = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					S[i][j] = 0.0;
				}
				S[i][i] = this.s[i];
			}
			return X;
		}

		/** Two norm
			@return     max(S)
			*/

		public function norm2 ():Number {
			return s[0];
		}

		/** Two norm condition number
			@return     max(S)/min(S)
			*/

		public function cond ():Number {
			return s[0]/s[Math.min(m,n)-1];
		}

		/** Effective numerical matrix rank
			@return     Number of nonnegligible singular values.
			*/

		public function rank ():int {
			var eps:Number = Math.pow(2.0,-52.0);
			var tol:Number = Math.max(m,n)*s[0]*eps;
			var r:int = 0;
			for (var i:int = 0; i < s.length; i++) {
				if (s[i] > tol) {
					r++;
				}
			}
			return r;
		}
	}

}
