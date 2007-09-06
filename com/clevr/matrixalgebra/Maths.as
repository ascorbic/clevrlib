/* AS3
	Copyright 2007 Sphex LLP
	Ported from Jama: http://math.nist.gov/javanumerics/jama/
	
	This work is licensed under the Creative Commons Attribution 2.0 UK: England & Wales License. 
	To view a copy of this license, visit http://creativecommons.org/licenses/by/2.0/uk/ or send a 
	letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

*/
package com.clevr.matrixalgebra {
	
	public class Maths {

	   /** sqrt(a^2 + b^2) without under/overflow. **/

	   public static function hypot(a:Number, b:Number):Number {
	      var r:Number;
	      if (Math.abs(a) > Math.abs(b)) {
	         r = b/a;
	         r = Math.abs(a)*Math.sqrt(1+r*r);
	      } else if (b != 0) {
	         r = a/b;
	         r = Math.abs(b)*Math.sqrt(1+r*r);
	      } else {
	         r = 0.0;
	      }
	      return r;
	   }
	}
	
}
