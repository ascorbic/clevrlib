/* AS3
	Copyright 2007 Sphex LLP.
	
	This work is licensed under the Creative Commons Attribution 2.0 UK: England & Wales License. 
	To view a copy of this license, visit http://creativecommons.org/licenses/by/2.0/uk/ or send a 
	letter to Creative Commons, 171 Second Street, Suite 300, San Francisco, California, 94105, USA.

*/
package com.clevr.graphics {
	
	import flash.display.BitmapData;
	import flash.utils.ByteArray;
	import flash.geom.Rectangle;
	
	/*
	*	Adds interpolation to bitmapdata
	*
	*	@langversion ActionScript 3.0
	*	@playerversion Flash 9.0
	*
	*	@author Matt Kane
	*	@since  26.05.2007
	*/
	public class InterpolatedBitmapData extends BitmapData {

		public function InterpolatedBitmapData(width:int, height:int, transparent:Boolean = true, fillColor:uint = 0xFFFFFFFF) {
			super(width, height, transparent, fillColor);
		}
		
	   /**
	      Computes the value of a pixel that is not on a pixel boundary.

	      @param x The sub-pixel precision x-coordinate.
	      @param y The sub-pixel precision y-coordinate.
	   */
	   public function getPixelBicubic(x:Number, y:Number): Number {
	      var i:int = (int)(x);
	      var j:int = (int)(y);

	      if(((i - 1) < 0) || ((j - 1) < 0))
	         return getPixel(i, j);
	      else if(((i + 1) >= width) || ((i + 2) >= width) ||
	         ((j + 1) >= height) || ((j + 2) >= height))
	         return getPixel(i, j);

	      var dx:Number = x - i;
	      var dy:Number = y - j;
	      var rsum:Number = 0;
	      var gsum:Number = 0;
	      var bsum:Number = 0;

	      var Rmdx:Array = new Array(4);
	      var Rdyn:Array = new Array(4);
	      for(var m:int = -1; m <= 2; ++m)
	         Rmdx[m + 1] = A(m - dx);
	      for(var n:int = -1; n <= 2; ++n)
	         Rdyn[n + 1] = A(dy - n);

	      for(m = -1; m <= 2; ++m)
	      {
	         for(n = -1; n <= 2; ++n)
	         {
	            var rgb:int = getPixel(i + m, j + n);
	            var r:int = (rgb >> 16) & 0x0FF;
	            var g:int = (rgb >>  8) & 0x0FF;
	            var b:int = (rgb >>  0) & 0x0FF;

	            var Rres:Number = Rmdx[m + 1] * Rdyn[n + 1];
	            rsum += r * Rres;
	            gsum += g * Rres;
	            bsum += b * Rres;
	         }
	      }

	      var red:int = (int)(rsum + 0.5);
	      if(red < 0)
	         red = 0;
	      else if(red > 255)
	         red = 255;
	      var green:int = (int)(gsum + 0.5);
	      if(green < 0)
	         green = 0;
	      else if(green > 255)
	         green = 255;
	      var blue:int = (int)(bsum + 0.5);
	      if(blue < 0)
	         blue = 0;
	      else if(blue > 255)
	         blue = 255;

	      return (red << 16) | (green << 8) | (blue << 0);
	   }
		
		
		
	   /**
	      Computes the value of a pixel that is not on a pixel boundary.

	      @param theX The sub-pixel precision x-coordinate.
	      @param theY The sub-pixel precision y-coordinate.
	   */
		public function getPixelBilinear(theX:Number, theY:Number): Number {
			
			var x:int;
			var y:int;
			var x_ratio:Number;
			var y_ratio:Number;
			var y_opposite:Number;
			var x_opposite:Number;
			var a:int;
			var be:int;
			var c:int;
			var d:int;
			var red:int;
			var green:int;
			var blue:int;
			
			
			x = Math.floor(theX);
			y = Math.floor(theY);
			
			if((x < 1) || (y < 1) || ((x + 2) >= width) || ((y + 2) >= height))
				return getPixel(x, y);
			
			x_ratio = theX - x;
			y_ratio = theY - y;
			x_opposite = 1 - x_ratio;
			y_opposite = 1 - y_ratio;
						
			a = getPixel(x, y);
			be =getPixel(x + 1, y);
			c = getPixel(x, y + 1);
			d = getPixel(x + 1, y + 1);
			red 	= (r(a)  * x_opposite  + r(be)   * x_ratio) * y_opposite + (r(c) * x_opposite  + r(d) * x_ratio) * y_ratio;
			green 	= (g(a)  * x_opposite  + g(be)   * x_ratio) * y_opposite + (g(c) * x_opposite  + g(d) * x_ratio) * y_ratio;
			blue 	= (b(a)  * x_opposite  + b(be)   * x_ratio) * y_opposite + (b(c) * x_opposite  + b(d) * x_ratio) * y_ratio;
				/*red = r(a);
				green = g(a);
				blue = b(a);*/
				
			if(red < 0)
				red = 0;
			else if(red > 255)
				red = 255;
			if(green < 0)
				green = 0;
			else if(green > 255)
				green = 255;
			if(blue < 0)
				blue = 0;
			else if(blue > 255)
				blue = 255;

			return (red << 16) | (green << 8) | (blue << 0);
		}
		
		
		/**
	      Retrieve the intensity of a pixel;

	      @param x The x coordinate of the pixel
	      @param y The y coordinate of the pixel
	      @return The intensity at pixel (x, y)
	   */
	   public function getIntensity (x:int, y:int):int {
	      var rgb:int = getPixel(x, y);
	      var r:int = (rgb >> 16) & 0x0FF;
	      var g:int = (rgb >>  8) & 0x0FF;
	      var b:int = (rgb >>  0) & 0x0FF;

	      return Math.ceil((r + g + b) / 3);
	   }
	  
		

	   /**
	      Support function for bicubic interpolation.
	   */
	   private static function A(x:Number):Number
	   {
	      var p0:Number = ((x + 2) > 0) ? (x + 2) : 0
	      var p1:Number = ((x + 1) > 0) ? (x + 1) : 0
	      var p2:Number = (x > 0) ? x : 0
	      var p3:Number = ((x - 1) > 0) ? (x - 1) : 0

	      return (1 / 6) * (p0 * p0 * p0 - 4 * (p1 * p1 * p1) + 6 * (p2 * p2 * p2) -
	         4 * (p3 * p3 * p3));
	   }

	   /**
	      Support function for bicubic interpolation.
	   */
	   private static function P(x:Number):Number
	   {
	      return (x > 0) ? x : 0;
	   }
	
		/**
		*	RGB convenience methods	
		*/
	
		private static function r(rgb:int):int {
			return (rgb >> 16) & 0x0FF;
		}

		private static function g(rgb:int):int {
			return (rgb >> 8) & 0x0FF;
		}

		private static function b(rgb:int):int {
			return (rgb >> 0) & 0x0FF;
		}

		
	}
	
}
