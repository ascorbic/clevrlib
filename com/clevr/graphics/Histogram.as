/* AS3
	Copyright 2007 Sphex LLP.
*/
package com.clevr.graphics {
	import flash.filters.ColorMatrixFilter;
	import flash.display.*;
	import flash.geom.*;
	
	
	/*
	*	Image histogram utilities.
	*
	*	@langversion ActionScript 3.0
	*	@playerversion Flash 9.0
	*
	*	@author Matt Kane
	*	@since  24.07.2007
	*/
	public class Histogram extends Object {
		
		public function Histogram(im:BitmapData){
			super();
			image = im;
			analyse();
		}
		
		public function analyse():void {
			hist = new Array(256);
			for (var i:int=0; i < 256; i++) {
				hist[i] = [0, 0, 0];
			}
			for (var x:int = 0; x < image.width; x++) {
				for (var y:int = 0; y < image.height; y++) {
					var pix:int	 = image.getPixel(x, y);
					var r:int = r(pix);
					var g:int = g(pix);
					var b:int = b(pix);
					hist[r][R]++;
					hist[g][G]++;
					hist[b][B]++;
				}
			}
		}
		
		public function stretchFilter(threshold:int, percent:Number=100):ColorMatrixFilter {
			var max:Array = getMax(threshold);
			var min:Array = getMin(threshold);
			
			var rMult:Number = 255 / (max[R] - min[R]) * (percent / 100);
			var gMult:Number = 255 / (max[G] - min[G]) * (percent / 100);
			var bMult:Number = 255 / (max[B] - min[B]) * (percent / 100);
			
			var mat:Array = [
				[rMult,	0,		0,		0, 	-min[R]],
				[0, 	gMult,	0,		0, 	-min[G]],
				[0,		0,		bMult,	0, 	-min[B]],
				[0,		0,		0,		1, 	0]
			];
			return new ColorMatrixFilter(mat);
		}
		
		public function getMin(threshold:int):Array {
			var ret:Array = [null, null, null];
			for (var i:int =0; i < 256; i++ ) {
				if (ret[R] == null && hist[i][R] > threshold) {
					ret[R] = i;
				}
				if (ret[G] == null && hist[i][G] > threshold) {
					ret[G] = i;
				}
				if (ret[B] == null && hist[i][B] > threshold) {
					ret[B] = i;
				}
				if (ret[R] != null && ret[G] != null && ret[B] != null) {
					break;
				}
			}
			return ret;
		}
		
		public function getMax(threshold:int):Array {
			var ret:Array = [null, null, null];
			for (var i:int = 255; i > -1; i--) {
				if (ret[R] == null && hist[i][R] > threshold) {
					ret[R] = i;
				}
				if (ret[G] == null && hist[i][G] > threshold) {
					ret[G] = i;
				}
				if (ret[B] == null && hist[i][B] > threshold) {
					ret[B] = i;
				}
				if (ret[R] != null && ret[G] != null && ret[B] != null) {
					break;
				}
			}
			return ret;
		}
		
		public function stretchTransform(threshold:int, percent:Number=100):ColorTransform {
			var max:Array = getMax(threshold);
			var min:Array = getMin(threshold);
			
			var rMult:Number = 255 / (max[R] - min[R]) * (percent / 100);
			var gMult:Number = 255 / (max[G] - min[G]) * (percent / 100);
			var bMult:Number = 255 / (max[B] - min[B]) * (percent / 100);
			
			return new ColorTransform(rMult, gMult, bMult, 1, -min[R], -min[G], -min[B], 0);
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
		
		private static var R:int = 0;
		private static var G:int = 1;
		private static var B:int = 2;
		private var image:BitmapData;
		private var hist:Array;
		
	}
	
}
