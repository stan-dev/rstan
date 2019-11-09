/**
 * Bar Chart plotter is adapted from http://dygraphs.com/tests/plotters.html
 */
(function() {
  "use strict";

	function stackedBarChartPlotter(e) {
    // We need to handle all the series simultaneously.
    if (e.seriesIndex !== 0) return;

    var g = e.dygraph;
    var ctx = e.drawingContext;
    var sets = e.allSeriesPoints;
    var y_bottom = e.dygraph.toDomYCoord(0);
  
					
		//extracting and reducing the Dygraph.stackPoints_ function
  	function stackPoints(points, cumulativeYval, seriesExtremes, fillMethod) {
  	  var lastXval = null;
  	  var prevPoint = null;
  	  var nextPoint = null;
  	  var nextPointIdx = -1;
  	
  	  // Find the next stackable point starting from the given index.
  	  var updateNextPoint = function(idx) {
  	    // If we've previously found a non-NaN point and haven't gone past it yet,
  	    // just use that.
  	    if (nextPointIdx >= idx) return;
  	
  	    // We haven't found a non-NaN point yet or have moved past it,
  	    // look towards the right to find a non-NaN point.
  	    for (var j = idx; j < points.length; ++j) {
  	      // Clear out a previously-found point (if any) since it's no longer
  	      // valid, we shouldn't use it for interpolation anymore.
  	      nextPoint = null;
  	      if (!isNaN(points[j].yval) && points[j].yval !== null) {
  	        nextPointIdx = j;
  	        nextPoint = points[j];
  	        break;
  	      }
  	    }
  	  };
  	
  	  for (var i = 0; i < points.length; ++i) {
  	    var point = points[i];
  	    var xval = point.xval;
  	    if (cumulativeYval[xval] === undefined) {
  	      cumulativeYval[xval] = 0;
  	    }
  	
  	    var actualYval = point.yval;
  	    if (isNaN(actualYval) || actualYval === null) {
  	      if(fillMethod == 'none') {
  	        actualYval = 0;
  	      } else {
  	        // Interpolate/extend for stacking purposes if possible.
  	        updateNextPoint(i);
  	        if (prevPoint && nextPoint && fillMethod != 'none') {
  	          // Use linear interpolation between prevPoint and nextPoint.
  	          actualYval = prevPoint.yval + (nextPoint.yval - prevPoint.yval) *
  	              ((xval - prevPoint.xval) / (nextPoint.xval - prevPoint.xval));
  	        } else if (prevPoint && fillMethod == 'all') {
  	          actualYval = prevPoint.yval;
  	        } else if (nextPoint && fillMethod == 'all') {
  	          actualYval = nextPoint.yval;
  	        } else {
  	          actualYval = 0;
  	        }
  	      }
  	    } else {
  	      prevPoint = point;
  	    }
  	
  	    var stackedYval = cumulativeYval[xval];
  	    if (lastXval != xval) {
  	      // If an x-value is repeated, we ignore the duplicates.
  	      stackedYval += actualYval;
  	      cumulativeYval[xval] = stackedYval;
  	    }
  	    lastXval = xval;
  	
  	    point.yval_stacked = stackedYval;
  	    
  	    if (stackedYval > seriesExtremes[1]) {
  	      seriesExtremes[1] = stackedYval;
  	    }
  	    if (stackedYval < seriesExtremes[0]) {
  	      seriesExtremes[0] = stackedYval;
  	    }
  	
  	  }
  	};
  	
  	var setNames = g.getLabels().slice(1);  // remove x-axis
  	
		var points = e.points;
  	var sets = e.allSeriesPoints;
  	var minIdx = Infinity;
    
		var fillColors = [];
    var strokeColors = g.getColors();
    for (var i = 0; i < strokeColors.length; i++) {
      fillColors.push(strokeColors[i]);
    }
  	
		var seriesExtremes = [];
    
		var tmpExtremes = [];
    tmpExtremes[0] = Infinity;
    tmpExtremes[1] = -Infinity;
  	for (var j = 0; j < sets.length; j++) {
    	seriesExtremes.push(tmpExtremes);
		}
     
  	// Find the minimum separation between x-values.
  	// This determines the bar width.
  	var points = sets[0];
		
		var min_sep = Infinity;
  	for (var i = 1; i < points.length; i++) {
  	  var sep = points[i].canvasx - points[i - 1].canvasx;
  	  if (sep < min_sep) min_sep = sep;
  	}
  	var bar_width = Math.floor(2.0 / 3 * min_sep);

  	// set up cumulative records
  	var cumulativeYval = [];
  	var packed = g.gatherDatasets_(g.rolledSeries_, null);
  	var extremes = packed.extremes;
  	var seriesName;
  	
  	for (var j = sets.length - 1; j >= 0; j--) {
  	  
			points = sets[j];
  	  seriesName = setNames[j];
  	  
  	  //  stack the data 
  	  stackPoints(points, cumulativeYval, seriesExtremes[j],
  	        g.getBooleanOption("stackedGraphNaNFill"));
  	  
  	  extremes[seriesName] = seriesExtremes[j];
  	}
  	 
  	// There is currently no way to update the axes height from inside the plotter...
  	// Will have to wait until update can be made to underlying dygraphs lib
  	// Preferring to do issue or pull request to main library on github instead of modifying here
  	// g.computeYAxisRanges_(extremes);
  	// g.layout_.setYAxes(g.axes_);
  	var axis;
  	var logscale;
  	var connectSeparated;
  	  
  	// Do the actual plotting.
  	for (var j = 0; j < sets.length; j++) {
  	  seriesName = setNames[j];
  	  connectSeparated = g.getOption('connectSeparatedPoints', seriesName);
  	  logscale = g.attributes_.getForSeries("logscale", seriesName);
  	  
  	  axis = g.axisPropertiesForSeries(seriesName);
  	  
  	  points = sets[j];
  	  
  	  for (var i = 0; i < points.length; i++) {
  	    var point = points[i];
  	    
  	    var yval = point.yval;
  	    
  	    point.y_stacked = DygraphLayout.calcYNormal_(
  	        axis, point.yval_stacked, logscale);
  	        
  	    if (yval !== null && !isNaN(yval)) {
  	      yval = point.yval_stacked;
  	    }
  	    if (yval === null) {
  	      yval = NaN;
  	      if (!connectSeparated) {
  	        point.yval = NaN;
  	      }
  	    }
  	    point.y = DygraphLayout.calcYNormal_(axis, yval, logscale);
  	  
  	    point.canvasx = g.plotter_.area.w * point.x + g.plotter_.area.x;
  	    point.canvasy = g.plotter_.area.h * point.y + g.plotter_.area.y;
  	    
  	    var center_x = point.canvasx;
  	    
  	    ctx.fillStyle = fillColors[j];
  	    ctx.strokeStyle = fillColors[j];
  	  
  	    ctx.fillRect(center_x - bar_width / 2, point.canvasy,
  	      bar_width, y_bottom - point.canvasy);
  	  
  	    ctx.strokeRect(center_x - bar_width / 2, point.canvasy,
  	      bar_width, y_bottom - point.canvasy);
  	  }
  	}
	}

  Dygraph.Plotters.StackedBarChart = stackedBarChartPlotter;

})();
