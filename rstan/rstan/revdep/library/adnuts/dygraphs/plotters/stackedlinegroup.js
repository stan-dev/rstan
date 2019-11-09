// the _lineplotter function extracted from dygraphs-combined-dev.js, available for use in conjunction with other per-series plotters and group plotters

function linePlotter(e) {
  //extracting and reducing the Dygraph.stackPoints_ function
  stackPoints = function(points, cumulativeYval, seriesExtremes, fillMethod) {
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
  
  // BEGIN HEADER BLOCK
  // This first block can be copied to other plotters to capture the group 
  var g = e.dygraph;
  
  var group;
  var groupIdx = [];
  var sets = [];
  var allSets = e.allSeriesPoints;
  var minIdx = Infinity;
  var setName = e.setName;
  var setNames = g.getLabels().slice(1);
  var groupSetNames = [];
  var fillColors = [];
  var strokeColors = g.getColors();
  // this next one we use further down, but will be populated in a decreasing loop,
  // so we'll establish the size in this forward loop so it has the structure to accept
  // later on.
  var seriesExtremes = [];
  
  var currGroup = g.attr_("group", setName);
  
  for (var setIdx = 0; setIdx < allSets.length; setIdx++) {
    // get the name and group of the current setIdx
    setName = setNames[setIdx];
    group = g.attr_("group", setName);

    if (group === currGroup) {
      //save the indv index and the points
      groupIdx.push(setIdx);
      sets.push(allSets[setIdx]);
      groupSetNames.push(setName);
      fillColors.push(strokeColors[setIdx]);
     
      // the aforementioned stuff for later on 
      var tmpExtremes = [];
      tmpExtremes[0] = Infinity;
      tmpExtremes[1] = -Infinity;
      
      seriesExtremes.push(tmpExtremes);
      
      // capturing the min indx helps to ensure we don't render the plotter
      // multiple times
      if (setIdx < minIdx) minIdx = setIdx;
    }
  }
  // END HEADER BLOCK
	

	//Stack the points
  // set up cumulative records
  var cumulativeYval = [];
  var packed = g.gatherDatasets_(g.rolledSeries_, null);
  var extremes = packed.extremes;
  var seriesName;
  
  for (var j = sets.length - 1; j >= 0; j--) {
    points = sets[j];
    seriesName = groupSetNames[j]; 
    
    //  stack the data 
    stackPoints(points, cumulativeYval, seriesExtremes[j],
          g.getBooleanOption("stackedGraphNaNFill"));
    
    extremes[seriesName] = seriesExtremes[j];
  }
  
	// Do the actual plotting.
	for (var j = 0; j < sets.length; j++) {
    setName = groupSetNames[j];
		if(setName !== e.setName) continue;
    
		var connectSeparated = g.getOption('connectSeparatedPoints', setName);
    var logscale = g.attributes_.getForSeries("logscale", setName);
    
    axis = g.axisPropertiesForSeries(setName);
    
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
		}

		e.points = points;

	  var strokeWidth = e.strokeWidth;
	
	  var borderWidth = g.getNumericOption("strokeBorderWidth", setName);
	  var drawPointCallback = g.getOption("drawPointCallback", setName) ||
	      Dygraph.Circles.DEFAULT;
	  var strokePattern = g.getOption("strokePattern", setName);
	  var drawPoints = g.getBooleanOption("drawPoints", setName);
	  var pointSize = g.getNumericOption("pointSize", setName);
	
	  if (borderWidth && strokeWidth) {
	    DygraphCanvasRenderer._drawStyledLine(e,
	        g.getOption("strokeBorderColor", setName),
	        strokeWidth + 2 * borderWidth,
	        strokePattern,
	        drawPoints,
	        drawPointCallback,
	        pointSize
	        );
	  }
	
	  DygraphCanvasRenderer._drawStyledLine(e,
	      e.color,
	      strokeWidth,
	      strokePattern,
	      drawPoints,
	      drawPointCallback,
	      pointSize
	  );
	}
}
