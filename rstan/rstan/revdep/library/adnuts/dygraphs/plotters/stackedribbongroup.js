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
  
	// Helper function to trace a line back along the baseline.
  var traceBackPath = function(ctx, baselineX, baselineY, pathBack) {
    ctx.lineTo(baselineX, baselineY);
    for (var i = pathBack.length - 1; i >= 0; i--) {
      var pt = pathBack[i];
      ctx.lineTo(pt[0], pt[1]);
    }
  };
  
	// Do the actual plotting.
				// First, we'll plot the line for this series, then...
				// Second, we'll add the fills
				// In contrast to stackedlinegroup, we do this in reverse
				// order to align with the fillplotter
	var area = e.plotArea;
  var fillAlpha = g.getNumericOption('fillAlpha');
	
	var baseline = {};
  var currBaseline;
  var prevStepPlot;  // for different line drawing modes (line/step) per series

  var ctx = e.drawingContext;
    
	// For filled charts, we draw points from left to right, then back along
  // the x-axis to complete a shape for filling.
  // For stacked plots, this "back path" is a more complex shape. This array
  // stores the [x, y] values needed to trace that shape.
  var pathBack = [];

	//We'll save the group indices of the current set,
				// so as to test later and hopefully skip 
				// past needless iterations of the loops
	var currSetIdx;
	//Now a quick FOR loop to capture the group indices
  for (var j = sets.length - 1; j >= 0; j--) {
    seriesName = groupSetNames[j]; 
		if (seriesName === e.setName) currSetIdx = j;
  }
  
	for (var j = sets.length - 1; j >= 0; j--) {
		// If we're not dealing with the immediate plotted series or it's
					// immediate predecesor, skip all this stuff below
		if (j > (currSetIdx + 1) || j < currSetIdx) continue;

    setName = groupSetNames[j];

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
		
		// We'll only draw the line plotter if the current series matches
					// the one in the top-level plotter call
		if (setName === e.setName) {
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
		//END OF THE LINE PLOTTER
					//
		//BEGIN THE FILL
    var stepPlot = g.getBooleanOption('stepPlot', setName);
		var color = e.color;
		var axis = g.axisPropertiesForSeries(setName);
    var axisY = 1.0 + axis.minyval * axis.yscale;
    if (axisY < 0.0) axisY = 0.0;
    else if (axisY > 1.0) axisY = 1.0;
    axisY = area.h * axisY + area.y;


		// setup graphics context
  	var prevX = NaN;
  	var prevYs = [-1, -1];
  	var newYs;
  	
		// should be same color as the lines but only 15% opaque.
  	var rgb = Dygraph.toRGB_(color);
  	var err_color =
  	    'rgba(' + rgb.r + ',' + rgb.g + ',' + rgb.b + ',' + fillAlpha + ')';
    

    var iter = Dygraph.createIterator(points, 0, points.length,
        DygraphCanvasRenderer._getIteratorPredicate(
            g.getBooleanOption("connectSeparatedPoints", setName)));

    ctx.fillStyle = err_color;
    ctx.beginPath();
    var last_x, is_first = true;

    // If the point density is high enough, dropping segments on their way to
    // the canvas justifies the overhead of doing so.
    if (points.length > 2 * g.width_ || Dygraph.FORCE_FAST_PROXY) {
      ctx = DygraphCanvasRenderer._fastCanvasProxy(ctx);
    }

    // TODO(danvk): there are a lot of options at play in this loop.
    //     The logic would be much clearer if some (e.g. stackGraph and
    //     stepPlot) were split off into separate sub-plotters.
    var point;

		// Throughout this loop, we test to see if the top-level called series
					// matches the one current in the parent FOR loop.  If not, we skip
					// the parts that draw the line but leave the others so the pathback
					// and prevYs and prevXs still get properly populated
    while (iter.hasNext) {
      point = iter.next();
      if (!Dygraph.isOK(point.y) && !stepPlot) {
				//
        if (e.setName === setName) traceBackPath(
								ctx, prevX, prevYs[1], pathBack);
        pathBack = [];
        prevX = NaN;
        if (point.y_stacked !== null && !isNaN(point.y_stacked)) {
          baseline[point.canvasx] = area.h * point.y_stacked + area.y;
        }
        continue;
      }
      if (!is_first && last_x == point.xval) {
        continue;
      } else {
        is_first = false;
        last_x = point.xval;
      }

      currBaseline = baseline[point.canvasx];
      var lastY;
      if (currBaseline === undefined) {
        lastY = axisY;
      } else {
        if(prevStepPlot) {
          lastY = currBaseline[0];
        } else {
          lastY = currBaseline;
        }
      }
      newYs = [ point.canvasy, lastY ];

      if (stepPlot) {
        // Step plots must keep track of the top and bottom of
        // the baseline at each point.
        if (prevYs[0] === -1) {
          baseline[point.canvasx] = [ point.canvasy, axisY ];
        } else {
          baseline[point.canvasx] = [ point.canvasy, prevYs[0] ];
        }
      } else {
        baseline[point.canvasx] = point.canvasy;
      }

      if (!isNaN(prevX)) {
				if (e.setName === setName) {
        	// Move to top fill point
        	if (stepPlot) {
        	  ctx.lineTo(point.canvasx, prevYs[0]);
        	  ctx.lineTo(point.canvasx, newYs[0]);
        	} else {
        	  ctx.lineTo(point.canvasx, newYs[0]);
        	}
        }

        // Record the baseline for the reverse path.
        pathBack.push([prevX, prevYs[1]]);
        if (prevStepPlot && currBaseline) {
          // Draw to the bottom of the baseline
          pathBack.push([point.canvasx, currBaseline[1]]);
        } else {
          pathBack.push([point.canvasx, newYs[1]]);
        }
      } else if (e.setName === setName) {
        ctx.moveTo(point.canvasx, newYs[1]);
        ctx.lineTo(point.canvasx, newYs[0]);
      }
      prevYs = newYs;
      prevX = point.canvasx;
		}
    prevStepPlot = stepPlot;
    if (newYs && point) {
      if (e.setName === setName) traceBackPath(
							ctx, point.canvasx, newYs[1], pathBack);
      pathBack = [];
    }
    ctx.fill();
	}
}
