// the _fillplotter combined with the _lineplotter function extracted from dygraphs-combined-dev.js, available for use in conjunction with other per-series plotters and group plotters

function filledlineplotter(e) {
  var g = e.dygraph;
  var setName = e.setName;

  var area = e.plotArea;
  var sets = e.allSeriesPoints;
  var setCount = sets.length;

  var fillAlpha = g.getNumericOption('fillAlpha');
  var stackedGraph = g.getBooleanOption("stackedGraph");
  var colors = g.getColors();

  // For stacked graphs, track the baseline for filling.
  //
  // The filled areas below graph lines are trapezoids with two
  // vertical edges. The top edge is the line segment being drawn, and
  // the baseline is the bottom edge. Each baseline corresponds to the
  // top line segment from the previous stacked line. In the case of
  // step plots, the trapezoids are rectangles.
  var baseline = {};
  var currBaseline;
  var prevStepPlot;  // for different line drawing modes (line/step) per series

  // Helper function to trace a line back along the baseline.
  var traceBackPath = function(ctx, baselineX, baselineY, pathBack) {
    ctx.lineTo(baselineX, baselineY);
    if (stackedGraph) {
      for (var i = pathBack.length - 1; i >= 0; i--) {
        var pt = pathBack[i];
        ctx.lineTo(pt[0], pt[1]);
      }
    }
  };

  var ctx = e.drawingContext;

  var stepPlot = g.getBooleanOption('stepPlot', setName);
  var color = e.color
  var axis = g.axisPropertiesForSeries(setName);
  var axisY = 1.0 + axis.minyval * axis.yscale;
  if (axisY < 0.0) axisY = 0.0;
  else if (axisY > 1.0) axisY = 1.0;
  axisY = area.h * axisY + area.y;

  var points = e.points;
  var iter = Dygraph.createIterator(points, 0, points.length,
      DygraphCanvasRenderer._getIteratorPredicate(
          g.getBooleanOption("connectSeparatedPoints", setName)));

  // setup graphics context
  var prevX = NaN;
  var prevYs = [-1, -1];
  var newYs;
  // should be same color as the lines but only 15% opaque.
  var rgb = Dygraph.toRGB_(color);
  var err_color =
      'rgba(' + rgb.r + ',' + rgb.g + ',' + rgb.b + ',' + fillAlpha + ')';
  ctx.fillStyle = err_color;
  ctx.beginPath();
  var last_x, is_first = true;

  // If the point density is high enough, dropping segments on their way to
  // the canvas justifies the overhead of doing so.
  if (points.length > 2 * g.width_ || Dygraph.FORCE_FAST_PROXY) {
    ctx = DygraphCanvasRenderer._fastCanvasProxy(ctx);
  }

  // For filled charts, we draw points from left to right, then back along
  // the x-axis to complete a shape for filling.
  // For stacked plots, this "back path" is a more complex shape. This array
  // stores the [x, y] values needed to trace that shape.
  var pathBack = [];

  // TODO(danvk): there are a lot of options at play in this loop.
  //     The logic would be much clearer if some (e.g. stackGraph and
  //     stepPlot) were split off into separate sub-plotters.
  var point;
  while (iter.hasNext) {
    point = iter.next();
    if (!Dygraph.isOK(point.y) && !stepPlot) {
      traceBackPath(ctx, prevX, prevYs[1], pathBack);
      pathBack = [];
      prevX = NaN;
      if (point.y_stacked !== null && !isNaN(point.y_stacked)) {
        baseline[point.canvasx] = area.h * point.y_stacked + area.y;
      }
      continue;
    }
    
      if (isNaN(point.canvasy) && stepPlot) {
        newYs = [ area.y + area.h, axisY ];
      } else {
        newYs = [ point.canvasy, axisY ];
      }
      
      if (!isNaN(prevX)) {
        // Move to top fill point
        if (stepPlot) {
          ctx.lineTo(point.canvasx, prevYs[0]);
          ctx.lineTo(point.canvasx, newYs[0]);
        } else {
          ctx.lineTo(point.canvasx, newYs[0]);
        }

      } else {
        ctx.moveTo(point.canvasx, newYs[1]);
        ctx.lineTo(point.canvasx, newYs[0]);
      }
      
      prevYs = newYs;
      prevX = point.canvasx;
      
    }
    
    prevStepPlot = stepPlot;
    
    if (newYs && point) {
      traceBackPath(ctx, point.canvasx, newYs[1], pathBack);
      pathBack = [];
    }
    
    ctx.fill();
    
}
