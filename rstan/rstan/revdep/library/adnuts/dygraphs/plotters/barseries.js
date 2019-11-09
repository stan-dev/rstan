/**
 * Bar Chart plotter is adapted from http://dygraphs.com/tests/plotters.html
 */
function barSeriesPlotter(e) {
  var g = e.dygraph;
  var ctx = e.drawingContext;
  var points = e.points;
  var axis = g.attr_("axis", e.setName);
  var y_bottom = g.toDomYCoord(0, axis == "y2" ? 1 : 0);

  ctx.fillStyle = e.color;
  ctx.strokeStyle = e.color;

  // Find the minimum separation between x-values.
  // This determines the bar width.
  var min_sep = Infinity;
  for (var i = 1; i < points.length; i++) {
    var sep = points[i].canvasx - points[i - 1].canvasx;
    if (sep < min_sep) min_sep = sep;
  }
  var bar_width = Math.floor(2.0 / 3 * min_sep);

  // Do the actual plotting.
  for (var i = 0; i < points.length; i++) {
    var p = points[i];
    var center_x = p.canvasx;

    ctx.fillRect(center_x - bar_width / 2, p.canvasy,
      bar_width, y_bottom - p.canvasy);

    ctx.strokeRect(center_x - bar_width / 2, p.canvasy,
      bar_width, y_bottom - p.canvasy);
  }
}

