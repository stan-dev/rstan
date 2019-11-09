/**
 * Multi-column Bar Chart plotter is adapted from http://dygraphs.com/tests/plotters.html
 */

(function() {
  "use strict";

  // Multiple column bar chart
  function multiColumnBarPlotter(e) {
    // We need to handle all the series simultaneously.
    if (e.seriesIndex !== 0) return;

    var g = e.dygraph;
    var ctx = e.drawingContext;
    var sets = e.allSeriesPoints;
    var y_bottom = e.dygraph.toDomYCoord(0);

    // Find the minimum separation between x-values.
    // This determines the bar width.
    var min_sep = Infinity;
    for (var j = 0; j < sets.length; j++) {
      var points = sets[j];
      for (var i = 1; i < points.length; i++) {
        var sep = points[i].canvasx - points[i - 1].canvasx;
        if (sep < min_sep) min_sep = sep;
      }
    }
    var bar_width = Math.floor(2.0 / 3 * min_sep);

    var fillColors = [];
    var strokeColors = g.getColors();
    for (var i = 0; i < strokeColors.length; i++) {
      fillColors.push(strokeColors[i]);
    }

    for (var j = 0; j < sets.length; j++) {
      ctx.fillStyle = fillColors[j];
      ctx.strokeStyle = strokeColors[j];
      for (var i = 0; i < sets[j].length; i++) {
        var p = sets[j][i];
        var center_x = p.canvasx;
        var x_left = center_x - (bar_width / 2) * (1 - j/(sets.length-1));

        ctx.fillRect(x_left, p.canvasy,
          bar_width/sets.length, y_bottom - p.canvasy);

        ctx.strokeRect(x_left, p.canvasy,
          bar_width/sets.length, y_bottom - p.canvasy);
      }
    }
  }

  Dygraph.Plotters.MultiColumn = multiColumnBarPlotter;
})();
