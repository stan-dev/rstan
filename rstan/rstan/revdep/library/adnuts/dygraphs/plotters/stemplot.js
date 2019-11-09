/**
 * stem plotter extracted from series.R
 */
function stemPlotter(e) { 
   var ctx = e.drawingContext; 
   var points = e.points; 
   var y_bottom = e.dygraph.toDomYCoord(0);
   ctx.fillStyle = e.color; 
   for (var i = 0; i < points.length; i++) { 
      var p = points[i]; 
      var center_x = p.canvasx;
      var center_y = p.canvasy; 
      ctx.beginPath(); 
      ctx.moveTo(center_x, y_bottom); 
      ctx.lineTo(center_x, center_y); 
      ctx.stroke();
      ctx.beginPath(); 
      ctx.arc(center_x, center_y, 3, 0, 2*Math.PI); 
      ctx.stroke();
   }
}
