/**
 * The Candle chart plotter is adapted from code written by
 * Zhenlei Cai (jpenguin@gmail.com)
 * https://github.com/danvk/dygraphs/pull/141/files
 * 
 * Adapted for use by dyGroup in the dygraphs R package
 * 
 */
  function candlestickgroupPlotter(e) {
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
    
    var currGroup = g.attr_("group", setName);
    
    for (var setIdx = 0; setIdx < allSets.length; setIdx++) {
      // get the name and group of the current setIdx
      setName = setNames[setIdx];
      group = g.attr_("group", setName);
  
      if (group === currGroup) {
        //save the indv index and the points
        groupIdx.push(setIdx);
        sets.push(allSets[setIdx]);
        
        // capturing the min indx helps to ensure we don't render the plotter
        // multiple times
        if (setIdx < minIdx) minIdx = setIdx;
      }
    }
    
    // We'll employ the plotter only on the first of the group
    if (e.seriesIndex !== minIdx) return;
    // END HEADER BLOCK
   
    // If the group doesn't have four series then revert to prior 
    if (groupIdx.length < 4) return;
  
    function getPrices(sets) {
      var prices = [];
      var price;
      for (var p = 0 ; p < sets[0].length; p++) {
        price = {
          open : sets[0][p].yval,
          high : sets[1][p].yval,
          low : sets[2][p].yval,
          close : sets[3][p].yval,
          openY : sets[0][p].y,
          highY : sets[1][p].y,
          lowY : sets[2][p].y,
          closeY : sets[3][p].y
        };
        prices.push(price);
      }
      return prices;
    }

    var prices = getPrices(sets);
    var area = e.plotArea;
    var ctx = e.drawingContext;
    ctx.strokeStyle = '#202020';
    ctx.lineWidth = 0.6;

    var minBarWidth = 2;
    var numBars = prices.length + 1; // To compensate the probably removed first "incomplete" bar
    var barWidth = Math.round((area.w / numBars) / 2);
    if (barWidth % 2 !== 0) {
      barWidth++;
    }
    barWidth = Math.max(barWidth, minBarWidth);

    var price;
    for (var p = 0 ; p < prices.length; p++) {
      ctx.beginPath();

      price = prices[p];
      var topY = area.h * price.highY + area.y;
      var bottomY = area.h * price.lowY + area.y;
      var centerX = Math.floor(area.x + sets[0][p].x * area.w) + 0.5; // crisper rendering
      ctx.moveTo(centerX, topY);
      ctx.lineTo(centerX, bottomY);
      ctx.closePath();
      ctx.stroke();
      var bodyY;
      if (price.open > price.close) {
        ctx.fillStyle ='#d9534f';
        bodyY = area.h * price.openY + area.y;
      }
      else {
        ctx.fillStyle ='#5cb85c';
        bodyY = area.h * price.closeY  + area.y;
      }
      var bodyHeight = area.h * Math.abs(price.openY - price.closeY);
      ctx.fillRect(centerX - barWidth / 2, bodyY, barWidth,  bodyHeight);
   }
}

