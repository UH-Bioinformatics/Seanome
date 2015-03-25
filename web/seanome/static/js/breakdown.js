function setcovVals(parent, child, values){
    var c = parseInt($(parent).val());
    if(c in values){$(child).val(values[c]);}else{$(child).val("???");}
}

function ondownload(url, cutoff, mincov, maxcov){
    var cval = parseInt($(cutoff).val());
    maxcov = parseInt(maxcov);
    mincov = parseInt(mincov);
    if (cval > maxcov || cval < mincov){alert("Coverage must be between: " + mincov + " and " + maxcov);return;}
    url += "?cutoff=" + $(cutoff).val();window.location.assign(url);
}




function createChart(datajson, parent, selectbx){
    var margin = {top: 20, right: 20, bottom: 50, left: 60},
    width = 900 - margin.left - margin.right,
    height = 300 - margin.top - margin.bottom,
    div = d3.select("#toolbox"),
    xe = d3.extent(datajson, function(d, i ) {return d[0]; }),
    ye = d3.extent(datajson, function(d) { return (d[1] == 0)?1:d[1]; }),
    x = "", y = "";
    
    if (xe[1] - xe[0] > 100)
	x = d3.scale.log().range([0, width]).domain(xe);
    else
	x = d3.scale.linear().range([0, width]).domain(xe);
    //if(ye[1] - ye[0] > 100)
    //y = d3.scale.log().range([height, 0]).domain(d3.extent(datajson, function(d) { return (d[1] == 0)?1:d[1]; }));
    //else
    y = d3.scale.linear().range([height, 0]).domain(d3.extent(datajson, function(d) { return (d[1] == 0)?1:d[1]; }));

    prettyprint = function(x) {return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");};
    
    var xAxis = d3.svg.axis()
	.scale(x)
	.orient("bottom"),
    yAxis = d3.svg.axis()
	.scale(y)
	.orient("left"),
    line = d3.svg.line()
	.x(function(d, i ) {return x(d[0]); })
	.y(function(d) { return y((d[1]==0)?1:d[1]); }),
    svg = d3.select(parent).attr("width", width + margin.left + margin.right)
	.attr("height", height + margin.top + margin.bottom)
	.append("g")
	.attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    svg.append("g")
	.attr("class", "x axis")
	.attr("transform", "translate(0," + height + ")")
	.call(xAxis);
    svg.append("g")
	.attr("class", "y axis")
	.call(yAxis);
    svg.append("path")
	.datum(datajson)
	.attr("class", "line")
	.attr("d", line)

    svg.append("text")
	.attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom - 8) + ")")
	.style("text-anchor", "middle")
	.text("Coverage");
/*
    svg.append("text")
	.attr("transform", "translate(" + (width / 2) + " ," + (margin.top) + ")")
	.style("text-anchor", "middle")
	.text( "Minimum average coverage");
*/
    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0 - margin.left)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Number of CSRs");

    svg.selectAll("dot")
        .data(datajson)
	.enter().append("circle")
        .attr("r", 2)
        .attr("cx", function(d) { return x(d[0]); })
        .attr("cy", function(d) { return y(d[1]); })
        .on("mouseover", function(d) {
            div.transition().duration(200).style("opacity", .9);
            div.html("Minimum coverage: " + prettyprint(d[0]) +"<br/>"  + "Total number of CSRs: " + prettyprint(d[1]))
                .style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
        })
	.on("click", function(d){
	    $(selectbx).select2('val',d[0]);
	    $(selectbx).trigger('change');
	})

        .on("mouseout", function(d) {div.transition().duration(500).style("opacity", 0);});

}
