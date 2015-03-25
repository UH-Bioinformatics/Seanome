function createChart(name, datajson, parent, dataid){
    var margin = {top: 20, right: 20, bottom: 30, left: 60},
    width = 900 - margin.left - margin.right,
    height = 300 - margin.top - margin.bottom,

    div = d3.select("#toolbox"),
    
    xt =  d3.extent(datajson, function(d, i ) {return d[0]; }),
    yt = d3.extent(datajson, function(d) { return (d[1] == 0)?1:d[1]; }),
    x= 0,
    y = 0,
    prettyprint = function(x) {return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");};

    if(xt[1] - xt[0] <= 100)
	x = d3.scale.linear().range([0, width]).domain(xt);
    else
	x = d3.scale.log().range([0, width]).domain(xt);

    if(yt[1] - yt[0] <= 100)
	y = d3.scale.linear().range([height, 0]).domain(yt);
    else
	y = d3.scale.log().range([height, 0]).domain(yt);

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
	.attr("transform", "translate(" + (width / 2) + " ," + (height + margin.bottom) + ")")
	.style("text-anchor", "middle")
	.text("Cluster size");

    svg.append("text")
	.attr("transform", "translate(" + (width / 2) + " ," + (margin.top) + ")")
	.style("text-anchor", "middle")
	.text( "Breakdown of cluster size for "+ name);

    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 0 - margin.left)
        .attr("x",0 - (height / 2))
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Number of clusters");

    svg.selectAll("dot")
        .data(datajson)
	.enter().append("circle")
        .attr("r", 2)
        .attr("cx", function(d) { return x(d[0]); })
        .attr("cy", function(d) { return y(d[1]); })
        .on("mouseover", function(d) {
            div.transition().duration(200).style("opacity", .9);
            div.html(prettyprint(d[2]) +" cluster(s) of size " + prettyprint(d[0]) +"<br/>"  + prettyprint(d[1]) + " cluster(s) of size >= "+ prettyprint(d[0]))
                .style("left", (d3.event.pageX) + "px").style("top", (d3.event.pageY - 28) + "px");
        })
        //.on("click", function(d){$(dataid).val(d[0]);})
        .on("mouseout", function(d) {div.transition().duration(500).style("opacity", 0);});
}
