var MSA;

function draw(baseurl){
    drawMSA(baseurl);
    drawTable(baseurl);
}


function drawMSA(baseurl){
    $("#msaParent").html("<div id='msa'></div>");
    if((!$("#aln").val()) || (!$("input[name=clean]:checked").val()) )
	return;
    var url = baseurl + "/ajax/msa/"+$("#aln").val() + $("input[name=clean]:checked").val() + "/";
    /* http://sniper.biojs-msa.org:9090/snippets/msa_show_menu */
    var opts = {};
    opts.el =  document.getElementById('msa');
    opts.el.textContent = "loading";
    opts.vis = {conserv: false, overviewbox: false, labelId : false};
    opts.zoomer = {boxRectHeight: 1,
                   boxRectWidth: 1,
                   labelWidth: 210,
                   labelIdLength: 0,
                   alignmentWidth: 940,
                   alignmentHeight: 330,
                   residueFont: "18px mono monospace",
                   labelFontsize: "9px",
                   labelLineHeight: "30px",
                   markerFontsize: "16px",
                   columnWidth: 30,
                   rowHeight: 20
                  };
    biojs.io.fasta.parse.read(url, function(seqs) {
	seqs.forEach(function(seq){ seq.name = $.trim(seq.name);
				    var idx = seq.name.indexOf(" ");
				    if(idx != -1)
					seq.name = seq.name.substring(0, idx);
	                            else if(seq.name.length > 22)
					seq.name = seq.name.substring(0,22);
				  });
	opts.seqs = seqs;
        opts.speed = true;
	opts.keyevents = true;
	opts.visibleElements = {ruler: false};
	MSA = new biojs.vis.msa.msa(opts);
	MSA.el.style.marginTop = "10px";
        MSA.render();
	MSA.g.off('residue:click');
        MSA.g.off('row:click');
        MSA.g.off('column:click');
    });
}


function gotopos(offset){MSA.g.zoomer.setLeftOffset(offset)}


function drawTable(baseurl){
    $("#snpParent").html("");
    if($("input[name=clean]:checked").val() != 1)
	return;
    url = baseurl + "/ajax/snp/"+$("#aln").val();
    $.get(url, function(data){
	if(data.length == 0){ return;}
	var keys = Object.keys(data[0].counts),
	hdr  = "<tr><th>Position</th><th>Reference</th><th>Alternate(s)</th>";
	html = "<table id='snptbl' class ='table table-striped table-bordered'>";
	keys.forEach(function(k){
	    hdr += "<th><div class='row' style='width:300px'><span class='col-md-12 text-center'>" + 
		k + "</span></div><div class='row'><div class='col-md-6 text-center'>Allele Freq.</div>" + 
		"<div class='col-md-6 text-center'>Counts</div></div></th>";
	});
	hdr += "</tr>";
	html += "<thead>" + hdr + "</thead><tfoot>"+ hdr +"</tfoot><tbody>";
	data.forEach(function(d){
	    html += "<tr><td style='cursor:pointer;' onclick='gotopos(" +'"' + (d.pos - 1) + '"' + ");'>" + d.pos + "</td><td>" + d.ref + "</td><td>" + d.alts.join(", ") + "</td>";
            MSA.g.selcol.add( new  biojs.vis.msa.selection.columnsel({xStart: d.pos - 1, xEnd: d.pos - 1}) );
	    keys.forEach(function(k){
		var tmpcnt = [];
		var total = 0
		/* sum the counts and create a string of the counts that we can display*/ 
		for(b in d.counts[k]){
		    tmpcnt.push( b + ": " + d.counts[k][b] );
		    total += d.counts[k][b]
		}
		var tmp = [];
		/*Using the total count, we can now compute a pct, for the frequency*/
		for(b in d.counts[k]){tmp.push(b + ": " + ( (total != 0) ? (d.counts[k][b] / total).toFixed(3) : "0.000" ) );}
		html += "<td><div class='row' style='width:300px'><div class='col-md-6 text-center'>" + 
		    tmp.join("<br />") + "</div><div class='col-md-6 text-center'>" + 
		    tmpcnt.join("<br />") + "</div></div></td>";
	    });
	    html += "</tr>";
	});
	html += "</tbody></table>";
	$("#snpParent").html(html);
	$('#snptbl').dataTable({"dom": 'Rlfrtip', "iDisplayLength": 10, "scrollX": true});
    }, 'json');


}