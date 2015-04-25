
var MSA;
var Fasta = require('biojs-io-fasta');


function buildAlignSelector(data, jid, prefix, baseurl){

    $("#tbldat").html('<table cellpadding="0" cellspacing="0" border="0" class="display" id="tabledat"></table>' );
    
    $("#tabledat").dataTable( { "order": [[2,'desc']], "pageLength": 5, "lengthChange": false, "data" : data,
				"columns":[{"title": "file id"},{"title": "File name"},{"title": "Has SNPs?"}],
				"columnDefs": [ { "targets": [0], "visible": false, "searchable": false}],
				"fnRowCallback": function( nRow, aData, iDisplayIndex, iDisplayIndexFull ) {
				    $('td:eq(0)', nRow).html(prefix + "_" + aData[1]);
				    $('td:eq(1)', nRow).html( ((aData[2] == 1)?"Yes":"No"));
				}
			      }
			    );
    $("#tabledat tbody").on('click', 'tr', function(e) {
	var aData = $("#tabledat").dataTable().fnGetData(this);
	key = jid + "/"+ aData[0] + "/";
	draw(key, baseurl);
    } );
}





function draw(part, baseurl){
    drawMSA(part, baseurl);
    drawTable(part, baseurl);
    $("#disinfo").val(part);
}


function compare(a, b){
    if(a.name == 'Consensus')
	return -1;
    return a.name.localeCompare(b.name);
}

function drawMSA(part, baseurl){
    $("#msaParent").html("<div id='msa'></div>");
    if( !(part) || ! $("input[name=clean]:checked").val() )
	return;
    var url = baseurl + "/ajax/msa/"+ part + $("input[name=clean]:checked").val() + "/";
    /* http://sniper.biojs-msa.org:9090/snippets/msa_show_menu */
     var opts = {};
    opts.el =  document.getElementById('msa');
    opts.el.textContent = "loading";
    opts.vis = {conserv: false, overviewbox: false, labelId : false};
    opts.zoomer = {boxRectHeight: 1,
                   boxRectWidth: 1,
                   labelWidth: 100,
                   labelIdLength: 0,
		   labelNameLength: 200,
                   alignmentWidth: 940,
                   alignmentHeight: 200,
                   residueFont: "16px mono monospace",
                   labelFontsize: "9px",
                   labelLineHeight: "30px",
                   markerFontsize: "16px",
                   columnWidth: 30,
                   rowHeight: 20
                  };
    
    Fasta.read(url, function(err, seqs) {
	seqs.forEach(function(seq){ seq.name = $.trim(seq.name);
				    //var idx = seq.name.indexOf(" ");
				    var idx = seq.name.indexOf("_");
				    if(idx != -1)
					seq.name = seq.name.substring(idx+1);
				    //seq.name = seq.name.substring(0, idx);
	                            else if(seq.name.length > 22)
					seq.name = seq.name.substring(0,22);
				  });
	//seqs.sort(compare);
	opts.zoomer['alignmentHeight'] = Math.min(20 * seqs.length, 200);
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


function drawTable(part, baseurl){
    $("#snpParent").html("");
    if($("input[name=clean]:checked").val() != 1)
	return;
    url = baseurl + "/ajax/snp/"+ part;
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
	    // DLS: race condition with the MSA creation.  Need to figure out how to create a dependency..
            //MSA.g.selcol.add( new  biojs.vis.msa.selection.columnsel({xStart: d.pos - 1, xEnd: d.pos - 1}) );
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