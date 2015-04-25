function processResult(data){
    if(data.status == 1){
	window.location.assign(data.page);
    } else {
	$("#hdr").html("Job " + data.id + " is being processed.");
	$("#statusdata").html(data.logdat);
    }
}

function pollStatus(url){$.PeriodicalUpdater(url, {method:'get', type: 'json', maxTimeout:60000}, processResult);}
