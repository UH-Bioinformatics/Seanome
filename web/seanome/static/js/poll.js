function processResult(data){
    if(data.status == 1){
	window.location.assign(data.page);
    } else {
	$("#statusdata").html(data.logdat);
    }
}

function pollStatus(url){$.PeriodicalUpdater(url, {method:'get', type: 'json', maxTimeout:60000}, processResult);}
