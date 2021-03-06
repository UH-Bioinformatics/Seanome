importScripts("readfq.min.js");


self.addEventListener('message', function(e){
    var VALIDBASE = {A:'',C:'',G:'',T:'',N:'', a:'',c:'',g:'',t:'',n:''};
    var errored = false;
    var  cnt = 0, handler;
    var granularity = e.data[1];
    if(e.data[2] == 'file') {
        if (!FileReaderSync) {
            postMessage(['txt', 'ERROR: FileReaderSync is not supported']);
	    self.close();
	    return;
        }
        var reader = new FileReaderSync();//FileReader();
        handler = new fastqaFileInput(e.data[4], reader);
    } else {
        handler = new fastqaTextInput(e.data[4]);
    }
    
    var pct = 0, view = 1;
    while (true) {
	try{
	    results = nextRecord(handler);
	} catch (e){
	    postMessage(['txt', "Error: Invalid file format or a sequence is too long"]);
	    postMessage(['result', false]);
            break;
	}
        if (results[0] == -1) { // eof
	    postMessage(['%', 100.00]);
	    if(cnt != 0){
		postMessage(['txt', "Input is valid"]);
		postMessage(['result', true]);
            } else {
                postMessage(['txt', "Error: No sequences found"]);
 		postMessage(['result', false]);
	    }
            break;
        } else if (results[0] == -2) {
            postMessage(['txt', "Error: Invalid file format"]);
	    postMessage(['%', handler.complete()]);
	    postMessage(['result', false]);
            break;
        } else {
            var fastq = false;
            if (results[1]['qual'])
                fastq = true;
            if (fastq && e.data[3]['fastq'] == false) {
                errored = true;
                postMessage(['txt', "Error: fastq format is not accepted"]);
            }
	    if (fastq == false && e.data[3]['fasta'] == false){
                errored = true;
                postMessage(['txt', "Error: fasta format is not accepted"]);
	    }
	    if(errored){
		postMessage(['%', handler.complete()]);
		postMessage(['result', false]);
		break;
	    }
	    
	    var valid = true;
	    
	    for(var i =0; i < results[1]['seq'].length && valid; ++i){
		var d = results[1]['seq'][i];
		valid = valid && d in VALIDBASE;
            }
	    if(valid == false || results[1]['seq'].length == 0){
                postMessage(['txt', "Error: Invalid file format"]);
		postMessage(['%', handler.complete()]);
		postMessage(['result', false]);
		break;
	    }
            ++cnt;
            pct = handler.complete();
            if( granularity == 0 || (pct / granularity) >= view ) {
                ++ view;
                postMessage(['%', pct.toFixed(2)]);
            }
        }
    }
    self.close();
}, false);
