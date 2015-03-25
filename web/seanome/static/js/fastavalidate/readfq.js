/*
This code loosely follows the C version found at:
https://github.com/lh3/readfq
*/


var SEP_SPACE = 0 // isspace(): \t, \n, \v, \f, \r
    ,SEP_TAB = 1 // isspace() && !' '
    //,SPACE_REGEXP = /\s/
    //,TAB_REGEXP =/[\f\n\r\t\v?\u00a0\u1680?\u180e\u2000?\u2001\u2002?\u2003\u2004?\u2005\u2006?\u2007\u2008?\u2009\u200a?\u2028\u2029??\u202f\u205f?\u3000]/;
    ,SPACE_REGEXP = /\s/
    ,TAB_REGEXP =/[\f\n\r\t\v]/
    ,BUFFER_MAX = 524288;//131072;//32768;//8192;//4096;//32;

// helper functions
function findregexp(str, delimiter , start) {
    var cnt = (start == 0) ? str.search(delimiter) : str.substr(start).search(delimiter);
    return (cnt != -1) ? start + cnt : cnt;
}

function findIdx(str, delimiter, start){
	return str.indexOf(delimiter, start);
}
/*
function  isspace(c){
    //return SPACE_REGEXP.test(c);
    return (c == ' ' || c == '\n' || c == '\t' ||  c == '\r' || c =='\v' ||  c == '\f');
}
*/
// end helper functions



// parent class for the parser.
function fastqaParse(){
    this.prev = 0;
    this.pos = 0;
    this.eof = false;
    this.last_char = 0;
    this.start = 0;
}


fastqaParse.prototype.read = function(count) {
    this.eof = true;
    this.start = 0;
    this.end = 0;
    return "";
}

// compute the % of the data that has been processed
fastqaParse.prototype.complete = function() {
    return ( (this.prev + this.start) / this.size) * 100.0;
}


fastqaParse.prototype.getc = function(){
    // emulate a getc call as seen in C, but with some additions
    if(this.start >= this.end){
        this.buffer = this.read(this.buffersize);
        if(this.end == 0)
            return "";
    }
    return this.buffer.charAt(this.start++);
};


// scan the buffer until we find a delimiter of interest. Once found, we concat
// all the data from the initial start point until the delimiter.
fastqaParse.prototype.getuntil = function(delimiter, buf) {
    if (this.start >= this.end && this.eof)
        return [false, "", ''];

    var finder, i;
    // remove the branch from the infinite loop
    if ((typeof delimiter != "number")) {
        finder = findIdx;
    } else {
        finder = findregexp;
        if (delimiter == SEP_SPACE)
            delimiter = SPACE_REGEXP;
        else
            delimiter = TAB_REGEXP;
    }

    for(;;){
		// can we somehow reduce the for loop body further?
   		if (this.start >= this.end) {
            this.buffer = this.read(this.buffersize);
            if (this.end == 0)
                break;
        }
        if ( (i = finder(this.buffer, delimiter, this.start)) != -1) {
            buf[buf.length] = this.buffer.substring(this.start, i);
            this.start = i + 1;
            if (i < this.end)
                break;
        } else {
            buf[buf.length] = (this.start == 0) ? this.buffer : this.buffer.substr(this.start);
            this.start = this.end;
        }
    }

    if( this.start < this.end){
        return [true, buf, this.buffer[this.start - 1]];
    } else {
        return [true, buf, ''];
    }
};


// skip over data in the file until we find the expression of interest
// Unlike getuntil, we do not save the region we skipped over, instead we
// return the character the expression matched
fastqaParse.prototype.findNext = function(expression) {

    var idx = -1;
	if (this.eof)
		return "";
	if(this.start != 0 && (this.start < this.end) ){
		// if start is not 0, we need to substr and offset the search string,
		// otherwise we can skip the substr operation and just work with the original string
		idx = this.buffer.substr(this.start).search(expression);
		this.start = (idx == -1) ? this.end : this.start + idx;
	}

	while(idx == -1){
		if (this.start >= this.end) {
			this.buffer = this.read(this.buffersize);
			if (this.end == 0)
				return "";
		}
		idx = this.buffer.search(expression);
		this.start = (idx == -1)? this.end : idx;
	}
    return this.buffer.charAt(this.start++);
}

// first child class.  this one handles fasta/q files if we are provided a file handle
fastqaFileInput.prototype = new fastqaParse();
fastqaFileInput.prototype.constructor = fastqaFileInput;
function fastqaFileInput(fHandle, reader){
	fastqaParse.call( this );

    this.fileHandle = fHandle;
    this.size = this.fileHandle.size;
    this.buffersize = Math.min(this.size, BUFFER_MAX);
    this.end = this.buffersize;

    this.reader = reader;
    this.readfunc = (this.reader.readAsBinaryString) ?
        function(reader, data){ return reader.readAsBinaryString(data);} :
        function(reader, data){return reader.readAsText(data);};
    this.buffer = this.read(this.buffersize);
}

// override the read function
fastqaFileInput.prototype.read = function(count){
    // read a chunk of the file into memeory so that we can process it
    if(this.pos >= this.size || this.end != count){
        this.start = 0;
        this.end = 0;
        this.eof = true;
        return "";
    }
    this.prev = this.pos;
    var data = this.readfunc(this.reader, this.fileHandle.slice(this.pos, this.pos + count));
    this.start = 0;
    this.pos += data.length;
    this.end = data.length;
    return data;
};


// first child class.  this one handles fasta/q input if we are provided text
fastqaTextInput.prototype = new fastqaParse();
fastqaTextInput.prototype.constructor = fastqaTextInput;
function fastqaTextInput(text){
	fastqaParse.call( this );
    this.size = text.length;
    this.buffersize = this.size
    this.end = this.buffersize;
    this.buffer = text;
}


// data processor
function nextRecord(handler) {
    var c,
        untildat = [false, [], ''],
        sdat = {name:'', comment:'', qual:'', seq:''},
        slen;

    if(handler.last_char == 0){
        handler.last_char = handler.findNext(/[>@]/);
        if(handler.last_char == '')
            return [-1, sdat];
    }

    untildat = handler.getuntil(SEP_SPACE, []);
    if (untildat[0] == false)
        return [-1, sdat];
    sdat['name'] = untildat[1].join("");

    if (untildat[2] != '\n'){
		sdat['comment'] = handler.getuntil('\n', [])[1].join("");
    }

    untildat[1] = [];
    while ((c = handler.getc()) != '' && c != '>' && c != '+' && c != '@') {
        //if(isspace(c))
        //    continue;
        untildat[1][untildat[1].length] = c;
        untildat = handler.getuntil('\n', untildat[1]);  //read the rest of the line
    }
    sdat['seq'] = untildat[1].join("");

	// fasta should stop here...
	handler.last_char = (c == '>' || c == '@')? c : 0;
    if (c != '+')
        return 	[sdat['seq'].length, sdat];

    if(handler.findNext('\n') == '') // skip the rest of '+' line
        return [-2, sdat];

    untildat[1] = [];
    slen = sdat['seq'].length;
    do{
        untildat = handler.getuntil('\n', untildat[1]);
        slen -= untildat[1].join("").length;
    }while(slen > 0 && untildat[0] == true);
	sdat['qual'] = untildat[1].join("");

	return [ (sdat['seq'].length == sdat['qual'].length) ? sdat['seq'].length : -2, sdat];
}
