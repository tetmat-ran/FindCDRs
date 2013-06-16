/*
  Many lines of codes are borrowed from:
  http://www.html5rocks.com/en/tutorials/file/dndfiles/
*/


// Global variables... yikes?
var seqFiles = [];

// Objects / Classes
function SeqFile(f, i) {
    this.i = i;
    this.file = f;
    this.filename = f.name;
    this.name = f.name.substr(0, f.name.length - 4);
    this.size = f.size;

    this.valid = 1; // Innocent unless proven guilty
    this.direction = ""; // "fwd" or "rev"

    lastModifiedMonth = f.lastModifiedDate.getMonth() + 1;
    lastModifiedHours = f.lastModifiedDate.getHours();
    lastModifiedMinutes = f.lastModifiedDate.getMinutes();
    lastModifiedSeconds = f.lastModifiedDate.getSeconds();

    this.lastModifiedDate = f.lastModifiedDate.getFullYear() + "-" + 
	(lastModifiedMonth < 10 ? "0" : "") + lastModifiedMonth + "-" + 
	f.lastModifiedDate.getDate() + " " +
	(lastModifiedHours < 10 ? "0" : "") + lastModifiedHours + ":" + 
	(lastModifiedMinutes < 10 ? "0" : "") + lastModifiedMinutes + ":" + 
	(lastModifiedSeconds < 10 ? "0" : "") + lastModifiedSeconds;

    this.addToTable();
    this.getSeq();
}

SeqFile.prototype.getSeq = function () {
    var reader = new FileReader();
    var seqFile = this;
    reader.onload = function (e) { 
	seqFile.processContent(e.target.result);
	//document.getElementById("comments").innerHTML += e.target.result + "<br />";
    }

    reader.readAsText(this.file);
}

SeqFile.prototype.processContent = function (fileContent) {
    // Make sure the content is in FASTA format
    if (fileContent[0] != ">") {
	// FASTA format starts with ">" followed by sequence name in the first line
	this.invalidate();
	return 0;
    }
    
    firstNewLineIndex = fileContent.search(/[\n\r]/);
    if (firstNewLineIndex == -1) {
	// no new line found...
	this.invalidate();
	return;
    }

    fileContent = fileContent.replace(/[\r\n\t ]/g, '');

    this.sequenceName = fileContent.substr(1, firstNewLineIndex - 1);
    this.sequence = fileContent.slice(firstNewLineIndex, fileContent.length);

    if (this.sequence.search(/[^ACGTN]/) != -1) {
	// non-standard DNA sequence
	this.invalidate();
	return;
    }

    // Check forward fro sequence tags
    positions = this.findSequenceTags();
    if (positions == 0) {
	this.sequence.original = this.sequence;
	this.sequence = revcomp(this.sequence);
	positions = this.findSequenceTags();
	if (positions == 0) {
	    this.invalidate();
	} else {
	    this.assignReverse();
	}
    } else {
	this.assignForward();
    }
}

SeqFile.prototype.sequenceTagIDs = ["LC_5", "LC_3", "HC1_5", "HC1_3", "HC2_5", "HC2_3", "HC3_5", "HC3_3"];
SeqFile.prototype.sequenceTags = ["TACTGTCAGCAA", "ACGTTCGGACAG", "TCTGGCTTCAAC", "CACTGGGTGCGT", "GAATGGGTTGCA", "TATGCCGATAGC", "TATTGTGCTCGC", "GACTACTGGGGT"];

/*
LC  = ["TACTGTCAGCAA", "ACGTTCGGACAG"]
HC1 = ["TCTGGCTTCAAC", "CACTGGGTGCGT"]
HC2 = ["GAATGGGTTGCA", "TATGCCGATAGC"]
HC3 = ["TATTGTGCTCGC", "GACTACTGGGGT"]
*/

SeqFile.prototype.findSequenceTags = function () {
    var nPositives = 0;
    positions = [];
    for (tag in this.sequenceTags) {
	var pos = this.sequence.indexOf(this.sequenceTags[tag])
	if (pos != -1) {
	    nPositives++;
	    document.getElementById(this.i + "_" + this.sequenceTagIDs[tag]).style.color = "black";
	}
	positions.push(pos);
    }

    if (nPositives == 0) {
	return 0;
    }

    return positions;
}

SeqFile.prototype.addToTable = function () {
    this.tableRow = document.createElement('tr');
    this.tableRow.className = "processing";

    this.tableName = document.createElement('td');
    this.tableName.appendChild(document.createTextNode(escape(this.name)));
    this.tableRow.appendChild(this.tableName);

    this.tableDirection = document.createElement('td');
    this.tableRow.appendChild(this.tableDirection);

    this.tableSeqTags = document.createElement('td');
    this.tableSeqTags.innerHTML = '<seqtags>o <font id="' + this.i + '_LC_5">[</font> LC <font id="' + this.i + '_LC_3">]</font> o o o <font id="' + this.i + '_HC1_5">[</font> HC1 <font id="' + this.i + '_HC1_3">]</font> o <font id="' + this.i + '_HC2_5">[</font> HC2 <font id="' + this.i + '_HC2_3">]</font> o <font id="' + this.i + '_HC3_5">[</font> HC3 <font id="' + this.i + '_HC3_3">]</font> o o</seqtags>';
    this.tableRow.appendChild(this.tableSeqTags);

    document.getElementById('table_seqfiles').appendChild(this.tableRow);
}

SeqFile.prototype.assignForward = function () {
    this.direction = "fwd";
    this.tableRow.className = "fwd";
    this.tableDirection.innerHTML = "fwd";
}

SeqFile.prototype.assignReverse = function () {
    this.direction = "rev";
    this.tableRow.className = "rev";
    this.tableDirection.innerHTML = "rev";
}

SeqFile.prototype.invalidate = function () {
    this.valid = 0;
    this.tableRow.className = "bad";
}

function init() {
    // We need these functionalities to read files locally
    if (!(window.File && window.FileReader && window.FileList && window.Blob)) {
	document.getElementById('choose_files').disabled = true;
	alert('Your browser does not support local file access.\n\nTry Chrome... if you distrust Google, try Chromium.');
	return;
    }

    document.getElementById('choose_files').addEventListener('change', handleFileSelect, false);

    var dropZone = document.getElementById('drop_zone');
    dropZone.addEventListener('dragover', handleDragOver, false);
    dropZone.addEventListener('drop', handleDrop, false);
}

// Drop zone
function handleDragOver(evt) {
    evt.stopPropagation();
    evt.preventDefault();
    evt.dataTransfer.dropEffect = 'copy'; // Explicitly show this is a copy.
}

function handleDrop(evt) {
    evt.stopPropagation();
    evt.preventDefault();

    var files = evt.dataTransfer.files; // FileList object.

    processFiles(files);
}

function handleFileSelect(evt) {
    var files = evt.target.files; // FileList object
    processFiles(files);
}

// Processing Files
function processFiles(files) {
    for (var i = 0, f; f = files[i]; i++) {

	// At this time only .seq files are used
	if (f.name.substr(-4, 4) == ".seq") {
	    seqFiles.push(new SeqFile(f, seqFiles.length));
	} else if (f.name.substr(-4, 4) == ".ab1") {
	    // TODO: Add support for .ab1 files
	} else {
	    document.getElementById("drop_zone_message_box").innerHTML += "Ignoring " + escape(f.name) + "<br />";
	}

    }
}

function revcomp(seq) {
    seq = seq.replace(/A/gi, "x");
    seq = seq.replace(/C/gi, "y");
    seq = seq.replace(/G/gi, "C");
    seq = seq.replace(/T/gi, "A");
    seq = seq.replace(/x/g, "T");
    seq = seq.replace(/y/g, "G");
    return seq.split("").reverse().join("");
}
