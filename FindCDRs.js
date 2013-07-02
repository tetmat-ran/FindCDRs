/*

FindCDRs.js
-Also require FindCDRs.html, FindCDRs.css

*only tested for Google Chrome v27.0.1453.110

Extracting CDR sequences from LC and HC sequencing results

Tet Matsuguchi <tet@alum.mit.edu>

Last modified: Jun 16, 2013 (v0.2)
 
Many lines on FileReader have been borrowed from:
http://www.html5rocks.com/en/tutorials/file/dndfiles/

http://thiscouldbebetter.wordpress.com/2012/12/18/loading-editing-and-saving-a-text-file-in-html5-using-javascrip/

*/

var version = "0.200";

// zip.js configuration
zip.useWebWorkers = false;

// Global variables... yikes?
var allowZipFiles = true;
var allowAB1Files = false;

var messageBox = 0;
var startDate = new Date();

var seqFiles = [];
var seqGroups = [];
var seqInGroups = [];
var seqNotInGroups = [];

var gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'};

// Objects / Classes
function SeqFile(f) {

    // Initialize
    this.file = f;
    this.valid = 1; // Innocent unless proven guilty
    this.direction = ""; // "fwd" or "rev"
    this.group = -1;
    this.CDRs_dna = ["", "", "", ""];

    // File information depends on compressed (zipped) vs uncompressed
    this.compressed = typeof f.compressedSize !== 'undefined';

    if (this.compressed) {
	this.filename = f.filename;
	this.size = f.uncompressedSize;
	this.lastModifiedDate = getTimestamp(f.lastModDate);
    } else {
	this.filename = f.name;
	this.size = f.size;
	this.lastModifiedDate = getTimestamp(f.lastModifiedDate);
    }

    this.name = this.filename.substr(0, this.filename.length - 4);

    this.tableRow = this.createRowElement();
    this.getSeq();
}

// SeqFile constants
SeqFile.prototype.sequenceTags = ["TACTGTCAGCAA", "", "ACGTTCGGACAG", "TCTGGCTTCAAC", "", "CACTGGGTGCGT", "GAATGGGTTGCA", "", "TATGCCGATAGC", "TATTGTGCTCGC", "", "GACTACTGGGGT"];
SeqFile.prototype.sequenceTagIDs = ["LC_5", "LC", "LC_3", "HC1_5", "HC1", "HC1_3", "HC2_5", "HC2", "HC2_3", "HC3_5", "HC3", "HC3_3"];
SeqFile.prototype.sequenceTagSymbols = ["[", "LC", "]", "[", "HC1", "]", "[", "HC2", "]", "[", "HC3", "]"];
SeqFile.prototype.CDRs_id = ["LC", "HC1", "HC2", "HC3"];


SeqFile.prototype.getSeq = function () {

    var seqFile = this;
    if (this.compressed) {
	// unzip contents using zip.js
	this.file.getData(new zip.TextWriter(), function(text) {
	    seqFile.processContent(text);
	});

    } else {
	// Uncompressed; use regular FileReader
	var reader = new FileReader();
	reader.onload = function (e) { 
	    seqFile.processContent(e.target.result);
	}
	
	reader.readAsText(this.file);
    }
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

    // Check forward from sequence tags
    this.seqTagPositions = this.findSequenceTags();
    if (this.seqTagPositions == 0) {
	this.sequence.original = this.sequence;
	this.sequence = revcomp(this.sequence);
	this.seqTagPositions = this.findSequenceTags();
	if (this.seqTagPositions == 0) {
	    this.invalidate();
	    return;
	} else {
	    this.assignReverse();
	}
    } else {
	this.assignForward();
    }

    this.findCDRs();
    this.translateCDRs();
}

SeqFile.prototype.findSequenceTags = function () {
    var nPositives = 0;
    var positions = [];
    for (var tag in this.sequenceTags) {
	var pos = -1;
	// Don't search for the "CDR's", they are not sequence tags
	if (this.sequenceTags[tag] != "") {
	    pos = this.sequence.indexOf(this.sequenceTags[tag])
	    if (pos != -1) {
		nPositives++;
		this.sequenceTagElements[tag].style.color = "black";
	    }
	}
	positions.push(pos);
    }

    if (nPositives == 0) {
	return 0;
    }

    return positions;
}

SeqFile.prototype.findCDRs = function () {
    this.CDRs_dna = [];
    
    for (var iCDR in this.CDRs_id) {
	var seq_dna = "";

	// Check to see whether boundary tags are found
	iSeqTag = this.sequenceTagIDs.indexOf(this.CDRs_id[iCDR]);
	if (this.seqTagPositions[iSeqTag - 1] != -1 &&
	    this.seqTagPositions[iSeqTag + 1] != -1) {

	    // If so, bounded region contains the CDR
	    seq_dna = this.sequence.slice(this.seqTagPositions[iSeqTag - 1] + this.sequenceTags[iSeqTag - 1].length,
						      this.seqTagPositions[iSeqTag + 1]);
	}

	this.CDRs_dna.push(seq_dna);
    }
}

SeqFile.prototype.translateCDRs = function () {
    for (var iCDR in this.CDRs_id) {
	// CDR sequence information is available
	if (this.CDRs_dna[iCDR] != "") {
	    // Does it code valid triplet aa's?
	    if (this.CDRs_dna[iCDR].length % 3 != 0) {
		this.sequenceTagElements[this.sequenceTagIDs.indexOf(this.CDRs_id[iCDR])].style.color = "red";
	    } else {
		this.sequenceTagElements[this.sequenceTagIDs.indexOf(this.CDRs_id[iCDR])].style.color = "black";
	    }
	}
    }
}

SeqFile.prototype.createRowElement = function () {
    var tableRow = document.createElement('tr');
    tableRow.className = "processing";

    this.tableName = document.createElement('td');
    this.tableName.appendChild(document.createTextNode(escape(this.name)));
    tableRow.appendChild(this.tableName);

    this.tableGroupName = document.createElement('td');
    tableRow.appendChild(this.tableGroupName);

    this.tableDirection = document.createElement('td');
    this.tableDirection.className = "direction";
    tableRow.appendChild(this.tableDirection);

    this.tableSeqTags = document.createElement('td');
    this.sequenceTagElement = document.createElement("seqtags");
    this.sequenceTagElements = [];
    for (sequenceTag in this.sequenceTags) {
	var tagElement = document.createElement("font");
	tagElement.innerHTML = this.sequenceTagSymbols[sequenceTag] + " ";

	this.sequenceTagElements.push(tagElement);
	this.sequenceTagElement.appendChild(tagElement);
    }
    this.tableSeqTags.appendChild(this.sequenceTagElement);
    tableRow.appendChild(this.tableSeqTags);

    return tableRow;
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
    this.tableDirection.innerHTML = "n/a";
}

SeqFile.prototype.getGroupTaggedName = function () {
    if (this.group == -1 || seqGroups[this.group] == "") {
	return this.name;
    }

    return this.name.slice(0, this.groupIdentifierAt) + 
	"<grouptag>" + seqGroups[this.group] + "</grouptag>"+ 
	this.name.slice(this.groupIdentifierAt + seqGroups[this.group].length, this.name.length);
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

    document.getElementById("comments").innerHTML += "Start time: " + getTimestamp(startDate) + "<br />\n";
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

    processInputFiles(files);
}

function handleFileSelect(evt) {
    var files = evt.target.files; // FileList object
    processInputFiles(files);
}

function processInputFiles(files) {
    processFiles(files);
    makeGroups();
    createTableByGroups();
    reportTime();
}

// Processing Files
function processFiles(files) {

    for (var i = 0, f; f = files[i]; i++) {
	if (typeof f.compressedSize !== 'undefined') {
	    // Compressed file has its filename defined differently
	    f.name = f.filename;
	}
	// At this time only .seq files are used
	if (f.name.substr(-4, 4) == ".seq") {
	    seqFiles.push(new SeqFile(f));
	} else if (f.name.substr(-4, 4) == ".ab1" && allowAB1Files) {
	    // TODO: Add support for .ab1 files
	} else if (f.name.substr(-4, 4) == ".zip" && allowZipFiles) {
	    zip.createReader(new zip.BlobReader(f), function(reader) {
		reader.getEntries(function (entries) {
		    processFiles(entries);
		    reader.close(function() {} );
		});
	    });
	} else {
	    document.getElementById("comments").innerHTML += "Ignoring " + escape(f.name) + "<br />";
	}
    }
}

function makeGroups() {

    resetGroups();

    var singleTable = document.getElementById("table_seqfiles");

    // Pull the RegExp from the dropdown menu to use for grouping sequences
    // Remember to escape \ with \\
    var re_grouping = new RegExp(makeRegExp(document.getElementById("group_regexp").value));

    for (seqFile in seqFiles) {
	var m = re_grouping.exec(seqFiles[seqFile].name);

	// If the sequence name does not have a match for the regexp,
	// this is a lone sequence with no group identifier
	if (!m) {
	    seqFiles[seqFile].group = -1;
	    seqNotInGroups.push(seqFiles[seqFile]);
	    seqFiles[seqFile].tableGroupName.innerHTML = "n/a";
	} else {

	    // The group identifier is either the full match to the regexp, or
	    // If there is one or more parenthical match, it is the last match
	    var group_identifier = m[m.length - 1];
	    // Note where in the name the identifier is (because sometimes the same identifier can occur more than once)
	    seqFiles[seqFile].groupIdentifierAt = 
		seqFiles[seqFile].name.slice(m.index, seqFiles[seqFile].name.length).indexOf(group_identifier) + m.index;

	    seqFiles[seqFile].tableGroupName.innerHTML = group_identifier;

	    var iGroup = seqGroups.indexOf(group_identifier); // is the identifier already in a group?

	    // If this is the first element of the group, then create a new group
	    if (iGroup == -1) {
		seqFiles[seqFile].group = iGroup = seqGroups.length; // new group
		seqGroups.push(group_identifier);
		seqInGroups.push([seqFiles[seqFile]]);
	    } else {
		seqFiles[seqFile].group = iGroup;
		seqInGroups[iGroup].push(seqFiles[seqFile]);
	    }
	} // end if (!m)

	// Add to table
	seqFiles[seqFile].tableName.innerHTML = seqFiles[seqFile].getGroupTaggedName(); // Highlight the seqname with the group identifier
    }
}

function resetGroups() {
    seqGroups = [];
    seqInGroups = [];
    seqNotInGroups = [];

    var singleTable = document.getElementById("table_seqfiles");

    while (singleTable.rows.length > 2) {
	singleTable.deleteRow(1);
    }
}

function changeGrouping() {
    regexpDOM = document.getElementById("group_regexp");
    if (regexpDOM.value == "Other...") {
	newOptionValue = prompt("Enter a new regular expression for grouping sequence names:");
	regexpDOM.add(new Option(newOptionValue));
	regexpDOM.selectedIndex = regexpDOM.length - 1;
    }

    makeGroups();
    createTableByGroups();
}

function createTableByGroups() {
    var singleTable = document.getElementById("table_seqfiles");
    for (iGroup in seqInGroups) {
	for (iSeq in seqInGroups[iGroup]) {
	    singleTable.appendChild(seqInGroups[iGroup][iSeq].tableRow);
	}
	var groupBorder = document.createElement('tr');
	groupBorder.className = "group_border";

	var groupBorderTD = document.createElement('td');
	groupBorderTD.colSpan = "100";

	groupBorder.appendChild(groupBorderTD);
	singleTable.appendChild(groupBorder);
    }

    for (iSeq in seqNotInGroups) {
	singleTable.appendChild(seqNotInGroups[iSeq].tableRow);
    }
}

function makeRegExp(str) {
    str = str.replace(/Row/g, "[A-H]");
    str = str.replace(/Column/g, "[01]?\\d");

    return str;
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

function translate(seq) {
    var aa = [];
    for (var i = 0, l = Math.floor(seq.length / 3); i < l ; i++) {
	aa.push(gencode[seq.slice(i * 3, i * 3 + 3)]);
    }
    return aa.join("") + new Array(seq.length % 3 + 1).join("^");
}

function getTimestamp(dateObject) {
    var month = dateObject.getMonth() + 1;
    var date = dateObject.getDate();
    var hours = dateObject.getHours();
    var minutes = dateObject.getMinutes();
    var seconds = dateObject.getSeconds();

    return dateObject.getFullYear() + "-" + 
	(month < 10 ? "0" : "") + month + "-" + 
	(date < 10 ? "0" : "") + date + "_" + 
	(hours < 10 ? "0" : "") + hours +
	(minutes < 10 ? "0" : "") + minutes +
	(seconds < 10 ? "0" : "") + seconds;
}

function saveTextAsFile()
{
    var textToWrite = "";
    var textLines = [["Group ID", "LC", "HC1", "HC2", "HC3", "LC", "HC1", "HC2", "HC3", "Sequence Names"].join("\t")];
    
    // Output ones in groups first
    for (var iGroup in seqGroups) {

	var seqs = [];
	var CDRs_dna = [[], [], [], []];
	var CDRs_aa = [[], [], [], []];

	for (seq_i in seqInGroups[iGroup]) {
	    var iSeq = seqInGroups[iGroup][seq_i];

	    seqs.push(seqFiles[iSeq].name);

	    for (iCDR in CDRs_dna) {
		if (seqFiles[iSeq].CDRs_dna[iCDR] != "") {
		    CDRs_dna[iCDR].push(seqFiles[iSeq].CDRs_dna[iCDR]);
		    CDRs_aa[iCDR].push(translate(seqFiles[iSeq].CDRs_dna[iCDR]));
		}
	    }
	}

	var lineElements = [seqGroups[iGroup],
			    CDRs_aa[0].join(","),
			    CDRs_aa[1].join(","),
			    CDRs_aa[2].join(","),
			    CDRs_aa[3].join(","),
			    CDRs_dna[0].join(","),
			    CDRs_dna[1].join(","),
			    CDRs_dna[2].join(","),
			    CDRs_dna[3].join(","),
			    seqs.join("\t")];

	textLines.push(lineElements.join("\t"));
    }

    // Those not in groups
    for (var iGroup in seqNotInGroups) {

	var seq = seqFiles[seqNotInGroups[iGroup]];
	var lineElements = [seq.name,
			    translate(seq.CDRs_dna[0]),
			    translate(seq.CDRs_dna[1]),
			    translate(seq.CDRs_dna[2]),
			    translate(seq.CDRs_dna[3]),
			    seq.CDRs_dna[0],
			    seq.CDRs_dna[1],
			    seq.CDRs_dna[2],
			    seq.CDRs_dna[3],
			    seq.name];

	textLines.push(lineElements.join("\t"));
    }

    textToWrite = textLines.join("\n");
    
    var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
    var fileNameToSaveAs = getTimestamp(new Date()) + "_FindCDRs_v" + version + ".txt";
    
    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";

    // Chrome allows the link to be clicked programmatically.
    downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
    downloadLink.click();
}

function describe(obj) {
    document.getElementById("comment_box").style.display = "block";

    var commentBox = document.getElementById("comments");
    for (elm in obj) {
	commentBox.innerHTML += elm + ": " + obj[elm] + "<br />\n";
    }
}

function message(msg) {
    document.getElementById('comments').innerHTML += msg + "<br />\n";
}

function toggleMessageBox() {
    messageBoxDom = document.getElementById("comment_box");
    messageBoxButton = document.getElementById("show_comments");

    if (messageBox) {
	messageBoxButton.className = "button";
	messageBoxDom.style.display = "none";
	messageBox = 0;
    } else {
	messageBoxButton.className = "pressed";
	messageBoxDom.style.display = "block";
	messageBox = 1;
    }
}

function reportTime() {
    var now = new Date();
    var elapsed = (now.getTime() - startDate.getTime()) / 1000;
    var timeMessage = elapsed + " seconds have elapsed.";
    document.getElementById('comments').innerHTML += timeMessage + "<br />\n";
}