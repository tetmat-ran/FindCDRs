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

// Global variables... yikes?
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
function SeqFile(f, i) {
    this.i = i;
    this.file = f;
    this.filename = f.name;
    this.name = f.name.substr(0, f.name.length - 4);
    this.size = f.size;

    this.valid = 1; // Innocent unless proven guilty
    this.direction = ""; // "fwd" or "rev"
    this.group = -1;
    this.CDRs_dna = ["", "", "", ""];

    this.lastModifiedDate = getTimestamp(f.lastModifiedDate);

    this.addToTable();
    this.getSeq();
}

// SeqFile constants
SeqFile.prototype.CDRs_id = ["LC", "HC1", "HC2", "HC3"];
SeqFile.prototype.sequenceTagIDs = ["LC_5", "LC_3", "HC1_5", "HC1_3", "HC2_5", "HC2_3", "HC3_5", "HC3_3"];
SeqFile.prototype.sequenceTags = ["TACTGTCAGCAA", "ACGTTCGGACAG", "TCTGGCTTCAAC", "CACTGGGTGCGT", "GAATGGGTTGCA", "TATGCCGATAGC", "TATTGTGCTCGC", "GACTACTGGGGT"];

SeqFile.prototype.getSeq = function () {
    var reader = new FileReader();
    var seqFile = this;
    reader.onload = function (e) { 
	seqFile.processContent(e.target.result);
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

SeqFile.prototype.findCDRs = function () {
    // CDRs_dna[0] = LC : seqTags[0], seqTags[1]
    // CDRs_dna[1] = HC1: seqTags[2], seqTags[3]
    // CDRs_dna[2] = HC2: seqTags[4], seqTags[5]
    // CDRs_dna[3] = HC3: seqTags[6], seqTags[7]

    this.CDRs_dna = ["", "", "", ""];

    for (iCDR = 0; iCDR < 4; iCDR++) {
	if (positions[2 * iCDR] != -1 && positions[2 * iCDR + 1] != -1) {
	    this.CDRs_dna[iCDR] = this.sequence.slice(positions[2 * iCDR] + this.sequenceTags[2 * iCDR].length,
						      positions[2 * iCDR + 1]);
	}
    }
}

SeqFile.prototype.translateCDRs = function () {
    for (iCDR = 0; iCDR < 4; iCDR++) {
	// CDR sequence information is available
	if (this.CDRs_dna[iCDR] != "") {
	    // Does it code valid triplet aa's?
	    if (this.CDRs_dna[iCDR].length % 3 != 0) {
		document.getElementById(this.i + "_" + this.CDRs_id[iCDR]).style.color = "red";
	    } else {
		document.getElementById(this.i + "_" + this.CDRs_id[iCDR]).style.color = "black";
	    }
	}
    }
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
    this.tableSeqTags.innerHTML = '<seqtags>o <font id="' + this.i + '_LC_5">[</font> <font id="' + this.i + '_LC">LC</font> <font id="' + this.i + '_LC_3">]</font> o o o <font id="' + this.i + '_HC1_5">[</font> <font id="' + this.i + '_HC1">HC1</font> <font id="' + this.i + '_HC1_3">]</font> o <font id="' + this.i + '_HC2_5">[</font> <font id="' + this.i + '_HC2">HC2</font> <font id="' + this.i + '_HC2_3">]</font> o <font id="' + this.i + '_HC3_5">[</font> <font id="' + this.i + '_HC3">HC3</font> <font id="' + this.i + '_HC3_3">]</font> o o</seqtags>';
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

    makeGroups();
}

function makeGroups() {

    resetGroups();

    var groupTable = document.getElementById("table_grouping");

    // Remember to escape \ with \\
    var re_grouping = new RegExp(makeRegExp(document.getElementById("group_regexp").value));

    for (seqFile in seqFiles) {
	var m = re_grouping.exec(seqFiles[seqFile].name);
	if (!m) {
	    seqFiles[seqFile].group = -1;
	    seqNotInGroups.push(seqFiles[seqFile].i);

	    // Add to grouping by itself as non-group
	    var newGroupSeq = document.createElement("td");
	    newGroupSeq.innerHTML = seqFiles[seqFile].name;

	    var newGroupSize = document.createElement("td");
	    newGroupSize.innerHTML = "0";

	    var newGroup = document.createElement("tr");
	    newGroup.className = "nongroup";
	    newGroup.appendChild(newGroupSeq);
	    newGroup.appendChild(newGroupSize);

	    groupTable.appendChild(newGroup);
	    
	    continue;
	}

	var group_identifier = m[m.length - 1];
	seqFiles[seqFile].groupIdentifierAt = 
	    seqFiles[seqFile].name.slice(m.index, seqFiles[seqFile].name.length).indexOf(group_identifier) + m.index;

	var iGroup = seqGroups.indexOf(group_identifier);
	if (iGroup == -1) {
	    seqFiles[seqFile].group = iGroup = seqGroups.length;
	    seqGroups.push(group_identifier);
	    seqInGroups.push([seqFiles[seqFile].i]);

	    // Add to grouping by itself as non-group
	    var newGroupSeqs = document.createElement("td");
	    newGroupSeqs.id = iGroup + "_groupSeqs";
	    newGroupSeqs.innerHTML = seqFiles[seqFile].getGroupTaggedName();

	    var newGroupSize = document.createElement("td");
	    newGroupSize.id = iGroup + "_groupSize";
	    newGroupSize.innerHTML = "1";

	    var newGroup = document.createElement("tr");
	    newGroup.appendChild(newGroupSeqs);
	    newGroup.appendChild(newGroupSize);

	    groupTable.appendChild(newGroup);

	} else {
	    seqFiles[seqFile].group = iGroup;
	    seqInGroups[iGroup].push(seqFiles[seqFile].i);

	    var thisGroupSeqs = document.getElementById(iGroup + "_groupSeqs");
	    thisGroupSeqs.innerHTML += "<br />" + seqFiles[seqFile].getGroupTaggedName();

	    var thisGroupSize = document.getElementById(iGroup + "_groupSize");
	    thisGroupSize.innerHTML = parseInt(thisGroupSize.innerHTML) + 1;
	}
    }
}

function resetGroups() {
    seqGroups = [];
    seqInGroups = [];
    seqNotInGroups = [];

    var groupTable = document.getElementById("table_grouping");

    while (groupTable.rows.length > 2) {
	groupTable.deleteRow(1);
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