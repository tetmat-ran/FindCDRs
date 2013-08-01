/*

  FindCDRs.js
  -Also require FindCDRs.html, FindCDRs.css

  *only tested for Google Chrome v27.0.1453.110

  Extracting CDR sequences from LC and HC sequencing results

  Tet Matsuguchi <tet@alum.mit.edu>

  Last modified: Aug 01 2013 (v0.610)
  
  Many lines on FileReader have been borrowed from:
  http://www.html5rocks.com/en/tutorials/file/dndfiles/

  http://thiscouldbebetter.wordpress.com/2012/12/18/loading-editing-and-saving-a-text-file-in-html5-using-javascrip/

*/

var version = "0.610";

// zip.js configuration
zip.useWebWorkers = false;

// Global variables... yikes?
var allowZipFiles = true;
var allowAB1Files = true;

var messageBox = 0;
var startDate = new Date();

var nSeqFilesPending = 0;
var nSeqFilesProcessed = 0;

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
SeqFile.prototype.sequenceTags = ["ATCAGCAGTCTG", "CAGCCGGAAGAC", "TTCGCAACTTAT",
				  "TACTGTCAGCAA", "", "ACGTTCGGACAG",
				  "GGTACCAAGGTG", "GAGATCAAACGA", "ACTGTGGCTGCA",
				  "CAGCCAGGGGGC", "TCACTCCGTTTG", "TCCTGTGCAGCT",
				  "TCTGGCTTCAAC", "", "CACTGGGTGCGT",
				  "CCCGGGTAAGGG",
				  "GAATGGGTTGCA", "", "TATGCCGATAGC",
				  "TCAAGGGCCGTT", "CACTATAAGCGC", "GACACATCCAAA", "ACACAGCCTACC", "ACAAATGAACAG", "TTAAGAGCTGAG", "ACACTGCCGTCT",
				  "TATTGTGCTCGC", "", "GACTACTGGGGT",
				 "CAAGGAACCCTG", "GTCACCGTCTCC", "TCGGCCTCCACC"];
SeqFile.prototype.sequenceTagIDs = ["LC01_1", "LC01_2", "LC01_3",
				    "LC_5", "LC", "LC_3",
				    "LC12_1", "LC12_2", "LC12_3",
				    "HC01_1", "HC01_2", "HC01_3",
				    "HC1_5", "HC1", "HC1_3",
				    "HC12",
				    "HC2_5", "HC2", "HC2_3",
				    "HC23_1", "HC23_2", "HC23_3", "HC23_4", "HC23_5", "HC23_6", "HC23_7",
				    "HC3_5", "HC3", "HC3_3",
				    "HC34_1", "HC34_2", "HC34_3"];
SeqFile.prototype.sequenceTagSymbols = ["o", "o", "o",
					"[", "LC", "]",
					"o", "o", "o",
					"o", "o", "o",
					"[", "HC1", "]",
					"o",
					"[", "HC2", "]",
					"o", "o", "o", "o", "o", "o", "o",
					"[", "HC3", "]",
					"o", "o", "o"];
SeqFile.prototype.CDRs_id = ["LC", "HC1", "HC2", "HC3"];


SeqFile.prototype.getSeq = function () {
    nSeqFilesPending++;

    var seqFile = this;
    if (this.compressed) {
	// unzip contents using zip.js
	this.file.getData(new zip.TextWriter(), 
			  function(text) { // onend callback function from uncompressing a file
			      seqFile.processContent(text);
			      nSeqFilesPending--;
			      nSeqFilesProcessed++;
			  });
    } else {
	// Uncompressed; use regular FileReader
	var reader = new FileReader();
	reader.onload = function (e) { 
	    seqFile.processContent(e.target.result);
	    nSeqFilesPending--;
	    nSeqFilesProcessed++;
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

    document.getElementById("messages").innerHTML += "Start time: " + getTimestamp(startDate) + "<br />\n";

    document.getElementById("version_box").innerHTML = "v" + version;
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
    status("Reading files...");
    processFiles(files);

    waitForProcessingFiles();
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
	    status("Unzipping files...");
	    nSeqFilesPending++;
	    zip.createReader(new zip.BlobReader(f), 
			     function(reader) { // callback function
				 reader.getEntries(function (entries) {
				     processFiles(entries);
				     nSeqFilesPending--;
				     // nSeqFilesProcessed++; // don't count zipfile
				     reader.close(function() {} );
				 });
			     },
			     function (error) { // onerror during zipfile reading
				 nSeqFilesPending--;
				 // nSeqFilesProcessed++; // don't count zipfile
				 message("Unzipping " + f.name + " failed... Skpping.");
			     });
	} else {
	    document.getElementById("messages").innerHTML += "Ignoring " + escape(f.name) + "<br />";
	}
    }
}

function waitForProcessingFiles() {
    // document.getElementById("nFiles").innerHTML = "<br />" + nSeqFilesProcessed + "/" + seqFiles.length + "<br />";
    //document.getElementById("seqFilesProcessed").innerHTML = nSeqFilesProcessed;
    //document.getElementById("totalSeqFiles").innerHTML = seqFiles.length;

    // document.getElementById("statusBar").style.width = (100.0 * nSeqFilesProcessed / seqFiles.length) + "%";
    // document.getElementById("statusBar").innerHTML = nSeqFilesProcessed + " / " + seqFiles.length;
    if (seqFiles.length == 0) {
	status("Unzipping files...");
    } else {
	status("Reading files... " + nSeqFilesProcessed + " / " + seqFiles.length);
    }

    if (nSeqFilesPending) {
	window.setTimeout(function() { waitForProcessingFiles(); }, 100);
    } else {
	reportGroups();
    }
}

function reportGroups() {
    status("Grouping sequences...");
    if (nSeqFilesPending) {
	window.setTimeout(function() { reportGroups(); }, 100);
    } else {
	makeGroups();
	createTableByGroups();
	status("Done processing!");
    }
}

function makeGroups() {

    resetGroups();

    // Pull the RegExp from the dropdown menu to use for grouping sequences
    // Remember to escape \ with \\
    var re_grouping = new RegExp(makeRegExp(document.getElementById("group_regexp").value));

    naElement = document.createElement('font');
    naElement.className = 'na';
    naElement.innerHTML = 'n/a';


    for (seqFile in seqFiles) {
	var m = re_grouping.exec(seqFiles[seqFile].name);

	// If the sequence name does not have a match for the regexp,
	// this is a lone sequence with no group identifier
	if (!m) {
	    seqFiles[seqFile].tableGroupName.innerHTML = "";
	    seqFiles[seqFile].group = -1;
	    seqNotInGroups.push(seqFiles[seqFile]);
	    seqFiles[seqFile].tableGroupName.appendChild(naElement.cloneNode(true));
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
	regexpDOM.add(new Option(newOptionValue), regexpDOM.options[regexpDOM.options.length - 1]);
	regexpDOM.selectedIndex = regexpDOM.length - 2;
    }

    reportGroups();
}

function createTableByGroups() {
    var singleTable = document.createDocumentFragment(); // instead of adding directly to existing table <= faster

    // Group Border Element
    var groupBorder = document.createElement('tr');
    groupBorder.className = "group_border";
    var groupBorderTD = document.createElement('td');
    groupBorderTD.colSpan = "100";
    groupBorder.appendChild(groupBorderTD);

    // Add ungrouped sequences first
    for (iSeq in seqNotInGroups) {
	singleTable.appendChild(seqNotInGroups[iSeq].tableRow);
	singleTable.appendChild(groupBorder.cloneNode(true));
    }

    // Add grouped sequences
    for (iGroup in seqInGroups) {
	for (iSeq in seqInGroups[iGroup]) {
	    singleTable.appendChild(seqInGroups[iGroup][iSeq].tableRow);
	}
	singleTable.appendChild(groupBorder.cloneNode(true));
    }

    document.getElementById("table_seqfiles").appendChild(singleTable);
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

function getOutput() {
    // Outputting array (lines) or arrays (each line)
    // It can be outputted as text, excel file, table, etc

    var lines = [];
    lines.push(["Group ID", "LC", "HC1", "HC2", "HC3", "LC", "HC1", "HC2", "HC3", "Sequence Names"]);

    // Those not in groups first
    for (var iSeq in seqNotInGroups) {

	var seq = seqNotInGroups[iSeq];
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

	lines.push(lineElements);
    }

    // Those in groups
    for (var iGroup in seqGroups) {

	var seqs = [];
	var CDRs_dna = [[], [], [], []];
	var CDRs_aa = [[], [], [], []];

	for (var iSeq in seqInGroups[iGroup]) {
	    seqs.push(seqInGroups[iGroup][iSeq].name);

	    for (iCDR in CDRs_dna) {
		if (seqInGroups[iGroup][iSeq].CDRs_dna[iCDR] != "") {
		    CDRs_dna[iCDR].push(seqInGroups[iGroup][iSeq].CDRs_dna[iCDR]);
		    CDRs_aa[iCDR].push(translate(seqInGroups[iGroup][iSeq].CDRs_dna[iCDR]));
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
			    CDRs_dna[3].join(",")];
	lineElements = lineElements.concat(seqs);

	lines.push(lineElements);
    }

    return lines;
}

function saveTextAsFile()
{
    var textToWrite = "";
    var textLines = [];
    var output = getOutput();

    for (var i in output) {
	textLines.push(output[i].join("\t"));
    }
    textToWrite = textLines.join("\n");

    var textFileAsBlob = new Blob([textToWrite], {type:'text/plain'});
    var fileNameToSaveAs =  getOutputFilename() + ".txt";
    
    var downloadLink = document.createElement("a");
    downloadLink.download = fileNameToSaveAs;
    downloadLink.innerHTML = "Download File";

    // Chrome allows the link to be clicked programmatically.
    downloadLink.href = window.webkitURL.createObjectURL(textFileAsBlob);
    downloadLink.click();
}

function export2excel() {
    var output = getOutput();
    
    var outputTable = '';
    for (var i in output) {
	outputTable += '<tr><td>' + output[i].join('</td><td>') + '</td></tr>';
    }
    
    var data_type = 'data:application/vnd.ms-excel';
    var downloadLink = document.createElement("a");
    downloadLink.download = getOutputFilename() + ".xls";
    downloadLink.innerHTML = "Download Excel File";

    var excelHeader = '<html><head><meta name=ProgId content=Excel.Sheet><xml><x:ExcelWorkbook><x:ExcelWorksheets><x:ExcelWorksheet><x:Name>CDRs</x:Name></x:ExcelWorksheet></x:ExcelWorksheets></x:ExcelWorkbook></xml><style>td { white-space:nowrap; }</style></head><body><table>'
    var excelFooter = '</table></body></html>';

    // Chrome allows the link to be clicked programmatically.
    downloadLink.href = data_type + ', ' + (excelHeader + outputTable + excelFooter).replace(/ /g, '%20');
    downloadLink.click();
}

function getOutputFilename() {
    return getTimestamp(new Date()) + "_FindCDRs_v" + version;
}

function describe(obj) {
    document.getElementById("message_box").style.display = "block";

    var messageBox = document.getElementById("messages");
    for (elm in obj) {
	messageBox.innerHTML += elm + ": " + obj[elm] + "<br />\n";
    }
}

function message(msg) {
    document.getElementById('messages').innerHTML += msg + "<br />\n";
}

function toggleMessageBox() {
    messageBoxDom = document.getElementById("message_box");
    messageBoxButton = document.getElementById("show_messages");

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

function status(msg) {
    document.getElementById("statusBox").innerHTML = msg;
}

function reportTime() {
    var now = new Date();
    var elapsed = (now.getTime() - startDate.getTime()) / 1000;
    var timeMessage = elapsed + " seconds have elapsed.";
    document.getElementById('messages').innerHTML += timeMessage + "<br />\n";
}