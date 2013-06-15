/*
  Many lines of codes are borrowed from:
  http://www.html5rocks.com/en/tutorials/file/dndfiles/
*/


// Global variables... yikes?
var seqFiles = [];

// Objects / Classes
function SeqFile(f) {
    this.file = f;
    this.filename = f.name;
    this.size = f.size;

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

    this.getSeq();
}

SeqFile.prototype.getSeq = function () {
    var reader = new FileReader();
    reader.onload = function (e) { 
	document.getElementById("comments").innerHTML += e.target.result + "<br />";
    }



    reader.readAsText(this.file);
}

SeqFile.prototype.setSeq = function (evt) {
    alert('read!');
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
	    seqFiles.push(new SeqFile(f));
	} else if (f.name.substr(-4, 4) == ".ab1") {
	    // TODO: Add support for .ab1 files
	} else {
	    document.getElementById("drop_zone_message_box").innerHTML += "Ignoring " + escape(f.name) + "<br />";
	}

    }
}
