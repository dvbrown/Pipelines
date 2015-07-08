/*
This macro expects files in a slide specific directory to be labelled 1_OLIG.tif etc
The second field of this slide is to be named 2_OLIG.tif. The macros goes through and does a background correction
binarisation and merging of channels. The merging is done by stack and not channel
if one of the images has no signal.
The number of fields is currently HARD CODED. I should change this eventually. There is no quantification of results here.
*/

// Retrieve the parameters for the current iteration
dir1 = getDirectory("Choose the parental folder containing the images");
list = getFileList(dir1);

// Background correct, threshold and make binary
function maskDAPI(image) 
	{
	open(dir1 + image);
    run("Subtract Background...", "rolling=125");
	run("Convert to Mask");
    run("Fill Holes");
    run("Outline");
	}

//.................................................................
function merge(image1, image2, image3, mergedFilename) {
	//Merge channels and save a new file
	run("Merge Channels...", "c1=" + image1 + " c2=" + image2 + " c3=" + image3 + " create");
	
    run("Stack to RGB");
    run("Size...", "width=900 height=668 constrain average interpolation=Bilinear");
    saveAs("Tiff", dir1 + mergedFilename);

    close();
	}

//.................................................................


//.................................................................
// Loop through the list of file names and get pairs of Olig2 and CD44 together

// Intialise a new array containing the first letter of the filenames
fields = newArray("1","2","3","4","5");
fields = newArray("1","2");

for (i=0; i<fields.length; i++){
	sampleStart = fields[i];
    dapi = sampleStart + "_" + "DAPI.tif";
	olig = sampleStart + "_" + "OLIG.tif";
    open(dir1 + olig);
    run("Subtract Background...", "rolling=500");
    // This can be 0.001 or 0.001 or no contrast at all
    run("Enhance Contrast...", "saturated=0.001");

	CD44 = sampleStart + "_" + "CD44.tif";
    open(dir1 + CD44);
    run("Subtract Background...", "rolling=500");
    run("Enhance Contrast...", "saturated=0.001");

	mergedFilename = sampleStart + "_merged.tif";
	print("The script is running ", olig, " ", CD44, " ", mergedFilename,"\n");

	maskDAPI(dapi);
	merge(CD44, olig, dapi, mergedFilename);

	close();	
	}
