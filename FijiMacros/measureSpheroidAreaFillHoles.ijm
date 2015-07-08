/* This macro measures the area of spheroid and writes it to an excel file
There is a prompt for manual input to set a threshold and adjust a circle
Haven't managed to get the threshold window to appear in the foreground.
*/

// Get the directries to read and write files
dir1 = getDirectory("Choose the folder that contains the images to be analysed ");
dir2 = getDirectory("Make a NEW folder to contain your processed images and output data ");

list = getFileList(dir1);

// Open a file and set the calibration based on 1mm ruler
open(dir1 + list[1]);
runMacro("invasionAssay/setScaleInvasion.ijm");
close();

// Set to record the sample labels in the results table
run("Set Measurements...", "area mean min display redirect=None decimal=3");

// Iterate over the files in the directory
for (i=0; i<list.length; i++) {
	showProgress(i+1, list.length);
	filename = dir1 + list[i];
	if (endsWith(filename, "tif")) {
		open(filename);

// ###################################  START MACRO  ################################### 
// Transform the image to 8 bit then make the zoom smaller

	run("8-bit");
	run("Set... ", "zoom=33 x=1296 y=972");

// Ask the user to set a threshold
	setAutoThreshold("Default");
	run("Threshold...");
	makeOval(612, 78, 1743, 1695);
	setTool("oval");
	waitForUser("Set the threshold of the image", "Use the sliding bar to the desired threshold. \nAdjust the circle around the neurosphere\n Press 'OK' on this window to continue");
	
	selectWindow("Threshold");
	
// Measure the area
	run("Analyze Particles...", "size=0-Infinity pixel circularity=0.00-1.00 show=Nothing include summarize");
	
	// Label the circle (also known as ROI)
	run("Label");

// Save the masked output image with scale bar in specified directory
	run("Make Binary", "thresholded remaining black");
	run("Scale Bar...", "width=250 height=12 font=42 color=Black background=White location=[Lower Right] bold");
	saveAs("TIFF", dir2+list[i]);
	
	close();
	}
}
// ###################################  END MACRO  ################################### 

// Save the measurements

dataFile = getString("What would you like to call your data file? ", "output.xls");
outData = dir2 + "/" + dataFile;
//saveAs("summary", outData);

selectWindow("Summary"); 
saveAs("Text",  outData);
