dir1 = getDirectory("Choose Source Directory ");
list = getFileList(dir1);
dir2 = getDirectory("Make a NEW folder to contain your processed images and output data ");

setBatchMode(true);
for (i=0; i<list.length; i++) {
  showProgress(i+1, list.length);
  open(dir1+list[i]);
// Transform the image to 8 bit then make the zoom smaller

  run("8-bit");
  run("Set... ", "zoom=33 x=1296 y=972");
  run("Scale Bar...", "width=250 height=12 font=42 color=Black background=White location=[Lower Right] bold");
  saveAs("TIFF", dir2+list[i]);
  close();
  run("Open Next");
}
