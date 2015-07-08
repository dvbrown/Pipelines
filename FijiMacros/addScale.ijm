dir1 = getDirectory("Choose Source Directory ");
list = getFileList(dir1);

setBatchMode(true);
for (i=0; i<list.length; i++) {
  showProgress(i+1, list.length);
  open(dir1+list[i]);
  run("Scale Bar...", "width=224 height=12 font=42 color=Black background=White location=[Lower Right] bold");
  run("Save");
  run("Open Next");
}
