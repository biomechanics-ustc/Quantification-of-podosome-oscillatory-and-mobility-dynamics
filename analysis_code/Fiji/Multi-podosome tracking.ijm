roiManager("Reset");
run("Clear Results");
run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
image=getImageID;

selectImage(image);
run("Z Project...", "projection=[Average Intensity]");
rename("Maxima");
maxima=getImageID;

run("Duplicate...", "title=Duplicate2");
rename("Clustermask");
clustermask=getImageID;

//Process podosome image
selectImage(maxima);
run("Subtract Background...", "rolling=10");
run("Unsharp Mask...", "radius=4.5 mask=0.80 stack");
for (j=1; j<=36; j++) {
	run("Smooth", "stack");
}
run("Subtract Background...", "rolling=10");
run("Unsharp Mask...", "radius=4.5 mask=0.80 stack");
for (k=1; k<=36; k++) {
	run("Smooth", "stack");
}

//Make mask for cluster
selectImage(clustermask);
//run("Subtract Background...", "rolling=20");
run("Median...", "radius=20");
setAutoThreshold("Huang dark");
run("Convert to Mask");
run("Analyze Particles...", "size=10000-Infinity add");

	
//Find Podosome Maxima
selectImage(maxima);
roiManager("Select", 0);
run("Enlarge...", "enlarge=-20");
run("Find Maxima...", "prominence=5 exclude output=List");

close("Maxima");
close("Clustermask");

x=newArray(nResults);
y=newArray(nResults);

for (i=0; i<nResults; i++) {
	x[i]=getResult("X", i);
	y[i]=getResult("Y", i);
}

selectImage(image);
run("Clear Results")
for (i=0; i<x.length; i++) {
		setResult("Profile "+i+1, 0, x[i]);
		setResult("Profile "+i+1, 1, y[i]);
}
updateResults;

for (i=1; i<=nSlices; i++) {
	setSlice(i);
		for (j=0; j<x.length; j++) {
			makeOval(x[j]-25, y[j]-25, 50, 50);
getSelectionBounds(roiX, roiY, roiWidth, roiHeight);
roiManager("add");
sumX = 0;
sumY = 0;
sumIntensity = 0;
for (yy = roiY; yy < roiY + roiHeight; yy++) {
    for (xx = roiX; xx < roiX + roiWidth; xx++) {
        if (selectionContains(xx, yy)) { 
            value = getPixel(xx, yy); 
             power = 25; 
             weightedValue = Math.pow(value, power);
             sumX += xx * weightedValue;
             sumY += yy * weightedValue;
             sumIntensity += weightedValue;       
        }
    }
}
centroidX = sumX / sumIntensity;
centroidY = sumY / sumIntensity;
print("centroidX", centroidX, "centroidY", centroidY, "Slice", i, "Oval",j);
//makePoint(centroidX, centroidY);
//setForegroundColor(0,0,0);
//run("Draw","slice");
        smallSize = 8;
        x1 = centroidX - smallSize;
        y1 = centroidY - smallSize;
        width = 2*smallSize ;
        height = 2*smallSize ;
makeOval(x1, y1, width, height);
getRawStatistics(nPixels, Mean, min, max, std, histogram);
                                    roiManager("add"); 
			setResult("Profile "+j+1, i+1, Mean);
		}
	updateResults;
}
roiManager("Show All"); 






