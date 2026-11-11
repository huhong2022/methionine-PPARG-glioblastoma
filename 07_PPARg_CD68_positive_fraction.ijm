// PPARg-positive fraction among CD68-positive cells
// Fiji/ImageJ macro for three aligned single-plane TIFF images per field.
//
// Required file names:
//   *_ch00.tif[f]  DAPI
//   *_ch01.tif[f]  CD68
//   *_ch02.tif[f]  PPARg
//
// Endpoint:
//   PPARg-positive CD68-positive cells / all CD68-positive cells * 100
//
// Review every DAPI center-detection image before accepting the results.


// User parameters -----------------------------------------------------------

inputDir = getDirectory("Select the folder containing ch00, ch01, and ch02 TIFF files");
outputFolderName = "PPARg_in_CD68_Analysis";
outputDir = inputDir + outputFolderName + File.separator;
File.makeDirectory(outputDir);

dapiProminence = 40;
dapiBackgroundSigma = 50;

ppargRadius = 20;
cd68Radius = 10;

cd68PositiveThreshold = 10;
ppargPositiveThreshold = 10;

edgeMargin = 30;

if (dapiProminence < 0)
    exit("dapiProminence must be non-negative.");
if (ppargRadius <= 0 || cd68Radius <= 0)
    exit("Measurement radii must be greater than zero.");
if (edgeMargin < maxOf(ppargRadius, cd68Radius))
    exit("edgeMargin must be at least as large as both measurement radii.");


// Output files --------------------------------------------------------------

summaryFile = outputDir + "PPARg_in_CD68_summary.csv";
detailFile = outputDir + "PPARg_in_CD68_cell_details.csv";
statusFile = outputDir + "Processing_Status.csv";
settingsFile = outputDir + "Analysis_Settings.txt";

summaryCSV = "Image,Analyzed_DAPI_centers,CD68_positive_cells,";
summaryCSV += "PPARg_positive_CD68_positive_cells,";
summaryCSV += "PPARg_positive_fraction_percent,";
summaryCSV += "CD68_threshold,PPARg_threshold\n";

detailCSV = "Image,Nucleus_ID,X,Y,CD68_mean,PPARg_mean,";
detailCSV += "Is_CD68_positive,Is_PPARG_positive\n";

statusCSV = "Image,Status,Note\n";

settings = "PPARg-positive fraction among CD68-positive cells\n";
settings += "DAPI channel=ch00\n";
settings += "CD68 channel=ch01\n";
settings += "PPARg channel=ch02\n";
settings += "dapiProminence=" + dapiProminence + "\n";
settings += "dapiBackgroundSigma=" + dapiBackgroundSigma + "\n";
settings += "ppargRadius=" + ppargRadius + "\n";
settings += "cd68Radius=" + cd68Radius + "\n";
settings += "cd68PositiveThreshold=" + cd68PositiveThreshold + "\n";
settings += "ppargPositiveThreshold=" + ppargPositiveThreshold + "\n";
settings += "edgeMargin=" + edgeMargin + "\n";
settings += "Intensity measurements use 8-bit images.\n";
settings += "Inspect every DAPI_MAXIMA_CHECK image before analysis acceptance.\n";

File.saveString(summaryCSV, summaryFile);
File.saveString(detailCSV, detailFile);
File.saveString(statusCSV, statusFile);
File.saveString(settings, settingsFile);


// Locate and analyze image triplets ----------------------------------------

fileList = getFileList(inputDir);
Array.sort(fileList);
processedCount = 0;

setBatchMode(false);

for (i = 0; i < fileList.length; i++) {
    fileName = fileList[i];
    lowerName = toLowerCase(fileName);
    dapiSuffix = "";

    if (endsWith(lowerName, "_ch00.tif"))
        dapiSuffix = "_ch00.tif";
    if (endsWith(lowerName, "_ch00.tiff"))
        dapiSuffix = "_ch00.tiff";

    if (dapiSuffix == "")
        continue;

    baseName = substring(
        fileName,
        0,
        lengthOf(fileName) - lengthOf(dapiSuffix)
    );

    dapiPath = inputDir + fileName;
    cd68Path = findChannelFile(inputDir, fileList, baseName, "_ch01");
    ppargPath = findChannelFile(inputDir, fileList, baseName, "_ch02");

    if (cd68Path == "") {
        statusCSV += baseName + ",SKIPPED,Missing CD68 ch01 image\n";
        File.saveString(statusCSV, statusFile);
        continue;
    }

    if (ppargPath == "") {
        statusCSV += baseName + ",SKIPPED,Missing PPARg ch02 image\n";
        File.saveString(statusCSV, statusFile);
        continue;
    }

    statusCSV += baseName + ",STARTED,Opening image triplet\n";
    File.saveString(statusCSV, statusFile);


    // Open aligned intensity images.

    dapiID = openIntensityImage(dapiPath, 3, baseName + "_DAPI");
    cd68ID = openIntensityImage(cd68Path, 2, baseName + "_CD68");
    ppargID = openIntensityImage(ppargPath, 1, baseName + "_PPARg");

    selectImage(dapiID);
    imageWidth = getWidth();
    imageHeight = getHeight();

    selectImage(cd68ID);
    cd68SizeMatches = (
        getWidth() == imageWidth &&
        getHeight() == imageHeight
    );

    selectImage(ppargID);
    ppargSizeMatches = (
        getWidth() == imageWidth &&
        getHeight() == imageHeight
    );

    if (!cd68SizeMatches || !ppargSizeMatches) {
        statusCSV += baseName + ",SKIPPED,Channel image-size mismatch\n";
        File.saveString(statusCSV, statusFile);
        closeImageIfOpen(dapiID);
        closeImageIfOpen(cd68ID);
        closeImageIfOpen(ppargID);
        continue;
    }


    // Detect DAPI nuclear centers after subtraction of a smooth background.

    selectImage(dapiID);
    dapiForMaximaTitle = baseName + "_DAPI_for_maxima";
    run("Duplicate...", "title=" + dapiForMaximaTitle);
    dapiForMaximaID = getImageID();

    selectImage(dapiID);
    dapiBackgroundTitle = baseName + "_DAPI_background";
    run("Duplicate...", "title=" + dapiBackgroundTitle);
    dapiBackgroundID = getImageID();
    run("Gaussian Blur...", "sigma=" + dapiBackgroundSigma);

    imageCalculator(
        "Subtract create",
        dapiForMaximaTitle,
        dapiBackgroundTitle
    );
    dapiSubtractedID = getImageID();
    rename(baseName + "_DAPI_subtracted");

    closeImageIfOpen(dapiForMaximaID);
    closeImageIfOpen(dapiBackgroundID);

    selectImage(dapiSubtractedID);
    if (bitDepth() != 8)
        run("8-bit");
    run("Gaussian Blur...", "sigma=1");
    run(
        "Find Maxima...",
        "prominence=" + dapiProminence +
        " output=[Point Selection] exclude"
    );

    if (selectionType() == 10) {
        getSelectionCoordinates(xPoints, yPoints);
    } else {
        xPoints = newArray(0);
        yPoints = newArray(0);
    }
    detectedCenters = xPoints.length;


    // Save a DAPI center-detection quality-control image.

    selectImage(dapiID);
    run("Duplicate...", "title=" + baseName + "_DAPI_MAXIMA_CHECK");
    qcID = getImageID();
    run("RGB Color");
    qcID = getImageID();
    setLineWidth(2);
    setForegroundColor(255, 0, 0);

    for (p = 0; p < detectedCenters; p++) {
        makeOval(xPoints[p] - 5, yPoints[p] - 5, 10, 10);
        run("Draw");
    }

    run("Select None");
    saveAs("Tiff", outputDir + safeFileName(baseName) + "_DAPI_MAXIMA_CHECK.tif");
    closeImageIfOpen(qcID);


    // Classify CD68 and PPARg status at each accepted DAPI center.

    analyzedCenters = 0;
    cd68PositiveCount = 0;
    doublePositiveCount = 0;

    for (p = 0; p < detectedCenters; p++) {
        centerX = xPoints[p];
        centerY = yPoints[p];

        if (centerX < edgeMargin ||
            centerY < edgeMargin ||
            centerX > imageWidth - edgeMargin ||
            centerY > imageHeight - edgeMargin)
            continue;

        analyzedCenters++;

        selectImage(ppargID);
        makeOval(
            centerX - ppargRadius,
            centerY - ppargRadius,
            2 * ppargRadius,
            2 * ppargRadius
        );
        getStatistics(ppargArea, ppargMean);

        isPpargPositive = 0;
        if (ppargMean >= ppargPositiveThreshold)
            isPpargPositive = 1;

        selectImage(cd68ID);
        makeOval(
            centerX - cd68Radius,
            centerY - cd68Radius,
            2 * cd68Radius,
            2 * cd68Radius
        );
        getStatistics(cd68Area, cd68Mean);

        isCd68Positive = 0;
        if (cd68Mean >= cd68PositiveThreshold) {
            isCd68Positive = 1;
            cd68PositiveCount++;
        }

        if (isCd68Positive == 1 && isPpargPositive == 1)
            doublePositiveCount++;

        detailCSV += baseName + "," +
            analyzedCenters + "," +
            centerX + "," +
            centerY + "," +
            cd68Mean + "," +
            ppargMean + "," +
            isCd68Positive + "," +
            isPpargPositive + "\n";

        run("Select None");
    }

    positiveFraction = "NA";
    if (cd68PositiveCount > 0)
        positiveFraction = "" + (
            doublePositiveCount / cd68PositiveCount * 100
        );

    summaryCSV += baseName + "," +
        analyzedCenters + "," +
        cd68PositiveCount + "," +
        doublePositiveCount + "," +
        positiveFraction + "," +
        cd68PositiveThreshold + "," +
        ppargPositiveThreshold + "\n";

    if (analyzedCenters > 0) {
        statusCSV += baseName + ",DONE,Image triplet analyzed\n";
    } else {
        statusCSV += baseName + ",DONE_NO_CENTERS,Review DAPI center detection\n";
    }

    File.saveString(summaryCSV, summaryFile);
    File.saveString(detailCSV, detailFile);
    File.saveString(statusCSV, statusFile);

    closeImageIfOpen(dapiSubtractedID);
    closeImageIfOpen(dapiID);
    closeImageIfOpen(cd68ID);
    closeImageIfOpen(ppargID);
    processedCount++;
}


// Complete the batch --------------------------------------------------------

File.saveString(summaryCSV, summaryFile);
File.saveString(detailCSV, detailFile);
File.saveString(statusCSV, statusFile);

showMessage(
    "PPARg in CD68-positive cells",
    "Completed " + processedCount + " image triplets.\n" +
    "Results were saved to:\n" + outputDir
);


// Helper functions ----------------------------------------------------------

function findChannelFile(folder, entries, baseName, channelSuffix) {
    targetTif = toLowerCase(baseName + channelSuffix + ".tif");
    targetTiff = toLowerCase(baseName + channelSuffix + ".tiff");

    for (index = 0; index < entries.length; index++) {
        candidate = toLowerCase(entries[index]);
        if (candidate == targetTif || candidate == targetTiff)
            return folder + entries[index];
    }
    return "";
}


function openIntensityImage(path, rgbComponent, imageTitle) {
    open(path);
    sourceID = getImageID();

    if (bitDepth() == 24) {
        run("RGB Stack");
        setSlice(rgbComponent);
        run("Duplicate...", "title=" + imageTitle);
        resultID = getImageID();
        closeImageIfOpen(sourceID);
    } else {
        rename(imageTitle);
        resultID = getImageID();
    }

    selectImage(resultID);
    if (bitDepth() != 8)
        run("8-bit");
    run("Set Scale...", "distance=1 known=1 unit=pixel");
    return resultID;
}


function closeImageIfOpen(imageID) {
    if (imageID > 0) {
        selectImage(imageID);
        close();
    }
}


function safeFileName(value) {
    value = replace(value, "\\", "__");
    value = replace(value, "/", "__");
    value = replace(value, ":", "_");
    value = replace(value, ",", "_");
    return value;
}
