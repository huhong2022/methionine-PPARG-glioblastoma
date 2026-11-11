// PPARg nuclear fluorescence quantification
// Fiji/ImageJ macro for paired DAPI and PPARg TIFF images.
//
// Supported file names:
//   *_CH1.tif[f]  DAPI
//   *_CH2.tif[f]  PPARg
//
// Supported exported channel files:
//   *c2x*.tif[f]  DAPI
//   *c1x*.tif[f]  PPARg
//   *c0-3x*.tif[f] merged images are ignored
//
// The macro requires Fiji with Bio-Formats. Images must be single-plane 2D.
// Review the saved DAPI masks and nuclear ROI overlays before using results.


// User parameters -----------------------------------------------------------

inputDir = getDirectory("Select the parent folder containing DAPI and PPARg TIFF files");
outputFolderName = "PPARg_Nuclear_Analysis";
outputDir = inputDir + outputFolderName + File.separator;
File.makeDirectory(outputDir);

maxFolderDepth = 4;
channelPaths = newArray(0);

useDapiRollingBackground = false;
rollingRadius = 35;
blurSigma = 3.0;

thresholdMethod = "Otsu dark";
dapiOtsuFactor = 0.70;
thresholdOffset = 0;

morphCloseRadius = 6;
useWatershed = false;

minNucleusArea = 800;
maxNucleusArea = 60000;
minCircularity = 0.05;

ringInnerPx = 3;
ringOuterPx = 12;
backgroundPercentile = 10;

if (ringOuterPx <= ringInnerPx)
    exit("ringOuterPx must be greater than ringInnerPx.");
if (minNucleusArea >= maxNucleusArea)
    exit("minNucleusArea must be smaller than maxNucleusArea.");
if (backgroundPercentile < 0 || backgroundPercentile > 100)
    exit("backgroundPercentile must be between 0 and 100.");


// Output files --------------------------------------------------------------

pairFile = outputDir + "Pairing_Report.csv";
statusFile = outputDir + "Processing_Status.csv";
cellFile = outputDir + "PPARg_Nuclear_Cell_Details.csv";
summaryFile = outputDir + "PPARg_Nuclear_Summary.csv";
settingsFile = outputDir + "Segmentation_Settings.txt";

pairCSV = "Sample,DAPI_File,PPARg_File,Pair_Status\n";
statusCSV = "Sample,Status,Note\n";

cellCSV = "Sample,ROI_ID,Nuclear_Area_px2,";
cellCSV += "Nuclear_PPARg_MFI,Perinuclear_PPARg_MFI,";
cellCSV += "Nuclear_Perinc_Ratio,";
cellCSV += "Nuclear_PPARg_Integrated,Background_PPARg\n";

summaryCSV = "Sample,Method,Detected_Nuclei,";
summaryCSV += "Measured_Nuclei,Background_PPARg,";
summaryCSV += "Mean_Nuclear_PPARg_MFI,";
summaryCSV += "Mean_Perinuclear_PPARg_MFI,";
summaryCSV += "Mean_Nuclear_Perinc_Ratio,";
summaryCSV += "Whole_Field_PPARg_MFI,";
summaryCSV += "DAPI_Threshold_8bit\n";

File.saveString(pairCSV, pairFile);
File.saveString(statusCSV, statusFile);
File.saveString(cellCSV, cellFile);
File.saveString(summaryCSV, summaryFile);

settings = "PPARg nuclear fluorescence quantification\n";
settings += "useDapiRollingBackground=" + useDapiRollingBackground + "\n";
settings += "rollingRadius=" + rollingRadius + "\n";
settings += "blurSigma=" + blurSigma + "\n";
settings += "thresholdMethod=" + thresholdMethod + "\n";
settings += "dapiOtsuFactor=" + dapiOtsuFactor + "\n";
settings += "thresholdOffset=" + thresholdOffset + "\n";
settings += "morphCloseRadius=" + morphCloseRadius + "\n";
settings += "useWatershed=" + useWatershed + "\n";
settings += "minNucleusArea=" + minNucleusArea + "\n";
settings += "maxNucleusArea=" + maxNucleusArea + "\n";
settings += "minCircularity=" + minCircularity + "\n";
settings += "ringInnerPx=" + ringInnerPx + "\n";
settings += "ringOuterPx=" + ringOuterPx + "\n";
settings += "backgroundPercentile=" + backgroundPercentile + "\n";
settings += "Nuclear and perinuclear values are background corrected.\n";
settings += "Inspect DAPI_MASK_CHECK and NUCLEAR_ROI_CHECK for every image.\n";
File.saveString(settings, settingsFile);


// Pair DAPI and PPARg files -------------------------------------------------

scanPairFiles(inputDir, 0);

readySamples = newArray(channelPaths.length);
readyDapi = newArray(channelPaths.length);
readyPparg = newArray(channelPaths.length);
readyTags = newArray(channelPaths.length);
readyCount = 0;

for (i = 0; i < channelPaths.length; i++) {
    dapiCandidate = channelPaths[i];
    lowerPath = toLowerCase(dapiCandidate);
    legacySuffix = "";

    if (endsWith(lowerPath, "_ch1.tif"))
        legacySuffix = "_ch1.tif";
    if (endsWith(lowerPath, "_ch1.tiff"))
        legacySuffix = "_ch1.tiff";

    isLegacy = (legacySuffix != "");
    isExport = false;
    targetLower = "";
    targetLowerTiff = "";

    if (isLegacy) {
        stemPath = substring(
            dapiCandidate,
            0,
            lengthOf(dapiCandidate) - lengthOf(legacySuffix)
        );
        sample = substring(stemPath, lengthOf(inputDir));
        targetLower = toLowerCase(stemPath + "_CH2.tif");
        targetLowerTiff = toLowerCase(stemPath + "_CH2.tiff");
    } else {
        markerPos = lastIndexOf(lowerPath, "c2x");
        separatorPos = lastIndexOf(lowerPath, File.separator);
        validExtension = endsWith(lowerPath, ".tif") ||
            endsWith(lowerPath, ".tiff");

        if (markerPos > separatorPos && validExtension &&
            indexOf(lowerPath, "c0-3x") < 0) {
            isExport = true;
            targetLower = substring(lowerPath, 0, markerPos) +
                "c1x" + substring(lowerPath, markerPos + 3);

            extensionLength = 4;
            if (endsWith(lowerPath, ".tiff"))
                extensionLength = 5;

            sampleStem = substring(
                dapiCandidate,
                0,
                lengthOf(dapiCandidate) - extensionLength
            );
            relativeStem = substring(sampleStem, lengthOf(inputDir));
            relativeMarker = lastIndexOf(toLowerCase(relativeStem), "c2x");
            sample = substring(relativeStem, 0, relativeMarker) +
                "x" + substring(relativeStem, relativeMarker + 3);
        }
    }

    if (!isLegacy && !isExport)
        continue;

    ppargMatch = "";
    for (j = 0; j < channelPaths.length; j++) {
        candidateLower = toLowerCase(channelPaths[j]);
        if (candidateLower == targetLower ||
            (isLegacy && candidateLower == targetLowerTiff)) {
            ppargMatch = channelPaths[j];
            break;
        }
    }

    if (ppargMatch == "") {
        pairCSV += sample + "," + dapiCandidate + ",,MISSING_PPARG\n";
        statusCSV += sample + ",SKIPPED,Missing PPARg channel\n";
    } else {
        readySamples[readyCount] = sample;
        readyDapi[readyCount] = dapiCandidate;
        readyPparg[readyCount] = ppargMatch;
        readyTags[readyCount] = safeSampleTag(sample);
        readyCount++;
        pairCSV += sample + "," + dapiCandidate + "," +
            ppargMatch + ",READY\n";
    }
}

File.saveString(pairCSV, pairFile);
File.saveString(statusCSV, statusFile);

if (readyCount == 0)
    exit("No matched DAPI and PPARg image pairs were found.");


// Batch analysis ------------------------------------------------------------

setBatchMode(true);
processed = 0;

for (s = 0; s < readyCount; s++) {
    sample = readySamples[s];
    sampleTag = readyTags[s];
    dapiPath = readyDapi[s];
    ppargPath = readyPparg[s];

    statusCSV += sample + ",STARTED,Opening images\n";
    File.saveString(statusCSV, statusFile);
    roiManager("Reset");

    dapiID = importIntensity(dapiPath, 3, "DAPI_Intensity");
    ppargID = importIntensity(ppargPath, 1, "PPARg_Intensity");

    selectImage(dapiID);
    imageWidth = getWidth();
    imageHeight = getHeight();
    run("Select None");
    getStatistics(tempArea, dapiMean);

    selectImage(ppargID);
    if (getWidth() != imageWidth || getHeight() != imageHeight) {
        statusCSV += sample + ",SKIPPED,DAPI and PPARg image-size mismatch\n";
        File.saveString(statusCSV, statusFile);
        selectImage(ppargID);
        close();
        selectImage(dapiID);
        close();
        continue;
    }

    run("Select None");
    getStatistics(tempArea, ppargMean);
    if (dapiMean <= 0 || ppargMean <= 0) {
        statusCSV += sample + ",SKIPPED,Empty DAPI or PPARg channel\n";
        File.saveString(statusCSV, statusFile);
        selectImage(ppargID);
        close();
        selectImage(dapiID);
        close();
        continue;
    }


    // Estimate the PPARg background from the selected image percentile.

    selectImage(ppargID);
    run("Select None");
    getHistogram(histogramValues, histogramCounts, 256);

    pixelTotal = 0;
    for (k = 0; k < histogramCounts.length; k++)
        pixelTotal += histogramCounts[k];

    percentileCutoff = pixelTotal * backgroundPercentile / 100;
    cumulativeCount = 0;
    backgroundPparg = histogramValues[0];

    for (k = 0; k < histogramCounts.length; k++) {
        cumulativeCount += histogramCounts[k];
        if (cumulativeCount >= percentileCutoff) {
            backgroundPparg = histogramValues[k];
            break;
        }
    }

    getStatistics(fieldArea, fieldMean);
    correctedFieldMean = fieldMean - backgroundPparg;
    if (correctedFieldMean < 0)
        correctedFieldMean = 0;


    // Segment nuclei on a copy of the DAPI image.

    selectImage(dapiID);
    run("Duplicate...", "title=DAPI_Mask_Work");
    maskID = getImageID();
    run("Select None");

    if (useDapiRollingBackground)
        run("Subtract Background...", "rolling=" + rollingRadius);

    run("Gaussian Blur...", "sigma=" + blurSigma);
    if (bitDepth() != 8)
        run("8-bit");

    setAutoThreshold(thresholdMethod);
    getThreshold(autoLow, autoHigh);
    selectedThreshold = round(autoLow * dapiOtsuFactor + thresholdOffset);

    if (selectedThreshold < 1)
        selectedThreshold = 1;
    if (selectedThreshold > 254)
        selectedThreshold = 254;

    setThreshold(selectedThreshold, 255);
    setOption("BlackBackground", true);
    run("Convert to Mask");

    selectImage(maskID);
    run("Select None");
    saveAs("Tiff", outputDir + sampleTag + "_DAPI_MASK_INITIAL.tif");

    if (morphCloseRadius > 0) {
        selectImage(maskID);
        run("Maximum...", "radius=" + morphCloseRadius);
        run("Minimum...", "radius=" + morphCloseRadius);
    }

    selectImage(maskID);
    setOption("BlackBackground", true);
    run("Fill Holes");

    if (useWatershed) {
        selectImage(maskID);
        setOption("BlackBackground", true);
        run("Watershed");
    }

    selectImage(maskID);
    saveAs("Tiff", outputDir + sampleTag + "_DAPI_MASK_CHECK.tif");


    // Extract nuclear ROIs.

    roiManager("Reset");
    selectImage(maskID);
    run(
        "Analyze Particles...",
        "size=" + minNucleusArea + "-" + maxNucleusArea +
        " pixel circularity=" + minCircularity + "-1.00" +
        " show=Nothing clear add exclude"
    );
    nucleiCount = roiManager("count");

    newImage(
        "Accepted_Nuclear_ROIs_QC",
        "8-bit black",
        imageWidth,
        imageHeight,
        1
    );
    acceptedMaskID = getImageID();
    setForegroundColor(255, 255, 255);

    for (r = 0; r < nucleiCount; r++) {
        selectImage(acceptedMaskID);
        roiManager("Select", r);
        run("Fill");
    }

    selectImage(acceptedMaskID);
    run("Select None");
    saveAs("Tiff", outputDir + sampleTag + "_DAPI_SELECTED_ROIS_MASK.tif");
    close();

    if (nucleiCount > 0) {
        roiManager("Save", outputDir + sampleTag + "_Nuclear_ROIs.zip");

        selectImage(dapiID);
        run("Duplicate...", "title=Nuclear_QC_Work");
        qcID = getImageID();
        run("Blue");
        run("RGB Color");
        qcID = getImageID();
        setForegroundColor(0, 255, 0);
        setLineWidth(1);

        for (r = 0; r < nucleiCount; r++) {
            selectImage(qcID);
            roiManager("Select", r);
            run("Draw");
        }

        selectImage(qcID);
        run("Select None");
        saveAs("Tiff", outputDir + sampleTag + "_NUCLEAR_ROI_CHECK.tif");
        close();
    }


    // Measure nuclear and perinuclear PPARg fluorescence.

    sumNuclear = 0;
    sumPerinuclear = 0;
    sumRatio = 0;
    measuredCount = 0;
    validRatioCount = 0;

    for (r = 0; r < nucleiCount; r++) {
        selectImage(ppargID);
        roiManager("Select", r);
        getSelectionBounds(roiX, roiY, roiWidth, roiHeight);

        if (roiX < ringOuterPx ||
            roiY < ringOuterPx ||
            roiX + roiWidth + ringOuterPx > imageWidth ||
            roiY + roiHeight + ringOuterPx > imageHeight)
            continue;

        getStatistics(nuclearArea, nuclearMean);

        roiManager("Select", r);
        run("Enlarge...", "enlarge=" + ringOuterPx);
        getStatistics(outerArea, outerMean);

        roiManager("Select", r);
        run("Enlarge...", "enlarge=" + ringInnerPx);
        getStatistics(innerArea, innerMean);

        ringArea = outerArea - innerArea;
        if (ringArea <= 0)
            continue;

        ringMean = (
            outerArea * outerMean -
            innerArea * innerMean
        ) / ringArea;

        correctedNuclear = nuclearMean - backgroundPparg;
        correctedRing = ringMean - backgroundPparg;
        if (correctedNuclear < 0)
            correctedNuclear = 0;
        if (correctedRing < 0)
            correctedRing = 0;

        nuclearIntegrated = correctedNuclear * nuclearArea;
        ratioText = "NA";

        if (correctedRing > 0) {
            nuclearPerinuclearRatio = correctedNuclear / correctedRing;
            ratioText = "" + nuclearPerinuclearRatio;
            sumRatio += nuclearPerinuclearRatio;
            validRatioCount++;
        }

        measuredCount++;
        sumNuclear += correctedNuclear;
        sumPerinuclear += correctedRing;

        cellCSV += sample + "," +
            (r + 1) + "," +
            nuclearArea + "," +
            correctedNuclear + "," +
            correctedRing + "," +
            ratioText + "," +
            nuclearIntegrated + "," +
            backgroundPparg + "\n";

        run("Select None");
    }


    // Save the image-level summary.

    meanNuclear = "NA";
    meanPerinuclear = "NA";
    meanRatio = "NA";

    if (measuredCount > 0) {
        meanNuclear = "" + (sumNuclear / measuredCount);
        meanPerinuclear = "" + (sumPerinuclear / measuredCount);
    }
    if (validRatioCount > 0)
        meanRatio = "" + (sumRatio / validRatioCount);

    summaryCSV += sample + ",Smooth_DAPI_Contour," +
        nucleiCount + "," +
        measuredCount + "," +
        backgroundPparg + "," +
        meanNuclear + "," +
        meanPerinuclear + "," +
        meanRatio + "," +
        correctedFieldMean + "," +
        selectedThreshold + "\n";

    File.saveString(cellCSV, cellFile);
    File.saveString(summaryCSV, summaryFile);

    if (measuredCount > 0)
        statusCSV += sample + ",DONE,Nuclear ROIs measured\n";
    else
        statusCSV += sample + ",DONE_NO_MEASUREMENTS,Review the DAPI mask\n";

    File.saveString(statusCSV, statusFile);

    roiManager("Reset");
    selectImage(maskID);
    close();
    selectImage(ppargID);
    close();
    selectImage(dapiID);
    close();
    processed++;
}


// Complete the batch --------------------------------------------------------

setBatchMode(false);
File.saveString(pairCSV, pairFile);
File.saveString(statusCSV, statusFile);
File.saveString(cellCSV, cellFile);
File.saveString(summaryCSV, summaryFile);

showMessage(
    "PPARg nuclear analysis",
    "Completed " + processed + " image pairs.\n" +
    "Results were saved to:\n" + outputDir
);


// Helper functions ----------------------------------------------------------

function importIntensity(path, component, imageTitle) {
    imageCountBefore = nImages();

    run(
        "Bio-Formats Importer",
        "open=[" + path + "] " +
        "ungroupFiles=true group_files=false " +
        "autoscale=false color_mode=Default " +
        "view=Hyperstack stack_order=XYCZT"
    );

    if (nImages() != imageCountBefore + 1)
        exit("Expected one image from: " + path);

    sourceID = getImageID();
    getDimensions(width, height, channels, slices, frames);
    sourceBitDepth = bitDepth();

    if (slices != 1 || frames != 1)
        exit("Only single-plane 2D images are supported: " + path);

    if (sourceBitDepth == 24) {
        run("RGB Stack");
        setSlice(component);
        run("Duplicate...", "title=" + imageTitle);
        resultID = getImageID();
        selectImage(sourceID);
        close();
    } else if (channels > 1) {
        if (channels < component)
            exit("Required channel is missing from: " + path);

        Stack.setPosition(component, 1, 1);
        run(
            "Duplicate...",
            "title=" + imageTitle +
            " duplicate channels=" + component + "-" + component +
            " slices=1-1 frames=1-1"
        );
        resultID = getImageID();
        selectImage(sourceID);
        close();
    } else {
        rename(imageTitle);
        resultID = getImageID();
    }

    selectImage(resultID);
    getDimensions(width, height, channels, slices, frames);
    if (channels != 1 || slices != 1 || frames != 1)
        exit("Channel extraction failed for: " + path);

    run("Set Scale...", "distance=1 known=1 unit=pixel");
    return resultID;
}


function scanPairFiles(folder, depth) {
    if (depth > maxFolderDepth)
        return;

    entries = getFileList(folder);
    Array.sort(entries);

    for (u = 0; u < entries.length; u++) {
        entry = entries[u];
        path = folder + entry;
        lowerEntry = toLowerCase(entry);

        if (File.isDirectory(path)) {
            if (indexOf(lowerEntry, toLowerCase(outputFolderName)) >= 0)
                continue;
            scanPairFiles(path, depth + 1);
        } else {
            validExtension = endsWith(lowerEntry, ".tif") ||
                endsWith(lowerEntry, ".tiff");
            legacyChannel = endsWith(lowerEntry, "_ch1.tif") ||
                endsWith(lowerEntry, "_ch1.tiff") ||
                endsWith(lowerEntry, "_ch2.tif") ||
                endsWith(lowerEntry, "_ch2.tiff");
            exportedChannel = indexOf(lowerEntry, "c1x") >= 0 ||
                indexOf(lowerEntry, "c2x") >= 0;
            mergedImage = indexOf(lowerEntry, "c0-3x") >= 0;

            if (validExtension && !mergedImage &&
                (legacyChannel || exportedChannel)) {
                channelPaths = Array.concat(
                    channelPaths,
                    newArray(path)
                );
            }
        }
    }
}


function safeSampleTag(sampleName) {
    sampleName = replace(sampleName, "\\", "__");
    sampleName = replace(sampleName, "/", "__");
    sampleName = replace(sampleName, ":", "_");
    sampleName = replace(sampleName, ",", "_");
    return sampleName;
}
