// Batch CD68/CD206 colocalization analysis using JACoP.
// Required file names: Overlay001 through Overlay012, paired as ch01/ch02.
// Image A = ch01 = green = CD206.
// Image B = ch02 = red = CD68.
// The original JACoP log is retained for each pair, and all metrics are
// combined in one CSV file.

dir = getDirectory("Select the folder containing the original ch01/ch02 TIFF files");

if (endsWith(dir, "CD68_CD206_JACoP_Analysis" + File.separator)) {
    dir = substring(dir, 0, lengthOf(dir) - lengthOf("CD68_CD206_JACoP_Analysis" + File.separator));
    print("The output folder was selected; using its parent folder instead:");
    print(dir);
}

ext = ".tif";
nPairs = 12;

outDir = dir + "CD68_CD206_JACoP_Analysis" + File.separator;
File.makeDirectory(outDir);

greenThreshold = 0;
redThreshold = 0;

jacopOptions = " pearson overlap mm";

summaryFile = outDir + "JACoP_summary.csv";

header = "CD206_ImageA_green,CD68_ImageB_red,Pearson_r,Overlap_r_original,k1_original,k2_original,Overlap_r_threshold,k1_threshold,k2_threshold,M1_CD206_original,M2_CD68_original,M1_CD206_threshold,M2_CD68_threshold,CD206_threshold,CD68_threshold,JACoP_log_file\n";
File.saveString(header, summaryFile);

setBatchMode(false);

for (i = 1; i <= nPairs; i++) {

    id = IJ.pad(i, 2);

    // Overlay001_ch01.tif / Overlay001_ch02.tif
    greenName = "Overlay0" + id + "_ch01" + ext;
    redName   = "Overlay0" + id + "_ch02" + ext;

    greenFile = dir + greenName;
    redFile   = dir + redName;

    if (!File.exists(greenFile)) {
        print("Missing Image A / green image: " + greenFile);
        continue;
    }

    if (!File.exists(redFile)) {
        print("Missing Image B / red image: " + redFile);
        continue;
    }

    open(greenFile);
    run("8-bit");
    rename(greenName);

    open(redFile);
    run("8-bit");
    rename(redName);

    print("\\Clear");

    run("JACoP ",
        "imga=[" + greenName + "] " +
        "imgb=[" + redName + "] " +
        "thra=" + greenThreshold + " " +
        "thrb=" + redThreshold +
        jacopOptions
    );

    wait(1000);

    logFileName = "Overlay0" + id + "_JACoP_CD206-green_CD68-red.txt";
    logFile = outDir + logFileName;

    if (isOpen("Log")) {

        selectWindow("Log");
        saveAs("Text", logFile);

        // Read the JACoP results from the saved text file.
        logText = File.openAsString(logFile);

        pearson_r = "NA";
        overlap_r_original = "NA";
        k1_original = "NA";
        k2_original = "NA";
        overlap_r_threshold = "NA";
        k1_threshold = "NA";
        k2_threshold = "NA";
        m1_original = "NA";
        m2_original = "NA";
        m1_threshold = "NA";
        m2_threshold = "NA";

        lines = split(logText, "\n");

        overlapCount = 0;
        kCount = 0;
        inPearson = false;
        inMandersOriginal = false;
        inMandersThreshold = false;

        for (j = 0; j < lines.length; j++) {

            line = trim(lines[j]);

            // Pearson's Coefficient:
            if (startsWith(line, "Pearson")) {
                inPearson = true;
                continue;
            }

            if (inPearson && startsWith(line, "r=")) {
                pearson_r = substring(line, 2, lengthOf(line));
                inPearson = false;
                continue;
            }

            // Overlap Coefficient:
            if (startsWith(line, "Overlap Coefficient:")) {
                overlapCount++;

                if (j + 1 < lines.length) {
                    rLine = trim(lines[j + 1]);

                    if (startsWith(rLine, "r=")) {
                        val = substring(rLine, 2, lengthOf(rLine));

                        if (overlapCount == 1) {
                            overlap_r_original = val;
                        }

                        if (overlapCount == 2) {
                            overlap_r_threshold = val;
                        }
                    }
                }
            }

            // r^2=k1xk2:
            if (startsWith(line, "r^2=k1xk2:")) {
                kCount++;

                if (j + 1 < lines.length) {
                    k1Line = trim(lines[j + 1]);

                    if (startsWith(k1Line, "k1=")) {
                        val = substring(k1Line, 3, lengthOf(k1Line));

                        if (kCount == 1) {
                            k1_original = val;
                        }

                        if (kCount == 2) {
                            k1_threshold = val;
                        }
                    }
                }

                if (j + 2 < lines.length) {
                    k2Line = trim(lines[j + 2]);

                    if (startsWith(k2Line, "k2=")) {
                        val = substring(k2Line, 3, lengthOf(k2Line));

                        if (kCount == 1) {
                            k2_original = val;
                        }

                        if (kCount == 2) {
                            k2_threshold = val;
                        }
                    }
                }
            }

            // Manders original
            if (startsWith(line, "Manders' Coefficients (original):")) {
                inMandersOriginal = true;
                inMandersThreshold = false;
                continue;
            }

            // Manders threshold
            if (startsWith(line, "Manders' Coefficients (using threshold")) {
                inMandersOriginal = false;
                inMandersThreshold = true;
                continue;
            }

            if (startsWith(line, "M1=")) {
                val = substring(line, 3, lengthOf(line));
                spacePos = indexOf(val, " ");

                if (spacePos >= 0) {
                    val = substring(val, 0, spacePos);
                }

                if (inMandersOriginal) {
                    m1_original = val;
                }

                if (inMandersThreshold) {
                    m1_threshold = val;
                }
            }

            if (startsWith(line, "M2=")) {
                val = substring(line, 3, lengthOf(line));
                spacePos = indexOf(val, " ");

                if (spacePos >= 0) {
                    val = substring(val, 0, spacePos);
                }

                if (inMandersOriginal) {
                    m2_original = val;
                }

                if (inMandersThreshold) {
                    m2_threshold = val;
                }
            }
        }

        row = greenName + "," +
              redName + "," +
              pearson_r + "," +
              overlap_r_original + "," +
              k1_original + "," +
              k2_original + "," +
              overlap_r_threshold + "," +
              k1_threshold + "," +
              k2_threshold + "," +
              m1_original + "," +
              m2_original + "," +
              m1_threshold + "," +
              m2_threshold + "," +
              greenThreshold + "," +
              redThreshold + "," +
              logFileName + "\n";

        File.append(row, summaryFile);

    } else {
        print("No JACoP Log output was detected for pair " + id + ".");
    }

    if (isOpen(greenName)) {
        selectWindow(greenName);
        close();
    }

    if (isOpen(redName)) {
        selectWindow(redName);
        close();
    }
}

print("Batch JACoP analysis finished.");
print("Image A = ch01 = green = CD206.");
print("Image B = ch02 = red = CD68.");
print("Individual JACoP logs saved to: " + outDir);
print("Combined CSV summary saved to: " + summaryFile);
