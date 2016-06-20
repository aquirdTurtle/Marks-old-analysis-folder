# AutoanalysisFunctions.py
# This File contains all of the main functions used by my c++ code that controls the experiment.

def pairAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    from numpy import array, zeros
    import ctypes  # An included library with Python install.
    
    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\My Data Analysis")
    from astropy.io import fits
    import numpy
    numpy.set_printoptions(threshold=numpy.nan)
    from matplotlib.pyplot import  subplots, show, gcf
    from matplotlib.cm import get_cmap
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian, 
                                       calculateAtomThreshold, getAnalyzedTunnelingData);    
    #
    dataRepositoryPath = "\\\\andor\\share\\Data and documents\\Data repository\\";
    todaysDataPath = dataRepositoryPath + date + "\\data_" + str(runNumber) + ".fits";
    keyPath = dataRepositoryPath + date + "\\key_" + str(runNumber) + ".txt";
    # Load Key
    key = numpy.array([]);
    with open(keyPath) as keyFile:
        for line in keyFile:
            key = numpy.append(key, float(line.strip('\n')))    
    # Load Fits File & Get Dimensions
    # Get the array from the fits file. That's all I care about.
    fitsInfo = (fits.open(todaysDataPath, "append"));
    rawData = fitsInfo[0].data;
    # the .shape member of an array gives an array of the dimesnions of the array.
    numberOfPictures = rawData.shape[0];
    numberOfExperiments = int(numberOfPictures / picturesPerExperiment);
    #
    numberAtomsToAnalyze = (array(analysisLocations).shape)[0];
    for atomInc in range(0, int(numberAtomsToAnalyze/4)):
        location1 = array([analysisLocations[4 * atomInc], analysisLocations[4 * atomInc + 1]]);        
        location2 = array([analysisLocations[4 * atomInc + 2], analysisLocations[4 * atomInc + 3]]);        
        allAtomData = [[],[]];
        firstExperimentData = [[],[]];
        allAtomData[0], firstExperimentData[0] = normalizeData(picturesPerExperiment, rawData, location1);
        allAtomData[1], firstExperimentData[1] = normalizeData(picturesPerExperiment, rawData, location2);
        binCenters = [[],[]];
        binnedData = [[],[]];
        binCenters[0], binnedData[0] = binData(5, allAtomData[0]);
        binCenters[1], binnedData[1] = binData(5, allAtomData[1]);
        guessLocation1 = [[],[]];
        guessLocation2 = [[],[]];
        guessLocation1[0], guessLocation2[0] = guessGaussianPeaks(allAtomData[0], binCenters[0], binnedData[0]);
        guessLocation1[1], guessLocation2[1] = guessGaussianPeaks(allAtomData[0], binCenters[1], binnedData[1]);
        guess = [[],[]];
        gaussianFitVals = [[],[]];
        thresholds = [[],[]];
        thresholdFidelity = [[],[]];
        guess[0] = numpy.array([100, guessLocation1[0], 30, 200, guessLocation2[0], 10]);
        guess[1] = numpy.array([100, guessLocation1[1], 30, 200, guessLocation2[1], 10]);
        gaussianFitVals[0] = fitDoubleGaussian(binCenters[0], binnedData[0], guess[0]);
        gaussianFitVals[1] = fitDoubleGaussian(binCenters[1], binnedData[1], guess[1]);
        thresholds[0], thresholdFidelity[0] = calculateAtomThreshold(gaussianFitVals[0]);
        thresholds[1], thresholdFidelity[1] = calculateAtomThreshold(gaussianFitVals[1]);
        #
        atomCount1 = 0;
        atomCount2 = 0;
        for experimentInc in range(0, firstExperimentData[0].size):
            if firstExperimentData[0][experimentInc] > thresholds[0]:
                atomCount1 += 1;
        for experimentInc in range(0, firstExperimentData[0].size):
            if firstExperimentData[1][experimentInc] > thresholds[1]:
                atomCount2 += 1;
        firstToFirst, firstToSecond, bothToBoth, bothToOne, captureData1, captureData2, summedData \
            = getAnalyzedTunnelingData(allAtomData, thresholds, key, accumulations, numberOfExperiments);
        myFigure, ((plot11, plot12, plot13), (plot21, plot22, plot23)) = subplots(2, 3, figsize = (25,12))
        myFigure.suptitle("Data for locations {" + str(location1[0] + 1) + "," 
                          + str(location1[1] + 1) + "}, {" + str(location1[0] + 1) 
                          + "," + str(location1[1] + 1) + "}", fontsize = 24);
        figObject = gcf()
        figObject.canvas.set_window_title("{" + str(location1[0] + 1) + "," 
                          + str(location1[1] + 1) + "}, {" + str(location1[0] + 1) 
                          + "," + str(location1[1] + 1) + "}")
        # make an image
        plot11.imshow(rawData[10], interpolation='none', cmap = get_cmap("bone"));
        plot11.set_title("Example Raw Image")
        
        # First counts histogram
        plot12.hist(firstExperimentData[0], 50, color='r');
        plot12.hist(firstExperimentData[1], 50, color='b');
        plot12.set_title("First Picture of Experiment Pixel Count Histogram");
        plot12.set_ylabel("Occurance Count");
        plot12.set_xlabel("Pixel Counts");
        
        # plot loading probabilities
        plot13.plot(captureData1[:,0], captureData1[:,1], linestyle = "none", color='r');
        plot13.plot(captureData2[:,0], captureData2[:,1], linestyle = "none", color='b');
        plot13.set_title("Average Loading Efficiencies: " + str(atomCount1 / firstExperimentData[0].size * 100) 
                         + "% and " + str(atomCount2 / firstExperimentData[1].size * 100) + "%");
        plot13.set_ylabel("Occurance Count");
        plot13.set_xlabel("Pixel Counts");
        
        # plot the fit on top of the histogram
        plot21.plot(binCenters[0], binnedData[0], "o", markersize = 3, color='r');
        plot21.plot(binCenters[1], binnedData[1], "o", markersize = 3, color='b');
        
        fineXData1 = numpy.linspace(min(allAtomData[0]), max(allAtomData[0]), 500);
        fineXData2 = numpy.linspace(min(allAtomData[1]), max(allAtomData[1]), 500);

        plot21.plot(fineXData1, doubleGaussian(fineXData1, *gaussianFitVals[0]), color='r');
        plot21.plot(fineXData2, doubleGaussian(fineXData2, *gaussianFitVals[1]), color='b');
        
        plot21.set_title("Fits. Threshold (red) = " + str(thresholds[0]) + "Threshold (blue) = " + str(thresholds[1]));
        plot21.axvline(x=thresholds[0], color='r', ls='dashed');
        plot21.axvline(x=thresholds[1], color='b', ls='dashed');
        plot21.set_ylabel("Occurance Count");
        plot21.set_xlabel("Pixel Counts");
        plot22.plot(allAtomData[0], ".", markersize = 1, color='r');
        plot22.plot(allAtomData[1], ".", markersize = 1, color='b');
        plot22.set_title("Pixel Count Data Over Time");
        plot22.set_ylabel("Count on Pixel");
        plot22.set_xlabel("Picture Number");

        # Plot survival data.
        #ctypes.windll.user32.MessageBoxW(0, "Made it.", "", 1)
        plot23.errorbar(firstToFirst[:,0], firstToFirst[:,1], yerr=firstToFirst[:,2], linestyle = "none", color = "r", fmt = "o", label = "first-to-first");

        plot23.errorbar(firstToSecond[:,0], firstToSecond[:,1], yerr=firstToSecond[:,2], linestyle = "none", color = "b", fmt = "o", label = "first-to-second");
        
        plot23.errorbar(bothToBoth[:,0], bothToBoth[:,1], yerr=bothToBoth[:,2], linestyle = "none", color = "c", fmt = "o", label = "both-to-both");
        
        plot23.errorbar(bothToOne[:,0], bothToOne[:,1], yerr=bothToOne[:,2], linestyle = "none", color = "g", fmt = "o", label = "both-to-one");
        
        plot23.errorbar(summedData[:,0],summedData[:,1], yerr = summedData[:,2], linestyle = "none", color = "m", fmt = "o", label = "first-survival");
        
        plot23.set_title("Two-Particle Data")
        plot23.set_ylabel("% Occurred")
        plot23.set_xlabel("Variation Parameter Value")
        plot23.set_ylim([-0.05, 1.05])
        xRange = max(firstToFirst[:,0]) - min(firstToFirst[:,0]);
        plot23.set_xlim([min(firstToFirst[:,0]) - xRange / 10.0, max(firstToFirst[:,0]) + xRange / 10.0]);
        plot23.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        # Export data exactly how it used to be exported using mathematica, so that a mathematica file can read it. 
        fitsInfo.close()
        #collatedData = {accumulations, peakData, key, toPlot, captProbs, dataRaw, toPlot2, toPlotBoth, toPlotBoth11, 
            #Date[]};
        outputName = (dataRepositoryPath + date + "\\" + fileName + "_run" + str(runNumber) + "_pair"
                      + "(" + str(location1[0] + 1) + "," + str(location1[1] + 1) + "), (" + str(location1[0] + 1) 
                          + "," + str(location1[1] + 1) + ")" + ".tsv");
        with open(outputName, "w") as record_file:
            # accumulations is special since it's an integer.
            record_file.write(str(accumulations) + "\n");
            ### Version for mathematica compatibility
            # Peak Data
            for peakInc in range(0, allAtomData[0].size):
                record_file.write(str(allAtomData[0][peakInc]) + "\t");
            record_file.write("\n");
            # key
            #ctypes.windll.user32.MessageBoxW(0, "Made it3.", "", 1)
            for keyInc in range(0, key.size-1):
                record_file.write(str(key[keyInc]) + "\t");
            record_file.write(str(key[key.size-1]));
            record_file.write("\n");
            # firstToFirst
            firstToFirstDimensions = firstToFirst.shape;
            for firstToFirstPoints in range(0, firstToFirstDimensions[0]-1):
                record_file.write(str("{{"));
                record_file.write(str(firstToFirst[firstToFirstPoints][0]) + ", ");
                record_file.write(str(firstToFirst[firstToFirstPoints][1]) + "}, ");
                record_file.write("ErrorBar[" + str(firstToFirst[firstToFirstPoints][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(firstToFirst[firstToFirstDimensions[0]-1][0]) + ", ");
            record_file.write(str(firstToFirst[firstToFirstDimensions[0]-1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(firstToFirst[firstToFirstDimensions[0]-1][2]) + "]");
            record_file.write(str("}"));
            record_file.write("\n")
            # capture probabilities data
            #ctypes.windll.user32.MessageBoxW(0, "capture data size:" + str(captureData1.shape) + ",", "", 1)
            for captureInc in range(0, captureData1.shape[0]):
                #ctypes.windll.user32.MessageBoxW(0, str(captureInc), "", 1)
                record_file.write(str(captureData1[captureInc, 1]) + " ");
            record_file.write("\n");
            #ctypes.windll.user32.MessageBoxW(0, "Made it5.", "", 1)
            # raw data
            rawDataDimensions = rawData.shape;
            for pictureInc in range(0, rawDataDimensions[0]):
                record_file.write("{");
                for rowInc in range(0, rawDataDimensions[1]):
                    record_file.write("{");
                    for columnInc in range(0, rawDataDimensions[2]):
                        record_file.write(str(rawData[pictureInc][rowInc][columnInc]) + " ")
                    record_file.write("} ")
                record_file.write("}");
            record_file.write("\n");
            #ctypes.windll.user32.MessageBoxW(0, "Made it.6", "", 1)
            # firstToSecond
            firstToSecondDimensions = firstToSecond.shape;
            for firstToSecondPoints in range(0, firstToSecondDimensions[0]-1):
                record_file.write(str("{{"));
                record_file.write(str(firstToSecond[firstToSecondPoints][0]) + ", ");
                record_file.write(str(firstToSecond[firstToSecondPoints][1]) + "}, ");
                record_file.write("ErrorBar[" + str(firstToSecond[firstToSecondPoints][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(firstToSecond[firstToSecondDimensions[0]-1][0]) + ", ");
            record_file.write(str(firstToSecond[firstToSecondDimensions[0]-1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(firstToSecond[firstToSecondDimensions[0]-1][2]) + "]");
            record_file.write(str("}"));
            record_file.write("\n")
            #ctypes.windll.user32.MessageBoxW(0, "Made it7.", "", 1)
            # bothToBoth
            bothToBothDimensions = bothToBoth.shape;
            for bothToBothPoints in range(0, bothToBothDimensions[0]-1):
                record_file.write(str("{{"));
                record_file.write(str(bothToBoth[bothToBothPoints][0]) + ", ");
                record_file.write(str(bothToBoth[bothToBothPoints][1]) + "}, ");
                record_file.write("ErrorBar[" + str(bothToBoth[bothToBothPoints][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(bothToBoth[bothToBothDimensions[0]-1][0]) + ", ");
            record_file.write(str(bothToBoth[bothToBothDimensions[0]-1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(bothToBoth[bothToBothDimensions[0]-1][2]) + "]");
            record_file.write(str("}"));
            record_file.write("\n")
            #ctypes.windll.user32.MessageBoxW(0, "Made iti8.", "", 1)
            # bothToOne
            bothToOneDimensions = bothToOne.shape;
            for bothToOnePoints in range(0, bothToOneDimensions[0]-1):
                record_file.write(str("{{"));
                record_file.write(str(bothToOne[bothToOnePoints][0]) + ", ");
                record_file.write(str(bothToOne[bothToOnePoints][1]) + "}, ");
                record_file.write("ErrorBar[" + str(bothToOne[bothToOnePoints][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(bothToOne[bothToOneDimensions[0]-1][0]) + ", ");
            record_file.write(str(bothToOne[bothToOneDimensions[0]-1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(bothToOne[bothToOneDimensions[0]-1][2]) + "]");
            record_file.write(str("}"));
            record_file.write("\n")
            #ctypes.windll.user32.MessageBoxW(0, "Made it.9", "", 1)
            # sensible version
            # for dataInc in range(1,len(collatedData)):
                #record_file.write(str(collatedData[dataInc][0:len(collatedData[dataInc])]) + "\n");
    show();
    return ("Threshold Found: " + str(thresholds) + "\r\nThreshold Fidelity: " + str(thresholdFidelity[0] * 100)  + ", " + str(thresholdFidelity[1] * 100)  + "%")

# calling imports before calling this function
def singlePointAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    from numpy import array
    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\My Data Analysis")
    from astropy.io import fits
    import numpy
    import ctypes
    numpy.set_printoptions(threshold=numpy.nan)
    from matplotlib.pyplot import  subplots, show, gcf
    from matplotlib.cm import get_cmap
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian, 
                                       calculateAtomThreshold, getAnalyzedSurvivalData);
    #
    #dataRepositoryPath = "C:\\Users\\Mark\\Documents\\Quantum Gas Assembly Control\\Data\\Camera Data\\"
    dataRepositoryPath = "\\\\andor\\share\\Data and documents\\Data repository\\";
    todaysDataPath = dataRepositoryPath + date + "\\Raw Data\\data_" + str(runNumber) + ".fits";
    keyPath = dataRepositoryPath + date + "\\Raw Data\\key_" + str(runNumber) + ".txt";
    # Load Key
    key = numpy.array([]);
    with open(keyPath) as keyFile:
        for line in keyFile:
            key = numpy.append(key, float(line.strip('\n')))
    # Load Fits File & Get Dimensions
    # Get the array from the fits file. That's all I care about.
    fitsInfo = (fits.open(todaysDataPath, "append"));
    rawData = fitsInfo[0].data;
    # the .shape member of an array gives an array of the dimesnions of the array.
    numberOfPictures = rawData.shape[0];
    numberOfExperiments = int(numberOfPictures / picturesPerExperiment);
    # get accumulation image
    accumulationImage = numpy.zeros((rawData.shape[1], rawData.shape[2]));
    for imageInc in range(0, int(rawData.shape[0])):
        accumulationImage += rawData[imageInc];
    # Initial Data Analysis
    #
    ###########################################################################
    #
    #       Loop for each atom to analyze
    #
    threshold = 0;
    thresholdFidelity = 0;
    numberAtomsToAnalyze = (array(analysisLocations).shape)[0];       
    for atomInc in range(0, int(numberAtomsToAnalyze/2)):
        atomLocation = array([analysisLocations[2 * atomInc], analysisLocations[2 * atomInc + 1]]);
        peakData = [];
        firstExperimentData = [];    
        # my function here.
        peakData, firstExperimentData = normalizeData(picturesPerExperiment, rawData, atomLocation);
        #normalizeData(picturesPerExperiment, rawData, atomLocation);
        # ### Histogram
        # Plot histogram 
    
    
        #return "stuff";
        # Get Binned Data
        binCenters, binnedData = binData(5, peakData);
        # Make educated Guesses for Peaks
        guess1, guess2 = guessGaussianPeaks(peakData, binCenters, binnedData);
    
        # Calculate Atom Threshold
        #define the fitting function
        guess = numpy.array([100, guess1, 30, 200, guess2, 10]);
        gaussianFitVals = fitDoubleGaussian(binCenters, binnedData, guess);
        threshold, thresholdFidelity = calculateAtomThreshold(gaussianFitVals);
        #     
        atomCount = 0;
        for experimentInc in range(0, firstExperimentData.size):
            if firstExperimentData[experimentInc] > threshold:
                atomCount += 1;
        print("Average Loading Efficiency: " + str(atomCount / firstExperimentData.size * 100) + "%")
        
        # ### Coalate Data        
        # Get Data in final form for exporting
        survivalData, fullCaptureData, captureArray = getAnalyzedSurvivalData(peakData, threshold, key, accumulations, 
                                                                      numberOfExperiments);
        collatedData = [accumulations, peakData, key, survivalData, captureArray, rawData, fullCaptureData];
        # Survival Data
        # Plot survival data.
        ###########
        # 
        # Plot Stuff
        #
        myFigure, ((plot11, plot12, plot13), (plot21, plot22, plot23)) = subplots(2, 3, figsize = (25,12))
        myFigure.suptitle(fileName + "_run" + str(runNumber) + "Data for location {" + str(atomLocation[0] + 1) + "," + str(atomLocation[1] + 1) + "}", fontsize = 24);        
        figObject = gcf()
        figObject.canvas.set_window_title("{" + str(atomLocation[0] + 1) + "," + str(atomLocation[1] + 1) + "}")
        # make an image
        plot11.imshow(accumulationImage, interpolation='none', cmap = get_cmap("bone"));
        plot11.set_title("Entire Run Accumulation Image")
        # First counts histogram
        plot12.hist(firstExperimentData, 50);
        plot12.set_title("First Picture of Experiment Pixel Count Histogram");
        plot12.set_ylabel("Occurance Count");
        plot12.set_xlabel("Pixel Counts");
        # plot loading probabilities
        plot13.plot(fullCaptureData[:,0], fullCaptureData[:,1], "o");
        plot13.set_title("Average Loading Efficiency: " + str(atomCount / firstExperimentData.size * 100) + "%");

        plot13.set_ylabel("Occurance Count");
        plot13.set_xlabel("Pixel Counts");
        plot13.set_ylim([-0.05, 1.05]);
        xRange = max(fullCaptureData[:,0]) - min(fullCaptureData[:,0]);
        plot13.set_xlim([min(fullCaptureData[:,0]) - xRange / 10.0, max(fullCaptureData[:,0]) + xRange / 10.0]);
        # plot the fit on top of the histogram
        plot21.plot(binCenters, binnedData, "o", markersize = 3);
        fineXData = numpy.linspace(min(peakData),max(peakData),500);
        plot21.plot(fineXData, doubleGaussian(fineXData,*gaussianFitVals))
        plot21.axvline(x=threshold, color='r', ls='dashed');
        plot21.set_title("Fits. Threshold = " + str(threshold));
        plot21.set_ylabel("Occurance Count");
        plot21.set_xlabel("Pixel Counts");
        plot22.plot(peakData, ".", markersize = 1);
        plot22.set_title("Pixel Count Data Over Time");
        plot22.set_ylabel("Count on Pixel");
        plot22.set_xlabel("Picture Number");
        plot23.errorbar(survivalData[:,0], survivalData[:,1], 
                                  yerr=survivalData[:,2], linestyle = "none");
        plot23.set_title("Survival Data for Each Variation")
        plot23.set_ylabel("% Survived")
        plot23.set_xlabel("Variation Parameter Value")
        plot23.set_ylim([-0.05, 1.05])
        xRange = max(survivalData[:,0]) - min(survivalData[:,0]);
        plot23.set_xlim([min(survivalData[:,0]) - xRange / 10.0, max(survivalData[:,0]) + xRange / 10.0]);
        # Export Data And Close
        # 
        fitsInfo.close()
        #collatedData = [accumulations, peakData, key, survivalData, captureProbabilities, rawData, captureProbabilities];
        outputName = dataRepositoryPath + date + "\\" + fileName + "_run" + str(runNumber) + "_well" + str(atomLocation[0] + 1) + str(atomLocation[1] + 1)+ ".tsv";
        print("data outputted to file " + outputName)
        with open(outputName, "w") as record_file:
            # accumulations is special since it's an integer.
            record_file.write(str(collatedData[0]) + "\n");
            ### Version for mathematica compatibility
            # Peak Data
            for peakInc in range(0, peakData.size):
                record_file.write(str(peakData[peakInc]) + "\t");
            record_file.write("\n");
            # key
            for keyInc in range(0, key.size-1):
                record_file.write(str(key[keyInc]) + "\t");
            record_file.write(str(key[key.size-1]));
            record_file.write("\n");
            # survival data
            survivalDimensions = survivalData.shape;
            for survivalPointsInc in range(0, survivalDimensions[0]-1):
                record_file.write(str("{{"));
                record_file.write(str(survivalData[survivalPointsInc][0]) + ", ");
                record_file.write(str(survivalData[survivalPointsInc][1]) + "}, ");
                record_file.write("ErrorBar[" + str(survivalData[survivalPointsInc][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(survivalData[survivalDimensions[0]-1][0]) + ", ");
            record_file.write(str(survivalData[survivalDimensions[0]-1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(survivalData[survivalDimensions[0]-1][2]) + "]");
            record_file.write(str("}"));
            record_file.write("\n")
                # capture probabilities data
            for captureInc in range(0, captureArray.size):
                record_file.write(str(captureArray[captureInc]) + " ");
            record_file.write("\n");
            # raw data
            rawDataDimensions = rawData.shape;
            for pictureInc in range(0, rawDataDimensions[0]):
                record_file.write("{");
                for rowInc in range(0, rawDataDimensions[1]):
                    record_file.write("{");
                    for columnInc in range(0, rawDataDimensions[2]):
                        record_file.write(str(rawData[pictureInc][rowInc][columnInc]) + " ")
                    record_file.write("} ")
                record_file.write("}");
            record_file.write("\n");
            # full Capture Probabilitiy Data (capture probabilities with x values?)
            fullCaptureDataDimensions = fullCaptureData.shape;
            for variationInc in range(0, fullCaptureDataDimensions[0]):
                record_file.write("{" + str(fullCaptureData[variationInc][0]) + ", "
                                  + str(fullCaptureData[variationInc][1]) + "} ");
            # sensible version
            # for dataInc in range(1,len(collatedData)):
                #record_file.write(str(collatedData[dataInc][0:len(collatedData[dataInc])]) + "\n");
    show();
    return ("Threshold Found: " + str(threshold) + "\r\nThreshold Fidelity: " + str(thresholdFidelity * 100) + "%")

 # In[]:
#date = "160521";
#fileName = "CarrierCalibration";
#runNumber = 55;
#wellIndicator = 4;
#accumulations = 150;
#### Zero-indexed!!!
#
## Vertical
#atomLocation = [1,wellIndicator-1];
## Horizontal
## atomLocation = [wellIndicator-1, 1];
#
#picturesPerExperiment = 2;
#                 
#singlePointAnalysis(date, runNumber, 1, 3, picturesPerExperiment, accumulations, "stuff")


