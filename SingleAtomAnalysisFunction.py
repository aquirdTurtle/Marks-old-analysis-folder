# # Simple Atom Analysis Notebook
# #### This notebook is designed to replicate the Mathematica "Single_Atom" Notebook. This notebook makes heavy use of the "dataAnalysisFunctions" module that I created.

# In[0]:

# In[1]:
# ### Imports
### Uninteresting stuff.
# Data is saved as "data.fits" files

# In[2]:
# ### Initialize Constants


# calling imports before calling this function
def singlePointAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    from numpy import array

    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\My Data Analysis")
    from astropy.io import fits
    import numpy
    numpy.set_printoptions(threshold=numpy.nan)
    from matplotlib.pyplot import  subplots, show, gcf
    from matplotlib.cm import get_cmap
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian, 
                                       calculateAtomThreshold, getAnalyzedSurvivalData);

    
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
    fitsInfo.info();
    rawData = fitsInfo[0].data;
    # the .shape member of an array gives an array of the dimesnions of the array.
    numberOfPictures = rawData.shape[0];
    numberOfExperiments = int(numberOfPictures / picturesPerExperiment)
    # Initial Data Analysis
    #
        #return str(array(analysisLocations).shape);
    ###########################################################################
    #
    #       Loop for each atom to analyze
    #
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
        myFigure.suptitle("Data for location {" + str(atomLocation[0] + 1) + "," + str(atomLocation[1] + 1) + "}", fontsize = 24);        
        figObject = gcf()
        figObject.canvas.set_window_title("{" + str(atomLocation[0] + 1) + "," + str(atomLocation[1] + 1) + "}")
        # make an image
        plot11.imshow(rawData[10], interpolation='none', cmap = get_cmap("bone"));
        plot11.set_title("Example Raw Image")
        # First counts histogram
        plot12.hist(firstExperimentData, 50);
        plot12.set_title("First Picture of Experiment Pixel Count Histogram");
        plot12.set_ylabel("Occurance Count");
        plot12.set_xlabel("Pixel Counts");
        # make a histogram
        plot13.plot(fullCaptureData[:,0], fullCaptureData[:,1], linestyle = "none");
        plot13.set_title("Average Loading Efficiency: " + str(atomCount / firstExperimentData.size * 100) + "%");
        plot13.set_ylabel("Occurance Count");
        plot13.set_xlabel("Pixel Counts");
        # plot the fit on top of the histogram
        plot21.plot(binCenters, binnedData, "o", markersize = 3);
        fineXData = numpy.linspace(min(peakData),max(peakData),500);
        plot21.plot(fineXData, doubleGaussian(fineXData,*gaussianFitVals))
        plot21.set_title("Gaussian Fits over Histogram");
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
        # Export Data And Close
        # 
        fitsInfo.close()
        #collatedData = [accumulations, peakData, key, survivalData, captureProbabilities, rawData, captureProbabilities];
        outputName = dataRepositoryPath + date + "\\" + fileName + "_run" + str(runNumber) + "_well" + str(atomLocation[0]) + str(atomLocation[1])+ ".tsv";
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
# In[];
def hiFunc():
    print("hi!");
    return "hello.";

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


