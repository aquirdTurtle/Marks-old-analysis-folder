# AutoanalysisFunctions.py
# This File contains all of the main functions used by my c++ code that controls the experiment.


def pairAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    from numpy import array
    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\My Data Analysis")
    from astropy.io import fits
    import numpy
    numpy.set_printoptions(threshold=numpy.nan)
    import matplotlib.pyplot as plt
    from matplotlib.pyplot import subplots, show, gcf
    from matplotlib.cm import get_cmap
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian,
                                       calculateAtomThreshold, getAnalyzedTunnelingData, getPostSelectedTunnelingData);
    from matplotlib.font_manager import FontProperties
    # paths for files
    dataRepositoryPath = "\\\\andor\\share\\Data and documents\\Data repository\\";
    todaysDataPath = dataRepositoryPath + date + "\\Raw Data\\data_" + str(runNumber) + ".fits";
    keyPath = dataRepositoryPath + date + "\\Raw Data\\key_" + str(runNumber) + ".txt";
    # Load Key
    key = numpy.array([])
    with open(keyPath) as keyFile:
        for line in keyFile:
            key = numpy.append(key, float(line.strip('\n')))
            # Load Fits File & Get Dimensions
    # Get the array from the fits file. That's all I care about.
    fitsInfo = (fits.open(todaysDataPath, "append"))
    rawData = fitsInfo[0].data
    accumulationImage = numpy.zeros((rawData.shape[1], rawData.shape[2]))
    for imageInc in range(0, int(rawData.shape[0])):
        accumulationImage += rawData[imageInc]
    # the .shape member of an array gives an array of the dimesnions of the array.
    numberOfPictures = rawData.shape[0]
    numberOfExperiments = int(numberOfPictures / picturesPerExperiment)
    #
    numberAtomsToAnalyze = array(analysisLocations).size

    for atomInc in range(0, int(numberAtomsToAnalyze / 4)):
        location1 = array([analysisLocations[4 * atomInc], analysisLocations[4 * atomInc + 1]])
        location2 = array([analysisLocations[4 * atomInc + 2], analysisLocations[4 * atomInc + 3]])
        allAtomData = [[], []]
        firstExperimentData = [[], []]
        allAtomData[0], firstExperimentData[0] = normalizeData(picturesPerExperiment, rawData, location1)
        allAtomData[1], firstExperimentData[1] = normalizeData(picturesPerExperiment, rawData, location2)
        binCenters = [[], []]
        binnedData = [[], []]
        binCenters[0], binnedData[0] = binData(5, allAtomData[0])
        binCenters[1], binnedData[1] = binData(5, allAtomData[1])
        guessLocation1 = [[], []]
        guessLocation2 = [[], []]
        guessLocation1[0], guessLocation2[0] = guessGaussianPeaks(allAtomData[0], binCenters[0], binnedData[0])
        guessLocation1[1], guessLocation2[1] = guessGaussianPeaks(allAtomData[0], binCenters[1], binnedData[1])
        guess = [[], []]
        gaussianFitVals = [[], []]
        thresholds = [[], []]
        thresholdFidelity = [[], []]
        guess[0] = numpy.array([100, guessLocation1[0], 30, 200, guessLocation2[0], 10])
        guess[1] = numpy.array([100, guessLocation1[1], 30, 200, guessLocation2[1], 10])
        gaussianFitVals[0] = fitDoubleGaussian(binCenters[0], binnedData[0], guess[0])
        gaussianFitVals[1] = fitDoubleGaussian(binCenters[1], binnedData[1], guess[1])
        thresholds[0], thresholdFidelity[0] = calculateAtomThreshold(gaussianFitVals[0])
        thresholds[1], thresholdFidelity[1] = calculateAtomThreshold(gaussianFitVals[1])
        #
        atomCount1 = 0
        atomCount2 = 0
        for experimentInc in range(0, firstExperimentData[0].size):
            if firstExperimentData[0][experimentInc] > thresholds[0]:
                atomCount1 += 1
        for experimentInc in range(0, firstExperimentData[0].size):
            if firstExperimentData[1][experimentInc] > thresholds[1]:
                atomCount2 += 1
        (average_1to1, error_1to1, average_1to2, error_1to2, average_2to1, error_2to1, average_2to2, error_2to2,
        averageBothToBoth, error_BothToBoth, averageBothToOne, error_BothToOne, captProbs1, captProbs2,
        survival_1, error_survival_1, survival_2, error_survival_2) \
            = getAnalyzedTunnelingData(allAtomData, thresholds, key, accumulations, numberOfExperiments)

        # ps stands for post-selected.
        (ps_average_1to1, ps_error_1to1, ps_average_1to2, ps_error_1to2, ps_average_2to1, ps_error_2to1, ps_average_2to2,
         ps_error_2to2, ps_averageBothToBoth, ps_error_BothToBoth, ps_averageBothToOne, ps_error_BothToOne,
         captProbs1, captProbs2, no, no, no, no) \
            = getPostSelectedTunnelingData(allAtomData, thresholds, key, accumulations, numberOfExperiments)

        myFigure = plt.figure(1, facecolor="white", figsize=(25, 12))
        noTransferPlot = plt.subplot2grid((3, 3), (0, 1), colspan=2, rowspan=1)
        transferPlot = plt.subplot2grid((3, 3), (1, 1), colspan=2, rowspan=1)
        twoParticleDataPlot = plt.subplot2grid((3, 3), (2, 1), colspan=2, rowspan=1)
        accumulationPlot = plt.subplot2grid((3, 3), (0, 0))
        loadingPlot = plt.subplot2grid((3, 3), (1, 0))
        signalPlot = plt.subplot2grid((3, 3), (2, 0))


        figObject = gcf()
        figObject.canvas.set_window_title("{" + str(location1[0] + 1) + ","
                                          + str(location1[1] + 1) + "}, {" + str(location2[0] + 1)
                                          + "," + str(location2[1] + 1) + "}")
        # make an image
        #print(accumulationImage)
        accumulationPlot.imshow(accumulationImage, interpolation='none', cmap=get_cmap("bone"))
        accumulationPlot.set_title("All Pictures Accumulation Image")

        # all atom data over time
        signalPlot.plot(allAtomData[0], ".", markersize=1, color='r')
        signalPlot.plot(allAtomData[1], ".", markersize=1, color='b')
        signalPlot.set_title("Camera Signal Data Over Time")
        signalPlot.set_ylabel("Count on Pixel")
        signalPlot.set_xlabel("Picture Number")
        signalPlot.axhline(thresholds[0], color='r')
        signalPlot.axhline(thresholds[1], color='b')

        # plot loading probabilities
        loadingPlot.plot(key, captProbs1, linestyle="none", marker="o", color='r', label="atom 1")
        loadingPlot.plot(key, captProbs2, linestyle="none", marker="o", color='b', label="atom 2")
        loadingPlot.set_title("Average Loading Efficiencies: " + "{:.1f}".format(atomCount1 / firstExperimentData[0].size * 100)
                         + "% and " + "{:.1f}".format(atomCount2 / firstExperimentData[1].size * 100) + "%")
        loadingPlot.set_ylabel("Capture Probabilies")
        loadingPlot.set_xlabel("Experiment Parameter")
        loadingPlot.grid(True)
        loadingPlot.set_ylim([-0.05, 1.05])
        loadingPlot.axhline(numpy.average(captProbs1), color='r')
        loadingPlot.axhline(numpy.average(captProbs2), color='b')

        xRange = max(key) - min(key)
        # No Transfer probability.
        noTransferPlot.errorbar(key, ps_average_1to1, yerr=ps_error_1to1, linestyle="none", color="r",
                                fmt="o", label="1to1")
        noTransferPlot.errorbar(key, ps_average_2to2, yerr=ps_error_2to2, linestyle="none", color="b",
                                fmt="o", label="2to2")
        noTransferPlot.errorbar(key, survival_1, yerr=error_survival_1, linestyle="none", color="g", fmt="o",
                                label="1-survival")
        noTransferPlot.errorbar(key, survival_2, yerr=error_survival_2, linestyle="none", color="k", fmt="o",
                                label="2-survival")
        noTransferPlot.set_title("Non-Transfer Probability")
        noTransferPlot.set_ylabel("% Occurred")
        noTransferPlot.set_xlabel("Variation Parameter Value")
        noTransferPlot.set_ylim([-0.05, 1.05])
        noTransferPlot.set_xlim([min(key) - xRange / 10.0, max(key) + xRange / 10.0])
        font = FontProperties()
        font.set_size('small')
        noTransferPlot.grid(True)
        noTransferPlot.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., prop=font)

        # transfer probability
        transferPlot.errorbar(key, ps_average_1to2, yerr=ps_error_1to2, linestyle="none", color="r",
                              fmt="o", label="1to2")
        transferPlot.errorbar(key, ps_average_2to1, yerr=ps_error_2to1, linestyle="none", color="b",
                              fmt="o", label="2to1")
        transferPlot.errorbar(key, survival_1, yerr=error_survival_1, linestyle="none", color="g", fmt="o",
                              label="1-survival")
        transferPlot.errorbar(key, survival_2, yerr=error_survival_2, linestyle="none", color="k", fmt="o",
                              label="2-survival")
        transferPlot.set_title("Transfer Probability")
        transferPlot.set_ylabel("% Occurred")
        transferPlot.set_xlabel("Variation Parameter Value")
        transferPlot.set_ylim([-0.05, 1.05])
        transferPlot.grid(True)
        transferPlot.set_xlim([min(key) - xRange / 10.0, max(key) + xRange / 10.0])
        font = FontProperties()
        font.set_size('small')
        transferPlot.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., prop=font)

        # two-particle data.
        twoParticleDataPlot.errorbar(key, averageBothToBoth, yerr=error_BothToBoth, linestyle="none", color="c", fmt="o",
                                     label="both-to-both")
        twoParticleDataPlot.errorbar(key, averageBothToOne, yerr=error_BothToOne, linestyle="none", color="m", fmt="o",
                                     label="both-to-one")
        twoParticleDataPlot.errorbar(key, survival_1, yerr=error_survival_1, linestyle="none", color="g", fmt="o",
                                     label="1-survival")
        twoParticleDataPlot.errorbar(key, survival_2, yerr=error_survival_2, linestyle="none", color="k", fmt="o",
                                     label="2-survival")
        twoParticleDataPlot.set_title("Two-Particle Data")
        twoParticleDataPlot.set_ylabel("% Occurred")
        twoParticleDataPlot.set_xlabel("Variation Parameter Value")
        twoParticleDataPlot.set_ylim([-0.05, 1.05])
        twoParticleDataPlot.grid(True)
        xRange = max(key) - min(key)
        twoParticleDataPlot.set_xlim([min(key) - xRange / 10.0, max(key) + xRange / 10.0])
        font = FontProperties()
        font.set_size('small')
        twoParticleDataPlot.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0., prop=font)

        fitsInfo.close()
        # Export data
        outputName = (dataRepositoryPath + date + "\\" + fileName + "_run" + str(runNumber) + "_pair"
                      + "(" + str(location1[0] + 1) + "," + str(location1[1] + 1) + "), (" + str(location1[0] + 1)
                      + "," + str(location1[1] + 1) + ")" + ".tsv")
        myFigure.suptitle("\"" + fileName + "_run" + str(runNumber) + "_pair"
                      + "(" + str(location1[0] + 1) + "," + str(location1[1] + 1) + "), (" + str(location2[0] + 1)
                      + "," + str(location2[1] + 1) + ")" + ".tsv" + "\"Data for locations {" + str(location1[0] + 1) + ","
                          + str(location1[1] + 1) + "} and {" + str(location2[0] + 1)
                          + "," + str(location2[1] + 1) + "}", fontsize=24)
        with open(outputName, "w") as record_file:
            # ## 1 ###
            # accumulations is special since it's an integer.
            record_file.write(str(accumulations) + "\n")
            # ## 2 ###
            # key
            for keyInc in range(0, key.size):
                record_file.write(str(key[keyInc]) + " ")
            record_file.write("\n")
            # ## 3 ###
            # Run number
            record_file.write(str(runNumber) + "\n")
            # ## 4 ###
            # Date
            record_file.write(str(date) + "\n")
            # ## 5,6 ###
            # Peak Data
            for peakInc in range(0, allAtomData[0].size):
                record_file.write(str(allAtomData[0][peakInc]) + " ")
            record_file.write("\n")
            for peakInc in range(0, allAtomData[1].size):
                record_file.write(str(allAtomData[0][peakInc]) + " ")
            record_file.write("\n")

            # ## 7,8 ###
            # 1to1 data and error respectively.
            for inc1to1 in range(0, average_1to1.size):
                record_file.write(str(average_1to1[inc1to1]) + " ")
            record_file.write("\n")
            for inc1to1 in range(0, error_1to1.size):
                record_file.write(str(error_1to1[inc1to1]))
            record_file.write("\n")
            # ## 9, 10 ###
            # post selected 1to1 data an error respectively.
            for inc1to1 in range(0, ps_average_1to1.size):
                record_file.write(str(ps_average_1to1[inc1to1]) + " ")
            record_file.write("\n")
            for inc1to1 in range(0, ps_error_1to1.size):
                record_file.write(str(ps_error_1to1[inc1to1]))
            record_file.write("\n")
            # ## 11,12 ###
            # 1to2 data and error respectively.
            for inc1to2 in range(0, average_1to2.size):
                record_file.write(str(average_1to2[inc1to2]) + " ")
            record_file.write("\n")
            for inc1to2 in range(0, error_1to2.size):
                record_file.write(str(error_1to2[inc1to2]))
            record_file.write("\n")
            # ## 13, 14 ###
            # post selected 1to2 data an error respectively.
            for inc1to2 in range(0, ps_average_1to2.size):
                record_file.write(str(ps_average_1to2[inc1to2]) + " ")
            record_file.write("\n")
            for inc1to2 in range(0, ps_error_1to2.size):
                record_file.write(str(ps_error_1to2[inc1to2]))
            record_file.write("\n")
            # ## 14, 15 ###
            # 2to1 data and error respectively.
            for inc2to1 in range(0, average_2to1.size):
                record_file.write(str(average_2to1[inc2to1]) + " ")
            record_file.write("\n")
            for inc2to1 in range(0, error_2to1.size):
                record_file.write(str(error_2to1[inc2to1]))
            record_file.write("\n")
            # ## 16, 17 ###
            # post selected 2to1 data an error respectively.
            for inc2to1 in range(0, ps_average_2to1.size):
                record_file.write(str(ps_average_2to1[inc2to1]) + " ")
            record_file.write("\n")
            for inc2to1 in range(0, ps_error_2to1.size):
                record_file.write(str(ps_error_2to1[inc2to1]))
            record_file.write("\n")
            # ## 18, 19 ###
            # 2to2 data and error respectively.
            for inc2to2 in range(0, average_2to2.size):
                record_file.write(str(average_2to2[inc2to2]) + " ")
            record_file.write("\n")
            for inc2to2 in range(0, error_2to2.size):
                record_file.write(str(error_2to2[inc2to2]))
            record_file.write("\n")
            # ## 20, 21 ###
            # post selected 2to2 data an error respectively.
            for inc2to2 in range(0, ps_average_2to2.size):
                record_file.write(str(ps_average_2to2[inc2to2]) + " ")
            record_file.write("\n")
            for inc2to2 in range(0, ps_error_2to2.size):
                record_file.write(str(ps_error_2to2[inc2to2]))
            record_file.write("\n")

            # ## 22, 23 ###
            # bothToBoth data and error respectively
            for inc_BToB in range(0, averageBothToBoth.size):
                record_file.write(str(averageBothToBoth[inc_BToB]) + " ")
            record_file.write("\n")
            for inc_BToB in range(0, error_BothToBoth.size):
                record_file.write(str(error_BothToBoth[inc_BToB]))
            record_file.write("\n")

            # ## 24, 25 ###
            # bothToOne data and error respectively.
            for inc_BTo1 in range(0, averageBothToOne.size):
                record_file.write(str(averageBothToOne[inc_BTo1]) + " ")
            record_file.write("\n")
            for inc_BTo1 in range(0, error_BothToOne.size):
                record_file.write(str(error_BothToOne[inc_BTo1]))
            record_file.write("\n")

            # ## 26, 27 ###
            # capture probabilities data for wells 1 and 2 respectively.
            for captureInc in range(0, captProbs1.size):
                record_file.write(str(captProbs1[captureInc]) + " ")
            record_file.write("\n")
            for captureInc in range(0, captProbs2.size):
                record_file.write(str(captProbs2[captureInc]) + " ")
            record_file.write("\n")
            # ## 28 ###
            # raw data
            rawDataDimensions = rawData.shape
            for pictureInc in range(0, rawDataDimensions[0]):
                record_file.write("{")
                for rowInc in range(0, rawDataDimensions[1]):
                    record_file.write("{")
                    for columnInc in range(0, rawDataDimensions[2]):
                        record_file.write(str(rawData[pictureInc][rowInc][columnInc]) + " ")
                    record_file.write("} ")
                record_file.write("}")
            record_file.write("\n")
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, right=0.94, left=0.03)
    show()
    return ""


def singlePointAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    from numpy import array
    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\My Data Analysis")
    from astropy.io import fits
    import numpy
    import ctypes
    numpy.set_printoptions(threshold=numpy.nan)
    from matplotlib.pyplot import subplots, show, gcf
    from matplotlib.cm import get_cmap
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian,
                                       calculateAtomThreshold, getAnalyzedSurvivalData);
    #
    # dataRepositoryPath = "C:\\Users\\Mark\\Documents\\Quantum Gas Assembly Control\\Data\\Camera Data\\"
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
    for atomInc in range(0, int(numberAtomsToAnalyze / 2)):
        atomLocation = array([analysisLocations[2 * atomInc], analysisLocations[2 * atomInc + 1]]);
        peakData = [];
        firstExperimentData = [];
        # my function here.
        peakData, firstExperimentData = normalizeData(picturesPerExperiment, rawData, atomLocation);
        # normalizeData(picturesPerExperiment, rawData, atomLocation);
        # ### Histogram
        # Plot histogram 
        # Get Binned Data
        binCenters, binnedData = binData(5, peakData);
        # Make educated Guesses for Peaks
        guess1, guess2 = guessGaussianPeaks(peakData, binCenters, binnedData);

        # Calculate Atom Threshold
        # define the fitting function
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
        myFigure, ((plot11, plot12, plot13), (plot21, plot22, plot23)) = subplots(2, 3, figsize=(25, 12))
        myFigure.suptitle(
            "\"" + fileName + "_run" + str(runNumber) + "\" Data for location {" + str(atomLocation[0] + 1) + "," + str(
                atomLocation[1] + 1) + "}", fontsize=24);
        figObject = gcf()
        figObject.canvas.set_window_title("{" + str(atomLocation[0] + 1) + "," + str(atomLocation[1] + 1) + "}")
        # make an image
        plot11.imshow(accumulationImage, interpolation='none', cmap=get_cmap("bone"));
        plot11.set_title("Entire Run Accumulation Image")
        # First counts histogram
        plot12.hist(firstExperimentData, 50);
        plot12.set_title("First Picture of Experiment Pixel Count Histogram");
        plot12.set_ylabel("Occurance Count");
        plot12.set_xlabel("Pixel Counts");
        # plot loading probabilities
        plot13.plot(fullCaptureData[:, 0], fullCaptureData[:, 1], "o");
        plot13.set_title("Average Loading Efficiency: " + str(atomCount / firstExperimentData.size * 100) + "%");

        plot13.set_ylabel("Occurance Count");
        plot13.set_xlabel("Pixel Counts");
        plot13.set_ylim([-0.05, 1.05]);
        xRange = max(fullCaptureData[:, 0]) - min(fullCaptureData[:, 0]);
        plot13.set_xlim([min(fullCaptureData[:, 0]) - xRange / 10.0, max(fullCaptureData[:, 0]) + xRange / 10.0]);
        # plot the fit on top of the histogram
        plot21.plot(binCenters, binnedData, "o", markersize=3);
        fineXData = numpy.linspace(min(peakData), max(peakData), 500);
        plot21.plot(fineXData, doubleGaussian(fineXData, *gaussianFitVals))
        plot21.axvline(x=threshold, color='r', ls='dashed');
        plot21.set_title("Fits. Threshold = " + str(threshold));
        plot21.set_ylabel("Occurance Count");
        plot21.set_xlabel("Pixel Counts");
        plot22.plot(peakData, ".", markersize=1);
        plot22.set_title("Pixel Count Data Over Time");
        plot22.set_ylabel("Count on Pixel");
        plot22.set_xlabel("Picture Number");
        plot23.errorbar(survivalData[:, 0], survivalData[:, 1],
                        yerr=survivalData[:, 2], linestyle="none");
        plot23.set_title("Survival Data for Each Variation")
        plot23.set_ylabel("% Survived")
        plot23.set_xlabel("Variation Parameter Value")
        plot23.set_ylim([-0.05, 1.05])
        xRange = max(survivalData[:, 0]) - min(survivalData[:, 0]);
        plot23.set_xlim([min(survivalData[:, 0]) - xRange / 10.0, max(survivalData[:, 0]) + xRange / 10.0]);
        # Export Data And Close
        # 
        fitsInfo.close()
        # collatedData = [accumulations, peakData, key, survivalData, captureProbabilities, rawData, captureProbabilities];
        outputName = dataRepositoryPath + date + "\\" + fileName + "_run" + str(runNumber) + "_well" + str(
            atomLocation[0] + 1) + str(atomLocation[1] + 1) + ".tsv";
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
            for keyInc in range(0, key.size - 1):
                record_file.write(str(key[keyInc]) + "\t");
            record_file.write(str(key[key.size - 1]));
            record_file.write("\n");
            # survival data
            survivalDimensions = survivalData.shape;
            for survivalPointsInc in range(0, survivalDimensions[0] - 1):
                record_file.write(str("{{"));
                record_file.write(str(survivalData[survivalPointsInc][0]) + ", ");
                record_file.write(str(survivalData[survivalPointsInc][1]) + "}, ");
                record_file.write("ErrorBar[" + str(survivalData[survivalPointsInc][2]) + "]");
                record_file.write(str("}\t"));
            record_file.write(str("{{"));
            record_file.write(str(survivalData[survivalDimensions[0] - 1][0]) + ", ");
            record_file.write(str(survivalData[survivalDimensions[0] - 1][1]) + "}, ");
            record_file.write("ErrorBar[" + str(survivalData[survivalDimensions[0] - 1][2]) + "]");
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
                # record_file.write(str(collatedData[dataInc][0:len(collatedData[dataInc])]) + "\n");
    show()
    print("function.")
    return ("Threshold Found: " + str(threshold) + "\r\nThreshold Fidelity: " + str(thresholdFidelity * 100) + "%")

    # In[]:

# date = "160521";
# fileName = "CarrierCalibration";
# runNumber = 55;
# wellIndicator = 4;
# accumulations = 150;
#### Zero-indexed!!!
#
## Vertical
# atomLocations = {1, 3, 1, 5};
## Horizontal
## atomLocation = [wellIndicator-1, 1];
#
# picturesPerExperiment = 2;
#                 
# singlePointAnalysis(date, runNumber, 1, 3, picturesPerExperiment, accumulations, "stuff")

# def pairAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
pairAnalysis("160805", 19, [3, 1, 5, 1], 2, 150, "test")
