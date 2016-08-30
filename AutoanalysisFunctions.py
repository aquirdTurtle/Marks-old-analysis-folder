# AutoanalysisFunctions.py
# This File contains all of the main functions used by my c++ code that controls the experiment.


def pairAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
    import matplotlib as mpl
    mpl.rcParams['text.color'] = '#ffffff'
    mpl.rcParams['figure.edgecolor'] = '#ffffff'
    mpl.rcParams['xtick.color'] = '#ffffff'
    mpl.rcParams['ytick.color'] = '#ffffff'
    mpl.rcParams['figure.facecolor'] = '#000000'
    mpl.rcParams['axes.facecolor'] = '#0a0a0a'
    mpl.rcParams['figure.figsize'] = (18.0, 8.0)
    mpl.rcParams['axes.labelcolor'] = '#ffffff'
    mpl.rcParams['grid.color'] = '#aaaaff'
    mpl.rcParams['axes.edgecolor'] = '#ffffff'
    mpl.rcParams['legend.facecolor'] = '#00001f'
    mpl.rcParams['axes.grid'] = True
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
        loadingPlot.set_title("Average Loading Efficiencies: "
                              + "{:.1f}".format(atomCount1 / firstExperimentData[0].size * 100)
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
        twoParticleDataPlot.errorbar(key, averageBothToBoth, yerr=error_BothToBoth, linestyle="none", color="c",
                                     fmt="o", label="both-to-both")
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
                      + "," + str(location2[1] + 1) + ")" + ".tsv" + "\"Data for locations {" + str(location1[0] + 1)
                          + "," + str(location1[1] + 1) + "} and {" + str(location2[0] + 1) + ","
                          + str(location2[1] + 1) + "}", fontsize=24)
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
    # return not null...
    return "Finished"


def singlePointAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, repetitions, fileName):
    import matplotlib as mpl
    mpl.rcParams['text.color'] = '#ffffff'
    mpl.rcParams['figure.edgecolor'] = '#ffffff'
    mpl.rcParams['xtick.color'] = '#ffffff'
    mpl.rcParams['ytick.color'] = '#ffffff'
    mpl.rcParams['figure.facecolor'] = '#000000'
    mpl.rcParams['axes.facecolor'] = '#0a0a0a'
    mpl.rcParams['figure.figsize'] = (18.0, 8.0)
    mpl.rcParams['axes.labelcolor'] = '#ffffff'
    mpl.rcParams['grid.color'] = '#aaaaff'
    mpl.rcParams['axes.edgecolor'] = '#ffffff'
    mpl.rcParams['legend.facecolor'] = '#00001f'
    mpl.rcParams['axes.grid'] = True
    from numpy import array
    import sys
    sys.path.append("C:\\Users\\Mark\\Documents\\Data-Analysis")
    from astropy.io import fits
    import numpy as np
    np.set_printoptions(threshold=np.nan)
    import matplotlib.pyplot as plt
    from matplotlib.cm import get_cmap
    from collections import OrderedDict as dic
    from dataAnalysisFunctions import (normalizeData, binData, guessGaussianPeaks, doubleGaussian, fitDoubleGaussian,
                                       calculateAtomThreshold, getCorrelationData, getAtomData,
                                       getSingleParticleSurvivalData)
    baseData = dic()
    baseData['Key List'] = ''
    baseData['Date'] = date
    baseData['Run Number'] = runNumber
    baseData['Repetitions'] = repetitions
    baseData['Pictures Per Experiment'] = picturesPerExperiment

    #dataRepositoryPath = "C:\\Users\\Mark\\Documents\\Quantum Gas Assembly Control\\Data\\Camera Data\\"
    dataRepositoryPath = "\\\\andor\\share\\Data and documents\\Data repository\\"
    todaysDataPath = dataRepositoryPath + date + "\\Raw Data\\data_" + str(baseData['Run Number']) + ".fits"
    keyPath = dataRepositoryPath + date + "\\Raw Data\\key_" + str(baseData['Run Number']) + ".txt"
    # Load Key
    baseData['Key'] = np.array([])
    with open(keyPath) as keyFile:
        for line in keyFile:
            baseData['Key'] = np.append(baseData['Key'], float(line.strip('\n')))
    # Load Fits File & Get Dimensions
    # Get the array from the fits file. That's all I care about.
    fitsInfo = fits.open(todaysDataPath, "append")
    #baseData['Raw Data'] = fitsInfo[0].data
    rawData = fitsInfo[0].data
    numberOfPictures = rawData.shape[0]
    numberOfExperiments = int(numberOfPictures / picturesPerExperiment)
    # get accumulation image
    accumulationImage = np.zeros((rawData.shape[1], rawData.shape[2]))
    for imageInc in range(0, int(rawData.shape[0])):
        accumulationImage += rawData[imageInc]
    fitsInfo.close()

    ###########################################################################
    #
    #       Loop for each atom to analyze
    #
    numberAtomsToAnalyze = np.array(analysisLocations).shape[0]
    # array of lists.
    allAtomSurvivalData = np.array([[]])
    for atomInc in range(0, int(numberAtomsToAnalyze / 2)):
        print('Analyzing atom #' + str(atomInc))
        tempData = dic()
        tempData['Key List'] = ''
        tempData['Atom Location'] = array([analysisLocations[2 * atomInc], analysisLocations[2 * atomInc + 1]])
        # my function here.
        tempData['Data Counts'], firstExperimentData = normalizeData(picturesPerExperiment, rawData,
                                                                 tempData['Atom Location'])
        # normalizeData(picturesPerExperiment, rawData, atomLocation);
        # ### Histogram
        # Plot histogram 
        # Get Binned Data
        binCenters, binnedData = binData(5, tempData['Data Counts'])
        # Make educated Guesses for Peaks
        guess1, guess2 = guessGaussianPeaks(tempData['Data Counts'], binCenters, binnedData)

        # Calculate Atom Threshold
        # define the fitting function
        guess = np.array([100, guess1, 30, 200, guess2, 10])
        gaussianFitVals = fitDoubleGaussian(binCenters, binnedData, guess)
        tempData['Threshold'], tempData['Threshold Fidelity'] = calculateAtomThreshold(gaussianFitVals)
        #
        atomCount = 0
        for experimentInc in range(0, firstExperimentData.size):
            if firstExperimentData[experimentInc] > tempData['Threshold']:
                atomCount += 1

        # Get Data in final form for exporting
        tempData['Atoms Data'] = getAtomData(tempData['Data Counts'], tempData['Threshold'], numberOfExperiments)
        allAtomSurvivalData.resize((allAtomSurvivalData.shape[0], len(tempData['Atoms Data'])))
        if atomInc == 0:
            allAtomSurvivalData[0] = tempData['Atoms Data']
        else:
            allAtomSurvivalData = np.append(allAtomSurvivalData, np.array([tempData['Atoms Data']]), 0)

        tempData['Survival Averages'], tempData['Survival Errors'], tempData['Loading Probabilities'] \
            = getSingleParticleSurvivalData(tempData['Atoms Data'], baseData['Repetitions'])

        ###########
        # 
        # Plot Data
        #
        # Carrier Plot
        myFig = plt.figure(atomInc)
        mng = plt.get_current_fig_manager()
        mng.window.showMaximized()
        mainPlot = plt.subplot2grid((4, 4), (0, 0), colspan=3, rowspan=4)
        mainPlot.errorbar(baseData["Key"], tempData["Survival Averages"],
                          yerr=tempData["Survival Errors"], ls='', marker='o', label="well 6", color='b',
                          capsize=6, elinewidth=3)
        mainPlot.set_ylim({-0.02, 1.01})
        mainPlot.set_title(str(tempData['Atom Location']) + " Survival Probability Throughout Experiment", fontsize=30)
        mainPlot.set_ylabel("Survival Probability", fontsize=20)
        mainPlot.set_xlabel("Key Value", fontsize=20)
        #mainPlot.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), fancybox=True, ncol=4)
        mainPlot.grid("on")
        # Capture Probabilities Plot
        capturePlot = plt.subplot2grid((4, 4), (0, 3))
        capturePlot.plot(baseData['Key'], tempData["Loading Probabilities"], ls='', marker='o', color='b')
        capturePlot.set_ylim({0, 1})
        capturePlot.set_xlabel("Key Value")
        capturePlot.set_ylabel("Loading %")
        capturePlot.set_title("Loading Probabilities")
        capturePlot.grid("on")
        # Count Series Plot
        countDataPlot = plt.subplot2grid((4, 4), (1, 3))
        countDataPlot.plot(tempData["Data Counts"], 'b', ls='', marker='.', markersize=1)
        countDataPlot.set_xlabel("Picture #")
        countDataPlot.set_ylabel("Camera Signal")
        countDataPlot.set_title("Camera Signal Over Time")
        countDataPlot.grid("on")
        countDataPlot.axhline(tempData["Threshold"], color='b')
        # Accumulation Image
        accumulationImagePlot = plt.subplot2grid((4, 4), (2, 3))
        zeroedImage = np.array(accumulationImage) - np.amin(np.array(accumulationImage))
        # make it slightly darker to make the red stand out.
        normalizedImage = np.array(zeroedImage) / (1.5 * np.amax(np.array(zeroedImage)))
        coloredImage = get_cmap("bone")(normalizedImage)
        # make the relevant pixel slightly redder.
        coloredImage[tempData['Atom Location'][0]][tempData['Atom Location'][1]][0] \
            += (1 - coloredImage[tempData['Atom Location'][0]][tempData['Atom Location'][1]][0]) / 2
        accumulationImagePlot.imshow(coloredImage, interpolation='none')
        accumulationImagePlot.set_title("Entire Run Accumulation Image")
        # info plot
        infoPlot = plt.subplot2grid((4, 4), (3, 3))
        infoPlot.axis("off")
        infoPlot.text(0, 0.0, "Repetitions: " + str(baseData["Repetitions"]))
        infoPlot.text(0, 0.2, "Analysis Location: " + str(tempData["Atom Location"]))
        infoPlot.text(0, 1.0, "Date: " + str(baseData["Date"]))
        infoPlot.text(0, 0.8, "Run #: " + str(baseData["Run Number"]))
        infoPlot.text(0, 0.6, "Fit Threshold: " + str(tempData['Threshold']))
        infoPlot.text(0, 0.4, "Fit Threshold Fidelity: " + str(tempData['Threshold Fidelity']))
        plt.tight_layout()
        # add the data to the main data object.
        tempData['Key List'] = list(tempData.keys())
        baseData[str(analysisLocations[2 * atomInc]) + ", " + str(analysisLocations[2 * atomInc + 1])] = tempData

    print('Getting Correlation Data...')
    baseData['Correlation Averages'], baseData['Correlation Errors'] = getCorrelationData(allAtomSurvivalData, baseData['Repetitions'])
    # Plot correlation Data
    myFig = plt.figure(atomInc + 1)
    numberSurvivedPlot = plt.subplot2grid((4, 4), (0, 0), colspan=2, rowspan=2)
    for atomInc in range(1, int(numberAtomsToAnalyze / 2) + 1):
        name = 'Load ' + str(atomInc) + ', 1 atoms survived'
        numberSurvivedPlot.errorbar(baseData['Key'], baseData['Correlation Averages'][name],
                                    yerr=baseData['Correlation Errors'][name], label=str(atomInc) + '-1',
                                    linestyle='none', marker='o', markersize=1)
    numberSurvivedPlot.set_ylim({-0.05, 1.05})
    xrange = max(baseData['Key']) - min(baseData['Key'])
    pixelNum = allAtomSurvivalData.shape[0]
    numberSurvivedPlot.set_xlim([min(baseData['Key']) - xrange / baseData['Key'].size,
                                max(baseData['Key']) + xrange / baseData['Key'].size])
    numberSurvivedPlot.set_xlabel('Key Value')
    numberSurvivedPlot.set_ylabel('Survival %')
    numberSurvivedPlot.set_ylim({-0.05, 1.05})
    numberSurvivedPlot.grid(True)
    numberSurvivedPlot.legend()
    numberSurvivedPlot.set_title('Chance that 1 Survived Probabilities')

    averagePlot = plt.subplot2grid((4, 4), (0, 2), colspan=2, rowspan=2)
    for atomsLoadedInc in range(1, int(numberAtomsToAnalyze/2) + 1):
        name = 'Load ' + str(atomsLoadedInc) + ', average single atom survival'
        averagePlot.errorbar(baseData["Key"], baseData['Correlation Averages'][name],
                             yerr=baseData['Correlation Errors'][name], linestyle='none', marker='o',
                             label='Load ' + str(atomsLoadedInc))
    averagePlot.set_ylim({-0.05, 1.05})
    xrange = max(baseData['Key']) - min(baseData['Key'])
    averagePlot.set_xlim([min(baseData['Key']) - xrange / baseData['Key'].size,
                          max(baseData['Key']) + xrange / baseData['Key'].size])
    averagePlot.set_xlabel('Key Value')
    averagePlot.set_ylabel('Survival %')
    averagePlot.set_ylim({-0.05, 1.05})
    averagePlot.grid(True)
    averagePlot.legend()
    averagePlot.set_title('(Averaged over Wells) Single Atom Survival Given Differnet Loading')

    repNum = allAtomSurvivalData.shape[1]
    wellPlot = plt.subplot2grid((4, 4), (2, 0), colspan=3, rowspan=2)
    if int(repNum / baseData['Repetitions']) == 1:
        # plot each pixel as a number, don't look at key
        pixels = list(range(1, pixelNum + 1))
        fourToOneData = []
        fourToOneError = []
        oneToOneData = []
        oneToOneError = []
        for pixelInc in range(1, pixelNum + 1):
            name = 'Load ' + str(pixelNum) + ', atom ' + str(pixelInc-1) + ' survived'
            fourToOneData.append(baseData['Correlation Averages'][name][0])
            fourToOneError.append(baseData['Correlation Errors'][name][0])
            name = 'Load 1, atom ' + str(pixelInc-1) + ' survived'
            oneToOneData.append(baseData['Correlation Averages'][name][0])
            oneToOneError.append(baseData['Correlation Errors'][name][0])
        wellPlot.errorbar(pixels, fourToOneData, yerr=fourToOneError, linestyle="none", marker='o', color='blue',
                          label=str(pixelNum) + '-1', markersize=1)
        wellPlot.errorbar(pixels, oneToOneData, yerr=oneToOneError, linestyle="none", marker='o', color='red',
                          label='1-1', markersize=1)
        wellPlot.set_xlim([0, pixelNum + 1])
        wellPlot.set_xlabel('Well #')
        wellPlot.set_ylabel('Survival %')
        wellPlot.set_ylim({-0.05, 1.05})
        wellPlot.grid(True)
        wellPlot.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        wellPlot.set_title('Individual Well Survival Given Different Loading')
    else:
        pixelNum = allAtomSurvivalData.shape[0]
        for pixelInc in range(1, pixelNum):
            name = 'Load ' + str(pixelNum) + ', atom ' + str(pixelInc) + ' survived'
            fourToOneData = baseData['Correlation Averages'][name]
            fourToOneError = baseData['Correlation Errors'][name]
            name = 'Load 1, atom ' + str(pixelInc) + ' survived'
            oneToOneData = baseData['Correlation Averages'][name]
            oneToOneError = baseData['Correlation Errors'][name]
            wellPlot.errorbar(baseData['Key'], fourToOneData, yerr=fourToOneError, linestyle="none", marker='o',
                              label="Load " + str(pixelNum) + ', atom ' + str(pixelInc) + " survives")
            wellPlot.errorbar(baseData['Key'], oneToOneData, yerr=oneToOneError, linestyle="none", marker='o',
                              label="Load Only 1, atom " + str(pixelInc) + " survives")
        xrange = max(baseData['Key']) - min(baseData['Key'])
        if xrange == 0:
            xrange = 1

        wellPlot.set_xlim([min(baseData['Key']) - xrange / baseData['Key'].size,
                           max(baseData['Key']) + xrange/baseData['Key'].size, pixelNum + 1])
        wellPlot.set_xlabel('Key Value')
        wellPlot.set_ylabel('Survival %')
        wellPlot.set_ylim({-0.05, 1.05})
        wellPlot.grid(True)
        wellPlot.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        wellPlot.set_title('Individual Well Survival Given Different Loading')
    # info plot
    infoPlot = plt.subplot2grid((4, 4), (3, 3))
    infoPlot.axis("off")
    infoPlot.text(0, 1.0, "Date: " + str(baseData["Date"]))
    infoPlot.text(0, 0.8, "Run #: " + str(baseData["Run Number"]))
    infoPlot.text(0, 0.6, "Total Atom #: " + str(int(numberAtomsToAnalyze/2)))
    infoPlot.text(0, 0.4, "Repetitions: " + str(baseData["Repetitions"]))

    myFig.tight_layout()

    # ########################################
    #
    # Export Data
    #
    print('Prepping Data')
    outputName = dataRepositoryPath + baseData['Date'] + "\\" + fileName + "_run" + str(
        baseData['Run Number']) + ".csv"
    baseData['Key List'] = list(baseData.keys())
    csvText = ''
    for keyHeader, value in baseData.items():
        if isinstance(value, str):
            # don't iterate through the string, just add it.
            csvText += '\n:' + keyHeader + ': ' + str(value)
            continue
        if isinstance(value, dict):
            # iterate through that! Assume no nested dictionaries.
            csvText += '\n:[' + keyHeader + ']:'
            for subHeader, subValue in value.items():
                if subHeader == "Raw Data":
                    # want to put this on last.
                    continue
                if isinstance(subValue, str):
                    # don't iterate through the string, just add it.
                    csvText += '\n\t;' + subHeader + '; ' + str(subValue)
                    continue
                try:
                    csvText += '\n\t;' + subHeader + '; ' + ", ".join(str(x) for x in subValue)
                except TypeError:
                    # catch integers.
                    csvText += '\n\t;' + subHeader + '; ' + str(subValue)
            continue
        try:
            csvText += '\n:' + keyHeader + ': ' + ", ".join(str(x) for x in value)
        except TypeError:
            # catch integers.
            csvText += '\n:' + keyHeader + ': ' + str(value)
    print("Writing Data...")
    with open(outputName, "w") as record_file:
        record_file.write(csvText)
    print('Complete!')
    #plt.switch_backend('QT4Agg')  # Widgits
    mng = plt.get_current_fig_manager()
    mng.window.showMaximized()
    plt.show()
    return "Finished"

# date = "160521";
# fileName = "CarrierCalibration";
# runNumber = 55;
# wellIndicator = 4;
# accumulations = 150;
# ### Zero-indexed!!!
#
# # Vertical
# atomLocations = {1, 3, 1, 5};
# # Horizontal
# # atomLocation = [wellIndicator-1, 1];
#
# picturesPerExperiment = 2;
#
#singlePointAnalysis("160824", 21, [3, 1, 5, 1, 8, 1, 10, 1], 2, 10000, "testAnalysis")


# def pairAnalysis(date, runNumber, analysisLocations, picturesPerExperiment, accumulations, fileName):
# pairAnalysis("160805", 19, [3, 1, 5, 1], 2, 150, "test")
