# -*- coding: utf-8 -*-
"""
Created on Sat May 21 10:45:58 2016

@author: Mark
"""


"""
This function analyzes raw data at the location "atom location" and returns a normalized version of
(1) every picture
(2) an array containing only the first pictures
"""
def normalizeData(picsPerExperiment, rawData, atomLocation):

    from numpy import append, array
    firstData = array([]); 
    allData = array([]);
    dimensions = rawData.shape;
    for imageInc in range(0, dimensions[0]):
        averageBackground = 1/4*(rawData[imageInc][0][0] 
                                 + rawData[imageInc][dimensions[1]-1][dimensions[2]-1] 
                                 + rawData[imageInc][0][dimensions[2]-1] 
                                 + rawData[imageInc][dimensions[1]-1][0])
        allData = append(allData, rawData[imageInc][atomLocation[0]][atomLocation[1]] - averageBackground);
        if imageInc % picsPerExperiment == 0:
            firstData = append(firstData, rawData[imageInc][atomLocation[0]][atomLocation[1]]
                                               - averageBackground);
    return allData, firstData

def binData(binWidth, data):
    from numpy import append, array, min, max, histogram
    binBorderLocation = min(data);
    binsBorders = array([]);
    # get bin borders
    while binBorderLocation < max(data):
        binsBorders = append(binsBorders, binBorderLocation);
        binBorderLocation = binBorderLocation + binWidth;
    # trash gets set but is unused.
    binnedData, trash = histogram(data, binsBorders);
    binCenters = binsBorders[0:binsBorders.size-1];
    return binCenters, binnedData;

""" 
This code is written to make educated guesses for where to set the initial guesses for the locations of the two 
gaussians in the fit below. It finds the largest-binned data point in the set and sets that as one guess. It then
finds the distance from this guess to the maximum pixel count in the set as well as the distance to the lowest pixel
count. It takes whichever distance is smaller, and uses that information to infer which gaussian (atom or no atom)
the guess should be for. It takes this smaller distance and takes double that distance away from the closest edge to
define the region for this gaussian. It then finds the maximum in the remaining region, and uses this as the second 
guess.
"""
def guessGaussianPeaks(rawData, binCenters, binnedData):
    from numpy import absolute, argmax, min, max
    import ctypes;
    leftBorder = min(rawData);
    rightBorder = max(rawData);
    dataRange = rightBorder - leftBorder;
    # get index corresponding to global max
    guess1Index = argmax(binnedData);
    # get location of global max
    guess1Location = binCenters[guess1Index];
    # find closest side
    distToLeft = absolute(guess1Location - leftBorder);
    distToRight = absolute(guess1Location - rightBorder);
    subBinnedData = {};
    if (distToLeft < distToRight ):
        # find index of dividing point:
        found = False;
        locInc = 0;
        boundary = 0;
        # set the boundary that the code will use to search for the second peak.
        if (2 * distToLeft < dataRange / 3) :
            boundary = leftBorder + dataRange / 3;
        else:
            boundary = leftBorder + 2 * distToLeft;      
        while (found == False):
            if (binCenters[locInc] < boundary):
                locInc += 1;
                if locInc == binCenters.size:
                    found = True;
                    locInc -= 1;
            else:
                found = True;
        subBinnedData = binnedData[locInc:binnedData.size];
        # get index corresponding to second local maxima
        guess2Index = argmax(subBinnedData) + binnedData.size - subBinnedData.size;
    elif (distToLeft >= distToRight ):
        found = False;
        locInc = 0;
        boundary = 0;
        # set the boundary that the code will use to search for the second peak.
        if (2 * distToRight < dataRange / 3) :
            boundary = rightBorder - dataRange / 3;
        else:
            boundary = rightBorder - 2 * distToRight;      
 
        while (found == False):
            if (binCenters[locInc] < (rightBorder - 2 * distToRight)):
                locInc += 1;
                if locInc == binCenters.size:
                    found = True;
                    locInc -= 1;
            else:
                found = True;
        subBinnedData = binnedData[0:locInc];
        # get index corresponding to second local maxima
        guess2Index = argmax(subBinnedData)    
    #ctypes.windll.user32.MessageBoxW(0, "Made it.", "", 1)
    # get location of second local maxima
    guess2Location = binCenters[guess2Index];    
    return guess1Location, guess2Location;

def doubleGaussian(data, A1, x1, sig1, A2, x2, sig2):
    from numpy import absolute,  exp, sqrt
    return ((absolute(A1) * exp(-((data-x1)/(sqrt(2)*sig1))**2)) 
            + absolute(A2) * exp(-((data-x2)/(sqrt(2)*sig2))**2))

def fitDoubleGaussian(binCenters, binnedData, fitGuess):
    from scipy.optimize import curve_fit
    fitVals, trash = curve_fit(doubleGaussian, binCenters, binnedData, fitGuess);
    return fitVals
    
def calculateAtomThreshold(fitVals):
    from numpy import (sqrt, abs)
    from scipy.special import erf
    TCalc = (fitVals[4] - fitVals[1])/(abs(fitVals[5]) + abs(fitVals[2]));
    threshold = fitVals[1] + TCalc * fitVals[2];
    fidelity = 1/2 * (1 + erf(abs(TCalc)/sqrt(2)))
    return threshold, fidelity
    
"""
This function assumes 2 pictures.
It returns
(1) Survival Data W/ Errors
(2) Full capture probabilty
(3) Capture Probability Array (for mathematica export format consistency only)
"""
def getAnalyzedSurvivalData(data, threshold, key, accumulations, numberOfExperiments):
    from numpy import (append, array, average, column_stack, std, sqrt)
    survivalRawData = array([]);
    survivalRawData.astype(int);
    numberTransferred = 0;
    # this doesn't take into account loss!
    for experimentInc in range(0, numberOfExperiments - 1):
        if data[2 * experimentInc] > threshold and data[2 * experimentInc + 1] >= threshold:
            numberTransferred += 1;
            #atom survived
            survivalRawData = append(survivalRawData, 1);
        elif data[2 * experimentInc] > threshold and data[2 * experimentInc + 1] < threshold:
            #atom didn't survive            
            survivalRawData= append(survivalRawData, 0);
        else:
            # no atom in the first place
            survivalRawData= append(survivalRawData, -1);
    # fractionTransferred = numberTransferred / numberOfExperiments;
    ###
    averageFractionTransfered = array([]);
    captureProbabilities = array([]);
    standardDeviation = array([]);
    ctsList = array([]);
    lctsList = array([]);
    stdevC = array([]);
    for variationInc in range(0, int(survivalRawData.size / (accumulations - 1)) ):
        cts = array([]);
        for accumulationInc in range(0, accumulations - 1):
            if survivalRawData[variationInc * accumulations + accumulationInc] != -1:
                cts = append(cts, survivalRawData[variationInc * accumulations + accumulationInc]);
        # catch the case where there's no relevant data, typically if laser becomes unlocked.
        if cts.size == 0:
            #print(cts.size)
            standardDeviation = append(standardDeviation, 0);
            captureProbabilities = append(captureProbabilities, 0);
            averageFractionTransfered = append(averageFractionTransfered, 0);
            ctsList = append(ctsList, 0);
            lctsList = append(lctsList, 0);
            stdevC = append(stdevC, 0);
        else:
            #print(cts.size)
            standardDeviation = append(standardDeviation, std(cts)/sqrt(cts.size));
            captureProbabilities = append(captureProbabilities, cts.size / accumulations);
            averageFractionTransfered = append(averageFractionTransfered, average(cts));
            ctsList = append(ctsList, average(cts.size));
            lctsList = append(lctsList, sqrt(cts.size));
            stdevC = append(stdevC, std(cts));
    dataSpectra = column_stack((key, averageFractionTransfered));
    survivalData = column_stack((dataSpectra, standardDeviation));
    fullCaptureProbabilityData = column_stack((array(key), captureProbabilities));
    return survivalData, fullCaptureProbabilityData, captureProbabilities 

def getAnalyzedTunnelingData(data, thresholds, key, accumulations, numberOfExperiments):
    from numpy import (append, array, average, column_stack, std, sqrt)
    # all the posibilities that I care about.
    firstToFirst = array([])
    firstToSecond = array([])
    secondToFirst = array([])
    secondToSecond = array([])
    bothToBoth = array([])
    bothToOne = array([])
    
    for experimentInc in range(0, numberOfExperiments):
        # start with atom in first and not in second
        if data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] < thresholds[1]:
            # if in second picture atom in first and not in second
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                firstToFirst = append(firstToFirst, 1)
                firstToSecond = append(firstToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                firstToFirst = append(firstToFirst, 0)
                firstToSecond = append(firstToSecond, 1)
            else:
                firstToFirst = append(firstToFirst, 0)
                firstToSecond = append(firstToSecond, 0)

            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            bothToBoth = append(bothToBoth, -1)
            bothToOne = append(bothToOne, -1)
        #start with atom in second and not first.
        elif data[0][2 * experimentInc] < thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)

            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                secondToFirst = append(secondToFirst, 1)
                secondToSecond = append(secondToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                secondToFirst = append(secondToFirst, 0)
                secondToSecond = append(secondToSecond, 1)
            else:
                secondToFirst = append(secondToFirst, 0)
                secondToSecond = append(secondToSecond, 0)
            bothToBoth = append(bothToBoth, -1)
            bothToOne = append(bothToOne, -1)
        # start with two atoms
        elif data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)
            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                bothToOne = append(bothToOne, 1)
                bothToBoth = append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = append(bothToOne, 1)
                bothToBoth = append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = append(bothToOne, 0)
                bothToBoth = append(bothToBoth, 1)
            else:
                bothToOne = append(bothToOne, 0)
                bothToBoth = append(bothToBoth, 0)
        # start with no atoms
        else:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)
            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            bothToOne = append(bothToOne, -1)
            bothToBoth = append(bothToBoth, -1)

    average_1to1 = array([])
    average_1to2 = array([])
    average_2to1 = array([])
    average_2to2 = array([])
    averageBothToBoth = array([])
    averageBothToOne = array([])
    captProbs = [[], []]
    error_1to1 = array([])
    error_1to2 = array([])
    error_2to1 = array([])
    error_2to2 = array([])
    error_BothToBoth = array([])
    error_BothToOne = array([])
    ctsList = array([])
    lctsList = array([])
    for variationInc in range(0, int(firstToFirst.size/accumulations)):
        firstToFirstData = array([])
        firstToSecondData = array([])
        secondToFirstData = array([])
        secondToSecondData = array([])
        data_BothToBoth = array([])
        data_BothToOne = array([])
        for accumInc in range (0, accumulations):
            picNum = variationInc * accumulations + accumInc
            if firstToFirst[picNum] != -1:
                firstToFirstData = append(firstToFirstData, firstToFirst[picNum])
            if firstToSecond[picNum] != -1:
                firstToSecondData = append(firstToSecondData, firstToSecond[picNum])
            if bothToBoth[picNum] != -1:
                data_BothToBoth = append(data_BothToBoth, bothToBoth[picNum])
            if bothToOne[picNum] != -1:
                data_BothToOne = append(data_BothToOne, bothToOne[picNum])
            if secondToFirst[picNum] != -1:
                secondToFirstData = append(secondToFirstData, secondToFirst[picNum])
            if secondToSecond[picNum] != -1:
                secondToSecondData = append(secondToSecondData, secondToSecond[picNum])
        captProbs[0] = append(captProbs[0], (firstToFirstData.size + data_BothToBoth.size) / accumulations)
        captProbs[1] = append(captProbs[1], (firstToSecondData.size + data_BothToBoth.size) / accumulations)
        average_1to1 = append(average_1to1, average(firstToFirstData))
        average_1to2 = append(average_1to2, average(firstToSecondData))
        average_2to1 = append(average_2to1, average(secondToFirstData))
        average_2to2 = append(average_2to2, average(secondToSecondData))
        averageBothToBoth = append(averageBothToBoth, average(data_BothToBoth))
        averageBothToOne = append(averageBothToOne, average(data_BothToOne))
        error_1to1 = append(error_1to1, std(firstToFirstData) / sqrt(firstToFirstData.size))
        error_1to2 = append(error_1to2, std(firstToSecondData) / sqrt(firstToSecondData.size))
        error_2to1 = append(error_2to1, std(secondToFirstData) / sqrt(secondToFirstData.size))
        error_2to2 = append(error_2to2, std(secondToSecondData) / sqrt(secondToSecondData.size))
        error_BothToBoth = append(error_BothToBoth, std(data_BothToBoth) / sqrt(data_BothToBoth.size))
        error_BothToOne = append(error_BothToOne, std(data_BothToOne) / sqrt(data_BothToOne.size))
        ctsList = append(ctsList, average(firstToFirstData))
        lctsList = append(lctsList, sqrt(firstToFirstData.size))

    # condense data for export.
    survival_1 = average_1to1+average_1to2
    error_survival_1 = sqrt(error_1to1**2 + error_1to2**2)
    survival_2 = average_2to1+average_2to2
    error_survival_2 = sqrt(error_2to1**2 + error_2to2**2)
    # lots to return. Might be better way to do this.
    return (average_1to1, error_1to1, average_1to2, error_1to2, average_2to1, error_2to1, average_2to2, error_2to2,
            averageBothToBoth, error_BothToBoth, averageBothToOne, error_BothToOne, captProbs[0], captProbs[1],
            survival_1, error_survival_1, survival_2, error_survival_2)


def getPostSelectedTunnelingData(data, thresholds, key, accumulations, numberOfExperiments):
    from numpy import (append, array, average, column_stack, std, sqrt)
    # all the posibilities that I care about.
    firstToFirst = array([])
    firstToSecond = array([])
    secondToFirst = array([])
    secondToSecond = array([])
    bothToBoth = array([])
    bothToOne = array([])
    # analyze all data points, looking for atoms.
    for experimentInc in range(0, numberOfExperiments):
        # start with atom in first and not in second
        if data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] < thresholds[1]:
            # if in second picture atom in first and not in second
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                firstToFirst = append(firstToFirst, 1)
                firstToSecond = append(firstToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                firstToFirst = append(firstToFirst, 0)
                firstToSecond = append(firstToSecond, 1)
            else:
                # this is the difference between post-selection and not. Not post-selected means that these data points
                # go to zero and get counted.
                firstToFirst = append(firstToFirst, -1)
                firstToSecond = append(firstToSecond, -1)
            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            bothToBoth = append(bothToBoth, -1)
            bothToOne = append(bothToOne, -1)
        # start with atom in second and not first.
        elif data[0][2 * experimentInc] < thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                secondToFirst = append(secondToFirst, 1)
                secondToSecond = append(secondToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                secondToFirst = append(secondToFirst, 0)
                secondToSecond = append(secondToSecond, 1)
            else:
                # post-select.
                secondToFirst = append(secondToFirst, -1)
                secondToSecond = append(secondToSecond, -1)
            bothToBoth = append(bothToBoth, -1)
            bothToOne = append(bothToOne, -1)
        # start with two atoms
        elif data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)
            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                bothToOne = append(bothToOne, 1)
                bothToBoth = append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                # other case for both-to-one.
                bothToOne = append(bothToOne, 1)
                bothToBoth = append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = append(bothToOne, 0)
                bothToBoth = append(bothToBoth, 1)
            else:
                # there is no post-select on two-atom loss.
                bothToOne = append(bothToOne, 0)
                bothToBoth = append(bothToBoth, 0)
        # start with no atoms
        else:
            firstToFirst = append(firstToFirst, -1)
            firstToSecond = append(firstToSecond, -1)
            secondToFirst = append(secondToFirst, -1)
            secondToSecond = append(secondToSecond, -1)
            bothToOne = append(bothToOne, -1)
            bothToBoth = append(bothToBoth, -1)

    average_1to1 = array([])
    average_1to2 = array([])
    average_2to1 = array([])
    average_2to2 = array([])
    averageBothToBoth = array([])
    averageBothToOne = array([])
    captProbs = [[], []]
    error_1to1 = array([])
    error_1to2 = array([])
    error_2to1 = array([])
    error_2to2 = array([])
    error_BothToBoth = array([])
    error_BothToOne = array([])
    ctsList = array([])
    lctsList = array([])
    for variationInc in range(0, int(firstToFirst.size / accumulations)):
        firstToFirstData = array([])
        firstToSecondData = array([])
        secondToFirstData = array([])
        secondToSecondData = array([])
        data_BothToBoth = array([])
        data_BothToOne = array([])
        for accumInc in range(0, accumulations):
            picNum = variationInc * accumulations + accumInc
            if firstToFirst[picNum] != -1:
                firstToFirstData = append(firstToFirstData, firstToFirst[picNum])
            if firstToSecond[picNum] != -1:
                firstToSecondData = append(firstToSecondData, firstToSecond[picNum])
            if bothToBoth[picNum] != -1:
                data_BothToBoth = append(data_BothToBoth, bothToBoth[picNum])
            if bothToOne[picNum] != -1:
                data_BothToOne = append(data_BothToOne, bothToOne[picNum])
            if secondToFirst[picNum] != -1:
                secondToFirstData = append(secondToFirstData, secondToFirst[picNum])
            if secondToSecond[picNum] != -1:
                secondToSecondData = append(secondToSecondData, secondToSecond[picNum])

        captProbs[0] = append(captProbs[0], (firstToFirstData.size + data_BothToBoth.size) / accumulations)
        captProbs[1] = append(captProbs[1], (secondToSecondData.size + data_BothToBoth.size) / accumulations)
        average_1to1 = append(average_1to1, average(firstToFirstData))
        average_1to2 = append(average_1to2, average(firstToSecondData))
        average_2to1 = append(average_2to1, average(secondToFirstData))
        average_2to2 = append(average_2to2, average(secondToSecondData))
        averageBothToBoth = append(averageBothToBoth, average(data_BothToBoth))
        averageBothToOne = append(averageBothToOne, average(data_BothToOne))
        error_1to1 = append(error_1to1, std(firstToFirstData) / sqrt(firstToFirstData.size))
        error_1to2 = append(error_1to2, std(firstToSecondData) / sqrt(firstToSecondData.size))
        error_2to1 = append(error_2to1, std(secondToFirstData) / sqrt(secondToFirstData.size))
        error_2to2 = append(error_2to2, std(secondToSecondData) / sqrt(secondToSecondData.size))
        error_BothToBoth = append(error_BothToBoth, std(data_BothToBoth) / sqrt(data_BothToBoth.size))
        error_BothToOne = append(error_BothToOne, std(data_BothToOne) / sqrt(data_BothToOne.size))
        ctsList = append(ctsList, average(firstToFirstData))
        lctsList = append(lctsList, sqrt(firstToFirstData.size))
    # condense data for export.
    survival_1 = average_1to1 + average_1to2
    error_survival_1 = sqrt(error_1to1 ** 2 + error_1to2 ** 2)
    survival_2 = average_2to1 + average_2to2
    error_survival_2 = sqrt(error_2to1 ** 2 + error_2to2 ** 2)

    # lots to return. Might be better way to do this.
    return (average_1to1, error_1to1, average_1to2, error_1to2, average_2to1, error_2to1, average_2to2, error_2to2,
            averageBothToBoth, error_BothToBoth, averageBothToOne, error_BothToOne, captProbs[0], captProbs[1],
            survival_1, error_survival_1, survival_2, error_survival_2)

