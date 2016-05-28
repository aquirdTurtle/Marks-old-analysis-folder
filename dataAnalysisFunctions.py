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

This probably doesn't work for data sets with very long trails. Should be easy to modify for that case though.
get range
"""
def guessGaussianPeaks(rawData, binCenters, binnedData):
    from numpy import absolute, argmax, min, max
    leftBorder = min(rawData);
    rightBorder = max(rawData);
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
        while (found == False):
            if (binCenters[locInc] < (leftBorder + 2 * distToLeft)):
                locInc += 1;
            else:
                found = True;
        subBinnedData = binnedData[locInc:binnedData.size];
        # get index corresponding to second local maxima
        guess2Index = argmax(subBinnedData) + binnedData.size - subBinnedData.size;
    elif (distToLeft >= distToRight ):
        found = False;
        locInc = 0;
        while (found == False):
            if (binCenters[locInc] < (rightBorder - 2 * distToRight)):
                locInc += 1;
            else:
                found = True;
        subBinnedData = binnedData[0:locInc];
        # get index corresponding to second local maxima
        guess2Index = argmax(subBinnedData)
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
    dataSpectra = column_stack((key, averageFractionTransfered))
    survivalData = column_stack((dataSpectra, standardDeviation));
    fullCaptureProbabilityData = column_stack((array(key), captureProbabilities));
    return survivalData, fullCaptureProbabilityData, captureProbabilities 

def getAnalyzedTunnelingData(data, thresholds, key, accumulations, numberOfExperiments):
    from numpy import (append, array, average, column_stack, std, sqrt)
    # all the posibilities that I care about.
    firstToFirst = array([]);
    firstToSecond = array([]);
    secondToFirst = array([]);
    secondToSecond = array([]);
    bothToBoth = array([]);
    bothToOne = array([]);
    
    for experimentInc in range(0, numberOfExperiments):
        # start with atom in first and not in second
        if data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] < thresholds[1]:
            # if in second picture atom in first and not in second
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                firstToFirst = append(firstToFirst, 1);
                firstToSecond = append(firstToSecond, 0);
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                firstToFirst = append(firstToFirst, 0);
                firstToSecond = append(firstToSecond, 1);
            else:
                firstToFirst = append(firstToFirst, 0);
                firstToSecond = append(firstToSecond, 0);

            secondToFirst = append(secondToFirst, -1);
            secondToSecond = append(secondToSecond, -1);                
            bothToBoth = append(bothToBoth, -1);
            bothToOne = append(bothToOne, -1);
        #start with atom in second and not first.
        elif data[0][2 * experimentInc] < thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1);
            firstToSecond = append(firstToSecond, -1);

            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                secondToFirst = append(secondToFirst, 1);
                secondToSecond = append(secondToSecond, 0);
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                secondToFirst = append(secondToFirst, 0);
                secondToSecond = append(secondToSecond, 1);
            else:
                secondToFirst = append(secondToFirst, 0);
                secondToSecond = append(secondToSecond, 0);
            bothToBoth = append(bothToBoth, -1);
            bothToOne = append(bothToOne, -1);
        # start with two atoms
        elif data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = append(firstToFirst, -1);
            firstToSecond = append(firstToSecond, -1);
            secondToFirst = append(secondToFirst, -1);
            secondToSecond = append(secondToSecond, -1);      
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                bothToOne = append(bothToOne, 1);
                bothToBoth = append(bothToBoth, 0);
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = append(bothToOne, 1);
                bothToBoth = append(bothToBoth, 0);
            elif data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = append(bothToOne, 0);
                bothToBoth = append(bothToBoth, 1);                
            else:
                bothToOne = append(bothToOne, 0);
                bothToBoth = append(bothToBoth, 0);
        # start with no atoms
        else:
            firstToFirst = append(firstToFirst, -1);
            firstToSecond = append(firstToSecond, -1);
            secondToFirst = append(secondToFirst, -1);
            secondToSecond = append(secondToSecond, -1);                  
            bothToOne = append(bothToOne, -1);
            bothToBoth = append(bothToBoth, -1);
    averageFirstToFirst = array([]);
    averageFirstToSecond = array([]);
    averageBothToBoth = array([]);
    averageBothToOne = array([]);
    captProbs = array([]);
    stdev = array([]);
    stdev2 = array([]);
    stdevBoth  = array([]);
    stdevBoth11 = array([]);
    ctsList = array([]);
    lctsList = array([]);
    stdevC = array([]);
    for variationInc in range(0, int(firstToFirst.size/accumulations)):
        cts = array([]);
        cts2 = array([]);
        ctsBoth = array([]);
        ctsBoth11 = array([]);
        for accumInc in range (0, accumulations):
            picNum = variationInc * accumulations + accumInc;
            if (firstToFirst[picNum] != -1):
                cts = append(cts, firstToFirst[picNum]);
            if (firstToSecond[picNum] != -1):
                cts2 = append(cts2, firstToSecond[picNum]);
            if (bothToBoth[picNum] != -1):
                ctsBoth = append(ctsBoth, bothToBoth[picNum]);
            if (bothToOne[picNum] != -1):
                ctsBoth11 = append(ctsBoth11, bothToOne[picNum]);
        captProbs = append(captProbs, cts.size / accumulations);
        averageFirstToFirst = append(averageFirstToFirst, average(cts));
        averageFirstToSecond = append(averageFirstToSecond, average(cts2));
        averageBothToBoth = append(averageBothToBoth, average(ctsBoth));
        averageBothToOne = append(averageBothToOne, average(ctsBoth11));
        stdev = append(stdev, std(cts) / sqrt(cts.size));
        stdev2 = append(stdev2, std(cts2) / sqrt(cts2.size));
        stdevBoth = append(stdevBoth, std(ctsBoth) / sqrt(ctsBoth.size));
        stdevBoth11 = append(stdevBoth11, std(ctsBoth11) / sqrt(ctsBoth11.size));
        ctsList = append(ctsList, average(cts));
        lctsList = append(lctsList, sqrt(cts.size));
        stdevC = append(stdevC, std(cts));
    # condense data for export.
    dataSpectra = column_stack((key, averageFirstToFirst));
    dataSpectra2 = column_stack((key, averageFirstToSecond));
    dataSpectraBoth = column_stack((key, averageBothToBoth));
    dataSpectraBoth11 = column_stack((key, averageBothToOne));
    toPlot = column_stack((dataSpectra, stdev));
    toPlot2 = column_stack((dataSpectra2, stdev2));
    toPlotBoth = column_stack((dataSpectraBoth, stdevBoth));
    toPlotBoth11 = column_stack((dataSpectraBoth11, stdevBoth11));
    toPlotSum = column_stack((key, averageFirstToFirst + averageFirstToSecond, sqrt(stdev**2 + stdev2**2)));
    return toPlot, toPlot2, toPlotBoth, toPlotBoth11, captProbs, toPlotSum;
   