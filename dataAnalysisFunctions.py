# -*- coding: utf-8 -*-
"""
Created on Sat May 21 10:45:58 2016

@author: Mark
"""


def normalizeData(picsPerExperiment, rawData, atomLocation):
    """
    This function analyzes raw data at the location "atom location" and returns a normalized version of
    (1) every picture
    (2) an array containing only the first pictures
    :param picsPerExperiment:
    :param rawData:
    :param atomLocation:
    :return:
    """
    from numpy import append, array
    firstData = array([])
    allData = array([])
    dimensions = rawData.shape
    for imageInc in range(0, dimensions[0]):
        averageBackground = 1/4*(rawData[imageInc][0][0] 
                                 + rawData[imageInc][dimensions[1]-1][dimensions[2] - 1]
                                 + rawData[imageInc][0][dimensions[2] - 1]
                                 + rawData[imageInc][dimensions[1]-1][0])
        allData = append(allData, rawData[imageInc][atomLocation[0]][atomLocation[1]] - averageBackground)
        if imageInc % picsPerExperiment == 0:
            firstData = append(firstData, rawData[imageInc][atomLocation[0]][atomLocation[1]] - averageBackground)
    return allData, firstData


def binData(binWidth, data):
    from numpy import append, array, min, max, histogram
    binBorderLocation = min(data)
    binsBorders = array([])
    # get bin borders
    while binBorderLocation < max(data):
        binsBorders = append(binsBorders, binBorderLocation)
        binBorderLocation = binBorderLocation + binWidth
    # trash gets set but is unused.
    binnedData, trash = histogram(data, binsBorders)
    binCenters = binsBorders[0:binsBorders.size-1]
    return binCenters, binnedData


def poisson(x, k, weight):
    """
    This function calculates $p_k{x} = norm * e^(-k) * k^x / x!.
    :param x: argument of the poissonian
    :param k: order or (approximate) mean of the poissonian.
    :param weight: a weight factor, related to the maximum data this is supposed to be fitted to, but typically over-
    weighted for the purposes of this function.
    :return: the poissonian evaluated at x given the parametes.
    """
    import numpy as np
    term = 1
    # calculate the term k^x / x!. Can't do this directly, x! is too large.
    for n in range(0, int(x)):
        term *= k / (x - n)
    return np.exp(-k) * term * weight


def guessGaussianPeaks(rawData, binCenters, binnedData):
    """
    This function guesses where the gaussian peaks of the data are. It assumes one is near the maximum of the binned
    data. Then, from the binned data it subtracts an over-weighted (i.e. extra tall) poissonion distribution e^-k k^n/n!
    From the binned data. This should squelch the peak that it found. It then assumes that the second peak is near the
    maximum of the (data-poissonian) array.
    :param binCenters: The pixel-numbers corresponding to the binned data data points.
    :param binnedData: the binned data data points.
    :return: the two guesses.
    """
    import numpy as np
    binCenters += 500
    # get index corresponding to global max
    guess1Index = np.argmax(binnedData)
    # get location of global max
    guess1Location = binCenters[guess1Index]
    binnedDataWithoutPoissonian = []
    for binInc in range(0, len(binCenters)):
        binnedDataWithoutPoissonian.append(binnedData[binInc]
                                           - poisson(binCenters[binInc], guess1Location, 2 * max(binnedData)
                                                     / poisson(guess1Location, guess1Location, 1)))
    guess2Index = np.argmax(binnedDataWithoutPoissonian)
    guess2Location = binCenters[guess2Index]
    return guess1Location - 500, guess2Location - 500


def doubleGaussian(data, A1, x1, sig1, A2, x2, sig2):
    from numpy import absolute,  exp, sqrt
    return ((absolute(A1) * exp(-((data-x1)/(sqrt(2)*sig1))**2)) 
            + absolute(A2) * exp(-((data-x2)/(sqrt(2)*sig2))**2))


def fitDoubleGaussian(binCenters, binnedData, fitGuess):
    from scipy.optimize import curve_fit
    fitVals, trash = curve_fit(doubleGaussian, binCenters, binnedData, fitGuess)
    return fitVals


def calculateAtomThreshold(fitVals):
    from numpy import (sqrt, abs)
    from scipy.special import erf
    TCalc = (fitVals[4] - fitVals[1])/(abs(fitVals[5]) + abs(fitVals[2]))
    threshold = fitVals[1] + TCalc * fitVals[2]
    fidelity = 1/2 * (1 + erf(abs(TCalc)/sqrt(2)))
    return threshold, fidelity


def getAtomData(data, threshold, numberOfExperiments):
    """
    This function assumes 2 pictures.
    It returns a raw array that includes every survival data point, including points where the the atom doesn't get
    loaded at all.
    """
    import numpy as np
    # this will include entries for when there is no atom in the first picture.
    survivalData = np.array([])
    survivalData.astype(int)
    # this doesn't take into account loss, since these experiments are feeding-back on loss.
    for experimentInc in range(0, numberOfExperiments):
        if data[2 * experimentInc] > threshold and data[2 * experimentInc + 1] >= threshold:
            # atom survived
            survivalData = np.append(survivalData, 1)
        elif data[2 * experimentInc] > threshold > data[2 * experimentInc + 1]:
            # atom didn't survive 
            survivalData = np.append(survivalData, 0)
        else:
            # no atom in the first place
            survivalData = np.append(survivalData, -1)
    return survivalData


def getSingleParticleSurvivalData(survivalData, repetitionsPerVariation):
    import numpy as np
    # Take the previous data, which includes entries when there was no atom in the first picture, and convert it to
    # an array of just loaded and survived or loaded and died.    
    survivalAverages = np.array([])
    loadingProbability = np.array([])
    survivalErrors = np.array([])
    for variationInc in range(0, int(survivalData.size / repetitionsPerVariation)):
        survivalList = np.array([])
        for repetitionInc in range(0, repetitionsPerVariation):
            if survivalData[variationInc * repetitionsPerVariation + repetitionInc] != -1:
                survivalList = np.append(survivalList,
                                         survivalData[variationInc * repetitionsPerVariation + repetitionInc])
        if survivalList.size == 0:
            # catch the case where there's no relevant data, typically if laser becomes unlocked.
            survivalErrors = np.append(survivalErrors, 0)
            loadingProbability = np.append(loadingProbability, 0)
            survivalAverages = np.append(survivalAverages, 0)
        else:
            # normal case
            survivalErrors = np.append(survivalErrors, np.std(survivalList)/np.sqrt(survivalList.size))
            loadingProbability = np.append(loadingProbability, survivalList.size / repetitionsPerVariation)
            survivalAverages = np.append(survivalAverages, np.average(survivalList))

    return survivalAverages, survivalErrors, loadingProbability


def getCorrelationData(allAtomSurvivalData, repetitionsPerVariation):
    from collections import OrderedDict as dic
    import numpy as np
    import math
    atomNum = allAtomSurvivalData.shape[0]
    repNum = allAtomSurvivalData.shape[1]
    correlationErrors = dic()
    correlationAverages = dic()
    correlationErrors['Key List'] = ''
    correlationAverages['Key List'] = ''
    # initialize the average dicts
    for atomsLoadedInc in range(1, atomNum + 1):
        for atomSurvivedInc in range(0, atomNum):
            name = 'Load ' + str(atomsLoadedInc) + ', atom ' + str(atomSurvivedInc) + ' survived'
            correlationErrors[name] = []
            correlationAverages[name] = []
    # data that doesn't discriminate between locations.
    for atomsLoadedInc in range(1, atomNum + 1):
        # + 1 because all atoms could survive.
        for atomSurvivedInc in range(0, atomsLoadedInc + 1):
            name = 'Load ' + str(atomsLoadedInc) + ', ' + str(atomSurvivedInc) + ' atoms survived'
            correlationErrors[name] = []
            correlationAverages[name] = []
    # holds the averages over wells
    for atomsLoadedInc in range(1, atomNum + 1):
        name = 'Load ' + str(atomsLoadedInc) + ', average single atom survival'
        correlationErrors[name] = []
        correlationAverages[name] = []
    # holds the average over all wells and loading scenarios.
    name = 'Total average single atom survival'
    correlationErrors[name] = []
    correlationAverages[name] = []
    # Start sorting data.
    for variationInc in range(0, int(repNum / repetitionsPerVariation)):
        totalNumberLoadedList = []
        totalNumberSurvivedList = []
        for repInc in range(0, repetitionsPerVariation):
            totalAtomsLoaded = 0
            totalAtomsSurvived = 0
            holeFlag = False
            for atomInc in range(0, atomNum):
                if allAtomSurvivalData[atomInc][variationInc * repetitionsPerVariation + repInc] == 0:
                    if holeFlag:
                        totalAtomsSurvived = math.nan
                        totalAtomsLoaded = math.nan
                        break
                    totalAtomsLoaded += 1
                elif allAtomSurvivalData[atomInc][variationInc * repetitionsPerVariation + repInc] == 1:
                    if holeFlag:
                        totalAtomsSurvived = math.nan
                        totalAtomsLoaded = math.nan
                        break
                    totalAtomsSurvived += 1
                    totalAtomsLoaded += 1
                else:
                    # no atom loaded here.
                    if totalAtomsLoaded > 0:
                        holeFlag = True
            totalNumberLoadedList = np.append(totalNumberLoadedList, totalAtomsLoaded)
            totalNumberSurvivedList = np.append(totalNumberSurvivedList, totalAtomsSurvived)
        # initialize entries in temporary dictionary. Sanity, mostly, probably a way around this.
        tempCorrelationData = dic()
        for atomsLoadedInc in range(1, atomNum + 1):
            for atomSurvivedInc in range(0, atomNum):
                name = 'Load ' + str(atomsLoadedInc) + ', atom ' + str(atomSurvivedInc) + ' survived'
                tempCorrelationData[name] = []
        # data that doesn't discriminate between locations.
        for atomsLoadedInc in range(1, atomNum + 1):
            # + 1 because all atoms could survive.
            for atomSurvivedInc in range(0, atomsLoadedInc + 1):
                name = 'Load ' + str(atomsLoadedInc) + ', ' + str(atomSurvivedInc) + ' atoms survived'
                tempCorrelationData[name] = []
        # holds the averages over wells
        for atomsLoadedInc in range(1, atomNum + 1):
            name = 'Load ' + str(atomsLoadedInc) + ', average single atom survival'
            tempCorrelationData[name] = []
        # holds the average over all wells and loading scenarios.
        tempCorrelationData['Total average single atom survival'] = []
        # get data for specific particles surviving.
        for atomInc in range(0, atomNum):
            for repInc in range(0, repetitionsPerVariation):
                if math.isnan(totalNumberLoadedList[repInc]) or totalNumberLoadedList[repInc] == 0:
                    # no atoms loaded or hole so throw the data out.
                    continue
                name = 'Load ' + str(int(totalNumberLoadedList[repInc])) \
                       + ', atom ' + str(atomInc) + ' survived'
                name2 = 'Load ' + str(int(totalNumberLoadedList[repInc])) + ', average single atom survival'
                # make sure *THIS* atom was loaded originally.
                if allAtomSurvivalData[atomInc][variationInc * repetitionsPerVariation + repInc] != -1:
                    value = allAtomSurvivalData[atomInc][variationInc * repetitionsPerVariation + repInc]
                    tempCorrelationData[name] = np.append(tempCorrelationData[name], value)
                    tempCorrelationData[name2] = np.append(tempCorrelationData[name2], value)
                    tempCorrelationData['Total average single atom survival'] \
                        = np.append(tempCorrelationData['Total average single atom survival'], value)
        # get indescriminatory data.
        for repInc in range(0, repetitionsPerVariation):
            if math.isnan(totalNumberLoadedList[repInc]) or int(totalNumberLoadedList[repInc]) == 0:
                # throw away holes
                continue
            for atomsSurvivedInc in range(0, int(totalNumberLoadedList[repInc]) + 1):
                name = 'Load ' + str(int(totalNumberLoadedList[repInc])) + ', ' + str(atomsSurvivedInc) \
                       + ' atoms survived'
                if totalNumberSurvivedList[repInc] == atomsSurvivedInc:
                    tempCorrelationData[name] = np.append(tempCorrelationData[name], 1)
                else:
                    tempCorrelationData[name] = np.append(tempCorrelationData[name], 0)
        # calculate averages
        for atomsLoadedInc in range(1, atomNum + 1):
            for atomSurvivedInc in range(0, atomNum):
                name = 'Load ' + str(atomsLoadedInc) + ', atom ' + str(atomSurvivedInc) + ' survived'
                correlationAverages[name] = np.append(correlationAverages[name], np.mean(tempCorrelationData[name]))
                correlationErrors[name] = np.append(correlationErrors[name], np.std(tempCorrelationData[name])
                                                    / np.sqrt(len(tempCorrelationData[name])))
            name2 = 'Load ' + str(atomsLoadedInc) + ', average single atom survival'

            correlationAverages[name2] = np.append(correlationAverages[name2], np.mean(tempCorrelationData[name2]))
            correlationErrors[name2] = np.append(correlationErrors[name2],
                                                 np.std(tempCorrelationData[name2])
                                                 / np.sqrt(len(tempCorrelationData[name2])))
        totalName = 'Total average single atom survival'
        correlationAverages[totalName] = np.append(correlationAverages[totalName],
                                                   np.mean(tempCorrelationData[totalName]))
        correlationErrors[totalName] = np.append(correlationErrors[totalName],
                                                   np.std(tempCorrelationData[totalName])
                                                 / np.sqrt(len(tempCorrelationData[totalName])))
        # data that doesn't discriminate between locations.
        for atomsLoadedInc in range(1, atomNum + 1):
            # + 1 because all atoms could survive.
            for atomSurvivedInc in range(0, atomsLoadedInc + 1):
                name = 'Load ' + str(atomsLoadedInc) + ', ' + str(atomSurvivedInc) + ' atoms survived'
                correlationAverages[name] = np.append(correlationAverages[name], np.mean(tempCorrelationData[name]))
                correlationErrors[name] = np.append(correlationErrors[name], np.std(tempCorrelationData[name])
                                                    / np.sqrt(len(tempCorrelationData[name])))
    correlationErrors['Key List'] = list(correlationErrors.keys())
    correlationAverages['Key List'] = list(correlationAverages.keys())
    return correlationAverages, correlationErrors


def getAnalyzedTunnelingData(data, thresholds, key, accumulations, numberOfExperiments):
    import numpy as np
    # all the possibilities that I care about.
    firstToFirst = np.array([])
    firstToSecond = np.array([])
    secondToFirst = np.array([])
    secondToSecond = np.array([])
    bothToBoth = np.array([])
    bothToOne = np.array([])
    
    for experimentInc in range(0, numberOfExperiments):
        # start with atom in first and not in second
        if data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] < thresholds[1]:
            # if in second picture atom in first and not in second
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                firstToFirst = np.append(firstToFirst, 1)
                firstToSecond = np.append(firstToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                firstToFirst = np.append(firstToFirst, 0)
                firstToSecond = np.append(firstToSecond, 1)
            else:
                firstToFirst = np.append(firstToFirst, 0)
                firstToSecond = np.append(firstToSecond, 0)

            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            bothToBoth = np.append(bothToBoth, -1)
            bothToOne = np.append(bothToOne, -1)
        # start with atom in second and not first.
        elif data[0][2 * experimentInc] < thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)

            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                secondToFirst = np.append(secondToFirst, 1)
                secondToSecond = np.append(secondToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                secondToFirst = np.append(secondToFirst, 0)
                secondToSecond = np.append(secondToSecond, 1)
            else:
                secondToFirst = np.append(secondToFirst, 0)
                secondToSecond = np.append(secondToSecond, 0)
            bothToBoth = np.append(bothToBoth, -1)
            bothToOne = np.append(bothToOne, -1)
        # start with two atoms
        elif data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)
            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                bothToOne = np.append(bothToOne, 1)
                bothToBoth = np.append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = np.append(bothToOne, 1)
                bothToBoth = np.append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = np.append(bothToOne, 0)
                bothToBoth = np.append(bothToBoth, 1)
            else:
                bothToOne = np.append(bothToOne, 0)
                bothToBoth = np.append(bothToBoth, 0)
        # start with no atoms
        else:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)
            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            bothToOne = np.append(bothToOne, -1)
            bothToBoth = np.append(bothToBoth, -1)

    average_1to1 = np.array([])
    average_1to2 = np.array([])
    average_2to1 = np.array([])
    average_2to2 = np.array([])
    averageBothToBoth = np.array([])
    averageBothToOne = np.array([])
    captProbs = [[], []]
    error_1to1 = np.array([])
    error_1to2 = np.array([])
    error_2to1 = np.array([])
    error_2to2 = np.array([])
    error_BothToBoth = np.array([])
    error_BothToOne = np.array([])
    ctsList = np.array([])
    lctsList = np.array([])
    for variationInc in range(0, int(firstToFirst.size/accumulations)):
        firstToFirstData = np.array([])
        firstToSecondData = np.array([])
        secondToFirstData = np.array([])
        secondToSecondData = np.array([])
        data_BothToBoth = np.array([])
        data_BothToOne = np.array([])
        for accumInc in range(0, accumulations):
            picNum = variationInc * accumulations + accumInc
            if firstToFirst[picNum] != -1:
                firstToFirstData = np.append(firstToFirstData, firstToFirst[picNum])
            if firstToSecond[picNum] != -1:
                firstToSecondData = np.append(firstToSecondData, firstToSecond[picNum])
            if bothToBoth[picNum] != -1:
                data_BothToBoth = np.append(data_BothToBoth, bothToBoth[picNum])
            if bothToOne[picNum] != -1:
                data_BothToOne = np.append(data_BothToOne, bothToOne[picNum])
            if secondToFirst[picNum] != -1:
                secondToFirstData = np.append(secondToFirstData, secondToFirst[picNum])
            if secondToSecond[picNum] != -1:
                secondToSecondData = np.append(secondToSecondData, secondToSecond[picNum])
        captProbs[0] = np.append(captProbs[0], (firstToFirstData.size + data_BothToBoth.size) / accumulations)
        captProbs[1] = np.append(captProbs[1], (firstToSecondData.size + data_BothToBoth.size) / accumulations)
        average_1to1 = np.append(average_1to1, np.average(firstToFirstData))
        average_1to2 = np.append(average_1to2, np.average(firstToSecondData))
        average_2to1 = np.append(average_2to1, np.average(secondToFirstData))
        average_2to2 = np.append(average_2to2, np.average(secondToSecondData))
        averageBothToBoth = np.append(averageBothToBoth, np.average(data_BothToBoth))
        averageBothToOne = np.append(averageBothToOne, np.average(data_BothToOne))
        error_1to1 = np.append(error_1to1, np.std(firstToFirstData) / np.sqrt(firstToFirstData.size))
        error_1to2 = np.append(error_1to2, np.std(firstToSecondData) / np.sqrt(firstToSecondData.size))
        error_2to1 = np.append(error_2to1, np.std(secondToFirstData) / np.sqrt(secondToFirstData.size))
        error_2to2 = np.append(error_2to2, np.std(secondToSecondData) / np.sqrt(secondToSecondData.size))
        error_BothToBoth = np.append(error_BothToBoth, np.std(data_BothToBoth) / np.sqrt(data_BothToBoth.size))
        error_BothToOne = np.append(error_BothToOne, np.std(data_BothToOne) / np.sqrt(data_BothToOne.size))
        ctsList = np.append(ctsList, np.average(firstToFirstData))
        lctsList = np.append(lctsList, np.sqrt(firstToFirstData.size))

    # condense data for export.
    survival_1 = average_1to1+average_1to2
    error_survival_1 = np.sqrt(error_1to1**2 + error_1to2**2)
    survival_2 = average_2to1+average_2to2
    error_survival_2 = np.sqrt(error_2to1**2 + error_2to2**2)
    # lots to return. Might be better way to do this.
    return (average_1to1, error_1to1, average_1to2, error_1to2, average_2to1, error_2to1, average_2to2, error_2to2,
            averageBothToBoth, error_BothToBoth, averageBothToOne, error_BothToOne, captProbs[0], captProbs[1],
            survival_1, error_survival_1, survival_2, error_survival_2)


def getPostSelectedTunnelingData(data, thresholds, key, accumulations, numberOfExperiments):
    """

    :param data:
    :param thresholds:
    :param key:
    :param accumulations:
    :param numberOfExperiments:
    :return:
    """
    import numpy as np
    # all the possibilities that I care about.
    firstToFirst = np.array([])
    firstToSecond = np.array([])
    secondToFirst = np.array([])
    secondToSecond = np.array([])
    bothToBoth = np.array([])
    bothToOne = np.array([])
    # analyze all data points, looking for atoms.
    for experimentInc in range(0, numberOfExperiments):
        # start with atom in first and not in second
        if data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] < thresholds[1]:
            # if in second picture atom in first and not in second
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                firstToFirst = np.append(firstToFirst, 1)
                firstToSecond = np.append(firstToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                firstToFirst = np.append(firstToFirst, 0)
                firstToSecond = np.append(firstToSecond, 1)
            else:
                # this is the difference between post-selection and not. Not post-selected means that these data points
                # go to zero and get counted.
                firstToFirst = np.append(firstToFirst, -1)
                firstToSecond = np.append(firstToSecond, -1)
            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            bothToBoth = np.append(bothToBoth, -1)
            bothToOne = np.append(bothToOne, -1)
        # start with atom in second and not first.
        elif data[0][2 * experimentInc] < thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                secondToFirst = np.append(secondToFirst, 1)
                secondToSecond = np.append(secondToSecond, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                secondToFirst = np.append(secondToFirst, 0)
                secondToSecond = np.append(secondToSecond, 1)
            else:
                # post-select.
                secondToFirst = np.append(secondToFirst, -1)
                secondToSecond = np.append(secondToSecond, -1)
            bothToBoth = np.append(bothToBoth, -1)
            bothToOne = np.append(bothToOne, -1)
        # start with two atoms
        elif data[0][2 * experimentInc] > thresholds[0] and data[1][2 * experimentInc] > thresholds[1]:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)
            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            if data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] < thresholds[1]:
                bothToOne = np.append(bothToOne, 1)
                bothToBoth = np.append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] < thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                # other case for both-to-one.
                bothToOne = np.append(bothToOne, 1)
                bothToBoth = np.append(bothToBoth, 0)
            elif data[0][2 * experimentInc + 1] > thresholds[0] and data[1][2 * experimentInc + 1] > thresholds[1]:
                bothToOne = np.append(bothToOne, 0)
                bothToBoth = np.append(bothToBoth, 1)
            else:
                # there is no post-select on two-atom loss.
                bothToOne = np.append(bothToOne, 0)
                bothToBoth = np.append(bothToBoth, 0)
        # start with no atoms
        else:
            firstToFirst = np.append(firstToFirst, -1)
            firstToSecond = np.append(firstToSecond, -1)
            secondToFirst = np.append(secondToFirst, -1)
            secondToSecond = np.append(secondToSecond, -1)
            bothToOne = np.append(bothToOne, -1)
            bothToBoth = np.append(bothToBoth, -1)

    average_1to1 = np.array([])
    average_1to2 = np.array([])
    average_2to1 = np.array([])
    average_2to2 = np.array([])
    averageBothToBoth = np.array([])
    averageBothToOne = np.array([])
    captProbs = [[], []]
    error_1to1 = np.array([])
    error_1to2 = np.array([])
    error_2to1 = np.array([])
    error_2to2 = np.array([])
    error_BothToBoth = np.array([])
    error_BothToOne = np.array([])
    ctsList = np.array([])
    lctsList = np.array([])
    for variationInc in range(0, int(firstToFirst.size / accumulations)):
        firstToFirstData = np.array([])
        firstToSecondData = np.array([])
        secondToFirstData = np.array([])
        secondToSecondData = np.array([])
        data_BothToBoth = np.array([])
        data_BothToOne = np.array([])
        for accumInc in range(0, accumulations):
            picNum = variationInc * accumulations + accumInc
            if firstToFirst[picNum] != -1:
                firstToFirstData = np.append(firstToFirstData, firstToFirst[picNum])
            if firstToSecond[picNum] != -1:
                firstToSecondData = np.append(firstToSecondData, firstToSecond[picNum])
            if bothToBoth[picNum] != -1:
                data_BothToBoth = np.append(data_BothToBoth, bothToBoth[picNum])
            if bothToOne[picNum] != -1:
                data_BothToOne = np.append(data_BothToOne, bothToOne[picNum])
            if secondToFirst[picNum] != -1:
                secondToFirstData = np.append(secondToFirstData, secondToFirst[picNum])
            if secondToSecond[picNum] != -1:
                secondToSecondData = np.append(secondToSecondData, secondToSecond[picNum])

        captProbs[0] = np.append(captProbs[0], (firstToFirstData.size + data_BothToBoth.size) / accumulations)
        captProbs[1] = np.append(captProbs[1], (secondToSecondData.size + data_BothToBoth.size) / accumulations)
        average_1to1 = np.append(average_1to1, np.average(firstToFirstData))
        average_1to2 = np.append(average_1to2, np.average(firstToSecondData))
        average_2to1 = np.append(average_2to1, np.average(secondToFirstData))
        average_2to2 = np.append(average_2to2, np.average(secondToSecondData))
        averageBothToBoth = np.append(averageBothToBoth, np.average(data_BothToBoth))
        averageBothToOne = np.append(averageBothToOne, np.average(data_BothToOne))
        error_1to1 = np.append(error_1to1, np.std(firstToFirstData) / np.sqrt(firstToFirstData.size))
        error_1to2 = np.append(error_1to2, np.std(firstToSecondData) / np.sqrt(firstToSecondData.size))
        error_2to1 = np.append(error_2to1, np.std(secondToFirstData) / np.sqrt(secondToFirstData.size))
        error_2to2 = np.append(error_2to2, np.std(secondToSecondData) / np.sqrt(secondToSecondData.size))
        error_BothToBoth = np.append(error_BothToBoth, np.std(data_BothToBoth) / np.sqrt(data_BothToBoth.size))
        error_BothToOne = np.append(error_BothToOne, np.std(data_BothToOne) / np.sqrt(data_BothToOne.size))
        ctsList = np.append(ctsList, np.average(firstToFirstData))
        lctsList = np.append(lctsList, np.sqrt(firstToFirstData.size))
    # condense data for export.
    survival_1 = average_1to1 + average_1to2
    error_survival_1 = np.sqrt(error_1to1 ** 2 + error_1to2 ** 2)
    survival_2 = average_2to1 + average_2to2
    error_survival_2 = np.sqrt(error_2to1 ** 2 + error_2to2 ** 2)

    # lots to return. Might be better way to do this.
    return (average_1to1, error_1to1, average_1to2, error_1to2, average_2to1, error_2to1, average_2to2, error_2to2,
            averageBothToBoth, error_BothToBoth, averageBothToOne, error_BothToOne, captProbs[0], captProbs[1],
            survival_1, error_survival_1, survival_2, error_survival_2)

