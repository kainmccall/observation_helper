# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:50:37 2019

@author: Kain
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargleFast
from PyAstronomy import pyasl

def calcMag(tCounts, listRefCounts, refMags, isInC): #TODO: MUST back-edit all functions which call this to include isInC! also, doesn't subtract mag from bright if isInC... issue ifthis is ever called with this as true and doesn't just pass it on to calcFlux
    numPoints = len(tCounts)
    numRefs = len(listRefCounts)
    relFlux = calcFlux(tCounts, listRefCounts, isInC)
    bright = 0
    for i in range(numRefs):
        bright = bright + 10**(-.4 * refMags[i])
    m_all = -2.5 * np.log10(bright)
    tMag = np.zeros(numPoints)
    for i in range(numPoints):
        tMag[i] = m_all - (2.5 * np.log10(relFlux[i]))
    return tMag

def calcFlux(tCounts, listRefCounts, isInC):
    numPoints = len(tCounts)
    numRefs = len(listRefCounts)
    tot_C_counts = np.zeros(numPoints)
    for i in range(numPoints):
        for j in range(numRefs):
            tot_C_counts[i] = tot_C_counts[i] + listRefCounts[j][i]
    if isInC: #TODO: need to make sure that if a list is specified of reference stars, and this reference star isn't in that list, that it doesn't still subtract these counts!!!
        for i in range(numPoints):
            tot_C_counts[i] = tot_C_counts[i] - tCounts[i]
    return tCounts / tot_C_counts

def getCharacteristicError(tCounts, listRefCounts):
    relFlux = calcFlux(tCounts, listRefCounts, True)
    dF = np.std(relFlux)
    F = np.average(relFlux)
    dM = (2.5 / np.log(10)) * (dF / F)
    return dM

def JD_to_HJD(JD, tCoords): #Converts a regular JD to an HJD (can be a float, or a list/array)
    ra, dec = tCoords
    if isinstance(JD, float):
        if JD > 2400000.5:
            JD1 = JD - 2400000.5
            return pyasl.helio_jd(JD1, ra, dec) + 2400000.5
        else:
            JD1 = JD
            return pyasl.helio_jd(JD1, ra, dec)
    else:
        HJD = np.zeros(len(JD))
        for i in range(len(JD)):
            if JD[i] > 2400000.5:
                mJD = JD[i] - 2400000.5
                HJD[i] = pyasl.helio_jd(mJD, ra, dec) + 2400000.5
            else:
                HJD[i] = pyasl.helio_jd(JD[i], ra, dec)
        return HJD

def getColumn(filename, label): #retrieves data (not including label) from "filename" in column labeled "label"
    data = np.genfromtxt(filename, skip_header=1)
    #print(filename)
    #print(data.shape)
    check = len(data.shape)
    labels = np.genfromtxt(filename, max_rows=1, dtype=str)
    numLabels = len(labels)
    index = 0
    for i in range(numLabels):
        realIndex = i + 1 #needed because there is a space in the output table in the first row!
        if str(labels[i])==str(label):
            index = realIndex
    #print(data[:,index])
    #print(label + ": " + str(index))
    if check == 1:
        print("Possible empty file! If this throws an error, it's because you gave an empty file somehow! Oops!")
        print("\n")
        return np.full(1, data[index])
    return data[:,index]

def phaseData(JD, period, **kwargs): #transforms a JD (or HJD) into a phase over a period. Assumes JD = 0 for the "initial point" for simplicity.
    if "hjd_offset" in kwargs:
        offset = float(kwargs["hjd_offset"])
        #print(offset)
        #print(JD)
        JD1 = np.zeros(len(JD))
        for i in range(len(JD)):
            JD1[i] = JD[i] - offset
    else:
        JD1 = JD
    phase = np.zeros(len(JD1))
    for i in range(len(JD1)):
        phase[i] = (JD1[i] % period) / period
    return phase

def getClippedNightFile(cCounts, listRefCounts, refStarMags, sigma, inFile, outFile): #creates a new night file with only clipped data using a sigma clipping proceedure.
    shouldClip = getClippingList(cCounts, listRefCounts, refStarMags, sigma)
    labels = np.genfromtxt(inFile, max_rows=1, dtype=str)
    numLabels = len(labels)
    labString = str(labels[0])
    for i in range(numLabels - 1):
        labString = labString + "   " + str(labels[i + 1])
    retFile = open(outFile, "w")
    retFile.write(labString)
    retFile.close()
    data = np.genfromtxt(inFile, skip_header=1, dtype=str)
    numR, numC = data.shape
    if numC != numLabels + 1:
        print("ERROR! Number of columns in files are not the same!")
        return
    retFile = open(outFile, "a")
    for j in range(numR):
        rowString = ""
        for k in range(numC - 1):
            rowString = rowString + str(data[j, k]) + "   "
        rowString = rowString + str(data[j, numC - 1])
        if shouldClip[j] == False:
            retFile.write("\n")
            retFile.write(rowString)
    retFile.close()
    return outFile

def getClippingList(cCounts, listRefCounts, refStarMags, sigma):
    cMags = calcMag(cCounts, listRefCounts, refStarMags, True)
    shouldClip = np.zeros(len(cMags), dtype=bool)
    avgM = np.nanmean(cMags)
    minM = avgM - (sigma * np.nanstd(cMags))
    maxM = avgM + (sigma * np.nanstd(cMags))
    for i in range(len(cMags)):
        if cMags[i] < minM or cMags[i] > maxM or np.isnan(cMags[i]):
            shouldClip[i] = True
        else:
            shouldClip[i]=False
    return shouldClip

def fourier(x, *a):
    f = a[0]
    for n in range(int((len(a) - 1)/2)):
        a_i = n + 1
        p_i = int(a_i + ((len(a) - 1)/2))
        f = f + (a[a_i] * np.cos((2 * a_i * np.pi * x) + a[p_i]))
    return f

def fourierFit(phases, allMags, alldMags, **kwargs):
    knownN = False
    if "forceN" in kwargs:
        topN = kwargs["forceN"]
        knownN = True
    else:
        topN = 1
    notDone = True
    while notDone:
        lenA = (2 * topN) + 1
        popt, pcov = curve_fit(fourier, phases, allMags, [1.0] * lenA, sigma=alldMags)
        resids = allMags - fourier(phases, *popt) 
        avResid = np.mean(resids)
        numer = 0
        denom = 0
        def takeFirst(elem):
            return elem[0]
        resids_xy = []
        for i in range(len(resids)):
            resids_xy.append((float(phases[i]), float(resids[i])))
        resids_xy.sort(key=takeFirst)
        phases1 = np.zeros(len(phases))
        resids1 = np.zeros(len(phases))
        for i in range(len(resids)):
            phases1[i], resids1[i] = resids_xy[i]
        for i in range(len(resids) - 1):
            numer = numer + ((resids1[i] - avResid) * (resids1[i + 1] - avResid))
            denom = denom + ((resids1[i] - avResid)**2)
        p = numer / denom
        p_c = (2 * (len(phases) - 1))**(-1/2)
        if p < p_c or knownN:
            notDone = False
        else:
            topN = topN + 1
    print("N: " + str(topN))
    print("p: " + str(p))
    print("p_c: " + str(p_c))
    if topN > 1 and use_f2 == False:
        print("R_21: " + str(popt[2] / popt[1]))
    return popt, pcov

def getPeriod(HJD, TMag, TdMag, **kwargs):
    p_min = .01
    p_max = 2
    giveP = False
    makeGraph = False
    if "p_range" in kwargs:
        p_min, p_max = kwargs["p_range"]
    if "graphP" in kwargs:
        makeGraph = kwargs["graphP"]
    if "printP" in kwargs:
        giveP = kwargs["printP"]
    lsf = LombScargleFast()
    lsf.optimizer.period_range = (p_min, p_max)
    lsf.fit(HJD, TMag, TdMag)
    period = lsf.best_period
    if giveP:
        print("\n")
        print("Period: " + str(period))
        print("\n")
    if makeGraph:
        periods = np.linspace(p_min, p_max, 1000)
        scores = lsf.score(periods)
        plt.plot(periods, scores)
        plt.xlabel("Period")
        plt.ylabel("Probability")
        plt.show()
    return period

def fluxFromMag(tMags, refStarMags): #takes in an array or list of magnitudes and turns them into relative fluxes
    if isinstance(refStarMags, float):
        bright = 10**(-.4 * refStarMags)
    else:
        numRefStars = len(refStarMags)
        bright = 0  
        for i in range(numRefStars):
            bright = bright + 10**(-.4 * refStarMags[i])
    m_all = -2.5 * np.log10(bright)
    tFlux = np.zeros(len(tMags))
    for i in range(len(tMags)):
        tFlux[i] = 10**((m_all - tMags[i]) / 2.5)
    return tFlux


def vertexParabolaFunc(x, *p): #Point-slope form
    a = p[0]
    h = p[1]
    k = p[2]
    return a*((x-h)**2) + k

def parabolaFit(xdata, ydata, dy, guessBeta):
    popt, pcov = curve_fit(vertexParabolaFunc, xdata, ydata, guessBeta, sigma=dy)
    err = np.sqrt(np.diag(pcov))
    plt.errorbar(xdata, ydata, yerr=dy, linestyle="None", fmt="b.")
    x = np.linspace(min(xdata), max(xdata), 1000)
    y = vertexParabolaFunc(x, *popt)
    plt.plot(x, y)
    plt.gca().invert_yaxis()
    plt.show()
    return popt, err