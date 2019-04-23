# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:50:04 2019

@author: Kain
"""

import matplotlib.pyplot as plt
import basic_functions_oh as bf
import numpy as np
import upsilon as up

class AIJFile: #TODO: rename to AIJFile???
    def __init__(self, filename, filter_type, target_star_dict, ref_star_dict, **kwargs):
        self.isClipped = False #initially, when made, assumes data is not clipped, unless otherwise specified in kwargs
        if "clipped" in kwargs: #optionally tell whether it is clipped
            self.isClipped = kwargs["clipped"]
        if "point_format" in kwargs:
            self.pltFmt = kwargs["point_format"]
        else:
            self.pltFmt = "b."
        self.raw_filename = filename #string containing the name of the file #TODO: this will be messed up if it's already clipped...
        self.filter_type = filter_type #string representing the type of filter
        self.target_star_dict = target_star_dict #dictionary with the name of the target star as the key, and the number of the target star as the value
        self.ref_star_dict = ref_star_dict #dictionary, with the name of the reference star as the keys, and the value as a tuple with (refStarNum, magInFilter)
        if filename.endswith(".txt") or filename.endswith(".dat"): #create a filename for the clipped file version (may not ever use)
            self.clip_filename = filename[:-4] + "_clipped" + filename[-4:]
        else:
            self.clip_filename = filename + "_clipped"
        if self.isClipped: #set initial filename as either clipped or unclipped version, depending on input
            self.filename = self.clip_filename
        else:
            self.filename = self.raw_filename 
        
    def clipFile(self, ref_star, sigma, **kwargs):
        if self.isClipped:
            print("This data file is already clipped!")
            return
        refStarNum, refStarMag = self.ref_star_dict[ref_star]
        refStarNums, refStarMags = self.getRefStarInfo(**kwargs)
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        cCounts = self.getColumn("Source-Sky_C" + str(refStarNum))
        bf.getClippedNightFile(cCounts, listRefCounts, refStarMags, sigma, self.filename, self.clip_filename) #TODO: edit this call-- may need to calculate more stuff first!
        self.isClipped = True
        self.filename = self.clip_filename
    
    def unclipFile(self): #reverts filename to the unclipped version
        if not self.isClipped:
            print("This data file is already not clipped!")
            return
        self.filename = self.raw_filename
        self.isClipped = False     
    
    def plotFile(self, xLabel, yLabel, **kwargs):
        showPlot = True
        xData = self.getColumn(xLabel)
        yData = self.getColumn(yLabel)
        if "showPlot" in kwargs:
            showPlot = kwargs["showPlot"]
        if showPlot:
            plt.plot(xData, yData, self.pltFmt)
            plt.xlabel(xLabel)
            plt.ylabel(yLabel)
            plt.show()
        return xData, yData
    
    def getColumn(self, label):
        data = bf.getColumn(self.filename, label)
        return data
    
    def getRefStarInfo(self, **kwargs): #returns a list of the reference star numbers and their magnitudes for the file; you can specify specific reference stars to use if you don't want to use all of them.
        refStarNums = []
        refStarMags = []
        kw_refNameList = []
        refStarNames = list(self.ref_star_dict.keys())
        useAllRefs = True
        if "use_ref_stars" in kwargs:
            useAllRefs = False
            kw_refNameList = kwargs["use_ref_stars"]
        for i in range(len(self.ref_star_dict)):
            sName = refStarNames[i]
            sNum, sMag = self.ref_star_dict[sName]
            if useAllRefs or (sName in kw_refNameList):
                refStarNums.append(sNum)
                refStarMags.append(sMag)
        return refStarNums, refStarMags
    
    def getMagData(self, targetName, targetCoords, ref_star, **kwargs):
        refStarNums, refStarMags = self.getRefStarInfo(**kwargs)
        targetNum = self.target_star_dict[targetName]
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        tCounts = self.getColumn(str("Source-Sky_T" + str(targetNum)))
        tMag = bf.calcMag(tCounts, listRefCounts, refStarMags, False) #TODO: is this EVER a comparison star, rather than a target? if so, CHANGE!!!
        JD = self.getColumn("JD_UTC")
        HJD = bf.JD_to_HJD(JD, targetCoords)
        ref_star_num, ref_star_mag = self.ref_star_dict[ref_star]
        cCounts = self.getColumn(str("Source-Sky_C" + str(ref_star_num)))
        dMag = bf.getCharacteristicError(cCounts, listRefCounts)
        return HJD, tMag, dMag
    
    def getFluxData(self, targetName, targetCoords, ref_star, **kwargs):
        refStarNums, refStarMags = self.getRefStarInfo(**kwargs)
        targetNum = self.target_star_dict[targetName]
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        tCounts = self.getColumn(str("Source-Sky_T" + str(targetNum)))
        tFlux = bf.calcFlux(tCounts, listRefCounts, False)
        JD = self.getColumn("JD_UTC")
        HJD = bf.JD_to_HJD(JD, targetCoords)
        ref_star_num, ref_star_mag = self.ref_star_dict[ref_star]
        cCounts = self.getColumn(str("Source-Sky_C" + str(ref_star_num)))
        dFlux = np.std(bf.calcFlux(cCounts, listRefCounts, True)) #TODO: you may need to subtract the reference star counts from each value in listRefCounts... or not? Depends on if you're fully treating it like a target.
        return HJD, tFlux, dFlux
    """
    def calcComparisonMag(): #TODO: implement later
        pass
    
    def calcComparisonFlux(self, ref_star, **kwargs): #TODO: implement later
        pass
    """    
        

class OtherDataFile: #TODO: decide about clipping for these types of files...!
    def __init__(self, filename, filter_type, column_positions, **kwargs):
        self.isClipped = False #initially, when made, assumes data is not clipped, unless otherwise specified in kwargs
        if "clipped" in kwargs: #optionally tell whether it is clipped
            self.isClipped = kwargs["clipped"]
        if "point_format" in kwargs:
            self.pltFmt = kwargs["point_format"]
        else:
            self.pltFmt = "b."
        if "has_header" in kwargs:
            self.has_header = kwargs["has_header"]
        else:
            self.has_header = 0 #assumes by default that the file contains labels/a header before the data
        self.raw_filename = filename #string containing the name of the file
        self.filter_type = filter_type #string representing the type of filter
        self.column_positions = column_positions
        if filename.endswith(".txt") or filename.endswith(".dat"): #create a filename for the clipped file version (may not ever use)
            self.clip_filename = filename[:-4] + "_clipped" + filename[-4:]
        else:
            self.clip_filename = filename + "_clipped"
        if self.isClipped: #set initial filename as either clipped or unclipped version, depending on input
            self.filename = self.clip_filename
        else:
            self.filename = self.raw_filename
            
    def getMagData(self, tName, tCoords, refStars, **kwargs):
        keepNAN = False
        if "keepNAN" in kwargs:
            keepNAN = kwargs["keepNAN"]
        data = np.genfromtxt(self.filename, skip_header=self.has_header)
        HJD1 = data[:,self.column_positions[0]]
        mag1 = data[:,self.column_positions[1]]
        dMag1 = data[:,self.column_positions[2]]
        HJD = []
        mag = []
        dMag = []
        #print(len(HJD1), len(mag1), len(dMag1))
        for i in range(len(HJD1)):
            if np.isnan(HJD1[i]) or np.isnan(mag1[i]) or np.isnan(dMag1[i]):
                print(str(self.filename) + " contains a non-numerical or infinite value in a column. Ignoring this point, unless otherwise specified...")
            else:
                HJD.append(HJD1[i])
                mag.append(mag1[i])
                dMag.append(dMag1[i])
        if keepNAN:
            return HJD1, mag1, dMag1
        else:
            return HJD, mag, dMag
    
    def getFluxData(self, tName, tCoords, refStars, **kwargs): #TODO: verify that there is no 'nan' data, as was done above!!!
        if len(self.column_positions) == 5:
            data = np.genfromtxt(self.filename, skip_header=self.has_header)
            HJD = data[:,self.column_positions[0]]
            flux = data[:,self.column_positions[3]]
            dFlux = data[:,self.column_positions[4]]
            return HJD, flux, dFlux
        else:
            print("No flux data in this file!")
            return
    
    def clipFile(self, sigma, iterations, **kwargs):
        tName = ""
        tCoords = ""
        refStars = ""
        if self.isClipped:
            print("This data file is already clipped!")
            return
        HJD1, mag1, dMag1 = self.getMagData(tName, tCoords, refStars, keepNAN=False, **kwargs)
        #print(HJD1)
        HJDa = np.asarray(HJD1)
        maga = np.asarray(mag1)
        dMaga = np.asarray(dMag1)
        HJD, mag, dMag = up.utils.sigma_clipping(HJDa, maga, dMaga, threshold=sigma, iteration=iterations)
        #print(HJD1)
        #print(HJD)
        origFile = open(self.raw_filename, "r")
        lines = origFile.readlines()
        origFile.close()
        #print(lines)
        data = np.genfromtxt(self.raw_filename, skip_header=self.has_header)
        numR, numC = data.shape
        retFile = open(self.clip_filename, "w")
        for i in range(self.has_header):
            retFile.write(lines[i])
        for j in range(numR):
            isInHJD1 = False
            for k in range(len(HJD)):
                if data[j, self.column_positions[0]] == HJD[k]:
                    isInHJD1 = True
            if isInHJD1:
                dataLine = ""
                for l in range(numC - 1):
                    dataLine = dataLine + str(data[j, l]) + "   "
                dataLine = dataLine + str(data[j, numC - 1])
                retFile.write(dataLine)
                retFile.write("\n")
        retFile.close()
        self.isClipped = True
        self.filename = self.clip_filename
        print("Clipped " + str(len(HJD1) - len(HJD)) + " points from " + str(self.raw_filename))
    
    def unclipFile(self): #reverts filename to the unclipped version
        if not self.isClipped:
            print("This data file is already not clipped!")
            return
        self.filename = self.raw_filename
        self.isClipped = False 

class Target:
    def __init__(self, name, coords, given_period, **kwargs):
        self.name = name
        self.coords = coords
        self.files = []
        self.refStars = []
        self.periods = {"given_period": given_period}
        self.period = self.periods["given_period"]
        if "past_hjd" in kwargs: #TODO: decide if this is how you want to handle the past HJD stuff...
            self.hjd = kwargs["past_hjd"]
        else:
            self.hjd = ""
        if "phase_hjd" in kwargs:
            self.phase_hjd = kwargs["phase_hjd"]
        else:
            self.phase_hjd = 0
    
    def getMagData(self, filter_type, **kwargs): #returns HJD, mag, dmag, and plot color lists. Maybe I should make this a thing in the AIJFile class also so I can just call that from this function for each file?
        return_1D_lists = False #by default, returns lists of lists (less computation)
        if "return_1D_lists" in kwargs:
            return_1D_lists = kwargs["return_1D_lists"]
        use_files, use_ref_stars = self.getFilesInFilter(filter_type, **kwargs)
        listPltColor = []
        listHJD = []
        listMag = []
        listDMag = []
        for i in range(len(use_files)):
            currentFile = use_files[i]
            #print(currentFile.filename)
            fileHJD, fileMag, fileDMag = currentFile.getMagData(self.name, self.coords, use_ref_stars[i], **kwargs)
            listPltColor.append(currentFile.pltFmt) #Returns a list of lists for HJD and Mag, where each list is the data from a file! Then the other 2 are just lists of values which correspond to each file. Or, if ASAS-SN, dmagvalue for that file will also be a list.
            listHJD.append(fileHJD)
            listMag.append(fileMag)
            listDMag.append(fileDMag)
        if return_1D_lists:
            HJD = []
            mag = []
            dMag = []
            pltColor = []
            isAIJ = []
            for i in range(len(listHJD)):
                fileHJD = listHJD[i] #always a list
                fileMag = listMag[i] #always a list
                fileDMag = listDMag[i] #sometimes a list
                fileColor = listPltColor[i] #never a list
                
                if isinstance(use_files[i], AIJFile):
                    fromAIJFile = True
                else:
                    fromAIJFile = False
                
                for j in range(len(fileHJD)):
                    HJD.append(fileHJD[j])
                    mag.append(fileMag[j])
                    if isinstance(fileDMag, float):
                        dMag.append(fileDMag)
                    else:
                        dMag.append(fileDMag[j])
                    pltColor.append(fileColor)
                    isAIJ.append(fromAIJFile)
            #print(mag)
            #return HJD, mag, dMag, pltColor #returns 4 lists, all equal in length
            return HJD, mag, dMag, pltColor, isAIJ
        else:
            #print(listMag)
            return listHJD, listMag, listDMag, listPltColor #returns lists of lists
    
    def getFluxData(self, filter_type, **kwargs): #returns HJD, flux, and plot color lists (and maybe flux error...?)
        return_1D_lists = False #by default, returns lists of lists (less computation)
        if "return_1D_lists" in kwargs:
            return_1D_lists = kwargs["return_1D_lists"]
        use_files, use_ref_stars = self.getFilesInFilter(filter_type, **kwargs)
        listPltColor = []
        listHJD = []
        listFlux = []
        listDFlux = []
        for i in range(len(use_files)):
            currentFile = use_files[i]
            fileHJD, fileFlux, fileDFlux = currentFile.getFluxData(self.name, self.coords, use_ref_stars[i], **kwargs)
            listPltColor.append(currentFile.pltFmt) #Returns a list of lists for HJD and Mag, where each list is the data from a file! Then the other 2 are just lists of values which correspond to each file. Or, if ASAS-SN, dmagvalue for that file will also be a list.
            listHJD.append(fileHJD)
            listFlux.append(fileFlux)
            listDFlux.append(fileDFlux)
        if return_1D_lists:
            HJD = []
            flux = []
            dFlux = []
            pltColor = []
            for i in range(len(listHJD)):
                fileHJD = listHJD[i] #always a list
                fileFlux = listFlux[i] #always a list
                fileDFlux = listFlux[i] #sometimes a list
                fileColor = listPltColor[i] #never a list
                for j in range(len(fileHJD)):
                    HJD.append(fileHJD[j])
                    flux.append(fileFlux[j])
                    if isinstance(fileDFlux, float):
                        dFlux.append(fileDFlux)
                    else:
                        dFlux.append(fileDFlux[j])
                    pltColor.append(fileColor)
            return HJD, flux, dFlux, pltColor #returns 4 lists, all equal in length
        else:
            return listHJD, listFlux, listDFlux, listPltColor #returns lists of lists
    
    def addData(self, file, comparison_star):
        self.files.append(file)
        self.refStars.append(comparison_star)
        
    def getFilesInFilter(self, filter_type, **kwargs): #returns lists of the files and their reference stars which match a filter type given by user (and are in list of files given by user, if speficied)
        matchingFiles = []
        matchingRefStars = []
        if "use_files" in kwargs:
            use_files = kwargs["use_files"]
        else:
            use_files = self.files
        include_non_AIJ = True
        if "include_non_AIJ" in kwargs:
            include_non_AIJ = kwargs["include_non_AIJ"]
        numFiles = len(self.files)
        for i in range(numFiles):
            currentFile = self.files[i]
            currentStar = self.refStars[i]
            if currentFile.filter_type == filter_type and currentFile in use_files:
                if include_non_AIJ or isinstance(currentFile, AIJFile):
                    matchingFiles.append(currentFile)
                    matchingRefStars.append(currentStar)
        return matchingFiles, matchingRefStars #TODO: if you use this on a set of files which includes an ASAS-SN file, what is returned for matchingRefStars?
    
    def makeLightCurve(self, filter_type, **kwargs):
        yIsFlux = False #Default is to plot in magnitudes on y-axis, but can change this so it plots in fluxes instead
        xIsHJD = False
        if "runShowPlot" in kwargs:
            runShowPlot = kwargs["runShowPlot"]
        else:
            runShowPlot = True
        if "plotFlux" in kwargs:
            yIsFlux = kwargs["plotFlux"]
        if "plotHJD" in kwargs:
            xIsHJD = kwargs["plotHJD"]
        if yIsFlux:
            HJDs, yData, dYData, pltFmts = self.getFluxData(filter_type, return_1D_lists=False, **kwargs)
            ylabel = "Flux in " + str(filter_type) + " Filter"
            yData_maxs = []
            for i in range(len(yData)):
                yData_maxs.append(max(yData[i]))
            errYpos = max(yData_maxs)
        else:
            HJDs, yData, dYData, pltFmts = self.getMagData(filter_type, return_1D_lists=False, **kwargs)
            ylabel = str(filter_type) + " Mag"
            yData_mins = []
            for i in range(len(yData)):
                yData_mins.append(min(yData[i]))
            errYpos = min(yData_mins)
        if xIsHJD:
            xData = HJDs
            xlabel="HJD"
            HJD_maxs = []
            HJD_mins = []
            for i in range(len(HJDs)):
                HJD_maxs.append(max(HJDs[i]))
                HJD_mins.append(min(HJDs[i]))
            errXpos = min(HJD_mins)
            errXdelta = (max(HJD_maxs) - min(HJD_mins)) / 20
        else:
            xlabel = "Phase"
            errXpos = 0
            errXdelta = .05
            xData = []
            for i in range(len(HJDs)):
                xData.append(bf.phaseData(HJDs[i], self.period, hjd_offset=self.phase_hjd))
        if "bar_position" in kwargs:
            errYpos = kwargs["bar_position"]
        for i in range(len(HJDs)):
            #print(dYData)
            if isinstance(dYData[i], float):
                plt.plot(xData[i], yData[i], pltFmts[i])
                plt.errorbar(errXpos, errYpos, yerr=dYData[i], linestyle="None", fmt=pltFmts[i])
                errXpos = errXpos + errXdelta
            else:
                plt.errorbar(xData[i], yData[i], yerr=dYData[i], linestyle="None", fmt=pltFmts[i])
        if yIsFlux == False:
            plt.gca().invert_yaxis()
        if "set_title" in kwargs:
            plt.title(str(kwargs["set_title"]))
        else:
            plt.title(str(self.name))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if runShowPlot:
            plt.show()
        
    def getUpsilonType(self, filter_type, **kwargs): # TODO: add in kwargs for the upsilon sigma clipping routine, if desired?
        HJD, Mags, dMags, pltColors, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        rf_model = up.load_rf_model()
        e_features = up.ExtractFeatures(HJD, Mags, dMags)
        e_features.run()
        features = e_features.get_features()
        label, probability, flag = up.predict(rf_model, features)
        print(label, probability, flag)
        feats = list(features.items())
        for i in range(len(feats)):
            title, value = feats[i]
            print(str(title) + ": " + str(value))
        p_label, up_period = feats[6]
        self.periods["upsilon_period"] = up_period
    
    def getFourierFit(self, filter_type, **kwargs):
        if "lineFmt" in kwargs:
            lineColor = kwargs["lineFmt"]
        else:
            lineColor = "k-"
        if "showPlot" in kwargs:
            showPlot = kwargs["showPlot"]
        else:
            showPlot = True
        HJD, Mag, dMag, pltColor, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        phase = bf.phaseData(HJD, self.period, hjd_offset=self.phase_hjd)
        labs = ("Phase", str(filter_type) + " Mag")
        popt, pcov = bf.fourierFit(phase, Mag, dMag, plot_title=self.name, plot_labels=labs, **kwargs)
        if showPlot:
            self.makeLightCurve(filter_type, runShowPlot=False, **kwargs)
            xAxis = np.linspace(0, 1, 100)
            plt.plot(xAxis, bf.fourier(xAxis, *popt), lineColor)
            plt.show()
        return popt, pcov
    
    def calcPeriod(self, filter_type, **kwargs):
        HJD, Mag, dMag, pltColor, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        calcP = bf.getPeriod(HJD, Mag, dMag, **kwargs)
        self.periods["calculated_period"] = calcP
        return calcP
    
    def getColorFit(self, filter1, filter2, **kwargs): # TODO: add kwargs to make graph not print, so that when this is called inside another function it doesn't make a graph too???
        if "forceN1" in kwargs:
            popt_c1, pcov_c1 = self.getFourierFit(filter1, forceN=kwargs["forceN1"], **kwargs)
        else:
            popt_c1, pcov_c1 = self.getFourierFit(filter1, **kwargs)
        if "forceN2" in kwargs:
            popt_c2, pcov_c2 = self.getFourierFit(filter2, forceN=kwargs["forceN2"], **kwargs)
        else:
            popt_c2, pcov_c2 = self.getFourierFit(filter2, **kwargs)
        xPhase = np.linspace(0, 1, 1000)
        mag_c1 = bf.fourier(xPhase, *popt_c1)
        mag_c2 = bf.fourier(xPhase, *popt_c2)
        color = mag_c1 - mag_c2
        y_label = filter1 + " - " + filter2
        plt.plot(xPhase, color, "b.")
        plt.xlabel("Phase")
        plt.ylabel(y_label)
        plt.title(self.name)
        plt.show()
        return xPhase, color
    
    def getColorTemperature(self, filter1, filter2, tempModelFit, **kwargs):
        a, b, c = tempModelFit
        xPhase, color = self.getColorFit(filter1, filter2, **kwargs)
        temp = a + (b*color) + (c*(color**2))
        plt.xlabel("Phase")
        plt.ylabel("Color Temperature")
        plt.plot(xPhase, temp, "b.")
        plt.title(self.name)
        plt.show()
        return xPhase, temp

    def plotFiles(self, xLabel, yLabel, **kwargs): # TODO: take in multiple labels (one for each file) if desired, and plot each on same axes ---> issue since the order of the files is unclear to the user... maybe a dictionary?
        if "use_files" in kwargs:
            plot_files = kwargs["use_files"]
        else:
            plot_files = self.files
        for i in range(len(plot_files)):
            file_x, file_y = plot_files[i].plotFile(xLabel, yLabel, showPlot=False)
            plt.plot(file_x, file_y, plot_files[i].pltFmt)
        plt.xlabel(xLabel)
        plt.ylabel(yLabel)
        plt.show()
    
    def setPeriod(self, period_label):
        self.period = self.periods[period_label]
    
    def multiFourierFit(self, filter_type, period_list, topN, **kwargs):
        if "showPlot" in kwargs:
            showPlot = kwargs["showPlot"]
        else:
            showPlot = True
        HJD, Mag, dMag, pltColor, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        popt, pcov = bf.multiFourierFit(HJD, Mag, dMag, period_list, forceN=topN, **kwargs)
        if showPlot:
            for i in range(len(HJD)):
                plt.plot(HJD[i], Mag[i], pltColor[i])
            xAxis = np.linspace(min(HJD), max(HJD), 1000)
            plt.plot(xAxis, bf.make_multiFourier(period_list, topN)(xAxis, *popt), "k-")
            plt.gca().invert_yaxis()
            plt.xlabel("HJD")
            plt.ylabel(str(filter_type) + " Mag")
            plt.title(self.name)
            if "xLims" in kwargs:
                plt.xlim(kwargs["xLims"])
            plt.show()
        return popt, pcov
    
    def addPeriod(self, label, period):
        self.periods[label] = period
    """
    def timeOnTarget(self, time_threshold, **kwargs): #TODO: implement! also, this WILL return ALL files added to target, soo... don't include multiple versions of the same file!
        allFiles = self.files
        numFiles = len(allFiles)
        timeSum = 0
        numImages = 0
        for i in range(numFiles): #TODO: implement a checker to skip any files that aren't AIJ files...!
            JDs = allFiles[i].getColumn("JD_UTC")
            numImages = numImages + len(JDs)
            JD_sets = []
            JD_sets.append([JDs[0]])
            listNum = 0
            j = 0
            while j < len(JDs) - 1:
                if JDs[j+1] - JDs[j] <= time_threshold:
                    JD_sets[listNum].append(JDs[j+1])
                    j = j + 1
                else:
                    JD_sets.append([JDs[j+1]])
                    listNum = listNum + 1
                    j = j + 1
            file_timeSum = 0
            for k in range(len(JD_sets)):
                file_timeSum = file_timeSum + (max(JD_sets[k]) - min(JD_sets[k]))
            timeSum = timeSum + file_timeSum
        numHours = int(timeSum * 24)
        numMins = int(((timeSum * 24) - numHours)*60)
        numSecs = ((((timeSum * 24) - numHours)*60) - numMins)*60
        totalHours = timeSum * 24
        print("Target: " + str(self.name))
        print("Number of Images Taken: " + str(numImages))
        print("Time on Target: " + str(numHours) + "h:" + str(numMins) + "m:" + str(numSecs) + "s")
        print("Total Hours on Target: " + str(totalHours))
        print("\n")
        """