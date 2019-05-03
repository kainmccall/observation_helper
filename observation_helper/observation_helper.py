# -*- coding: utf-8 -*-
#Created on Thu Apr 25 17:53:18 2019
#@author: Kain

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from gatspy.periodic import LombScargleFast
from PyAstronomy import pyasl

#from utils import calcMag, calcFlux, getCharacteristicError, JD_to_HJD, getCol, phaseData, getClippedNightFile, getClippingList, fourier, fourierFit, getPeriod, fluxFromMag, vertexParabolaFunc, parabolaFit

class AIJFile:
    """This is a class of objects which represent the output text files from AstroImageJ. Each object must be created with properties encoded to tell the program what filter the data were taken in, which stars belong to which aperture, etc.
    
    :param filename: name of output text file
    :type filename: str
    :param filter_type: filter in which data were taken
    :type filter_type: str
    :param target_star_dict: dictionary containing the name of the target star in each target aperture (designated with a 'T' in the label by AstroImageJ) as keys, and the integer aperture number (i.e., 1 for 'T1', 5 for 'T5') as the corresponding value
    :type target_star_dict: dict{str:int}
    :param ref_star_dict: dictionary, similar to above, containing name of each reference star as key, and a tuple as value. Tuple should be of the format (aperture number, magnitude).
    :type ref_star_dict: dict{str:(int, float)}
    :param point_format: matplotlib color/style for plotting anything from this file (e.g., "b." or "k-"). Default is "b."
    :type point_format: str
    """
    
    def __init__(self, filename, filter_type, target_star_dict, ref_star_dict, point_format=None): #removed kwargs since I'm not using isClipped anymore as an option!
        """AIJFile objects represent the output text files from AstroImageJ. Each object must be created with properties encoded to tell the program what filter the data were taken in, which stars belong to which aperture, etc.
        
        
        :param filename: name of output text file
        :type filename: str
        :param filter_type: filter in which data were taken
        :type filter_type: str
        :param target_star_dict: dictionary containing the name of the target star in each target aperture (designated with a 'T' in the label by AstroImageJ) as keys, and the integer aperture number (i.e., 1 for 'T1', 5 for 'T5') as the corresponding value
        :type target_star_dict: dict{str:int}
        """
        self.isClipped = False #initially, when made, assumes data is not clipped, unless otherwise specified in kwargs
        #if "clipped" in kwargs: #optionally tell whether it is clipped #TODO: remove the ability to mark a file as already clipped? It's rather annoying...
            #self.isClipped = kwargs["clipped"]
        if "point_format" is None:
            self.pltFmt = "b."
        else:
            self.pltFmt = point_format
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
        """
        Clips data points from file if the magnitude of ref_star is more that sigma standard deviations away from its mean magnitude. Actually creates a new .txt file in the same directory with the same name, with _clipped appended to the end and uses this new file for all following calculations, etc.
        
        :param ref_star: name of the reference star whose magnitude to calculate and use to determine which data points should be clipped. MUST match name in ref_star_dict.
        :type ref_star: str
        :param sigma: limit on number of standard deviations away a value can be without being eliminated
        :type sigma: int or float
        
        Kwargs:
        :param use_ref_stars: a list of reference star names (see ``ref_star_dict``) to use for calculation of magnitude, if use of only some reference stars is desired
        :type use_ref_stars: list
        """
        if self.isClipped:
            print("This data file is already clipped!")
            return
        refStarNum, refStarMag = self.ref_star_dict[ref_star]
        refStarNums, refStarMags = self._getRefStarInfo(**kwargs)
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        cCounts = self.getColumn("Source-Sky_C" + str(refStarNum))
        getClippedNightFile(cCounts, listRefCounts, refStarMags, sigma, self.filename, self.clip_filename) #TODO: edit this call-- may need to calculate more stuff first!
        self.isClipped = True
        self.filename = self.clip_filename
    
    def unclipFile(self): #reverts filename to the unclipped version
        """
        Un-clips the data file. In reality, goes back to using old, untouched file for calculations, etc.
        """
        if not self.isClipped:
            print("This data file is already not clipped!")
            return
        self.filename = self.raw_filename
        self.isClipped = False     
    
    def plotFile(self, xLabel, yLabel, showPlot=None):
        """
        Creates a plot of column ``yLabel`` vs. column ``xLabel``, where ``xLabel`` and ``yLabel`` are the labels of each column in the AstroImageJ output file.
        
        :param xLabel: label of column in AstroImageJ output file whose data will be considered the x-values of plot
        :type xLabel: str
        :param yLabel: label of column in AstroImageJ output file whose data will be considered the y-values of plot
        :type yLabel: str
        :param showPlot: determines whether to print plot to console (default is True)
        :type showPlot: bool
        :returns: column in file with label ``xLabel``, column in file with label ``yLabel``
        :rtype: ndarray, ndarray
        """
        if showPlot is None:
            showPlot = True
        xData = self.getColumn(xLabel)
        yData = self.getColumn(yLabel)
        if showPlot:
            plt.plot(xData, yData, self.pltFmt)
            plt.xlabel(xLabel)
            plt.ylabel(yLabel)
            plt.show()
        return xData, yData
    
    def getColumn(self, label):
        """
        Returns all values from the column in the AstroImageJ file with the label label.
        
        :param label: label of the column whose values are returned
        :type label: str
        :returns: column in file with label ``label``
        :rtype: ndarray
        """
        data = getCol(self.filename, label)
        return data
    
    def _getRefStarInfo(self, **kwargs): #returns a list of the reference star numbers and their magnitudes for the file; you can specify specific reference stars to use if you don't want to use all of them.
        """
        Returns the aperture number and magnitudes for reference stars in self.ref_star_dict
        
        :returns: refStarNums, refStarMags -- list of aperture numbers for reference stars, list of magnitudes for reference stars
        :rtype: list, list
        
        kwargs here
        """
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
    
    def _getMagData(self, targetName, targetCoords, ref_star, **kwargs):
        """
        AIJFile getMagData docstring here
        """
        refStarNums, refStarMags = self._getRefStarInfo(**kwargs)
        targetNum = self.target_star_dict[targetName]
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        tCounts = self.getColumn(str("Source-Sky_T" + str(targetNum)))
        tMag = calcMag(tCounts, listRefCounts, refStarMags, False) #TODO: is this EVER a comparison star, rather than a target? if so, CHANGE!!!
        JD = self.getColumn("JD_UTC")
        HJD = JD_to_HJD(JD, targetCoords)
        ref_star_num, ref_star_mag = self.ref_star_dict[ref_star]
        cCounts = self.getColumn(str("Source-Sky_C" + str(ref_star_num)))
        dMag = getCharacteristicError(cCounts, listRefCounts)
        return HJD, tMag, dMag
    
    def _getFluxData(self, targetName, targetCoords, ref_star, **kwargs):
        """
        AIJFile getFluxData docstring here
        """
        refStarNums, refStarMags = self._getRefStarInfo(**kwargs)
        targetNum = self.target_star_dict[targetName]
        listRefCounts = []
        for i in range(len(refStarMags)):
            refCounts = self.getColumn(str("Source-Sky_C" + str(refStarNums[i])))
            listRefCounts.append(refCounts)
        tCounts = self.getColumn(str("Source-Sky_T" + str(targetNum)))
        tFlux = calcFlux(tCounts, listRefCounts, False)
        JD = self.getColumn("JD_UTC")
        HJD = JD_to_HJD(JD, targetCoords)
        ref_star_num, ref_star_mag = self.ref_star_dict[ref_star]
        cCounts = self.getColumn(str("Source-Sky_C" + str(ref_star_num)))
        dFlux = np.std(calcFlux(cCounts, listRefCounts, True)) #TODO: you may need to subtract the reference star counts from each value in listRefCounts... or not? Depends on if you're fully treating it like a target.
        return HJD, tFlux, dFlux
    #TODO: calcComparisonMag, calcComparisonFlux methods should be made eventually to plot individual reference stars...   
        

class OtherDataFile:
    """
    Class of objects which represent data files not from AstroImageJ. These could be, for instance, .txt files of data from ASAS-SN or CSS, etc. Text files must contain at least columns for HJD, magnitude, and error in magnitude.
    
    :param filename: name of the file which contains the data
    :type filename: str
    :param filter_type: filter in which data were taken
    :type filter_type: str
    :param column_positions: list of zero-indexed columns corresponding to HJD, magnitude, and error in magnitude, in that order. (e.g. [0, 2, 3] would mean that HJD data is in first column, magnitude is in the third column, and error in magnitude is in 4th column)
    :type column_position: list
    :param point_format: matplotlib color/style for plotting this file (e.g., "b." or "k-"). Default is "b."
    :type point_format: str
    :param has_header: number of rows to ignore in beginning of file (in case of headers, column labels, etc.). Default is 0
    :type has_header: int
    """
    
    def __init__(self, filename, filter_type, column_positions, point_format=None, has_header=None):
        """
        OtherDataFile init docstring here (won't be shown unless specified somehow...)
        """
        if point_format is None:
            self.pltFmt = "b."
        else:
            self.pltFmt = point_format
        if has_header is None:
            self.has_header = 0
        else:
            self.has_header = has_header #assumes by default that the file contains labels/a header before the data
        self.filename = filename #string containing the name of the file
        self.filter_type = filter_type
        self.column_positions = column_positions
            
    def _getMagData(self, tName, tCoords, refStars, **kwargs): #TODO: edit so people don't need to include irrelevant stuff?
        """
        otherDataFile getMagData docstring here
        """
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
    
    def _getFluxData(self, tName, tCoords, refStars, **kwargs): #TODO: verify that there is no 'nan' data, as was done above!!!
        """
        otherDataFile getFluxData docstring here
        """
        if len(self.column_positions) == 5:
            data = np.genfromtxt(self.filename, skip_header=self.has_header)
            HJD = data[:,self.column_positions[0]]
            flux = data[:,self.column_positions[3]]
            dFlux = data[:,self.column_positions[4]]
            return HJD, flux, dFlux
        else:
            print("No flux data in this file!")
            return

class Target:
    """
    Class of objects which represent target stars. Each has properties like a name, coordinates, a period, etc.
    
    :param name: name of the target
    :type name: str
    :param coords: tuple containing the J2000 RA and Dec of the target (RA, Dec), in decimal format (NOT sexagesimal)
    :type coords: tuple (float, float)
    :param given_period: period of target to use as default, initially (can be changed-- see :func:`calcPeriod`, :func:`addPeriod`, and :func:`setPeriod`)
    :type given_period: float
    :param phase0_hjd: when phasing, this will be the HJD with the initial phase, or phase=0 (default is 0)
    :type phase0_hjd: float
    """
    
    def __init__(self, name, coords, given_period, phase0_hjd=None):
        """
        class Target init docstring here (won't be shown unless specified somehow...)
        """
        self.name = name
        self.coords = coords
        self.files = []
        self.refStars = []
        self.periods = {"given_period": given_period}
        self.period = self.periods["given_period"]
        if phase0_hjd is None:
            self.phase_hjd = 0
        else:
            self.phase_hjd = phase0_hjd
    
    def getMagData(self, filter_type, return_1D_lists=None, **kwargs): #returns HJD, mag, dmag, and plot color lists. Maybe I should make this a thing in the AIJFile class also so I can just call that from this function for each file?
        """
        Calculates magnitude of Target object for all data points added and returns HJD, magnitude, error in magnitude, and plot color
        
        :param filter_type: Filter in which data to be returned was taken
        :type filter_type: str
        :param return_1D_lists: determines whether to return lists of equal length with the data from all files combined into the same list, or whether to return lists of lists, where each returned list contains lists of values which correspond to a particular file (default is False)
        :type return_1D_lists: bool
        :returns: HJD, mag, dMag, pltColor
        :rtype: list, list, list, list
        """
        if return_1D_lists is None:
            return_1D_lists = False
        use_files, use_ref_stars = self._getFilesInFilter(filter_type, **kwargs)
        listPltColor = []
        listHJD = []
        listMag = []
        listDMag = []
        for i in range(len(use_files)):
            currentFile = use_files[i]
            #print(currentFile.filename)
            fileHJD, fileMag, fileDMag = currentFile._getMagData(self.name, self.coords, use_ref_stars[i], **kwargs)
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
            return HJD, mag, dMag, pltColor, isAIJ #TODO: am I using isAIJ for something? If not, get rid of it...!
        else:
            #print(listMag)
            return listHJD, listMag, listDMag, listPltColor #returns lists of lists
    
    def getFluxData(self, filter_type, return_1D_lists=None, **kwargs): #returns HJD, flux, and plot color lists (and maybe flux error...?)
        """
        Calculates the relative flux of the target star
        
        :param filter_type: Filter in which data to be returned was taken
        :type filter_type: str
        :param return_1D_lists: determines whether to return lists of equal length with the data from all files combined into the same list, or whether to return lists of lists, where each returned list contains lists of values which correspond to a particular file (default is False)
        :type return_1D_lists: bool
        :returns: HJD, flux, dFlux, pltColor
        :rtype: list, list, list, list
        """
        if return_1D_lists is None:
            return_1D_lists = False
        use_files, use_ref_stars = self._getFilesInFilter(filter_type, **kwargs)
        listPltColor = []
        listHJD = []
        listFlux = []
        listDFlux = []
        for i in range(len(use_files)):
            currentFile = use_files[i]
            fileHJD, fileFlux, fileDFlux = currentFile._getFluxData(self.name, self.coords, use_ref_stars[i], **kwargs)
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
    
    def addData(self, file, comparison_star): #TODO: make comparison_star optional...? But, don't want it to be optional for AIJFiles... hmmm...
        """
        Adds an AIJFile object or an OtherDataFile object to the target object.
        
        :param file: File to be added to the target
        :type file: AIJFile or OtherDataFile
        :param comparison_star: name of the comparison star to use to calculate characteristic error
        :type comparison_star: str
        """
        self.files.append(file)
        self.refStars.append(comparison_star)
        
    def _getFilesInFilter(self, filter_type, **kwargs): #returns lists of the files and their reference stars which match a filter type given by user (and are in list of files given by user, if speficied)
        """
        Returns all files added to target with filter_type.
        
        :param filter_type: filter in which data in files to be returned was taken
        :type filter_type: str
        
        kwargs
        """
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
    
    def makeLightCurve(self, filter_type, plotHJD=None, plotFlux=None, runShowPlot=None, bar_position=None, set_title=None, **kwargs):
        """
        Creates a light curve from target data added in filter_type
        
        :param filter_type: filter in which data to be plotted was taken
        :type filter_type: str
        :param plotHJD: indicates whether x-axis should be HJD or phase. Default is phase (False)
        :type plotHJD: bool
        :param plotFlux: indicates whether y-axis should be magnitude or flux. Default is magnitude (False)
        :type plotFlux: bool
        :param runShowPlot: indicates whether to run show() command at end of sequence-- if False, can add/edit plot (i.e., plot additional items on light curve, like a fit line, etc.). Default is True
        :type runShowPlot: bool
        :param bar_position: determines the y-axis position of the characteristic error bars. Calculates automatically by default, but to set manually, this argument can be used.
        :type bar_position: float
        :param set_title: title to be used for light curve. Default is the target name.
        :type set_title: str
        """
        if runShowPlot is None:
            runShowPlot = True
        if plotFlux is None:
            yIsFlux = False
        else:
            yIsFlux = plotFlux
        if plotHJD is None:
            xIsHJD = False
        else:
            xIsHJD = plotHJD
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
                xData.append(phaseData(HJDs[i], self.period, hjd_offset=self.phase_hjd))
        if bar_position is not None:
            errYpos = bar_position
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
        if set_title is None:
            plt.title(str(self.name))
        else:
            plt.title(str(set_title))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if runShowPlot:
            plt.show()
    
    def getFourierFit(self, filter_type, lineFmt=None, showPlot=None, forceN=None, **kwargs):
        """
        Fits a Fourier series to data added to target in filter_type
        
        :param filter_type: filter in which data to be fitted was taken
        :type filter_type: str
        :param lineFmt: matplotlib parameter for style and color of line indicating Fourier Fit. Default is "k-"
        :type lineFmt: str
        :param showPlot: determines whether or not to show the graph of the resulting fit.
        :type showPlot: bool
        :param forceN: highest N out to which to calculate fit (see definition of Fourier series). By default, uses unit-lag auto-correlation to determine best N for fit, but this parameter overries that process.
        :type forceN: int
        :returns: popt, pcov (fitted parameters and their covariance-- see curve_fit in SciPy package)
        :rtype: array, 2d array
        """
        if lineFmt is None:
            lineColor = "k-"
        else:
            lineColor = lineFmt
        if showPlot is None:
            showPlot = True
        HJD, Mag, dMag, pltColor, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        phase = phaseData(HJD, self.period, hjd_offset=self.phase_hjd)
        labs = ("Phase", str(filter_type) + " Mag")
        popt, pcov = fourierFit(phase, Mag, dMag, plot_title=self.name, plot_labels=labs, forceN=forceN)
        if showPlot:
            self.makeLightCurve(filter_type, runShowPlot=False, **kwargs)
            xAxis = np.linspace(0, 1, 100)
            plt.plot(xAxis, fourier(xAxis, *popt), lineColor)
            plt.show()
        return popt, pcov
    
    def calcPeriod(self, filter_type, period_range=None, graphP=None, printP=None, **kwargs):
        """
        Calculates a period for the target based upon data added in filter_type.
        
        :param filter_type: filter in which data to be used for calculation of period was taken
        :type filter_type: str
        :param period_range: renge, in days, in which to search for periods
        :type period_range: tuple
        :param graphP: determines whether or not to print periodogram to the console. Default is False.
        :type graphP: bool
        :param printP: determines whether or not to print the value of the best period to the console. Default is False.
        :type printP: bool
        :returns: calcP -- calculated period
        :rtype: float 
        """
        HJD, Mag, dMag, pltColor, isCharError = self.getMagData(filter_type, return_1D_lists=True, **kwargs)
        calcP = getPeriod(HJD, Mag, dMag, p_range=period_range, makeGraph=graphP, giveP=printP)
        self.periods["calculated_period"] = calcP
        return calcP
    
    def getColorFit(self, filter1, filter2, forceN1=None, forceN2=None, **kwargs): # TODO: add kwargs to make graph not print, so that when this is called inside another function it doesn't make a graph too???
        """
        Fits a Fourier series each to data in filter1 and filter2, then uses these fits to find color index (filter1 - filter2) over phase for a Target.
        
        :param filter1: filter in which data to be fitted was taken (e.g. V in V - R)
        :type filter1: str
        :param filter2: filter in which data to be fitted was taken (e.g. R in V - R)
        :type filter2: str
        :param forceN1: highest N out to which to calculate fit for filter 1 (see definition of Fourier series). By default, uses unit-lag auto-correlation to determine best N for fit, but this parameter overrides the process
        :type forceN1: int
        :param forceN2: highest N out to which to calculate fit for filter 2 (see definition of Fourier series). By default, uses unit-lag auto-correlation to determine best N for fit, but this parameter overrides the process
        :type forceN2: int
        :returns: xPhase, color-- lists of the phases and color index at each of the phases, respectively
        :rtype: list, list
        """
        if forceN1 is None:
            popt_c1, pcov_c1 = self.getFourierFit(filter1, **kwargs)
        else:
            popt_c1, pcov_c1 = self.getFourierFit(filter1, forceN=forceN1, **kwargs)
        if forceN2 is None:
            popt_c2, pcov_c2 = self.getFourierFit(filter2, **kwargs)
        else:
            popt_c2, pcov_c2 = self.getFourierFit(filter2, forceN=forceN2, **kwargs)
        xPhase = np.linspace(0, 1, 1000)
        mag_c1 = fourier(xPhase, *popt_c1)
        mag_c2 = fourier(xPhase, *popt_c2)
        color = mag_c1 - mag_c2
        y_label = filter1 + " - " + filter2
        plt.plot(xPhase, color, "b.")
        plt.xlabel("Phase")
        plt.ylabel(y_label)
        plt.title(self.name)
        plt.show()
        return xPhase, color
    
    def getColorTemperature(self, filter1, filter2, tempModelFit, forceN1=None, forceN2=None, **kwargs):
        """
        Uses :func:`getColorFit` and a quadratic model relating temperature to color index to calculate the color temperature of the Target.
        
        :param filter1: filter in which data to be fitted was taken (e.g. V in V - R)
        :type filter1: str
        :param filter2: filter in which data to be fitted was taken (e.g. R in V - R)
        :type filter2: str
        :param tempModelFit: (a, b, c) where color temperature = a + b(filter1 - filter2) + c(filter1 - filter2)^2
        :type tempModelFit: tuple of length 3
        :param forceN1: highest N out to which to calculate fit for filter 1 (see definition of Fourier series). By default, uses unit-lag auto-correlation to determine best N for fit, but this parameter overrides the process
        :type forceN1: int
        :param forceN2: highest N out to which to calculate fit for filter 2 (see definition of Fourier series). By default, uses unit-lag auto-correlation to determine best N for fit, but this parameter overrides the process
        :type forceN2: int
        :returns: xPhase, temp -- lists of the phases and the calculated temperature at each of those phases, respectively
        :rtype: list, list
        """
        a, b, c = tempModelFit
        xPhase, color = self.getColorFit(filter1, filter2, forceN1=forceN1, forceN2=forceN2, **kwargs)
        temp = a + (b*color) + (c*(color**2))
        plt.xlabel("Phase")
        plt.ylabel("Color Temperature")
        plt.plot(xPhase, temp, "b.")
        plt.title(self.name)
        plt.show()
        return xPhase, temp

    def plotFiles(self, xLabel, yLabel, **kwargs): # TODO: take in multiple labels (one for each file) if desired, and plot each on same axes ---> issue since the order of the files is unclear to the user... maybe a dictionary?
        """
        Plots column with label ``xLabel`` vs. column with label ``yLabel`` for files added to Target.
        
        :param xLabel: label of column in AstroImageJ output files whose data will be considered the x-values of plot
        :type xLabel: str
        :param yLabel: label of column in AstroImageJ output files whose data will be considered the y-values of plot
        :type yLabel: str
        """
        if "use_files" in kwargs: #Still will plot non-AIJ files, right...? Could this be an issue...?
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
        """
        Sets the default period to the period with label ``period_label``.
        
        :param period_label: label of period to be set as default
        :type period_label: str
        """
        self.period = self.periods[period_label]
    
    def addPeriod(self, label, period):
        """
        Adds a period to the list of periods for Target (does NOT set this period as the default period. See :func:`setPeriod`)
        
        :param label: label of period to be added
        :type label: str
        :param period: value of period to be added
        :type period: float
        """
        self.periods[label] = period
    
    def timeOnTarget(self, time_threshold):
        """
        Prints a summary of time spent observing Target based upon files added.
        
        [note that this WILL NOT work for things like ASAS-SN, since images are taken every few days, unless threshold is set to be larger than the time difference between images being taken]
        
        :param time_threshold: value (in days) of how long of a break between images being taken to not include it in the calculation of time on target.
        :type time_threshold: float
        """
        allFiles = self.files
        numFiles = len(allFiles)
        timeSum = 0
        numImages = 0
        for i in range(numFiles): #TODO: implement a checker to skip any files that aren't AIJ files, since those don't have a column called JD_UTC!
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

def getCol(filename, label): #retrieves data (not including label) from "filename" in column labeled "label"
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

def fourierFit(phases, allMags, alldMags, forceN=None):
    if forceN is None:
        topN = 1
        knownN = False
    else:
        topN = forceN
        knownN = True
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
    if topN > 1:
        print("R_21: " + str(popt[2] / popt[1]))
    return popt, pcov

def getPeriod(HJD, TMag, TdMag, p_range=None, makeGraph=None, giveP=None):
    if p_range is None:
        p_min = .01
        p_max = 2
    else:
        p_min, p_max = p_range
    if makeGraph is None:
        makeGraph = False
    if giveP is None:
        giveP = False
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