# -*- coding: utf-8 -*-
import observation_helper.observation_helper as oh
import os.path
"""
Welcome to the observation_helper tutorial file! This should help you get a basic understanding of how this code might be most helpful to you, and how to use some of the basic functions laid out below.

First, we must use a dictionary to define which apertures from AstroImageJ were used for which target stars. Here, GSC was T1, Variable A was T2, Variable B was T3, etc. This must be included, even if you only had one target aperture.
"""
target_dictionary = {"GSC 02950-01137":1, "Variable A":2, "Variable B":3, "Variable C":4, "Variable D":5, "Variable E":6}
"""
At first, we'll just focus on the first target star, GSC 02950-01137 (if you're interested in any of the targets or why they were chosen in the included data, please feel free to read my capstone paper)

Similar to above, the following line defines, for each comparison/reference star, the aperture it used (first number in parentheses) and the known magnitude of that star in a particular filter from the literature (second number in parentheses)
"""
ref_dictionary_V = {"Ref7":(7,15.163), "Ref8":(8,14.387), "Ref9":(9,15.517), "Ref10":(10,13.788), "Ref11":(11,14.112), "Ref12":(12,16.082), "Ref13":(13,15.060), "Ref14":(14,16.137), "Ref15":(15,15.986), "Ref16":(16,15.606)}

"""
Here, we're creating an object for each text file output by AIJ. I have included my own data files so that we have something to work with for this tutorial.
The first argument is the name of the actual text file; next is the filter (whatever you want to call it, as long as you use that name consistently throughout), then the target star dictionary, then the reference star dictionary for the apropriate filter
The last argument is the color/style the points on the light curve will show up as for that file; see the documentation on matplotlib color format for more information
"""
jan07_V = oh.AIJFile(os.path.join("test_data", "allStars_jan07_V_data.txt"), "V", target_dictionary, ref_dictionary_V, point_format="r.") 

"""
Now, I'm creating an object to represent my target star. The first argument is the name of the target (MUST be the same as its name in the target dictionary in line 10)
The second argument is a tuple with RA and Dec in J2000 (decimal format, RA first and then Dec), and third is a "given" period of the object (perhaps one known from literature, or a guess-- I'll deminstrate how to calculate a period from the data later in this tutorial)
A couple of targets have an additional optional argument, phase0_hjd, which is just an HJD to offset the phased light curve by (i.e. what HJD corresponds to a phase of 0). If phase0_hjd is not specified, it defaults to 0. Basically, this argument shifts your phased light curve left or right a bit, and you can change it as you'd like.
"""
gsc = oh.Target("GSC 02950-01137", (102.7425089, 41.2720900), .25207, phase0_hjd=2458480.41)

"""
Now, I want to add the data file to my target. This must be done for EACH data file we want to work with.
For addData, the first argument is the file to add, and the second argument is the reference star to use in order to calculate the standard error for that particular target star (see my Capstone paper for more details on that process)
"""
gsc.addData(jan07_V, "Ref7") 

"""
Let's plot a light curve of GSC 02950-01137.
"""
print("Light Curve of GSC in V filter from Jan. 10th, phased using given period")
gsc.makeLightCurve("V")
"""
Hooray, we did it! Let's take a look at some more features.

The data look a little bit scattered in that light curve-- let's see of we can clean it up a bit. We can use a sigma-clipping algorithm to adjust the scatter of the data in this file. For a detailed explanation of this particular sigma-clipping proceedure, see my Capstone paper(s).
The clipFile function has the name of the reference star to use to clip as the first argument, and the sigma over which to clip as the second argument (i.e., in this example, points more than 2 standard deviations from the mean are discarded)
"""
print("------------")
jan07_V.clipFile("Ref14", 2)
print("Clipping " + str(jan07_V.filename))

"""
Let's see if that did anything! We'll plot the light curve again to see.
"""
print("------------")
print("Light Curve of GSC in V filter from Jan. 7th, phased using given period, post sigma-clipping")
gsc.makeLightCurve("V")
"""
Looks like that got rid of a few data points-- we'll leave it like this for now. Beware-- from this point on in the code, this data file will always be using the clipped data (unless you unclip it)! 

What if we want to see this light curve as a function of HJD instead of phase, and flux instead of magnitude? Let's do that. 
"""
print("------------")
print("Light Curve of GSC in V filter from Jan. 7th, plotted with HJD on the x-axis and flux on the y-axis, post sigma-clipping")
gsc.makeLightCurve("V", plotHJD=True, plotFlux=True)

"""
Now, what if we want to look at some data from another file? All we need to do is add the file, just like we did the first one. Let's add some more V-filter data from a different night. The same target stars and comparison stars were used in the same apertures for this file.
"""
jan04_V = oh.AIJFile(os.path.join("test_data", "allStars_jan04_V_data.txt"), "V", target_dictionary, ref_dictionary_V, point_format="b.") #Creating the file object
jan04_V.clipFile("Ref14", 2) #clipping the file
gsc.addData(jan04_V, "Ref7") #adding the data to the target object

#
"""
Now, let's take a look at the light curve with the new data added
"""
print("------------")
print("Light Curve of GSC in V filter with data from Jan. 7th and Jan. 4th")
gsc.makeLightCurve("V")

"""
What if we wanted to look at some of the raw data from the AIJ text files, like the flux for different comparison stars? We can do that pretty easily!
If we just want to look at the relative flux values for comparison star 14 from the January 4th file, we can do the following:
"""
print("------------")
print("Relative flux for C14 from January 4th file in V filter")
jan04_V.plotFile("JD_UTC", "rel_flux_C14")

"""
We can look at the sky background per pixel values for our target star from the 7th (and compare to the 4th) all in one plot by using plotFiles on our target object
"""
print("------------")
print("Sky background per pixel for GSC, for January 4th and 7th, in V filter")
gsc.plotFiles("JD_UTC", "Sky/Pixel_T1")


"""
Up until now, we've been using a value from the literature for our period for GSC 02950-01137. What if we want to calculate our own using the data we have so far? Let's do that!
"""
print("------------")
print("Calculating a period for GSC 02950-01137...")
newly_calculated_period = gsc.calcPeriod("V")
print("Calculated period for GSC 02950-01137: " + str(newly_calculated_period))
"""
Check out the API documentation for options for showing the periodogram, setting a range in which to search for a period, etc.,

Now, if we want to use this period for calculations, graphs, etc, we need to set this period as the one to use for our target object. We can do that as follows:
"""
gsc.setPeriod("calculated_period")

"""
Let's take a look at our light curve, phased with this new period.
"""
print("------------")
print("Light curve for GSC, phased using calculated period")
gsc.makeLightCurve("V")

"""
What if we want to add some data we took in the R filter? We can do that by adding a new file, like what we did above with the V-filter data.
We need to define the magnitudes of the reference stars in the R filter, so we'll do that with a new dictionary.
"""
ref_dictionary_R = {"Ref7":(7,14.898), "Ref8":(8,14.164), "Ref9":(9,15.401), "Ref10":(10,13.610), "Ref11":(11,14.168), "Ref12":(12,15.536), "Ref13":(13,14.708), "Ref14":(14,15.754), "Ref15":(15,15.760), "Ref16":(16,15.291)} 

"""
Now, we create file objects for the two R-filter data files we have:
"""
jan04_R = oh.AIJFile(os.path.join("test_data", "allStars_jan04_R_data.txt"), "R", target_dictionary, ref_dictionary_R, point_format="b.")
jan06_R = oh.AIJFile(os.path.join("test_data", "allStars_jan06_R_data.txt"), "R", target_dictionary, ref_dictionary_R, point_format="g.")

"""
Now, let's add this data file to our target object.
"""
gsc.addData(jan04_R, "Ref7")
gsc.addData(jan06_R, "Ref7")

"""
Now, let's take a look at that R-filter data!
"""
print("------------")
print("Light curve for GSC using data in R filter")
gsc.makeLightCurve("R")

"""
These AstroImageJ data files actually contain data on more than one target. So, what if we wanted to look at data from one of the other targets? To do that, we need to create a new target object. One target in particular, which I called Variable B, was pretty interesting to look at-- let's take a look!
"""
varB = oh.Target("Variable B", (102.9408086, 41.1729448), 0.121677522372, phase0_hjd=2458480.41)
varB.addData(jan04_V, "Ref12")
varB.addData(jan07_V, "Ref12")
varB.addData(jan04_R, "Ref12")
varB.addData(jan06_R, "Ref12")
print("------------")
print("Light Curve of Variable B in V filter, using data from Jan. 4th and Jan 7th")
varB.makeLightCurve("V")

"""
Say we wanted to only use specific reference stars to calculate the magnitude of the target star instead of all of them. We could do that for all of the graphs and calculations in this tutorial by going and deleting certain reference stars from our reference star dictionaries. Alternatively, we can use the keyword argument "use_ref_stars" in almost any of the functions in observation_helper. For example, we can do the following:
"""
print("------------")
print("Light Curve of Variable B in V filter, Jan. 4th & 7th, magnitude calculated using only comparison stars 7, 8, and 10")
varB.makeLightCurve("V", use_ref_stars=["Ref7", "Ref8", "Ref10"])

"""
Now that we've plotted some light curves in a couple of different ways, we can take a look at some of the other things built into observation_helper. First, we can look at how to fit a Fourier Series to these data. Let's do this for our GSC target.
"""
print("------------")
print("Fourier Fitted Light Curve for GSC using data in V filter")
gsc.getFourierFit("V")
"""
The resulting plot shows the light curve with the Fourier Fit as a black line fitted over it. This also prints certain parameters, like N and p. N is the top N the fourier series stopped at when attempting to fit (i.e., the series didn't go to infinity, it just went to N). p is related to the unit-lag auto-correlation process, which finds the best N. If N>1, this function will also print out R_21 (a comparison of the fitted values for N=2 and N=1), which can be a helpful piece of information when attempting to classify variable stars. You can override the unit-lag auto-correlation process by including forceN as an integer argument in the function call (see API documentation). The function returns the actual Fourier Series parameters of best fit from the resulting N as an array, as well as the covariance matrix associated with the parameters as a second returned argument.
"""

"""
There are a few functions included which use a Fourier Series fit to look at some interesting properties of a target star. Let's look at the color of GSC using a color fitting function.
"""
print("------------")
print("Color of GSC by phase, calculated by fitting a Fourier series to the V filter data and the R filter data, and calculating V - R for each phase using the Fourier series fitting results")
gsc.getColorFit("V", "R")

"""
This graph, on its own, might look a bit strange-- the color of the star might not vary quite so wildly. If we take a look at the individual Fourier Series fits which gave rise to this plot, we can see that one has N=1 (in V), and the other has N=3 (in R). Likely, this is because there is not enough data in the R light curve for the algorithm to correctly fit to, so the algorithm over-fit to an N=3 in the R filter. We can fix this by forcing N=1 in the R filter, like so:
"""
print("------------")
print("Color of GSC by phase, calculated by fitting a Fourier series to the V filter data and the R filter data, and calculating V - R for each phase using the Fourier series fitting results, but forcing N=1 for the R filter Fourier series")
gsc.getColorFit("V", "R", forceN2=1)

"""
Be sure to check out the documentation (included PDF, or online at kainmccall.github.io/observation_helper/) to see all functions!
"""