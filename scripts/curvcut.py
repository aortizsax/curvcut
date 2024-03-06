# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 20:56:13 2020

@author: aorti

This code is meant to help the user define how many zeros to keep in their count
table with mathematical backing. 

Load values min to max.

Input is a histogram of the number of features (yaxis)
with a certain number of zeros(xaxis). 

The histogram is then tansformed with log(hist + c)
y = log(hist + c)

A curve representing the transfromed histogram data is be made with cubic-splines. 
ys = cubicspline of y

By maximixing curvature (minizing the radius of the curve), the user can see 
when the histogram columns start to ramp up. That point where it ramps up will 
the cut off.

curvature: kappa = 1/Radius = |ys''| * (1 + ys'^2)^(-3/2)

Usage:
    python 2022-09-13_autoCutoff.py filename
    
Where:
    filename: is a csv or tsv formatted count table 

"""


import os
import pathlib
import sys
import argparse

import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
import time
from jenkspy import JenksNaturalBreaks

startTime = time.time()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument(
        "-ct",
        "--counttable",
        action="store",
        metavar="FILE",
        help="Path to source file.",
    )
    parser.add_argument(
        "-af",
        "--minorallelefreq",
        action="store",
        metavar="FILE",
        help="Path to source file.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        action="store",
        default="output",
        help="Prefix for output files [default=%(default)s].",
    )
    parser.add_argument(
        "-u",
        "--user-cutoff",
        action="store",
        default=-1,
        help="Prefix for output files [default=%(default)s].",
    )

    args = parser.parse_args()
    print("Parsing")

    print("")

    if args.minorallelefreq == None:
        # Get path to make new directory
        # filename = sys.argv
        filename = args.counttable
#        print(filename)
        ID = filename[:-3].split("/")[-1].replace(".", "")
        # set output path
        # path = "/".join(sys.argv[1].split("/")[:-1]) + "/" + zerofiltered
        path = "./" + args.output_prefix
        print("saving to:", path)

        # Check whether the specified path exists or not
        isExist = os.path.exists(path)

        # unmask
        os.umask(0)

        # check output folder
        if isExist == False:
            # Create a new directory because it does not exist
            os.makedirs(path, mode=0o755, exist_ok=False)
            os.chmod(path, 0o755)  # just in case mode above doesn't process
            print("The new directory is created! path:", path)

        # load data
        # eliminate rows/features that are below cutoff
        if filename[-3:] == "tsv":
            data = pd.read_csv(filename, sep="\t")  # index_col=0,header = 0)
            if data.columns.to_list()[0] == "# Constructed from biom file":
                data = pd.read_csv(
                    filename, sep="\t", skiprows=1
                )  # index_col=0,header = 0)
        elif filename[-3:] == "csv":
            data = pd.read_csv(filename, sep=",")  # index_col=0,header = 0)
            if data.columns.to_list()[0] == "# Constructed from biom file":
                data = pd.read_csv(filename, sep=",", skiprows=1)
        dimen = data.shape

        if args.user_cutoff != -1:
            cutoff = int(args.user_cutoff)
            print("Trimming features at user cutoff of", cutoff)
            # vectors of number of presence and absence in each row
            ispres = data != 0.0

            # filter table
            OTUzerofiltered = data[ispres.sum(axis=1) > cutoff]

            # save filtered table
            OTUzerofiltered.to_csv(path + "/" + ID + "table.zerofiltered.csv")

            print(
                "Zero Filtered OTU table saved:",
                path + "/" + ID + "table.zerofiltered.csv",
            )
            exit()

        # vectors of number of presence and absence in each row
        iszero = data == 0.0

        # counts of zeros per row
        zero_sum = iszero.sum(axis=1)

        # bin count of zeros per row
        hist_zero = np.bincount(zero_sum)
#        print(len(hist_zero),hist_zero)
        dimen = [len(hist_zero), len(hist_zero)]

    else:
        # Get path to make new directory
        # filename = sys.argv
        filename = args.minorallelefreq
        ID = filename[:-3].split("/")[-1].replace(".", "")
        print(filename)
        # set output path
        path = "/".join(sys.argv[1].split("/")[:-1]) + "/filtered-sfs"
        path = "./" + args.output_prefix
        print("saving to:", path)

        # Check whether the specified path exists or not
        isExist = os.path.exists(path)

        # unmask
        os.umask(0)

        # check output folder
        if isExist == False:
            # Create a new directory because it does not exist
            os.makedirs(path, mode=0o755, exist_ok=False)
            os.chmod(path, 0o755)  # just in case mode above doesn't process
            print("The new directory is created! path:", path)

        # load data
        # eliminate rows/features that are below cutoff
        hist_zero = pd.read_csv(filename, delimiter=",", header=0, index_col=0)
        hist_zero = hist_zero.iloc[1:, 0].tolist()
        hist_zero.reverse()
 #       print("hist_zero", len(hist_zero), hist_zero)

        dimen = [len(hist_zero), len(hist_zero)]

    # ID = filename[:-3].replace('/','').replace('.','')

    bact_zero_table = {}
    # histogram information as dictionary
    for i, zero_i in enumerate(hist_zero):
        bact_zero_table[i] = int(zero_i)

    # log transformation constant
    constant = 100  # con[0]

    # yvalues are hieght from histograms
    y_axis = list(hist_zero)[:]
    

#    jnb = JenksNaturalBreaks(2)
#    jnb.fit(y_axis)

    # initialize postion from histogram
    x_axis = list(bact_zero_table.keys())
    n = len(y_axis)

    # save orginal
    y_original = y_axis.copy()
    # intialize array for reverse original
    y_reverse = np.zeros((1, len(y_axis)))[0]

    acc = 0  # acculumation of yvalues
    i = 0  # index

    # Sorting by acculuating
    for i, y in enumerate(y_axis):
        # Accultion of histogram for monotonically decreasing funciton
        acc += y
        y_axis[i] = acc

    # initialize yaxis 2 become
    y_axis2 = np.zeros((1, len(y_axis)))[0]
    for i in range(len(y_axis)):  ####reverse!!!!
        index = i  # len(y_axis) - i -1
        y_reverse[index] = y_original[i]

    # yaxix 2 TRANSFORMATION
    y_axis2 = []
    for y in y_axis:
        trans = np.log(y + constant)
        y_axis2.append(trans)
    y_axis = y_axis2.copy()
    y_avg = sum(y_axis) / len(y_axis)

    # Cubic slpine
    cs = CubicSpline(x_axis, y_axis)  # yaxis spline and derivatives
    print("Cubic Spline Finshed")

    # xaxis spline
    step = 0.05
    xs = np.arange(x_axis[0], x_axis[-1], step)

    # initalize curvature array
    kappa = []
    for i in range(len(xs)):
        kTemp = abs(cs(xs, 2)[i])  # numerator
        kTemp /= (1 + (cs(xs, 1)[i]) ** 2) ** (3 / 2)  # denominator
        kappa.append(kTemp)

    # find the last max local maximum
    maxCurv = 0
    maxPoints = []
    for k in range(0, len(kappa) - 1):
        # print(kappa[k])
        if kappa[k] - kappa[k - 1] > 0:
            if kappa[k + 1] - kappa[k] < 0:
                if kappa[k] > maxCurv:
                    maxCurv = kappa[k]
                    maxPoints.append(kappa[k])

    # maxCurv = max(kappa[20:]) old find max     cutoff = dimen[1] - cutoff
    index = kappa.index(maxCurv)
    cutoff = xs[index] - 1
   # print('CUTOFF',cutoff)
    cutoffpres = dimen[1] - cutoff
    cutoffspres = []
    for point in maxPoints:
        index = kappa.index(point)
        cutoffspres.append(dimen[1] - xs[index])

    executionTime = time.time() - startTime
    print("Execution time in seconds: " + str(executionTime))

    if args.minorallelefreq == None:
        # cutoff calculated for absent species
        cutoffabsent = cutoff
        
        print(
            "Recommended  Cutoff: Trim Features that are present in less than "
            + str(round(cutoffpres, 3)) #+ #' or ' + str(len(jnb.groups_[-1]))
            + " samples. Close graph when you have made a decision on the cuttoff."
        )
  
 
        # calculate cutoff for number of present species/features
        # cutoff = dimen[1] - cutoff
        cutoff = dimen[1] - cutoff
        presentheight = y_original.copy()
        presentheight.reverse()
        presentx = x_axis.copy()
        presentx.reverse()
        presentx = np.subtract(dimen[1], presentx)
        start = 0
        end = -1

        plt.bar(x=presentx, height=presentheight)
        plt.xlabel("# samples a feature is present in")
        plt.ylabel("# of features")
        plt.title(
            "Cutoffs to consider: Curvcut algo: " 
            + str(round(cutoff, 3)) 
            + "; Fisher-Jenks: "
   #         + str(len(jnb.groups_[-1])) 
            + " "
        )
        plt.axvline(x=cutoff, color="red",ls='--',
                    label="Curvcut: "+str(cutoff))
   #     plt.axvline(x=len(jnb.groups_[-1]), color="black",ls='-.',
    #                label="Fisher-Jenks: "+str(len(jnb.groups_[-1])))
        plt.xlim([0, dimen[1]])
        plt.legend(loc='upper right')
        #        plt.savefig(path + "/"+ID+"processingcutoff.png")
        plt.savefig(
            path + "/" + ID + "processingcutofffinal.svg", format="svg", dpi=1200
        )
        plt.show()
        
        
        # calculate cutoff for number of present species/features

        cutoff = int(input('After visual inspection of the graph, what cutoff would you like to use?\nPlease enter integer:'))

        # vectors of number of presence and absence in each row
        ispres = data != 0.0

        # filter table
        OTUzerofiltered = data[ispres.sum(axis=1) > cutoff]

        # save filtered table
        OTUzerofiltered.to_csv(path + "/" + ID + "table.zerofiltered.csv")
        print(
            "Zero Filtered OTU table saved:", path + "/" + ID + "table.zerofiltered.csv"
        )
     
        
    else:
        print(
            "Recommended  Cutoff: Trim alleles that are present in less than "
            + str(round(dimen[1] - cutoff, 3))
            + " sequences"
        )
        print("")
        print("")
        print("")

        # plot Curvature
        plt.plot(xs, kappa)  # ,color=colmapBIO[clustNum])
        plt.gca().xaxis.grid(True)
        plt.title("Maximizing Curvature")
        plt.xlabel("# samples species is absent")
        plt.ylabel("# of species")
        plt.savefig(path + "/" + ID + "maxcurvplot.png")
        plt.show()

        # cutoff calculated for absent species
        cutoffabsent = cutoff
        plt.bar(x=x_axis, height=y_original)
        plt.xlabel("# samples feature is absent in")
        plt.ylabel("# of features")
        plt.title(
            "Trim Features that are absent in more than "
            + str(round(cutoffabsent, 3))
            + " samples"
        )
        plt.axvline(x=cutoffabsent, color="red")
        plt.xlim([dimen[1] - 50, dimen[1]])
        #        plt.savefig(path + "/"+ID+"processingcutoff.png")
        plt.savefig(path + "/" + ID + "processingcutoffA.svg", format="svg", dpi=1200)
        plt.show()


        plt.clf()        
        
        # calculate cutoff for number of present species/features
        cutoff = dimen[1] - cutoff

        # calculate cutoff for number of present species/features
        # cutoff = dimen[1] - cutoff
        presentheight = y_original.copy()
        presentheight.reverse()
        presentx = x_axis.copy()
        presentx.reverse()
        presentx = np.subtract(dimen[1], presentx)
        start = 0
        end = -1

        plt.bar(x=presentx, height=presentheight)
        plt.xlabel("# samples a feature is present in")
        plt.ylabel("# of features")
        plt.title(
            "Cutoffs to consider: curvature: " 
            + str(round(cutoff, 3)) 
            + "; Fisher-Jenks: "
   #         + str(len(jnb.groups_[-1])) 
            + " "
        )
        plt.axvline(x=cutoff, color="red",ls='--')
    #    plt.axvline(x=len(jnb.groups_[-1]), color="black",ls='-.')
        plt.xlim([0, 50])
        #        plt.savefig(path + "/"+ID+"processingcutoff.png")
        plt.savefig(path + "/" + ID + "processingcutoffB.svg", format="svg", dpi=1200)
        plt.show()
        
        
        
