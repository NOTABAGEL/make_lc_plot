import numpy as np
import warnings
import matplotlib.pyplot as plt
import time
import sys
import math
import numpy.polynomial.polynomial as ply

MJDREF = 51910.0 + 7.428703703703703E-4

spacer = "\n--------------------------------------------------\n"

#by Jonathan Katz and Stack Overflow I know this is terrible OOP sorry

debug = False  #you know what this is for

fname = sys.argv[1]

fit_dat = sys.argv[-1]

#fname = raw_input("Input path to .dat file: ") #Alternative way to get file path, IDK why you would ever do this though

if not fname.endswith(".dat"):
    warnings.warn("Is the file a .dat?")

def met_mjd(np_arr):
    return np.divide(np_arr, 86400) #calculates the proportion of the met/mjd, not absolute time.  this is nescessary to calc the x error bars.  Add MJDREF for absolute time hey you are still reading this

dat = np.genfromtxt(fname)

xcoords = dat[: ,0]

xcoords = met_mjd(xcoords) + MJDREF

xerrors = dat[: ,1]

xerrors = met_mjd(xerrors)

ycoords = dat[: ,2]

yerrors = dat[: ,3]

xcoords_no_ul = []
ycoords_no_ul = []
yerrors_no_ul = []

for i in range(len(yerrors)):
    if yerrors[i] != 0:
        xcoords_no_ul.append(xcoords[i])
        ycoords_no_ul.append(ycoords[i])
        yerrors_no_ul.append(yerrors[i])


#yerrors_no_ul, ycoords_no_ul, xcoords_no_ul  = filter(lambda a: a != 0, yerrors)

print yerrors_no_ul, spacer, ycoords_no_ul, spacer, xcoords_no_ul, spacer

if debug:
    print yerrors

uls = np.invert(dat[: ,3].astype(bool))*1e-7

if debug:
    print uls

def fit(xcoords_fit, ycoords_fit, dat_fit): #does fit heavy lifting
    x_fit = np.linspace(xcoords_fit[0], xcoords_fit[-1], 50)
    np.concatenate([xcoords_fit, x_fit])
    dat_fit = int(dat_fit)
    coefs = ply.polyfit(xcoords_fit, ycoords_fit, dat_fit, w=np.divide(1, yerrors_no_ul))
    y_fit = ply.polyval(x_fit, coefs)
    return (x_fit, y_fit, coefs)

def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c


def get_chi_squared(xcoords_chi, ycoords_chi, yerrors_chi):
    chi_2 = 0
    for i in range(len(ycoords_chi)):
        expect = fit(xcoords_chi, ycoords_chi, fit_dat)[1][i]
        if yerrors_chi[i] != 0:
            chi_2 += ((ycoords_chi[i]-expect)**2)/(yerrors_chi[i]**2)
    return chi_2

print spacer, "Chi-squared: ", get_chi_squared(xcoords_no_ul, ycoords_no_ul, yerrors_no_ul), spacer

plt.plot(fit(xcoords_no_ul, ycoords_no_ul, fit_dat)[0], fit(xcoords_no_ul, ycoords_no_ul, fit_dat)[1])

plt.xlabel("Time (MJD)")

plt.ylabel("Flux (photon cm^-2 s^-1)")

print fit(xcoords_no_ul, ycoords_no_ul, fit_dat)[2]

if debug:
    print spacer, xcoords, spacer, xerrors, spacer, ycoords, spacer, yerrors, spacer, uls

plt.errorbar(xcoords, ycoords, yerr = yerrors, xerr = xerrors, uplims = uls, ls = "")

plt.show()
