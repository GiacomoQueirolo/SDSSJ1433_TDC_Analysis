# This script takes in the output sci and weight images and generates the error image.
"""
Hi Giacomo, here is the python script:

If you just run this script it will open file dialogues, so it’s interactive.
The first file that is read in is the original whole _sci.fits file. This is just to read out the CCD gain from the header, so you can just read it out yourself and plug in the value in the code.
# ACTUALLY all my images are already in e-/sec 
Next, it reads in the _sci.fits cutout you use for the modeling, as well as the _wht.fits cutout.

Then you need to provide a DS9 region file where you mark with a box an empty region in your cutout. It might be more accurate to use an empty region in the full _sci.fits file, not so close to the lens, so you don’t bias your background estimation with the lens light. This region file is used to calculate the standard deviation inside this box which will be used as background noise.

Then we add in quadrature Poisson noise to the background noise. Please check if your data is in units counts or counts per seconds (usually it’s the latter), because you need to use a different line in the code depending on the units. Currently it’s set to counts per second.

You can read more about the error map construction and the use of this code in my TDCOSMO X paper. And please don’t hesitate to ask me if you have more questions or something doesn’t work.

Cheers,
Sebastian

"""

from astropy.io import fits
import numpy
import sys, math, os
import astropy.io.fits as pyfits
import re
import tkinter, tkinter.filedialog
root = tkinter.Tk()
root.withdraw()

print("Choose sci file path.")
scifile_orig = tkinter.filedialog.askopenfilename()

#Obtain CCD gain
#hdu_number = 0
#header=fits.getheader(scifile_orig, hdu_number)
#ccdgain=header['CCDGAIN']
ccdgain=1
print("Choose sci file cutout path.")
scifile = tkinter.filedialog.askopenfilename()
print("Choose wht file cutout path.")
whtfile = tkinter.filedialog.askopenfilename()
outfile = input("Choose output file name for error map: ")

# read in the sci and EXP weight image
sciIm = pyfits.open(scifile)[0].data
whtIm = pyfits.open(whtfile)[0].data

# pixel of lower left corner of box to compute sigma (1-based)
print("Open region file for errormap.")
err_region_path = tkinter.filedialog.askopenfilename()
err_region = open(err_region_path,'r')
for line in err_region:
    if 'box' in line:
        err_box = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)

print(err_box)
llx = int(round(float(err_box[0])-(float(err_box[2])/2)))
lly = int(round(float(err_box[1])-(float(err_box[2])/2))) # note: assumes a symmetric box!

bsize = int(round(float(err_box[2])))

errForZeroWhtPix = 100. # a very high number to effectively discard pixel

blankIm = sciIm[(lly-1):(lly+bsize),(llx-1):(llx+bsize)]

# standard deviation of blankIm to use as pedestal
sped = blankIm.std()
print("sigma_pedestal = ", sped)

# Poisson noise from astrophysical sources
#for cps
poissonIm = numpy.where (whtIm>0., numpy.sqrt(numpy.abs(sciIm) / (whtIm * ccdgain)), float('nan'))
#for counts
#poissonIm = numpy.where (whtIm>0., numpy.sqrt((numpy.abs(sciIm) * whtIm) / ccdgain), float('nan'))  for cps instead of counts

# If Poisson > 2*Pedestal, add Poisson and Pedestal in quadrature.  Else use Pedestal.
errIm = numpy.where (poissonIm > 2*sped, numpy.sqrt(poissonIm**2 + sped**2), sped)
errIm = numpy.where (whtIm > 0, errIm, errForZeroWhtPix) #for lanczos3


errFits=pyfits.PrimaryHDU(errIm)
errFits.writeto(outfile)
