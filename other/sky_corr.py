from astropy.io import fits
import numpy as np
import argparse

def load_fits(image_path,HDU=0):
    #load the image and read it as numpy array
    with fits.open(image_path) as hdulist:
        image   = hdulist[HDU].data
    return image
def get_header_fits(image_path):
    with fits.open(image_path,ignore_missing_end=True) as target:
        scihdr  = target[0].header
    return scihdr
    
    
def exposure_time(image_path):
    with fits.open(image_path,ignore_missing_end=True) as hdulist:
        hdr = hdulist[0].header
    return hdr["EXPTIME"]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="correct for sky value")
    parser.add_argument("-s", "--sky", type=float, dest="sky",help="Sky value")
    parser.add_argument("-t", "--time", type=float, dest="time",help="Exposure time")
    parser.add_argument("-p", "--prefix", type=str,default="sc", dest="prefix",help="Output prefix")
    parser.add_argument("FITS",help='fits file to analyse')
    args  = parser.parse_args()
    data =  str(args.FITS) 
    sky  =  float(args.sky) 
    time  = float(args.time) 
    prefix = str(args.prefix)
    if prefix =="":
        print("Warning! You are going to overwrite the input data with the corrected one")
        if input("Continue? [Y/N]")!="Y":
            exit()
        
    overwrite = True
    WD = "./"+"/".join(data.split("/")[:-1])+"/"
    WD = WD.replace("//","/")
    print(data)
    data_image = load_fits(data)
    scihdr  = get_header_fits(data)
    scihdr["HISTORY"]  = "used sky_corr.py to subtract the sky value: "+str(sky) 
    corrected_data_image = data_image - sky
    hdu = fits.PrimaryHDU(data=corrected_data_image,header=scihdr)
    output = WD +prefix+data.split("/")[-1]
    print("Saving output as ",output)
    hdu.writeto(output, overwrite=overwrite)
    errdata = WD+"e."+data.split("/")[-1]
    print("Correcting the error frame ",errdata)
    err_image = load_fits(errdata)
    errhdr    = get_header_fits(errdata)
    scihdr["HISTORY"]  = "used sky_corr.py to correct the error frame considering the sky value: "+str(sky)+ " and the exposure time :"+str(time) 
    corr_err_image = np.sqrt(err_image**2 +  (sky/time) ) 
    errhdu = fits.PrimaryHDU(data=corr_err_image,header=errhdr)
    output_err = WD +"e."+prefix+data.split("/")[-1]
    print("Saving error output as ",output_err)
    errhdu.writeto(output_err, overwrite=overwrite)
    print("Success")
