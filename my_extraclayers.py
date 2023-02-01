from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser(description="Extract lens light model LLM as model_* ")
parser.add_argument("-n", dest="LLM", help="Name of LLM")
LLM = argparse.LLM
with fits.open(LLM) as f:
    for i,fi in enumerate(f[:-1]):
        if i==0:
            image = fi.data
            hdr = fi.header
        elif  i==len(f[:-1])-1:
            model = fi.data
hdu = fits.PrimaryHDU(data=image,header=hdr)
name_new_model = "model_"+LLM
hdu.writeto(name_new_model)

