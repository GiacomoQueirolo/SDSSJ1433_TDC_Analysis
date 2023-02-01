from astropy.io import fits

def load_fits(image_path,HDU=0):
    #load the image and read it as numpy array
    with fits.open(image_path) as hdulist:
        image   = hdulist[HDU].data
    return image

def psf_correction(psf,setting):
    corr_psf = psf
    dl_mean0 = np.mean(psf[0])
    dl_mean1 = np.mean(psf[:,0])
    corr_psf = np.ones((len(psf)+2,len(psf[0])+2))
    corr_psf[0]*=dl_mean0
    corr_psf[-1]*=dl_mean0
    corr_psf[:,0]*=dl_mean1
    corr_psf[:,-1]*=dl_mean1
    corr_psf[0][0]=dl_mean0
    corr_psf[0][-1]=dl_mean0
    corr_psf[-1][-1]=dl_mean1
    corr_psf[-1][0]=dl_mean1
    for i in range(1,len(corr_psf)-1):
        for j in range(1,len(corr_psf[0])-1):
            corr_psf[i][j]=psf[i-1][j-1]
    return corr_psf


def get_kwargs_psf(psf_file,err_psf_file):
    
    psf_image = load_fits(psf_file)
    err_psf_image = load_fits(err_psf_file)
    err_psf_image = psf_correction(err_psf_image)

    kwargs_psf = {'psf_type': "PIXEL", 
          'kernel_point_source':psf_image,
          'point_source_supersampling_factor': 5,
          'psf_error_map':err_psf_image}
    return kwargs_psf