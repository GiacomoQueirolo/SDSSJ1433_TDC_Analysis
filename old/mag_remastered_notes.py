def flux_ratio(setting,kwargs_result,kwargs_numerics=None,kwargs_data=None,
               lens_model_list=None,light_model_list=None,point_source_list=['LENSED_POSITION'],outnames=False):
    setting=get_setting_module(setting,1)
    if setting.fixed_mag[0]:
        raise RuntimeError("To implement better")
        from lenstronomy.Data.imaging_data import ImageData
        from lenstronomy.Data.psf import PSF
        from lenstronomy.ImSim.image_model import ImageModel
        from lenstronomy.PointSource.point_source import PointSource
        from lenstronomy.LensModel.lens_model import LensModel
        from lenstronomy.LightModel.light_model import LightModel  
        
        lens_model_class = LensModel(lens_model_list=lens_model_list)

        #In order to do that we create the unconvolved modelled image of the ps alone, 
        # and consider the ratio of flux considering the "NONE" psf

        kwargs_lens = kwargs_result["kwargs_lens"]
        data_class = ImageData(**kwargs_data)
        psf_class  = PSF("NONE")
        kwargs_lens_light = kwargs_result["kwargs_lens_light"]
        lens_light_model_class = LightModel(light_model_list=light_model_list)

        kwargs_ps = kwargs_result["kwargs_ps"]
        # note: the relative magnification of point sources is not used as constraints in the fitting in the default settings of lenstronomy.
        # you can set this constraint with the keyword 'fixed_magnification_list' (see next block). 
        # The images are treated otherwise as separate linear amplitudes that are constraint independently of each other.
        point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[True])
        imageModel = ImageModel(data_class, psf_class, lens_model_class,
                                lens_light_model_class,point_source_class,
                                kwargs_numerics=kwargs_numerics)

        for i in range(len(kwargs_lens_light)):
            kwargs_lens_light[i]["amp"]=0

        image_ps = imageModel.image(kwargs_lens, kwargs_lens_light, kwargs_ps,unconvolved=False)

        
        cent_x = kwargs_result["kwargs_ps"][0]["ra_image"]/setting.pix_scale
        cent_y = kwargs_result["kwargs_ps"][0]["dec_image"]/setting.pix_scale

        rng_rd = 1#arcsec of the 1/2 square we consider
        rng    = rng_rd/setting.pix_scale
        numPix = kwargs_data["image_data"].shape[0] 
        amp_i  = []
        for j in range(len(cent_x)):
            x_i = cent_x[j] - rng +(numPix/2.) 
            x_f = cent_x[j] + rng +(numPix/2.) 
            y_i = cent_y[j] - rng +(numPix/2.) 
            y_f = cent_y[j] + rng +(numPix/2.) 
            im_i = []
            for i in range(len(image_ps)):
                im_line =[]
                if i>y_i and i<y_f:
                    for j in range(len(image_ps[i])):
                        if j>x_i and j<x_f:
                            im_line.append(image_ps[i][j])
                    im_i.append(im_line)
            im_i = np.array(im_i)
            amp_i.append(np.sum(im_i))
    else:
        kwargs_ps = kwargs_result["kwargs_ps"][0]
        amp_i=kwargs_ps["point_amp"]
    # Let's do it wrt image A
    #amp_max = np.max(amp_i)
    #FR = np.array(amp_i)/amp_max
    ####################
    # reorder them so that it is in alphabetical order
    new_order    = get_new_image_order(setting,starting_from_A=True)
    amp_i        = [amp_i[i] for i in new_order] 
    ####################
    amp_A = amp_i[0]
    FR    = np.array(amp_i[1:])/amp_A
    if outnames:
        return FR,fr_nms
    else:
        return FR