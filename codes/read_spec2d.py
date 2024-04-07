from astropy.io import fits
import numpy as np

def rebin(arr2d, new_shape):
    shape = (new_shape[0], arr2d.shape[0] // new_shape[0],
             new_shape[1], arr2d.shape[1] // new_shape[1])
    return arr2d.reshape(shape).mean(-1).mean(1)

def read_spec2d(spec2dfilepath):
    # ------------------ 2D Spectra ---------------------
    hdul = fits.open(spec2dfilepath)
    # hdul.info()
    #coldef = hdul[1].columns
    # print(hdul[0].header)
    #
    fdata = hdul[1].data['FLUX']
    idata = hdul[1].data['IVAR']
    wdata = hdul[1].data['LAMBDA']
    farr = fdata[0,:,:]
    iarr = idata[0,:,:]
    warr = wdata[0,:,:]
    #
    f_xlen = len(farr[0,:])
    f_ylen = len(farr[:,0])
    i_xlen = len(iarr[0,:])
    i_ylen = len(iarr[:,0])
    w_xlen = len(warr[0,:])
    w_ylen = len(warr[:,0])
    #
    xbinfac, ybinfac = 1, 1
    f_xlen_rb = f_xlen//xbinfac # lengths after rebin
    i_xlen_rb = i_xlen//xbinfac
    w_xlen_rb = w_xlen//xbinfac
    f_ylen_rb = f_ylen//ybinfac
    i_ylen_rb = i_ylen//ybinfac
    w_ylen_rb = w_ylen//ybinfac
    #
    fxlen = f_xlen_rb*xbinfac # corresponding lengths before
    fylen = f_ylen_rb*ybinfac
    ixlen = i_xlen_rb*xbinfac
    iylen = i_ylen_rb*ybinfac
    wxlen = w_xlen_rb*xbinfac
    wylen = w_ylen_rb*ybinfac
    #
    farr = farr[:fylen,:fxlen] # delete extra pixels
    iarr = iarr[:iylen,:ixlen]
    warr = warr[:wylen,:wxlen]
    #
    farr = rebin(farr, (f_ylen_rb,f_xlen_rb))
    iarr = rebin(iarr, (i_ylen_rb,i_xlen_rb))
    warr = rebin(warr, (w_ylen_rb,w_xlen_rb))
    #
    farr = np.nan_to_num(farr, nan=-1000) # replace nan with -1000
    iarr = np.nan_to_num(iarr, nan=-1000)
    warr = np.nan_to_num(warr, nan=-1000)
    #
    skylv = np.median(farr [farr > -1000]) # skylv must be under median
    rms   =    np.std(farr [farr > -1000])
    vmin  = skylv -0.5*rms
    vmax  = skylv +2.0*rms
    #
    ystrfac = 1 # strech factor (vertical)
    ffxlen = len(farr[0,:])
    ffylen = len(farr[:,0])
    iixlen = len(iarr[0,:])
    iiylen = len(iarr[:,0])
    wwxlen = len(warr[0,:])
    wwylen = len(warr[:,0])
    farr =   np.tile(farr, ystrfac) # copying horizontally
    farr = np.resize(farr, (ffylen*ystrfac,ffxlen)) # move extra parts down
    iarr =   np.tile(iarr, ystrfac)
    iarr = np.resize(iarr, (iiylen*ystrfac,iixlen))
    warr =   np.tile(warr, ystrfac)
    warr = np.resize(warr, (wwylen*ystrfac,wwxlen))
    return farr,vmin,vmax,warr