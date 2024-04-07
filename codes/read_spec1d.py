from astropy.io import fits
from astropy    import units as u
from codes.binning_1d import binning_1d

def read_spec1d(spec1dfilepath, binfactor):
    # **observed** spectrum from Keck/DEIMOS
    # Reduced by PypeIt and pypeittospecpro.py
    hdulist = fits.open(spec1dfilepath)
    print('Opened...         ' + spec1dfilepath)
    hdr01data = hdulist[1].data
    w = hdr01data["LAMBDA"][0]
    f = hdr01data["FLUX"][0]
    e = hdr01data["IVAR"][0]
    arr_obs_wave = binning_1d(w, binfactor)
    arr_obs_flux = binning_1d(f, binfactor)
    arr_obs_erro = binning_1d(e, binfactor)
    
    return arr_obs_wave, arr_obs_flux, arr_obs_erro