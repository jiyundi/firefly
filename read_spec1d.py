from astropy.io import fits
from astropy    import units as u
from binning_1d import binning_1d

def read_spec1d(spec1dfilepath):
    # **observed** spectrum from Keck/DEIMOS
    # Reduced by PypeIt and pypeittospecpro.py
    binfactor = 24

    hdulist = fits.open(spec1dfilepath)
    print('Opened...         ' + spec1dfilepath)
    hdr01data = hdulist[1].data
    w = hdr01data["LAMBDA"][0]
    f = hdr01data["FLUX"][0]
    e = hdr01data["IVAR"][0]
    arr_obs_wave = binning_1d(w, binfactor) * u.AA
    arr_obs_flux = binning_1d(f, binfactor) * u.dimensionless_unscaled
    arr_obs_erro = binning_1d(e, binfactor) * u.dimensionless_unscaled
    
    return arr_obs_wave, arr_obs_flux, arr_obs_erro