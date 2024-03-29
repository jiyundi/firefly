import numpy as np
from astropy.io import fits

def read_templates():
    """"
        SpecPro has 19 templates in formats of fits and dat.
        The following are functions reading in both formats.
    """
    def read_fits_template(temfilename):
        hdulist=fits.open(temfilename)
        # print('Opened...'+temfilename)
        
        # Primary header keywords
        hdr00kwds = hdulist[0].header
        hdr00data = hdulist[0].data
        
        # Calculate wavelength array
        CRVAL = hdr00kwds['CRVAL1']
        CDELT = hdr00kwds['CDELT1']
        CRPIX = hdr00kwds['CRPIX1']
        wave_pixel    = np.arange(len(hdr00data))
        wave_angstrom = CRVAL + CDELT * (wave_pixel + 1 - CRPIX)
        flux_ergscm2A = hdr00data
        arr_template  = np.append([wave_angstrom],
                                  [flux_ergscm2A],
                                  axis=0)
        return arr_template
    
    def read_dat_template(temfilename):
        infile=open(temfilename)
        # print('Opened...'+temfilename)
        
        # 1st column: wave (in Angstrom)
        # 2nd column: flux (unknown unit)
        arr_template = np.zeros((1,2))
        for row in infile:
            if row[0]!='#': # Ignore first-line tag
                lk      = float( row[0:7]  )
                ukf_a0v = float( row[9:17] )
                arr_row = [lk, ukf_a0v]
                arr_template = np.append(arr_template,
                                         [arr_row], 
                                         axis=0)
        arr_template = np.delete(arr_template, (0), axis=0)
        arr_template = arr_template.T
        infile.close()
        return arr_template
    
    def read_sdss_template(temfilename):
        hdulist=fits.open(temfilename)
        # print('Opened...'+temfilename)
        
        # Primary header keywords
        hdr00kwds = hdulist[0].header
        hdr00data = hdulist[0].data
        
        # Calculate wavelength array
        logwavemin = hdr00kwds['CRVAL1']
        Nmaxpixel  = hdr00kwds['NAXIS1']
        log10disp  = hdr00kwds['CD1_1']
        
        logwave       = logwavemin + log10disp * np.arange(Nmaxpixel)
        
        wave_angstrom = 10**logwave
        flux_ergscm2A = hdr00data[0] * 1e-17
        erro_ergscm2A = hdr00data[2] * 1e-17
        arr_template  = np.concatenate(([wave_angstrom],
                                        [flux_ergscm2A],
                                        [erro_ergscm2A]),
                                        axis=0)
        return arr_template
    
    dic_templates = {}
    
    temname1      = 'vvds_lbg'
    temfilename1  = 'templates/template5_olf.fits'
    arr_template1 = read_fits_template(temfilename1)
    dic_templates[temname1] = arr_template1
    
    temname2      = 'vvds_elliptical'
    temfilename2  = 'templates/template0_olf.fits'
    arr_template2 = read_fits_template(temfilename2)
    dic_templates[temname2] = arr_template2
    
    temname3      = 'vvds_s0'
    temfilename3  = 'templates/template1_olf.fits'
    arr_template3 = read_fits_template(temfilename3)
    dic_templates[temname3] = arr_template3
    
    temname4      = 'vvds_early_spiral'
    temfilename4  = 'templates/template2_olf.fits'
    arr_template4 = read_fits_template(temfilename4)
    dic_templates[temname4] = arr_template4
    
    temname5      = 'vvds_spiral'
    temfilename5  = 'templates/template3_olf.fits'
    arr_template5 = read_fits_template(temfilename5)
    dic_templates[temname5] = arr_template5
    
    temname6      = 'vvds_starburst'
    temfilename6  = 'templates/template4_olf.fits'
    arr_template6 = read_fits_template(temfilename6)
    dic_templates[temname6] = arr_template6
    
    temname7      = 'vvds_qso'
    temfilename7  = 'templates/sdss_qso.fits'
    arr_template7 = read_fits_template(temfilename7)
    dic_templates[temname7] = arr_template7
    
    temname8      = 'red_galaxy'
    temfilename8  = 'templates/gal001vel0.fits'
    arr_template8 = read_fits_template(temfilename8)
    dic_templates[temname8] = arr_template8
    
    temname9      = 'green_galaxy'
    temfilename9  = 'templates/gal015vel0.fits'
    arr_template9 = read_fits_template(temfilename9)
    dic_templates[temname9] = arr_template9
    
    temname10      = 'blue_galaxy'
    temfilename10  = 'templates/gal025vel0.fits'
    arr_template10 = read_fits_template(temfilename10)
    dic_templates[temname10] = arr_template10
    
    temname11      = 'lbg_shapley'
    temfilename11  = 'templates/lbg.fits'
    arr_template11 = read_fits_template(temfilename11)
    dic_templates[temname11] = arr_template11
    
    temname12      = 'LoBAL'
    temfilename12  = 'templates/first_lobals.fits'
    arr_template12 = read_fits_template(temfilename12)
    dic_templates[temname12] = arr_template12
    
    temname13      = 'HiBAL'
    temfilename13  = 'templates/first_hibals.fits'
    arr_template13 = read_fits_template(temfilename13)
    dic_templates[temname13] = arr_template13
    
    temname14      = 'a0_star'
    temfilename14  = 'templates/uka0v.dat'
    arr_template14 = read_dat_template(temfilename14)
    dic_templates[temname14] = arr_template14
    
    temname15      = 'f0_star'
    temfilename15  = 'templates/ukf0v.dat'
    arr_template15 = read_dat_template(temfilename15)
    dic_templates[temname15] = arr_template15
    
    temname16      = 'g0_star'
    temfilename16  = 'templates/ukg0v.dat'
    arr_template16 = read_dat_template(temfilename16)
    dic_templates[temname16] = arr_template16
    
    temname17      = 'k0_star'
    temfilename17  = 'templates/ukk0v.dat'
    arr_template17 = read_dat_template(temfilename17)
    dic_templates[temname17] = arr_template17
    
    temname18      = 'm0_star'
    temfilename18  = 'templates/ukm0v.dat'
    arr_template18 = read_dat_template(temfilename18)
    dic_templates[temname18] = arr_template18
    
    temname19      = 'm6_star'
    temfilename19  = 'templates/ukm6v.dat'
    arr_template19 = read_dat_template(temfilename19)
    dic_templates[temname19] = arr_template19
    
    temname20      = 'sdss_early_type'
    temfilename20  = 'templates/sdss_early_type.fits'
    arr_template20 = read_sdss_template(temfilename20)
    dic_templates[temname20] = arr_template20
    
    temname21      = 'sdss_emi_high'
    temfilename21  = 'templates/sdss_emi_high.fits'
    arr_template21 = read_sdss_template(temfilename21)
    dic_templates[temname21] = arr_template21
    
    temname22      = 'sdss_emi_low'
    temfilename22  = 'templates/sdss_emi_low.fits'
    arr_template22 = read_sdss_template(temfilename22)
    dic_templates[temname22] = arr_template22
    
    temname23      = 'sdss_emi_mid'
    temfilename23  = 'templates/sdss_emi_mid.fits'
    arr_template23 = read_sdss_template(temfilename23)
    dic_templates[temname23] = arr_template23
    
    temname24      = 'sdss_late_type'
    temfilename24  = 'templates/sdss_late_type.fits'
    arr_template24 = read_sdss_template(temfilename24)
    dic_templates[temname24] = arr_template24
    
    temname25      = 'sdss_luminous_red'
    temfilename25  = 'templates/sdss_luminous_red.fits'
    arr_template25 = read_sdss_template(temfilename25)
    dic_templates[temname25] = arr_template25
    
    return dic_templates



