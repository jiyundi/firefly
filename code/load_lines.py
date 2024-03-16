import numpy as np

def read_lines():
    """"
        SpecPro has 8 sets of lines in the dat format.
        The following will load 2 of these that are most-commonly used.
    """
    def read_dat_line(linfilename, dic_lines):
        infile=open(linfilename)
        arr_line = np.zeros((1,2))
        for row in infile:
            if row[0]!='#': # Ignore first-line tag
                wave = row[0:7]
                name = row[8:-1]
                arr_row  = [wave, name]
                arr_line = np.append(arr_line, [arr_row], axis=0)
        arr_line = np.delete(arr_line, (0), axis=0)
        arr_line = arr_line.T
        infile.close()
        del name, wave, row, arr_row, linfilename
        for j in range(len(arr_line[0,:])): 
            name  =  str(arr_line[1,j])
            wave  =  str(arr_line[0,j])
            try:
                waves = dic_lines[name]
                wave_ = np.array([wave], dtype='f')
                waves = np.append(waves, wave_)
                dic_lines[name] = np.array(waves,  dtype='f')
            except KeyError:
                dic_lines[name] = np.array([wave], dtype='f')
        return dic_lines
        
    dic_emi_lines = {}
    linfilename1  = 'lines/emlines.dat'
    dic_emi_lines = read_dat_line(linfilename1, dic_emi_lines)
    
    dic_abs_lines = {}
    linfilename2  = 'lines/prominentellipticalabsorptionlines.dat'
    dic_abs_lines = read_dat_line(linfilename2, dic_abs_lines)
    
    return dic_emi_lines, dic_abs_lines
