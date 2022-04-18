import numpy as np
import scipy.linalg as la
import scipy.special



#################################
def diagonalize(H,S=None):
    if S is None:
         S = np.identity(len(H))

    E,C = la.eigh(H,S)

    return E,C

##################################

def printarray( array, filename='array.dat', long_fmt=False ):
#subroutine to print out an ndarry of 2,3 or 4 dimensions to be read by humans

    dim = len(array.shape)
    filehandle = open(filename,'w')
    comp_log = np.iscomplexobj( array )

    if( comp_log ):

        if( long_fmt ):
            fmt_str = '%25.14e%+.14ej'
        else:
            fmt_str = '%15.4f%+.4fj'
    else:

        if( long_fmt ):
            fmt_str = '%25.14e'
        else:
            fmt_str = '%15.4f'

    if ( dim == 1 ):

        Ncol = 1
        np.savetxt(filehandle, array, fmt_str*Ncol )

    elif ( dim == 2 ):

        Ncol = array.shape[1]
        np.savetxt(filehandle, array, fmt_str*Ncol )

    elif ( dim == 3 ):

        for dataslice in array:
            Ncol = dataslice.shape[1]
            np.savetxt(filehandle, dataslice, fmt_str*Ncol)
            filehandle.write('\n')

    elif ( dim == 4 ):

        for i in range( array.shape[0] ):
            for dataslice in array[i,:,:,:]:
                Ncol = dataslice.shape[1]
                np.savetxt(filehandle, dataslice, fmt_str*Ncol )
                filehandle.write('\n')
        filehandle.write('\n')

    else:
        print('ERROR: Input array for printing is not of dimension 2, 3, or 4')
        exit()

    filehandle.close()
