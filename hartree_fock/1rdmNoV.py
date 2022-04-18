import numpy as np

Nelec = 6
Nsites = Nelec
Vbias = 0.5
bias = False



def hubbard_1D( Nsites, Vbias, bias = True):
    t = 1.0
    NsitesArray = np.arange(1, Nsites)
    NA = NsitesArray[:len(NsitesArray)//2]
    NB = NsitesArray[len(NsitesArray)//2:]

    Tmat = np.zeros( ( Nsites, Nsites ) )

    for i in range(Nsites-1):
        Tmat [ i, i+1 ] = Tmat [ i+1, i ] = -t

    if (bias):
        for i in NA:
            Tmat [ i, i ] = -Vbias/2.0
        for i in NB:
            Tmat [ i, i ] = Vbias/2.0

    Tmat [ 0, Nsites-1 ] = Tmat [ Nsites - 1, 0 ] = -t

    return Tmat

Tmat = hubbard_1D(Nsites, Vbias, bias)
print(Tmat)
