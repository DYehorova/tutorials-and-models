import numpy as np
import math
import sys
import os
import utils
sys.path.append('/Users/user/GT')
import codes
class harm_osc():

    def __init__( self, m, k, step, Nprint, Ncycle, temp, xi):

        self.m = m
        self.k = k
        self.Ncycle = Ncycle
        self.Nprint = Nprint
        self.step = step
        self.beta = 1.0 / (temp /  3.1577465e5 )

        self.positions = np.zeros(self.Ncycle+1)
        self.positions[0] = xi
        self.xi = xi
        self.file_output = open( 'output.dat', 'w')

###########################################################

    def kernel( self ):
        print('****** RUNING MC CALCULATION ******')
        print('Running MC for 1D Harmonic Oscillator for ', self.Ncycle, 'step')
        print()

        print(self.xi)
        xo= self.xi
        print(xo)
        old_pe = self.pe( xo )
        print('old pE is', old_pe)

        for cycle in range( self.Ncycle ):
   #         if(np.mod( cycle, Nprint ) == 0):
            print( 'Writing data at cycle', cycle )

            print('old pE is', old_pe)
#random displacement
            random_1=np.random.rand()
            random_2=np.random.rand()
            print('random_1 is', random_1, 'random_2 is', random_2)

            xn = xo + ( random_1 - 0.5) * self.step  #determine step size for a one atom system (~20% acceptance?) (kbT/f(u))per unit of cpu time
            print('new position is', xn)

            new_pe = self.pe( xn )
            print('new pE is', new_pe)

            chk = np.exp( - self.beta * (new_pe - old_pe))
            print('chk is', chk)
            if  ( chk > random_2): #no need to include 'and != 1' because fcn=1 only if x_new=xi 
                xo = xn
                old_pe = new_pe
            else:
                xo = xo
                old_pe = old_pe
#calculate last cycle
            self.positions[ cycle+1 ] = xo

            if( np.mod( cycle, 100 ) == 0 ):
                print( 'Writing data at cycle', cycle )

        print()
        print('calculating last cycle')
        print()
        self.calc_data()

###########################################################

    def move( self, beta ):


        print('old pE is', old_pe)
#random displacement
        random_1=np.random.uniform(0,1)
        random_2=np.random.uniform(0,1)
        print('random_1 is', random_1, 'random_2 is', random_2)

        xn = xo + ( random_1 - 0.5) * self.step  #determine step size for a one atom system (~20% acceptance?) (kbT/f(u))per unit of cpu time
        print('new position is', xn)

        new_pe = self.get_pe(xn)
        print('new pE is', new_pe)

        chk = np.exp( - self.beta * (new_pe - old_pe))
        print('chk is', chk)
        if  ( chk > random_2): #no need to include 'and != 1' because fcn=1 only if x_new=xi 
            xo = xn

        else:
            xo = xi

###########################################################

 #   def sample( self ):
#output energy for each move or avarage of all moves? 
 #       self.get_pe(self.xxx)
 #       positions.append(self.xxx)#do 'positions' have to be a self.obj? 

###########################################################    

    def init( self ):
        positions = []
        positions.append(self.xxx)
#do i need it to be sample then move so I record the initial position in the list? 
###########################################################


    def pe(self, xxx):
        print(self.k, xxx)
        return 0.5 * self.k * xxx**2

###########################################################

    def calc_data( self ):

        #split data into chunks and calculate P(x) for each chunk
        Nchunk = 5
        Nbins  = 100
        histo = np.zeros( [Nbins,Nchunk+1] ) #dont know
        #positions_array = np.array(positions)
        split_pos_array = np.split( self.positions, Nchunk )

        for i in range(Nchunk):

            if( i == 0):
             #calculate midpoint of each bin
                bin_edges = np.histogram( split_pos_array[i], bins=Nbins, range=( self.positions.min(), self.positions.max() ), density=True )[1]
                for j in range( Nbins ):
                    histo[j,0] = ( bin_edges[j] + bin_edges[j+1] ) / 2.0
             #calculate normalized histogram for each chunk, note that all histograms have the same range
            histo[:,i+1]   = np.histogram( split_pos_array[i], bins=Nbins, range=( self.positions.min(), self.positions.max() ), density=True )[0]#why is it 0 or 1
         #calculate average and error over all instances of P(X)
        avg_histo = np.zeros( [Nbins,3] )
        avg_histo[:,0] = np.copy( histo[:,0] )
        avg_histo[:,1] = np.mean( histo[:,1:], axis=1 )
        avg_histo[:,2] = np.std( histo[:,1:], axis=1 ) / np.sqrt(Nchunk-1)

         #print out P(x)
        codes.printarray(avg_histo,'prob_x.dat')

############################################################

