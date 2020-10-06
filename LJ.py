import numpy as np
import sys
import os
sys.path.append('/Users/user/GT')
import codes




#        self.ngr = 0 # number of times the gr is sampled, (but if Im putting it into an md loop doesn't it already get resamled? or does it have to do resampling within resampling)



class LJ():

    def __init__(self, tequil, delt, a, Natoms, Nstep, Nprint, density, temp=0.728, Ntemp=100, resample=False):


        self.Natoms = Natoms
        self.Nstep = Nstep
        self.Nprint = Nprint
        self.temp = temp
        self.beta = 1.0/ ( temp )
        self.density = density
        self.resample = resample
        self.Ntemp = Ntemp
        self.a = a # (cut-off distance)
        self.tequil = tequil
        self.Ndim = 3
        self.rcut2 = a**2
        self.delt = delt

        #box size
        self.box = ( Natoms/density ) ** (1/self.Ndim)

        #GR
        self.nhis = 100
        self.ngr = 0
        self.hist = np.zeros([self.nhis, 2])
        self.delg = self.box / ( 2 * self.nhis )

        #initialize
        self.init()

        #output
        self.file_output = open( 'output.dat', 'w' )
        #################################################

    def kernel (self):
        print('******* RUNNING MD CALCULATION ********')
        print('Running MD for Lenard Jones fluid for ',self.Nstep,'steps')
        print('*****************************************')
        print()

        for step in range( self.Nstep ):
            currtime = step * self.delt

            if( self.resample ):
                if( np.mod(step, self.Nprint ) == 0):

                    self.sample_vel(self.beta )
            self.velocity_verlet()

            if( np.mod( step, self.Nprint ) == 0 and currtime >= self.tequil ):
                print( 'Writing data at step', step, 'an time', currtime )
                self.calc_data(step, currtime )

            #self.velocity_verlet()

        step += 1
        curtime = step*self.delt
        print( 'Writing data at step', step, 'and time', currtime )
        self.calc_data(step, currtime)

        self.calc_gr()
        ###################################################
    #def gr_sample( self,  ):
#
 #       ngr += 1
   #     for i in range(self.Natoms-1)
  #          for j in range( i+1, self.Natoms ):
    #            x

   # def gr_result( self ):

        ###################################################   

    def velocity_verlet( self ):

        self.get_velocities()
        self.get_positions()
        self.get_forces()
        self.get_velocities()

        ###################################################

    def get_velocities( self ):

        self.vvv += 0.5 * self.fff * self.delt
#        print('velocity')
        ###################################################

    def get_positions( self ):

        self.xxx += self.vvv * self.delt

 #       print('position')
        ###################################################
    def get_forces( self ):

        self.fff = np.zeros( [ self.Natoms, self.Ndim ] )

        for i in range(self.Natoms-1 ):
            for j in range( i+1, self.Natoms ):
                xr = self.xxx[i] - self.xxx[j]
                xr = xr - self.box * np.rint( xr / self.box )
                r2 =np.sum( xr**2 )
                r = np.sqrt(r2)
            #    if (r2 < self.rcut2) and (xr.all() != 0):
                if (r2 < self.rcut2):
                    r2i = 1/r2
                    r6i = r2i**3
                    ff = 48 * r2i * r6i * (r6i - 0.5)
                    self.fff[i, :] += ff*xr
                    self.fff[j, :] -= ff*xr

        ##################################################

    def calc_gr(self):
        for i in range(self.nhis):
            self.hist[i,0] = self.delg*(i+0.5)
            vb = ((i+1)**3 - i**3) * self.delg**3
            nid = (4/3) * np.pi * vb * self.density
            self.hist[i,1] = self.hist[i,1]/(self.ngr * nid * self.Natoms)
            codes.printarray(self.hist, 'gr.dat', True)
           # print('calc_gr works')
        ##################################################

    def sample_vel( self, beta ):

        #from boltzman distribution
        sigma = np.sqrt(1.0/(beta))
        self.vvv = np.random.normal(0.0, sigma, self.Natoms*self.Ndim).reshape( self.Natoms, self.Ndim )

        #eliminate center of mass motion 

        center_of_mass = np.sum( self.vvv, axis=0)/ self.Natoms #axis=0 sum over column
        self.vvv -= center_of_mass

        #re-scale velocities to proper temperatur


        scale = np.sqrt( self.temp / ( np.sum( self.vvv**2 ) / ( self.Ndim*self.Natoms ) ) )
        self.vvv *= scale

        #calculate temperature 
        #velocity at the center of mass is zero 
        ##################################################

    def get_pe( self ):
        enP = 0.0
        for  i in range(self.Natoms-1):
            for j in range( i+1, self.Natoms):
                xr = self.xxx[i] - self.xxx[j]
#                print('potential xr', xr)
                xr = xr - self.box * np.rint( xr/self.box)
                r2 = np.sum( xr**2 )
            #    if (r2 < self.rcut2) and (xr.all() != 0):
                if (r2 < self.rcut2):
                    r2i = 1/r2
                    r6i = r2i**3
                    enP += 4 * r6i * (r6i - 1) - self.ecut
#                    print('pe', enP)
 #               else:
  #                  enP += 0
        return enP

        ####q###############################################

    def get_ke( self ):
        ke = 0.5 * np.sum(self.vvv**2)
       # print ('kinetic E', ke)
        return ke
        ###################################################

    def calc_data( self, step, currtime ):
    #g(r) sampling
        self.ngr += 1
        for i in range(self.Natoms-1 ):
            for j in range( i+1, self.Natoms ):
                xr = self.xxx[i] - self.xxx[j]
                xr = xr - self.box * np.rint( xr / self.box )
                r2 =np.sum( xr**2 )
                r = np.sqrt(r2)

                if(r < self.box/2):
                    ig = int(np.floor( r / self.delg )) #why floor???
                    self.hist[ig,1] +=2
                   # print('gr completed')
    # output energies
        fmt_str = '%20.8e'
        engpe = self.get_pe()
        engke = self.get_ke()
        etot = engpe + engke
       # print('etot', etot)
        output = np.zeros(4)
        output[0] = currtime
        output[1] = etot
        output[2] = engpe
        output[3] = engke

        np.savetxt( self.file_output, output.reshape(1, output.shape[0]), fmt_str )
        self.file_output.flush()

        ####################################################

    def init( self ):

        self.xxx = np.zeros( [ self.Natoms, self.Ndim ] )

        N1d = self.Natoms**(1/self.Ndim)
        N1d = int(round(N1d))

        delL = self.box / N1d

        cnt = 0

        for x in  range(N1d):
            for y in range(N1d):
                for z in range(N1d):
                    self.xxx[ cnt, 0] = x*delL
                    self.xxx[ cnt, 1] = y*delL
                    self.xxx[ cnt, 2] = z*delL
                    cnt += 1
        #print('positions', self.xxx )
        self.sample_vel( self.beta )
        self.ecut = 4 * (1/self.rcut2**6 - 1/self.rcut2**3)
        self.get_forces()








