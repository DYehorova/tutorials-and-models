import numpy as np
import sys
import os
sys.path.append('/Users/user/GT')
from codes import diagonalize
#sys.path.append('/Users/user/GT/HF')
#import make_ham


class hartree_fock():

    def __init__(self, S, Hcore, G, itrmax, Nu, Nd, Ns, orbs, Eorb_u, Eorb_d, orb_u, orb_d, Enew, P, restrict = True):

        self.Ns = Ns #number of spacial orbitals
        self.Nu = Nu #number of up electrons 
        self.Nd = Nd #number of down electrons
        self.Hcore = Hcore
        self.S = S
        self.G = G
        self.restrict = restrict
        self.itrmax =itrmax
        self.Eorb_u = Eorb_u
        self.orbs = orbs
        self.Eorb_d = Eorb_d
        self.orb_u = orb_u
        self.orb_d = orb_d
        self.Enew = Enew
        self.P = P
       # self.init_P()
       # self.file_output =open( 'output_hf.dat', 'w')

    ################################

    def kernel (self):
        print('******RUNNING HARTREE_FOCK CALCULATION******')
        print('Running HF for', self.itrmax, 'steps')
        print('********************************************')


       #initiate P 
        eorb,self.orbs = diagonalize(self.Hcore, self.S) #orbitals=eigonvectors, orbital energies =eigenvalues

        Pu = self.rdm_1el(self.orbs, self.Nu )
        if( self.restrict ):
            Pd = self.rdm_1el(self.orbs, self.Nd )#density matrix for down electronsi
        else:
            Pd = (self.Nd*1.0/self.Ns)*np.identity(self.Ns)

        self.P = Pu + Pd
        print('p',self.P)
        print('G',self.G)
        itr = 0
        Enew = 9999.9
        Ediff = 10.0
        while itr<self.itrmax and Ediff > 1e-8:

            print('Iteration',  itr)

                # calculate G (2 electron contribution)
            h2el_u = self.make_h2el(self.G,self.P,Pu)
            h2el_d = self.make_h2el(self.G,self.P,Pd)
                #
            print('h2el_u',h2el_u)
            print('h2el_d',h2el_d)

            #fock operator
            fock_u = self.Hcore + h2el_u
            fock_d = self.Hcore + h2el_d
            print('focka', fock_u, 'fockb', fock_d)
             #or with an einsum
               # fock_u
                #solve fock operators
            self.Eorb_u,self.orbs_u = diagonalize(fock_u,self.S)
            self.Eorb_d,self.orbs_d = diagonalize(fock_d,self.S)

                #make a new density matrix
            Pu = self.rdm_1el(self.orbs_u, self.Nu)
            Pd = self.rdm_1el(self.orbs_d, self.Nd)
            P = Pu + Pd

                #get the energy, check convergence
            Eold = self.Enew
            self.Enew =  0.5*np.trace(np.dot(P,self.Hcore)+np.dot(Pu,fock_u)+np.dot(Pd,fock_d))#Change into eignsum!
            Ediff = np.fabs(Eold-self.Enew)

            itr += 1

        print( "Energy= ", self.Enew )
        print( "Evalsa= " )
        print( self.Eorb_u )
        print( "Evalsb= " )
        print( self.Eorb_d )
        print( "Orbsa= " )
        print( self.orbs_u )
        print( "Orbsb= " )
        print( self.orbs_d )
###########################################################
#1 electron density matrix, C =spatial  orbitals, Ne = # of el
    def rdm_1el(self,C, Ne):
        Coc = C[:,:Ne] #occupied orbitals
        P = np.dot(Coc, np.transpose(np.conjugate(Coc)))
        return P
###########################################################

    def make_h2el(self, G,P,Pu):
        #the two-electron contribution to the fock matrix

       # h2el = (np.tensordot(P,G,axes =([0,1],[3,2])) - np.tensordot(Pu, G, axes =([0,1],[1,2])))#how the axis work and change into eignsum 3.241,242
        h2el = (np.einsum('pq,rsqp->rs', P, G)-(np.einsum('pq,rsqp->rs', Pu, G))) # why subtracting a tensordot Pa?

        return h2el
