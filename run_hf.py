import make_ham
import hf

Hcore, S, G, Enuc = make_ham.make_ham_sto3g('h2', 1.4, 3) #type, radious and STO-3G

itrmax =900
Nu =1 #up electrons
Nd =1 #down electrons
Ns =2 #spatial orbitalsi
orbs = None
Eorb_u = None
Eorb_d = None
orb_u = None
orb_d = None
Enew = 0
#Eold = 0
P = None

hartree_f = hf.hartree_fock(S, Hcore, G, itrmax, Nu, Nd, Ns, orbs,Eorb_u, Eorb_d,orb_u,orb_d, Enew,P )
hartree_f.kernel()
