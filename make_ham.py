import numpy as np
from scipy.special import erf

def make_ham_sto3g(typ, R, ng):
#page 18, table 3.8
    zeta_he_sq = 2.0925**2
    zeta_h_sq = 1.24**2

#first atom
    if typ=="h2":
        z1 = 1.0
        zeta_1_sq = zeta_h_sq
    else:
        z1 = 2.0
        zeta_1_sq = zeta_he_sq

#second atom
    z2 = 1.0
    zeta_2_sq = zeta_h_sq

#nuclear repulsion 
    Enuc = (z1*z2)/R
#initialize  position vectors 
    Rvec = np.zeros(2)
    zvec = np.zeros(2)
    Rvec[0] = 0.0
    Rvec[1] = R
    zvec[0] = z1
    zvec[1] = z2

#contraction coefficients and exponents for the gausian sto3g
#(ф_GF(а)=(2a/п)^(3/4)e^(-ar^2) --> ф_CGF(r)=sum_n(d_n*ф_n_GF(a))
    alpha = np.zeros((ng, 2))
    d_coef = np.zeros((ng, 2))
#numbers from page 185  a, d for 1S orbital
    d_coef[0,0] = d_coef[0,1] = 0.444635
    d_coef[1,0] = d_coef[1,1] = 0.535328
    d_coef[2,0] = d_coef[2,1] = 0.154329

    alpha[0,0] = zeta_1_sq*0.109818
    alpha[0,1] = zeta_2_sq*0.109818
    alpha[1,0] = zeta_1_sq*0.405771
    alpha[1,1] = zeta_2_sq*0.405771
    alpha[2,0] = zeta_1_sq*2.227660
    alpha[2,1] = zeta_2_sq*2.227660

#normalization of d coef for gaussian basis fcns
    for i in range(ng):
        for j in range(2):
            d_coef[i,j] = d_coef[i,j]*(2*alpha[i,j]/np.pi)**(3./4.)

    K= 2 #number of basis K
    S= np.zeros((K,K))
    T= np.zeros((K,K))#kinetic E
    V= np.zeros((K,K))#nuc-e atraction
    Hcore= np.zeros((K,K))
    G= np.zeros((K,K,K,K))


    for i in range(K): #    why not range? (range - list, arange - array, but arent we choosing 0,1)
        for j in range(K):
            S[i,j] = get_S(i,j,Rvec,ng,alpha,d_coef)
            T[i,j] = get_T(i,j,Rvec,ng, alpha, d_coef)#since they are symmetric does it make sence to try making them triangular and then make equal around diagonal? 
            V[i,j] = get_V(i, j, Rvec, ng, alpha, d_coef, zvec)
            for k in range(K):
                for l in range(K):
                    G[i,j,k,l] = get_G(i,j,k,l,Rvec,ng,alpha,d_coef)
    Hcore = V+T
    return Hcore, S, G, Enuc
#    return Hcore, S, G, T, V
####################################
def get_T(mu, nu, R, ng, alpha, d_coef):
    Rmunu=R[mu]-R[nu]
    Tmunu=0.0
    for i in range(ng):
        for j in range(ng):
            a=alpha[i,mu]
            b=alpha[j,nu]

            T = a*b/(a+b)*(3-2*a*b/(a+b)*Rmunu**2)*(np.pi/(a+b))**(3/2)*np.exp(-a*b/(a+b)*Rmunu**2)
            Tmunu += d_coef[i,mu]*d_coef[j,nu]*T

    return Tmunu
#####################################
def get_V(mu, nu,R, ng, alpha, d_coef, zvec):

    Natoms = len(R)

    Rmunu = R[mu]-R[nu]
    Vmunu=0.0
    for i in range(ng):
        for j in range(ng):
            a=alpha[i,mu]
            b=alpha[j,nu]

            Rp=(a*R[mu]+b*R[nu])/(a+b)

            for k in range(Natoms):
                #position and charge of nuclei
                Z = zvec[k]
                Rnuc = R[k]

                V = -2*np.pi/(a+b)*Z*np.exp(-a*b/(a+b)*Rmunu**2)*F((a+b)*(Rp-Rnuc)**2)
                Vmunu += d_coef[i,mu]*d_coef[j,nu]*V
    return Vmunu

#####################################

def F(x):
    if x<1e-3:
        return 1.0
    else:
        return 0.5*np.sqrt(np.pi/x)*erf(np.sqrt(x))
#try with and without limit

#####################################

def get_S(mu, nu,R, ng, alpha, d_coef):
    Smunu=0.0
    Rmunu = R[mu]-R[nu]
    for i in range(ng):
        for j in range(ng):
            if mu == nu:
                Smunu = 1
            else:
                a = alpha[i,mu]
                b = alpha[j,nu]

                S = (np.pi/(a+b))**1.5*np.exp(-a*b/(a+b)*Rmunu**2)
                Smunu += d_coef[i,mu]*d_coef[j,nu]*S

    return Smunu

#####################################
def get_G(i,j,k,l, R,ng, alpha,d_coef):
    #p,q,r,s represent basis function indices and i,j,k,l represent primitive gaussian indices

    Ri = R[i]
    Rj = R[j]
    Rk = R[k]
    Rl = R[l]
    Gijkl = 0.0
    for m in range(ng):

        a = alpha[m,i]
        d_a = d_coef[m,i]

        for n in range(ng):

            b = alpha[n,j]
            d_b = d_coef[n,j]

            Rppp = (a*Ri+b*Rj)/(a+b)

            for o in range(ng):

                c = alpha[o,k]
                d_c = d_coef[o,k]

                for p in range(ng):

                    d = alpha[p,l]
                    d_d = d_coef[p,l]

                    Rqqq = (c*Rk+d*Rl)/(c+d)

                    Gmnop = 2*np.pi**(5/2)/((a+b)*(c+d)*np.sqrt(a+b+c+d))*(np.exp(-a*b/(a+b)*(Ri-Rj)**2-c*d/(c+d)*(Rk-Rl)**2)*F((a+b)*(c+d)/(a+b+c+d)*(Rppp-Rqqq)**2))
                    Gijkl += d_a*d_b*d_c*d_d*Gmnop
    print('G',Gmnop)
    print('Gijkl',Gijkl)
    return Gijkl
