import numpy as np
import matplotlib.pyplot as plt
import math


for N in np.arange(1, 100, 10):
    probs = []
    for n in range(N):
        omega=(math.factorial(int(N))/(math.factorial(int(n))*math.factorial(int(N-n))))
        P = (omega*(1/2)**n)/(3/2)**N
        probs.append(P)
    plt.plot (np.array(range(0,N))/N,np.array(list(probs))*N,label='%s particles'%N)
#plt.legend(loc='upper right')
    plt.xlabel('n/N')
    plt.ylabel('P(n/N)')
    plt.title('')
plt.legend(loc='upper right')
plt.show()
    

