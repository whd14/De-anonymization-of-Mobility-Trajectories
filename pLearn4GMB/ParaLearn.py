# -*- coding:utf-8 -*-


import os
import math
import numpy as np



Dmax=20;

def GaussPorb(dis, pi, sigma2):
    dis[dis>Dmax]=Dmax
    #res=1/np.sqrt(2*np.pi*sigma2.T)*pi.T*np.exp((dis**2)*(-1/2)/sigma2.T)
    #res=np.exp(-1*(dis**2)/(sigma2.T*2))/np.sqrt(2*np.pi*sigma2.T)
    res=pi.T*np.exp(-1*(dis**2)/(sigma2.T*2))/np.sqrt(2*np.pi*sigma2.T)
    return res

def main():
    fp = open('Mismatches')
    Mismatch = np.loadtxt(fp)
    N=Mismatch.shape[0]
    Mismatch = np.hstack((Mismatch,np.zeros([N,1])+Dmax))
    H=Mismatch.shape[1]
    pi=np.zeros([H,1])+1/float(H)
    sigma2=np.ones([H,1])
    gamma=np.zeros([N,H])
    likeli=np.zeros([N,H])
    

    Niter=200
    for iter in range(1,Niter):
        for s in range(0,N):
            likeli[s,:]=GaussPorb(Mismatch[s,:],pi,sigma2)
            gamma[s,:]=likeli[s,:]/np.sum(likeli[s,:])
        Ngamma=np.sum(gamma,axis=0)
        pi=Ngamma/np.sum(Ngamma)
        sigma2=np.sum(np.multiply(gamma,Mismatch**2),axis=0)/Ngamma
    
    fw=open('para','w')
    fw.write('pi=['+','.join([str(x) for x in pi])+']\n')
    fw.write('sigma2=['+','.join([str(x) for x in sigma2])+']\n')


if __name__ == '__main__':
    main()
