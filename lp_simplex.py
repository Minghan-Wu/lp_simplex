#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 00:11:02 2025

@author: mwu
"""

import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv

def simplex_iteration(A,b,C,m:int, n:int):
    Iteration = 0
    Z = 0
    X = np.zeros((n+m))
    XB = np.zeros((m))
    CB = np.zeros((m))
    XN = np.zeros((n))
    CN = np.zeros((n))
    RC = np.zeros((n+m))
    Basis:int = np.zeros((m))
    B = np.zeros((m,m))
    NB = np.zeros((m,n))
    Index_Enter = -1
    Index_Leave = -1
    eps = 1e-12
    
    for i in range(0,m):
        Basis[i] = n+i
        for j in range(0, m):
            B[i,j] = A[i, n+j]
        for j in range(0,n):
            NB[i, j] = A[i, j]
            
    for i in range(0,n):
        CN[i] = C[i]
        print("CN: ", CN[i])
    
    RC = C-np.dot(CB.transpose(),np.dot(inv(B), A))
    MaxRC = 0
    for i in range(0, n+m):
        if(MaxRC < RC[i]):
            MaxRC = RC[i]
            Index_Enter = i
    
    print("Basis", Basis)
    while(MaxRC > eps):
        Iteration = Iteration + 1
        print("=> Iteration: ", Iteration)
        
        print("Index_Enter: ", Index_Enter)
        Index_Leave = -1
        MinVal = 1000000
        print("Enter B: ", B)
        for i in range(0,m):
            if(np.dot(inv(B),A)[i, Index_Enter]> 0):
                bratio = np.dot(inv(B), b)[i]/np.dot(inv(B),A)[i, Index_Enter]
                print(" bratio: ", bratio)
                if(MinVal > bratio):
                    Index_Leave = i
                    print(" Index_Leave: ", Index_Leave)
                    MinVal = bratio
                    print(" MinVal: ", MinVal)
        if(Index_Leave == -1):
            print("Problem Unbounded. ")
            return Z,X,RC
        Basis[Index_Leave]=Index_Enter
        print("before updated Basis", Basis)
        print(" Index_Leave: ", Index_Leave)
        for i in range(m-1,0,-1):
            if(Basis[i] < Basis[i-1]):
                temp = Basis[i-1]
                Basis[i-1]=Basis[i]
                Basis[i] = temp
                
        print("updated Basis", Basis)
        
        for i in range(0, m):
            for j in range(0, n+m):
                if(j==Basis[i]):
                    B[:,j] = A[:, j]
                    CB[i] = C[j]
        
        print("Exit Basis", Basis)
        print("Exit B: ", B)
        
        RC = C-np.dot(CB.transpose(), np.dot(inv(B),A))
        MaxRC= 0
        for i in range(0, n+m):
            if(MaxRC < RC[i]):
                MaxRC = RC[i]
                Index_Enter = i
        print("MaxRC", MaxRC)
        X=np.dot(inv(B), b)
        Z=np.dot(CB, X)
    return Z,X,BC
    
    