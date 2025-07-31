#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 00:11:02 2025

@author: mwu
Code source: https://www.youtube.com/watch?v=6rc_9kPCzUU
"""

import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv

def simplex_iteration(A,b,C,m:int, n:int):
    '''
    A: coefficient matrix (constraints of the linear program) m * (n+m) array
    b: right side of the conditions
    C: objective function
    m: number of rows
    n: overall number of decision values
    This function is an implementation of simplex method, giving optimal solution to a linear program
    Max C^tX
    Subject to: AX = b
    X >= 0
    
    Returns:
        X: solution to AX = b
        RC: reduced costs of X and slack variables
        Z: optimal objective value
    
    Process:
        B: m * m Basis Matrix
        NB: n * m Non-basis Matrix
    
    '''
    # Initialization
    Iteration = 0
    Z = 0
    X = np.zeros((n+m))
    XB = np.zeros((m)) # Basic variables
    CB = np.zeros((m)) # Coefficients of basic variables
    XN = np.zeros((n)) # Non-basic variables
    CN = np.zeros((n)) # Coefficients of non-basic variables
    RC = np.zeros((n+m)) 
    Basis:int = np.zeros((m)) # Take the indexes of the basic variables
    B = np.zeros((m,m))
    NB = np.zeros((m,n))
    Index_Enter = -1 # Index of entering variable = -1
    Index_Leave = -1 # Index of leaving variable = -1
    eps = 1e-12 # epsilon (for numerical precision)
    
    # Suppose starting with standard form
    for i in range(0,m):
        Basis[i] = n+i # Take the indexes of all slack variables as basis
        for j in range(0, m):
            B[i,j] = A[i, n+j] 
        for j in range(0,n):
            NB[i, j] = A[i, j]
            
    for i in range(0,n):
        '''
        print("C =", C)
        print("type(C) =", type(C))
        print("C[0] =", C[0])
        print("type(C[0]) =", type(C[0]))
        '''
        CN[i] = C[i]
        print("CN: ", CN[i])
    
    RC = C-np.dot(CB.transpose(),np.dot(inv(B), A))
    MaxRC = 0
    for i in range(0, n+m):
        if(MaxRC < RC[i]):
            MaxRC = RC[i]
            Index_Enter = i
    
    print("Basis", Basis)
    
    # Main loop
    while(MaxRC > eps):
        Iteration = Iteration + 1
        print("=> Iteration: ", Iteration)
        
        print("Index_Enter: ", Index_Enter)
        Index_Leave = -1
        MinVal = 1000000 # Minimum ratio when dividing b's by the entries in the A's
        print("Enter B: ", B)
        for i in range(0,m):
            if(np.dot(inv(B),A)[i, Index_Enter]> 0): # Eliminating all negative or zero entries
                bratio = np.dot(inv(B), b)[i]/np.dot(inv(B),A)[i, Index_Enter] # the ratio to pick the minimum val and new basic variables
                print(" bratio: ", bratio)
                if(MinVal > bratio):
                    Index_Leave = i
                    print(" Index_Leave: ", Index_Leave)
                    MinVal = bratio
                    print(" MinVal: ", MinVal)
                    
        if(Index_Leave == -1):
            print("Problem Unbounded. ")
            return Z,X,RC
        
        Basis[Index_Leave]=Index_Enter # The entering variable takes the place of the leaving variable
        print("before updated Basis", Basis)
        print(" Index_Leave: ", Index_Leave)
        
        # Reordering cols of the basis to initialization
        for i in range(m-1,0,-1):
            if(Basis[i] < Basis[i-1]):
                temp = Basis[i-1]
                Basis[i-1]=Basis[i]
                Basis[i] = temp # Switch the positions
                
        print("updated Basis", Basis)
        
        for i in range(0, m):
            for j in range(0, n+m):
                if(j==Basis[i]):
                    B[:,i] = A[:, j]
                    CB[i] = C[j]
        
        print("Exit Basis", Basis)
        print("Exit B: ", B)
        
        RC = C-np.dot(CB.transpose(), np.dot(inv(B),A)) # Update RC
        MaxRC= 0
        for i in range(0, n+m):
            if(MaxRC < RC[i]):
                MaxRC = RC[i]
                Index_Enter = i
        print("MaxRC", MaxRC)
        X=np.dot(inv(B), b)
        Z=np.dot(CB, X)
    return Z,X,RC

# Example 1
A = np.array([[1,1,1,3,1,1,0,0],
              [1,4,1,3,1,0,1,0],
              [1,2,1,4,1,0,0,1]])
b = np.array([[1],[2],[3]])
C = np.array([5,3,4,2,3,0,0,0])
Z,X,RC = simplex_iteration(A, b, C, 3, 5)



# Example 2
A = np.array([[1, -2, 3, 1, 1, 0, 0, 0, 0],
              [1, 0, 1, 0, 0, 1, 0, 0, 0],
              [4, 9, 1, 4, 0, 0, 1, 0, 0],
              [2,  2, 1, 1, 0, 0, 0, 1, 0],
              [2, -1, 5, 0, 0, 0, 0, 0, 1]])
b = np.array([[99],[40],[106],[60],[170]])
C = np.array([20, 12, 15, 6,0,0,0,0,0])
Z,X,RC = simplex_iteration(A,b,C,5,4)


print("Z: ", Z)
print("X: ", X)
print("RC: ", RC)
    