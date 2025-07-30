#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 00:11:02 2025

@author: mwu
"""

import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv

def simplex_iteration(A,b,c,m:int, n:int):
    