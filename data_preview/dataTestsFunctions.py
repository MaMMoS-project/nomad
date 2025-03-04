#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:21:49 2025

@author: santapile
"""
import numpy as np


def read2coldata(filename):
    """
    Returns a tuple of two numpy-arrays where the first array
    contains all the data from the first column in the file,
    and the second array contains all the data in the second column
    of the data in the file.

    """
    alldata = np.loadtxt(filename)
    return alldata[:, 0], alldata[:, 1]


def f(x, a, b, c):
    return a * np.exp(-b * x) + c


def spontaneous_magnetization(
    T: np.ndarray, M_0: float, s: float, T_C: float, p: float, beta: float
) -> np.ndarray:
    """Formula 4 from
    https://link.springer.com/article/10.1140/epjp/s13360-020-00294-y
    """
    return M_0 * (1 - s * (T / T_C) ** (3 / 2) - (1 - s) * (T / T_C) ** p) ** beta


def findlinevaldict(fileName, valname):
    """
    Find all lines in lines (list) with valname (string) and
    return a dict of IDs and valuse coming after valname
    """
    with open(fileName, 'rt') as f:
        lines = f.read().splitlines()
    # print(lines)
    alllineswithvalname = {}
    for line in lines:
        if valname in line:
            # print(line)
            pos = line.find(valname) + len(valname)
            key, value = line.split()[0], [float(x) for x in line[pos:].split()]
            alllineswithvalname[key] = value
    return alllineswithvalname

def getEnergieFromFile(fileName): # reading file into lines
    with open(fileName, 'rt') as f:
        lines = f.read().splitlines()
        
    valname = "Moment"
    pos = 0
    for line in lines:
        pos += 1
        if valname in line:
            posE = pos + 1

    return float(lines[posE].split()[-1])