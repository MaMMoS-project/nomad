#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:21:49 2025

@author: santapile
"""
import numpy as np
import os


def read_2_col_data(filename):
    """
    Returns a tuple of two numpy-arrays where the first array
    contains all the data from the first column in the file,
    and the second array contains all the data in the second column
    of the data in the file.

    """
    alldata = np.loadtxt(filename)
    return alldata[:, 0], alldata[:, 1]


def spontaneous_magnetization(
    T: np.ndarray, M_0: float, s: float, T_C: float, p: float, beta: float
) -> np.ndarray:
    """Formula 4 from article "Exchange stiffness of ferromagnets" by Kuz’min et al.
    See: https://link.springer.com/article/10.1140/epjp/s13360-020-00294-y
    """
    return M_0 * (1 - s * (T / T_C) ** (3 / 2) - (1 - s) * (T / T_C) ** p) ** beta


def find_line_val_dict(fileName, valname, verbose=False):
    """
    Find all lines in lines (list) with valname (string) and
    return a dict of IDs and valuse coming after valname
    """
    # Check if file exists
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"The file '{fileName}' does not exist.")

    with open(fileName, 'rt') as f:
        lines = f.read().splitlines()

    if verbose:
        print(lines)

    alllineswithvalname = {}
    for line in lines:
        if valname in line:
            if verbose:
                print(line)
            pos = line.find(valname) + len(valname)
            key, value = line.split()[0], [float(x)
                                           for x in line[pos:].split()]
            alllineswithvalname[key] = value
    return alllineswithvalname


def get_energy_from_file(fileName):
    """
    Extracts energy value from a file.

    This function reads a file line by line, searches for a line containing the keyword "Moment",
    and then retrieves the energy value from the line immediately following it.

    Args:
        fileName (str): The path to the file from which to extract the energy value.

    Returns:
        float: The extracted energy value if it is found and otherwise None.

    Raises:
        ValueError: If the keyword "Moment" is not found in the file or if the energy value cannot be converted to a float.
        FileNotFoundError: If the file specified by fileName does not exist.
    """
    # Check if file exists
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"The file '{fileName}' does not exist.")

    with open(fileName, 'rt') as f:
        lines = f.read().splitlines()

    valname = "Moment"
    pos = 0
    posE = None

    for line in lines:
        pos += 1
        if valname in line:
            posE = pos + 1

    if posE is None:
        # raise ValueError(f"Keyword '{valname}' not found in file '{fileName}'.")
        return None
    else:
        energy = float(lines[posE].split()[-1])
        return energy
