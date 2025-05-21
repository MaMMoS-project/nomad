#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 10:21:49 2025

@author: santapile
"""

import numpy as np
import os
from functools import partial
from scipy.optimize import curve_fit

from scipy.constants import mu_0
from scipy.constants import physical_constants
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({"font.size": 14})


def read_2_col_data(filename):
    """
    Returns a tuple of two numpy-arrays where the first array
    contains all the data from the first column in the file,
    and the second array contains all the data in the second column
    of the data in the file.

    """

    try:
        alldata = np.loadtxt(filename)
    except ValueError:
        alldata = np.loadtxt(filename, skiprows=1)
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
    Find last line in lines (list) with valname (string) and
    return a dict of IDs and valuse coming after valname
    """
    # Check if file exists
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"The file '{fileName}' does not exist.")

    with open(fileName, "rt") as f:
        lines = f.read().splitlines()

    if verbose:
        print(lines)

    last_lines_with_valname = {}
    for line in lines:
        if valname in line:
            if verbose:
                print(line)
            pos = line.find(valname) + len(valname)
            # TODO: check if part of the line can be converted to float;
            # introduce boundaries in which the value should be
            key, value = line.split()[0], [float(x) for x in line[pos:].split()]
            last_lines_with_valname[key] = value
    return last_lines_with_valname


def get_energy_from_file(fileName):
    """
    Extracts energy value from a file.

    This function reads a file line by line, searches for a line containing
    the keyword "Moment",
    and then retrieves the energy value from the line immediately following it.

    Args:
        fileName (str): The path to the file from which to extract the energy value.

    Returns:
        float: The extracted energy value if it is found and otherwise None.

    Raises:
        ValueError: If the keyword "Moment" is not found in the file or if
        the energy value cannot be converted to a float.
        FileNotFoundError: If the file specified by fileName does not exist.
    """
    # Check if file exists
    if not os.path.isfile(fileName):
        raise FileNotFoundError(f"The file '{fileName}' does not exist.")

    with open(fileName, "rt") as f:
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


def round_to_significant_digits(value, sig_digits):
    """
    Rounds a value to a specified number of significant digits.
    :param value: The value to round.
    :param sig_digits: The number of significant digits to round to.
    :return: The rounded value.
    """
    if value == 0:
        raise ValueError("Cannot round a value of 0 to significant digits.")
    digits = -int(np.floor(np.log10(abs(value))) - (sig_digits - 1))
    value_rounded = round(value, digits)
    return value_rounded


def get_unit_cell_volume(file_name):
    """
    Extracts the unit cell volume from the given file.
    :param file_name: The name of the file to extract the unit cell volume from.
    :return: The unit cell volume in cubic angstroms (A^3).
    """
    ucv = find_line_val_dict(file_name, "unit cell volume:")
    ucvA = ucv[list(ucv.keys())[0]][0] / 1.8897259**3  # unit cell volume in A^3
    return ucvA


def compute_magnetization(tot_moments_D, dir_of_JD, file_Name_Ms):
    # Calculating total magnetic moment by summing all
    # (Total moment (of an orbital) * Direction of J
    # (takes the one with abs value > 0.9, should be +/-1))
    tot_magn_mom_C = 0
    for key in tot_moments_D.keys():
        tot_magn_mom_C += (
            tot_moments_D[key][0] * [x for x in dir_of_JD[key] if abs(x) > 0.9][0]
        )
        print(
            key,
            ": ",
            tot_moments_D[key][0],
            "*",
            [x for x in dir_of_JD[key] if abs(x) > 0.9][0],
            "=",
            tot_moments_D[key][0] * [x for x in dir_of_JD[key] if abs(x) > 0.9][0],
        )
        #    total magnetic moment J=L+S * orientation, wherever +/-1 is

    print(f"\nTotal magnetic moment = {tot_magn_mom_C}")

    # Getting unit cell volume in A^3 from the file
    ucvA = get_unit_cell_volume(file_Name_Ms)
    print(f"Unit cell volume: {ucvA} A\N{SUPERSCRIPT THREE}")

    # Calculating magnetization in Tesla
    magnetization_in_T = tot_magn_mom_C / ucvA * 11.654

    return magnetization_in_T, ucvA


def compute_anisotropy_constant(data_dir_GS, xyz_dirs, ucvA):
    energies = {}

    if f"out_MF_{xyz_dirs[0]}" in os.listdir(data_dir_GS + f"/{xyz_dirs[0]}"):
        for dirdir in xyz_dirs:
            fileName = data_dir_GS + f"/{dirdir}/out_MF_{dirdir}"
            eigenvalue_sum = find_line_val_dict(fileName, "Eigenvalue sum:")
            energies[dirdir] = eigenvalue_sum[list(eigenvalue_sum.keys())[0]][0]
    elif f"out_Etot_{xyz_dirs[0]}" in os.listdir(data_dir_GS + f"/{xyz_dirs[0]}"):
        for dirdir in xyz_dirs:
            fileName = data_dir_GS + f"/{dirdir}/out_Etot_{dirdir}"
            energies[dirdir] = get_energy_from_file(fileName)
    else:
        print("no files for anisotropy")

    allKs = list()
    if "z" in energies.keys():
        if "x" in energies.keys():
            Kxz = (energies["x"] - energies["z"]) / ucvA * 2179874
            allKs.append(Kxz)
        if "y" in energies.keys():
            Kyz = (energies["y"] - energies["z"]) / ucvA * 2179874
            allKs.append(Kyz)

    K1_in_JPerCubibm = (
        max(allKs) * 1e6
    )  # anisotropy J/m³; MagnetocrystallineAnisotropyConstantK1
    return K1_in_JPerCubibm


def extract_values_from_readme(data_dir):
    import re

    with open(data_dir + "/README", "rt") as f:
        readme_content = f.read()

    ms_match = re.search(r"Ms\s*=\s*([\d.]+)\s*T", readme_content)
    if ms_match:
        Ms_README_in_T = float(ms_match.group(1))
        print(f"Extracted Ms value from README: {Ms_README_in_T} T")
    else:
        print("Ms value not found in README")
        Ms_README_in_T = None

    mae_matches = re.findall(
        r"([xyz])\s*-\s*([xyz])\s*=\s*([-+]?\d*\.\d+|\d+)\s*MJ/m³", readme_content
    )
    if mae_matches:
        mae_values_in_MJPerCubicm = [float(match[2]) for match in mae_matches]
        max_MAE_README_in_MJPerCubicm = max(mae_values_in_MJPerCubicm)
        print(
            f"Extracted maximum MAE value from README: {max_MAE_README_in_MJPerCubicm} MJ/m³"
        )
    else:
        mae_match = re.search(r"MAE\s*=\s*([\d.]+)\s*MJ/m³", readme_content)
        if mae_match:
            max_MAE_README_in_MJPerCubicm = float(mae_match.group(1))
            print(
                f"Extracted MAE value from README: {max_MAE_README_in_MJPerCubicm} MJ/m³"
            )
        else:
            print("MAE values not found in README")
            max_MAE_README_in_MJPerCubicm = None

    return Ms_README_in_T, max_MAE_README_in_MJPerCubicm


def compute_exchange_and_anisotropy_constants(
    data_dir_MC, tot_moments_D, K1_in_JPerCubibm, ucvA, plot_Js=True
):
    # Define the form of the function you want to fit
    fn = data_dir_MC + "/M(T)"  # magnetic polarization
    TK, Js = read_2_col_data(fn)
    print(TK)

    mc_dirs = [
        dirdir for dirdir in os.listdir(data_dir_MC) if dirdir in ["momfile", "posfile"]
    ]

    if len(mc_dirs) > 0:
        with open(data_dir_MC + "/" + mc_dirs[0], "rt") as f:
            n_atoms = len(f.read().splitlines())
        print(f"Number of atoms: {n_atoms}")
    else:
        print(
            "There seem to be NO momfile and no posfile :( Getting n_atoms from GS/out_last"
        )
        atoms_Ns = set(
            [
                int(key[key.find(":") + 1 : key.find(":") + 3])
                for key in tot_moments_D.keys()
            ]
        )
        n_atoms = len(atoms_Ns)
    print(f"Number of atoms: {n_atoms}")

    Js = [item * n_atoms / ucvA * 11.654 for item in Js]

    poscut = np.argmin(np.diff(Js) / np.diff(TK)) + 2
    try:
        Tc = TK[poscut]
    except IndexError:
        Tc = TK[-1]
        poscut = len(TK) - 1
    print(f"Tc = {Tc} K")
    TKc = TK[:poscut].copy()
    Jsc = Js[:poscut].copy()

    xfine = np.linspace(0, Tc, 500)
    p = 5.0 / 2
    beta = 1.0 / 3
    m_s = partial(spontaneous_magnetization, p=p, beta=beta, T_C=Tc)

    popt, pcov = curve_fit(m_s, TKc, Jsc)
    Js_0, s = popt
    print(Js_0, s)
    # T_fit = np.linspace(min(TKc), max(TKc), 500)
    # Js_fit = m_s(T_fit, Js_0, s)

    g = 2
    k_b = physical_constants["Boltzmann constant"][0]
    mu_b = physical_constants["Bohr magneton"][0]

    M_0 = Js_0 / mu_0
    D = 0.1509 * ((g * mu_b) / (s * beta * M_0)) ** (2.0 / 3) * k_b * Tc
    print("Spin wave stiffness constant ", D)
    A_0 = M_0 * D / (2 * g * mu_b)
    print("Exchange constant A at T=0 (J/m) : ", A_0)

    # Magnetic polarization at 300 K
    Js_300 = m_s(300.0, Js_0, s)
    print("Js_300 (T) :", Js_300)

    A_300 = A_0 * (Js_300 / Js_0) ** 2
    print("A_300 (J/m) :", A_300)

    K_300 = K1_in_JPerCubibm * (Js_300 / Js_0) ** 3
    print("K_300 (MJ/m^3) :", K_300 / 1e6)

    # exchange length in nm
    print("Lex ", np.sqrt(mu_0 * A_300 / (Js_300 * Js_300)) / 1e-9)

    if plot_Js:
        fig, ax = plt.subplots()
        ax.scatter(
            TK,
            Js,
            marker="o",
            label="data points",
            facecolors="none",
            edgecolors="#4575b4",
        )
        label = "Kuz'min's fit"
        ax.plot(xfine, m_s(xfine, Js_0, s), label=label, color="#f46d43")
        ax.legend()
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Magnetic polarization J (T)")
        ax.grid()

    return A_0, A_300, K_300, Js_300
