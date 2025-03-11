import os
import numpy as np
import pandas as pd

import data_tests_functions as dtf


print('Importing done.')


data_dirs = []
positive_datasets = {}
structure_check_datasets = {}
# data_dir = '/Users/santapile/Santas/Projects/MaMMoS/NOMAD/dataExamples/UU/UU-sets/Fe3-xZnxY(x=0.22)'

data_dirs.append('Co2Fe2H4')
data_dirs.append('Co2Fe16Y6')
data_dirs.append('Fe3-xTaxY(x=0.22)')
data_dirs.append('Fe3-xVxY(x=0.22)')
data_dirs.append('Fe3-xZnxY(x=0.22)')


for data_dir in data_dirs:
    import check_file_and_folder_structure as cffs
    structure_check_datasets[data_dir] = cffs.check_structure(data_dir, check_README=True, verbose=True)

    data_dir_GS = data_dir + "/GS"
    data_dir_MC = data_dir + "/MC"

    xyz_dirs = [dirdir for dirdir in os.listdir(data_dir_GS) if len(dirdir) == 1]

    print(data_dir_GS)
    print(data_dir_MC)
    print(xyz_dirs)

    # Reading file into lines; all folders are equivalent according to UU-colleagues, so we can use the first one
    file_Name_Ms = data_dir_GS + f"/{xyz_dirs[0]}/out_last"
    print(file_Name_Ms)

    # Create dicts with orbital IDs and values of corresponding
    # Total moment and its direction (+/-1)
    tot_moments_D = dtf.find_line_val_dict(file_Name_Ms, 'Total moment [J=L+S] (mu_B):')
    dir_of_JD = dtf.find_line_val_dict(file_Name_Ms, 'Direction of J (Cartesian):')
    # print("tot_moments_D:"+str(tot_moments_D))
    # print("dir_of_JD:"+str(dir_of_JD))
    
    # params_merged = {}
    # for key in tot_moments_D.keys():
    #     params_merged[key] = tot_moments_D[key] + dir_of_JD[key]

    # # Pretty table for checking the entries for magnetic moments and directions
    # df = pd.DataFrame.from_dict(params_merged, orient='index',
    #                             columns = ['J=L+S (Cartesian)',
    #                                     'J=L+S (Spin axis)',
    #                                     'Direction of J (x)',
    #                                     'Direction of J (y)',
    #                                     'Direction of J (z)'])
    # df

    magnetization_in_T, ucvA = dtf.compute_magnetization(tot_moments_D, dir_of_JD, file_Name_Ms)
    #print(f'Magnetization Ms: {magnetization_in_T} T')

    K1_in_JPerCubibm = dtf.compute_anisotropy_constant(data_dir_GS, xyz_dirs, ucvA)
    print(f'Anisotropy constant (max of all): {K1_in_JPerCubibm} J/m\N{SUPERSCRIPT THREE}')

    A_0, A_300, K_300, Js_300 = dtf.compute_exchange_and_anisotropy_constants(data_dir_MC, tot_moments_D, K1_in_JPerCubibm, ucvA, plot_Js=True)

    Ms_README_in_T, max_MAE_README_in_MJPerCubicm = dtf.extract_values_from_readme(data_dir)

    MAE_in_MJPerCubicm = K1_in_JPerCubibm /1e6

    Ms_in_T = magnetization_in_T    
    # TODO: magnetization_in_T should be renamed polarisation_in_T or Js_in_T (SpontaneousMagneticPolarisation)

    print(" Ms_in_T = "+str(dtf.round_to_significant_digits(Ms_in_T, 4)) + " T")
    print(" MAE_in_MJPerCubicm = " + str(dtf.round_to_significant_digits(MAE_in_MJPerCubicm, 4)) + " MJ/m^3")

    deviation_Ms_in_percent = 100 - (Ms_in_T/ Ms_README_in_T) * 100
    print(" Deviation for magnetization in percent = "+str(dtf.round_to_significant_digits(deviation_Ms_in_percent, 4)) + " %")

    try:
        deviation_mae_percent = 100 - ((MAE_in_MJPerCubicm) / max_MAE_README_in_MJPerCubicm) * 100
        print(" Deviation for anisotropy energy/constant in percent = " + str(dtf.round_to_significant_digits(deviation_mae_percent, 4)) + " %")
    except NameError:
        print("MAE values not found in README")

    if (abs(deviation_Ms_in_percent)<2.) and (abs(deviation_mae_percent)<2.):
        print('The deviation for magnetization and anisotropy in dataset '+data_dir+' is\n **less** than 2% compared to the provided values.')
        positive_datasets[data_dir] = True
    else:
        print('The deviation for magnetization and anisotropy in dataset '+data_dir+' is\n **more** than 2% compared to the provided values.')
        positive_datasets[data_dir] = False

if all(structure_check_datasets.values()):
    print('\n# The structure of all datasets is correct.')
else:
    print('\n# The structure of some datasets is incorrect.')
    print('Datasets with incorrect structure:')
    for data_dir, is_correct in structure_check_datasets.items():
        if not is_correct:
            print(data_dir)
       
if all(positive_datasets.values()):
    print('\n# The results for all datasets are within the specified boundaries.')
else:
    print('\n# The results for some datasets are outside the specified boundaries.')
    print('Datasets with deviations more than 2%:')
    for data_dir, is_positive in positive_datasets.items():
        if not is_positive:
            print(data_dir)
