'''

Functions to create LaTex tables


'''

import pandas as pd
import numpy as np 
from math import floor
from global_settings import * 


def ParseNumber(number):
    '''
    Arguments:
        - number [float]
    Retrieves the power of 10 and the base of 
    a floating point number 

    Return the base and the exponent of the power of 10 
    '''
    exponent = floor(np.log10(number))
    base = number*10**(-1*exponent)
    return [base, exponent] 


def CompareParse(input_a,input_b):
    '''
    Arguments:
        - input_a [float]
        - input_b [float]

    This function expect b to be higher than a
    if b has a different order of magnitude this
    function scalets a in order to have both numbers
    multiplied by the same power of 10.

    Returns: The base of both numbers and the power of 10 
    '''
    number_a = ParseNumber(input_a)
    number_b = ParseNumber(input_b)

    base_a = number_a[0]
    exp_a = number_a[1]
    base_b = number_b[0]
    exp_b = number_b[1]

    if(base_b < base_a):
        base_a*=(10**(exp_a - exp_b))
        exp_a = exp_b
    else:
        pass
    return [base_a, base_b, exp_a]



def XSecTable(ratio):

    if(ratio == 0.6):
        masses = mass_dic["0.6"]
        pi0_keys_run1 = run1_xsec_uncert_a["pi0"]
        eta_keys_run1 = run1_xsec_uncert_a["eta"]
        pi0_keys_run3 = run3_xsec_uncert_a["pi0"]
        eta_keys_run3 = run3_xsec_uncert_a["eta"]
    else:
        masses = mass_dic["2.0"]
        pi0_keys_run1 = run1_xsec_uncert_b["pi0"]
        eta_keys_run1 = run1_xsec_uncert_b["eta"]
        pi0_keys_run3 = run3_xsec_uncert_b["pi0"]
        eta_keys_run3 = run3_xsec_uncert_b["eta"]



    buffer = "$\\boldsymbol{M_{A^\prime}} \,$ \\textbf{[GeV]} &  $\\boldsymbol{\pi^{0}}$ \\textbf{Run 1} & $\\boldsymbol{\eta}$ \\textbf{Run 1} & $\\boldsymbol{\pi^{0}}$ \\textbf{Run 3} & $\\boldsymbol{\eta}$ \\textbf{Run 3} \\\\ \n"

    # Add horizontal lines:
    buffer+= "\\hline \n"
    buffer+= "\\hline \n"
    
    for m in masses:


        if(ratio == 0.6):
            buffer += "{:.2f} & ".format(float(m))

        else: 
            buffer += "{:.3f} & ".format(float(m))

        if m in pi0_keys_run1:
            buffer += "{:.2f} & ".format(pi0_keys_run1[m])
        else:
            buffer += " -- & "

        if m in eta_keys_run1:
            buffer+= "{:.2f} & ".format(eta_keys_run1[m])
        else:
            buffer+= " -- & "

        if m in pi0_keys_run3:
            buffer += "{:.2f} & ".format(pi0_keys_run3[m])
        else:
            buffer += " -- & "
        
        if m in eta_keys_run3:
            buffer += "{:.2f} \\\\ \n".format(eta_keys_run3[m])
        else:
            buffer += " -- \\\\ \n"


    # Add horizontal lines:
    buffer+= "\\hline \n"
    buffer+= "\\hline \n"    

    return buffer 














def LimitsTable(input_df, ratio, exp_mode="times"):

    '''

    Arguments:
        - input_df [Pandas dataframe]
        - ratio: a float indicating the mass ratio
        - exp_mode: string used to check if \cdot or \times 

    This function expectes a data frame 
    containing the limits obtained from Pyhf
    or other limit setting framework.

    It returns an string that can be used
    to print the data frame into a 
    latex table 

    '''

    masses = input_df["mass"].to_numpy()
    observed = input_df["observed"].to_numpy()
    expected = input_df["epsilon_squared"].to_numpy()
    two_up = input_df["two_sig_up"].to_numpy()
    one_up = input_df["one_sig_up"].to_numpy()
    one_down = input_df["one_sig_down"].to_numpy()
    two_down = input_df["two_sig_down"].to_numpy()


    # Initialize buffer with:
    buffer = "$ \\boldsymbol{M_{A^\prime}} \\, \\boldsymbol{\\left[ \\textrm{MeV} \\right]}$ & \\textbf{Observed} & \\textbf{Expected} & \\textbf{Exp.} $\\boldsymbol{1\sigma}$ & \\textbf{Exp.} $\\boldsymbol{2\sigma}$ \\\\ \n"

    # Add horizontal lines:
    buffer+= "\\hline \n"
    buffer+= "\\hline \n"


    for i in range(len(observed)):
 
        obs_temp = ParseNumber(observed[i])
        exp_temp = ParseNumber(expected[i])
        one_band = CompareParse(one_down[i],one_up[i])
        two_band = CompareParse(two_down[i],two_up[i])

        if(exp_mode == "cdot"):
            obs_string = r"{:.2f}$\cdot 10^{{{exp}}}$".format(obs_temp[0], exp=str(obs_temp[1]))
            exp_string = r"{:.2f}$\cdot 10^{{{exp}}}$".format(exp_temp[0], exp=str(exp_temp[1]))
            two_sig_string = r"({:.2f} - {:.2f})$\cdot 10^{{{exp}}}$ ".format(two_band[0],two_band[1],exp=str(two_band[2]))
            one_sig_string = r"({:.2f} - {:.2f})$\cdot 10^{{{exp}}}$".format(one_band[0],one_band[1],exp=str(one_band[2]))
        else:
            obs_string = r"{:.2f}$\times 10^{{{exp}}}$".format(obs_temp[0], exp=str(obs_temp[1]))
            exp_string = r"{:.2f}$\times 10^{{{exp}}}$".format(exp_temp[0], exp=str(exp_temp[1]))
            two_sig_string = r"({:.2f} - {:.2f})$\times 10^{{{exp}}}$ ".format(two_band[0],two_band[1],exp=str(two_band[2]))
            one_sig_string = r"({:.2f} - {:.2f})$\times 10^{{{exp}}}$".format(one_band[0],one_band[1],exp=str(one_band[2]))



        if(ratio==0.6):
            buffer+= "{} & {} & {} & {} & {} \\\\ \n".format(int(masses[i]*1000), obs_string, exp_string, one_sig_string, two_sig_string)

        else:
            buffer+= "{} & {} & {} & {} & {} \\\\ \n".format(int(masses[i]*1000), obs_string, exp_string, one_sig_string, two_sig_string)

    # Add horizontal lines:
    buffer+= "\\hline \n"
    buffer+= "\\hline \n"

    return buffer




