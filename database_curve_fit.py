'''
The following script will access database information and fit to a curve for x number of points
'''
import json
import numpy as np
from numpy.random import beta
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def curve_func(x, a, b, f):

    return (a * x) + (b * x**2) + f


def func_curve_fit(f1, f2, n=100):

    popt, _ = curve_fit(curve_func, f1, f2)
    # summarize the parameter values
    a, b, f = popt
    # define a sequence of inputs between 0 and largest known inputs
    f1_new = np.arange(
        0.0001,
        max(f1),
        (max(f1) / n)
    )
    
    # calculate the output for the range
    f2_new = curve_func(f1_new, a, b, f)
    # create pandas series from lists to increase speed
    f1_new = pd.Series(f1_new)
    f2_new = pd.Series(f2_new)

    plt.figure()
    # plot original Q v H curve
    plt.scatter(f1, f2)
    # create a line plot for the mapping function for a sanity check
    plt.plot(f1_new, f2_new, '--', color='red')
    plt.title('Curve Fit Check Figure')
    plt.grid()
    plt.show(block=False)

    return f2_new, f1_new


def angle_curve_func(x, a, b, c, f):

    return (a * x) + (b * x**2) + (c * x**3) + f


def angle_curve_fit(f1, f2, n=100):

    popt, _ = curve_fit(angle_curve_func, f1, f2)
    # summarize the parameter values
    a, b, c, f = popt
    # define a sequence of inputs between the smallest and largest known inputs
    f1_new = np.arange(
        min(f1),
        max(f1),
        ((max(f1) - min(f1)) / n)
    )
    # calculate the output for the range
    f2_new = angle_curve_func(f1_new, a, b, c, f)
    # create pandas series from lists to increase speed
    f1_new = pd.Series(f1_new)
    f2_new = pd.Series(f2_new)

    plt.figure()
    # plot original Q v H curve
    plt.scatter(f1, f2)
    # create a line plot for the mapping function for a sanity check
    plt.plot(f1_new, f2_new, '--', color='red')
    plt.title('Curve Fit Check Figure')
    plt.grid()
    plt.show(block=False)

    return f2_new, f1_new


def thk_curve_func(x, a, b, f):

    return (a * x) + (b * x**2) + f


def thk_curve_fit(f1, f2, n=100):

    popt, _ = curve_fit(thk_curve_func, f1, f2)
    # summarize the parameter values
    a, b, f = popt
    # define a sequence of inputs between the smallest and largest known inputs
    f1_new = np.arange(
        min(f1),
        max(f1),
        ((max(f1) - min(f1)) / n)
    )
    # calculate the output for the range
    f2_new = thk_curve_func(f1_new, a, b, f)
    # create pandas series from lists to increase speed
    f1_new = pd.Series(f1_new)
    f2_new = pd.Series(f2_new)

    plt.figure()
    # plot original Q v H curve
    plt.scatter(f1, f2)
    # create a line plot for the mapping function for a sanity check
    plt.plot(f1_new, f2_new, '--', color='red')
    plt.title('Curve Fit Check Figure')
    plt.grid()
    plt.show(block=False)

    return f2_new, f1_new



if __name__ == "__main__":
    # inputs
    database_path = 'D:\\Scripts\\fishFriendly\\database.json'
    points = 30  # number of points to curve fit
    # load in the pump data    
    with open(database_path) as f:
        data = json.load(f)
    for pump in data['pumps']:
        # loop for head and flow curves
        if pump['type'] == 'af 120_05':
            # H, Q = func_curve_fit(pump['Q'], pump['H'], n=points)
            Q = [1868.75,
                1868.80,
                1868.88,
                2491.84,
                3114.81,
                3737.77,
                4360.73,
                4983.69,
                5606.65,
                6229.61,
                6852.57,
                7475.53,
                7787.02,
                7942.76,
                8098.50
            ]
            effy = [30.16,
                    30.16,
                    30.16,
                    40.21,
                    50.26,
                    60.31,
                    70.37,
                    78.41,
                    84.44,
                    86.45,
                    83.94,
                    68.36,
                    48.82,
                    34.58,
                    15.26
            ]
            effy, Q = func_curve_fit(Q, effy, n=points)
            # P, Q = func_curve_fit(pump['Q'], pump['P'], n=points)
            # pump['Q'] = pd.Series.tolist(Q) 
            # pump['H'] = pd.Series.tolist(H)
            pump['effy'] = pd.Series.tolist(effy)
            # pump['P'] = pd.Series.tolist(P)
            print(pump['effy'])
            
        # loop for beta and delta angles and radius
        # if type(pump['NEN_r_array']) == list:
        #     beta, r = angle_curve_fit(pump['NEN_r_array'], pump['NEN_beta_array'], n=points)
        #     delta, r = angle_curve_fit(pump['NEN_r_array'], pump['NEN_delta_array'], n=points)
        #     pump['NEN_r_array'] = pd.Series.tolist(r)
        #     pump['NEN_beta_array'] = pd.Series.tolist(beta)
        #     pump['NEN_delta_array'] = pd.Series.tolist(delta)
        # loop for blade thickness and radius
        # if type(pump['imp_thk']) == list:
        #     thk, r = thk_curve_fit(pump['r_imp_thk'], pump['imp_thk'], n=points)
        #     pump['imp_thk']= pd.Series.tolist(thk)
        #     pump['r_imp_thk'] = pd.Series.tolist(r)
    # save data to json file
    # with open(database_path, 'w') as f:
    #     json.dump(data, f)
    