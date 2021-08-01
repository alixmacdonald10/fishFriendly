'''
The following script will access database information and fit to a curve for x number of points
'''
import json
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def curve_func(x, a, b, c, d, e, f):

    return (a * x) + (b * x**2) + (c * x**3) + (d * x**4) + (e * x**5) + f


def func_curve_fit(f1, f2, n=100):

    popt, _ = curve_fit(curve_func, f1, f2)
    # summarize the parameter values
    a, b, c, d, e, f = popt
    # define a sequence of inputs between 0 and largest known inputs
    f1_new = np.arange(
        0.0001,
        max(f1),
        (max(f1) / n)
    )
    
    # calculate the output for the range
    f2_new = curve_func(f1_new, a, b, c, d, e, f)
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
    plt.show()

    return f2_new, f1_new



if __name__ == "__main__":
    # inputs
    database_path = 'D:\\Scripts\\fishFriendly\\database.json'
    points = 30
    # load in the pump data    
    with open(database_path) as f:
        data = json.load(f)
    for pump in data['pumps']:
        # loop for head and flow curves
        pump['Q']
        pump['H']
        pump['H'], pump['Q'] = func_curve_fit(pump['Q'], pump['H'], n=points)
        pump['effy'], pump['Q'] = func_curve_fit(pump['Q'], pump['effy'], n=points)
        pump['P'], pump['Q'] = func_curve_fit(pump['Q'], pump['P'], n=points)
        # loop for beta and delta angles and radius
        pump['NEN_beta_array'], pump['NEN_r_array'] = func_curve_fit(pump['NEN_r_array'], pump['NEN_beta_array'], n=points)
        pump['NEN_delta_array'], pump['NEN_r_array'] = func_curve_fit(pump['NEN_r_array'], pump['NEN_delta_array'], n=points)
        # loop for blade thickness and radius
        pump['imp_thk'], pump['r_imp_thk'] = func_curve_fit(pump['r_imp_thk'], pump['imp_thk'], n=points)
    

        