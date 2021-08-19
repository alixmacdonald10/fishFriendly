from numpy import dtype
import fishFriendly_funcs as main
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import timeit
from numba import jit, njit


def input():
    pump_des = 'CBF'
    pump_size = '140'
    if 'BF' in pump_des:
        pump_name = f'{pump_des} {pump_size}_12'
    else:
        pump_name = f'{pump_des} {pump_size}_05'
    # convert name to database naming convention
    database_name = pump_name.lower()[1:]
    # pump speed
    pump_speed = 196
    pump_speed = float(pump_speed)
    # fish type and size
    fish_type = 'fish'
    L_f = 0.5
    L_f = float(L_f)
    B_f = 0
    B_f = float(B_f)
    fish_db = {'fish_type': fish_type, 'L_f': L_f, 'B_f': B_f}
    # intake used?
    intake = 'Y'
    # duty points
    H_duty = 4.8
    H_duty = float(H_duty)
    Q_duty = 4
    Q_duty = float(Q_duty)
    # database info
    database_path = os.path.join(os.getcwd(), 'database.json')
    # load in the pump data    
    pump_db = main.load_pump(database_name, database_path)
    
    #scale duty to suit operating speed and overwrite database
    if pump_speed == 0:
        pump_db['Q'] = pd.Series(pump_db['Q'])
        pump_db['H'] = pd.Series(pump_db['H'])
        pump_db['effy'] = pd.Series(pump_db['effy'])
        pump_db['P'] = pd.Series(pump_db['P'])
    else:
        Q_scaled, H_scaled, P_scaled, N = main.scale_duty(
            pump_speed,
            pd.Series(pump_db['Q']),
            pd.Series(pump_db['H']),
            pd.Series(pump_db['P']),
            pump_db['N']
        )
        pump_db['Q'] = Q_scaled.to_numpy(dtype='float32')
        pump_db['H'] = H_scaled.to_numpy(dtype='float32')
        pump_db['effy'] = pd.Series(pump_db['effy']).to_numpy(dtype='float32')
        pump_db['P'] = P_scaled.to_numpy(dtype='float32')
        pump_db['N'] = N
    
    return pump_db, fish_db, intake, Q_duty, fish_type, pump_name, pump_speed, L_f


def run(pump_db, fish_db, intake, Q_duty):
    # analyse to standard
    L_f = fish_db['L_f']
    B_f = fish_db['B_f']
    fish_type = fish_db['fish_type']
    # pump inputs 
    Q = pump_db['Q']
    H = pump_db['H']
    N = pump_db['N']
    D_blade = pump_db['D_blade']
    d_blade = pump_db['d_blade']
    n_blade = pump_db['n_blade']
    NEN_r_array = np.array(pump_db['NEN_r_array'], dtype=np.float64)
    NEN_beta_array = np.array(pump_db['NEN_beta_array'], dtype=np.float64)
    NEN_delta_array = np.array(pump_db['NEN_delta_array'], dtype=np.float64)
    imp_thk = np.array(pump_db['imp_thk'], dtype=np.float64)
    r_imp_thk = np.array(pump_db['r_imp_thk'], dtype=np.float64)
    v_strike, P_th, f_MR, P_m = main.NEN_analyse(Q, H, N, D_blade, d_blade, n_blade, NEN_r_array, NEN_beta_array, NEN_delta_array, imp_thk, r_imp_thk, L_f, B_f, fish_type, intake, n_steps=30)
    if Q_duty is not None:
        idx_duty, _ = main.find_nearest(pump_db['Q'], Q_duty)
        duty_db = {'H': pump_db['H'][idx_duty], 'Q': pump_db['Q'][idx_duty]}
    else:
        BEP_effy = float(max(pump_db['effy']))
        idx_BEP, _ = main.find_nearest(pump_db['effy'], BEP_effy)
        duty_db = {'H': pump_db['H'][idx_BEP], 'Q': pump_db['Q'][idx_BEP]}   

    return P_m, duty_db


def test():
    
    pump_db, fish_db, intake, Q_duty, fish_type, pump_name, pump_speed, L_f = input()
    P_m, duty_db = run(pump_db, fish_db, intake, Q_duty)
    # plot
    title = f'NEN 8775 Mortality Probability\n\n Pump Type: {pump_name}\nPump Speed: {pump_speed} RPM\nFish Type: {fish_type}\nFish Length: {L_f} m'
    fig = main.plot_result(pump_db, P_m, duty_db, title)
    fig.show()  


if __name__ == '__main__':

    number = 1
    time = timeit.timeit(test, number=number)
    print(f'Total time = {time/number} s per run')
