import fishFriendly_funcs as main
import os
import pandas as pd
import matplotlib.pyplot as plt
import timeit



def test():
    
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
        pump_db['Q'] = Q_scaled
        pump_db['H'] = H_scaled
        pump_db['effy'] = pd.Series(pump_db['effy'])
        pump_db['P'] = P_scaled
        pump_db['N'] = N
    
    # analyse to standard
    v_strike, P_th, f_MR, P_m = main.NEN_analyse(pump_db, fish_db, intake, n_steps=30)
    if Q_duty is not None:
        idx_duty, _ = main.find_nearest(pump_db['Q'], Q_duty)
        duty_db = {'H': pump_db['H'][idx_duty], 'Q': pump_db['Q'][idx_duty]}
    else:
        BEP_effy = float(max(pump_db['effy']))
        idx_BEP, _ = main.find_nearest(pump_db['effy'], BEP_effy)
        duty_db = {'H': pump_db['H'][idx_BEP], 'Q': pump_db['Q'][idx_BEP]}
        
    Q_idx, _ = main.find_nearest(pump_db['Q'], duty_db['Q'])
    H_idx, _ = main.find_nearest(pump_db['H'], duty_db['H'])
    duty_mortality = format(P_m[Q_idx, H_idx] * 100, ".2f")
    print(f'Mortality probability: {(duty_mortality)}%')

    # plot
    title = f'NEN 8775 Mortality Probability\n\n Pump Type: {pump_name}\nPump Speed: {pump_speed} RPM\nFish Type: {fish_type}\nFish Length: {L_f} m'
    fig = main.plot_result(pump_db, P_m, duty_db, title)
    fig.show()
    


if __name__ == '__main__':

    number = 5
    time = timeit.timeit(test, number=number)
    print(f'Total time = {time/number} s per run')
