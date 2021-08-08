import streamlit as st
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt



def load_pump(database_name, database_path):
    # open database and return data 
    pump_db = None
    with open(database_path) as f:
        data = json.load(f)
    for pump in data['pumps']:
        if pump['type'] == database_name:
            pump_db = pump
            break

    return pump_db
    

def scale_duty(pump_speed, Q, H, P, N):
    Q_scaled = Q_scale(Q, N, pump_speed)
    H_scaled = H_scale(H, N, pump_speed)
    P_scaled = P_scale(P, N, pump_speed)
    N = pump_speed
    
    return Q_scaled, H_scaled, P_scaled, N     


def Q_scale(Q1, N1, N2, D1=1.0, D2=1.0):
    Q2 = ((N2 / N1) * (D2 / D1)) * Q1
    return Q2


def H_scale(H1, N1, N2, D1=1.0, D2=1.0):  
    H2 = ((N2**2 / N1**2) * (D1**2 / D2**2)) * H1
    return H2  


def P_scale(P1, N1, N2, D1=1.0, D2=1.0):
    P2 = ((N2**2 / N1**2) * (D1**2 / D2**2)) * P1
    return P2


def N_scale(N1, H1, H2, D1=1.0, D2=1.0):
    N2 = (H2**0.5 / H1**0.5) * (D1 / D2) * N1
    return N2


def find_nearest(array, value):
    idx = np.abs(array - value).argmin()
    val = array[idx]
    return idx, val


def NEN_analyse(pump_db, fish_db, intake, n_steps=30):
    
    # calculate effective length of fish
    L_eff = length_eff(fish_db['L_f'], fish_db['B_f'], fish_db['fish_type'], intake)
    # define max length of  head and flow curve series
    Q_points = len(pump_db['Q'])
    H_points = len(pump_db['H'])
    # change in radius per step
    R_o = pump_db['D_blade'] / 2
    R_i = pump_db['d_blade'] / 2
    d_r = (R_o - R_i) / n_steps
    #Initialise zero matrices as correct size for storing variables at the end of the loop
    v_strike_max_array= np.zeros([Q_points, H_points], dtype=np.float64, order='C')
    P_th = np.zeros([Q_points, H_points], dtype=np.float64, order='C')
    f_MR = np.zeros([Q_points, H_points], dtype=np.float64, order='C')
    P_m = np.zeros([Q_points, H_points], dtype=np.float64, order='C')
    # for each radial position along the blade and for each flowrate the
    # mortality rate will be calculated.
    for i in range(0, Q_points):
        # flow
        Q_temp = pump_db['Q'][i]
        # flow velocity
        area = np.pi * (R_o**2 - R_i**2)
        v_m = Q_temp / area
        # loop through head values
        for j in range(0, H_points):
            # head
            H_temp = pump_db['H'][j]
            # determine pump speed required for head and flow value
            N_pump_scale = N_scale(pump_db['N'], pump_db['H'][i], H_temp)
            # determine rotational speed
            omega = (2.0 * np.pi * N_pump_scale) / 60.0  # rad/s
            # initial area - reset at each duty position
            A = 0
            # step along blade at this duty and sum factors
            v_strike_max = 0
            for i_r in range(0, n_steps):
                # radial position
                r = R_i + (i_r) * d_r + (d_r / 2)  # removed - 1 due to index
                # change in area
                d_A = 2 * np.pi * r * d_r
                # new area - overwriting previous A value
                A = A + d_A
                # collision probabilty
                dP_th = collision_probability(L_eff, v_m, omega, r, pump_db['n_blade'])
                # determine angles at radial position
                r_angle_lookup = r / R_o
                r_angle_idx, _ = find_nearest(pd.Series(pump_db['NEN_r_array']), r_angle_lookup)
                beta_deg = pump_db['NEN_beta_array'][r_angle_idx]
                beta = np.radians(beta_deg)
                delta_deg = pump_db['NEN_delta_array'][r_angle_idx]
                delta = np.radians(delta_deg)
                # strike velocity
                v_strike = strike_velocity(v_m, omega, r, beta, delta)
                if v_strike > v_strike_max:
                    v_strike_max = v_strike
                # determine thickness at radial position
                r_thk_idx, _ = find_nearest(pd.Series(pump_db['r_imp_thk']), r * 1e3)
                d = pump_db['imp_thk'][r_thk_idx]
                # mortality factor
                df_MR = mortality_factor(L_eff, d, v_strike, fish_type)
                # mortality probabilty
                dP_m = dP_th * df_MR
                # determine cumulative probabilities and factors
                P_th[j][i] = ((P_th[j][i] + (dP_th * v_m * d_A)))
                f_MR[j][i] = ((f_MR[j][i] + (df_MR * v_m * d_A)))
                P_m[j][i] = ((P_m[j][i] + (dP_m * v_m * d_A)))
            v_strike_max_array[j][i] = v_strike_max
            # determine total mortality as per NEN 8775, section 9.9
            P_th[j][i] = P_th[j][i] / Q_temp
            f_MR[j][i] = f_MR[j][i] / Q_temp
            P_m[j][i] = P_m[j][i] / Q_temp

    return v_strike_max_array, P_th, f_MR, P_m


def length_eff(L_f, B_f, fish_type, intake):
    
    L_max = np.sqrt(L_f**2 + B_f**2)
    theeta = np.arctan(B_f / L_f)
    if intake.lower() == 'n':
        L_eff = L_max * np.cos(0 - theeta)
    else:
        if fish_type.lower() == 'eel':
            L_eff = 0.8 * L_f
        else:
            L_eff = ((2 * L_max) / np.pi) * (np.sin((np.pi / 2) - theeta) + np.sin(theeta))
            # L_eff_2 = ((2 * L_max) / np.pi) * (np.sin((np.pi / 2) + theeta) - np.sin(theeta))
            # L_eff_3 = ((2 * L_max) / np.pi) * (np.sin((np.pi / 2) - theeta) + np.sin(theeta))
            # L_eff_4 = ((2 * L_max) / np.pi) * (np.sin((np.pi / 2) - theeta) - np.sin(theeta))
            # L_eff = max(L_eff_1, L_eff_2, L_eff_3, L_eff_4)
    return L_eff


def collision_probability(L_eff, v_m, omega, r, n_blade, wf=0, alpha=0):

    # collision probability - assume no fish relative motion and no preswirl
    P_th = (L_eff * n_blade * omega * r) / (v_m  * 2 * np.pi * r)

    # probability can only be from 0 - 1 therefore:
    P_th = max(0, min(1, P_th))
    return P_th


def strike_velocity(v_m, omega, r, beta, delta):
    # assume no fish relative motion and no pre swirl
    v_strike = np.sqrt(
        (v_m * np.cos(beta))**2 + ((omega * r) * np.cos(delta))**2
    )
    return v_strike


def mortality_factor(L_f, t_blade, v_strike, fish_type):

    ratio = L_f / (t_blade * 1e-3)

    fish = fish_type.lower()
    if fish == 'fish':
        if ratio >= 0 and ratio < 2:
            a = 0.0531
            b = 0.0202
            v_crit = 4.8
        elif ratio >= 2 and ratio < 10:
            a = 0.0829
            b = -0.0021
            v_crit = 4.8
        elif ratio >= 10 and ratio < 25:
            a = 0.0327
            b = 0.1146
            v_crit = 4.8
        else:
            #warnings.warn(f'Lf/d ratio {ratio} exceeds model limits (25)')
            # place holder if data is revisted
            a = 0.0327
            b = 0.1146
            v_crit = 4.8
    elif fish == "eel":
        if ratio > 25:
            #warnings.warn(f'Lf/d ratio {ratio} exceeds model limits (25)')
            # place holder if data is revisted
            a = 0.0024
            b = 0
            v_crit = 8
        elif ratio < 1:
            #warnings.warn(f'Lf/d ratio {ratio} exceeds model lower limit (1)')
            # place holder if data is revisted
            a = 0.0024
            b = 0
            v_crit = 8
        else:
            a = 0.0024
            b = 0
            v_crit = 8

    # determine mortality factor
    if v_strike < v_crit:
        f_MR = 0
    else:
        if fish == "fish":
            f_MR = (a * np.log(ratio) + b) * (v_strike - v_crit)
        else:
            f_MR = (a * ratio) * (v_strike - v_crit)
    # bounded by probabilty between 0 and 1
    f_MR = max(0, min(1, f_MR))
    return f_MR


def plot_result(pump_db, result, duty_db, title):
    
    # set values
    X = pump_db['Q']
    Y = pump_db['H']
    Z = result
    
    # plot contour lines
    fig, ax = plt.subplots(figsize=(10, 4))
    CS = ax.contour(X, Y, Z, cmap='binary')
    ax.contourf(X, Y, Z,
                cmap='RdYlGn_r', extend='both', alpha=0.5
    )
    ax.clabel(CS, inline=True, fontsize=10)
    ax.set_xlim(left=0)
    ax.set_title(title)
    ax.set_xlabel('Flow ($m^3/s$)')
    ax.set_ylabel('Head (m)')
    # plot pump curve over contour plot
    plt.plot(pump_db['Q'], pump_db['H'])
    ax.annotate(f"{pump_db['N']} RPM", (pump_db['Q'][8], pump_db['H'][8]))
    # plot duty point
    plt.scatter(duty_db['Q'], duty_db['H'])
    Q_idx, _ = find_nearest(pump_db['Q'], duty_db['Q'])
    H_idx, _ = find_nearest(pump_db['H'], duty_db['H'])
    
    plt.text(duty_db['Q'], duty_db['H'], f'  {format(result[Q_idx, H_idx], ".2f")}')
    plt.show(block=False)
    
    return fig



# set title and subtitle
st.title('Fish Freindly Analysis')
subtitle = ('The following tool determines the mortality probability of a pump at a specific duty according to NEN 8775')
st.write(subtitle)

# sidebar inputs
st.sidebar.markdown('**Inputs**')
st.sidebar.markdown("""---""")
pump_des = st.sidebar.selectbox(
    'Pump designation',
     ['SAF', 'SBF', 'CAF', 'CBF']
)
pump_size = st.sidebar.selectbox(
    'Pump size',
     ['45', '60', '70', '90', '100', '120', '140']
)
pump_speed = st.sidebar.number_input(label="Pump speed (RPM)", value=0)
pump_speed = float(pump_speed)
fish_type = st.sidebar.selectbox('Fish type', ['eel', 'fish'])  # fish or eel
L_f = st.sidebar.number_input(label="Fish length (m)", min_value=0.0, max_value=1.5, step=1.,format="%.2f")
L_f = float(L_f)
B_f = st.sidebar.number_input(label="Fish height (m)", min_value=0.0, max_value=1.5, step=1.,format="%.2f")
B_f = float(B_f)
intake = st.sidebar.selectbox('Type 10 / Bedford intake?', ['Y', 'N']) 
H_duty = st.sidebar.number_input(label="Duty head (m)", min_value=0.0, max_value=18.0, step=1.,format="%.1f")
H_duty = float(H_duty)
Q_duty = st.sidebar.number_input(label="Duty flow (m\N{SUPERSCRIPT THREE}/s)", min_value=0.0, max_value=10.0, step=1.,format="%.3f")
Q_duty = float(Q_duty)

# database info
database_path = 'D:\\Scripts\\fishFriendly\\database.json'
st.markdown("""---""")
# run 
pressed = st.sidebar.button('Run Analysis')
if pressed:    
    if 'BF' in pump_des:
        pump_name = f'{pump_des} {pump_size}_12'
    else:
        pump_name = f'{pump_des} {pump_size}_05'
    # convert name to database naming convention
    database_name = pump_name.lower()[1:]
    # load in the pump data    
    pump_db = load_pump(database_name, database_path)
    placeholder = st.empty()
    if pump_db is None:
        placeholder.write('ERROR: Pump type not found!\nPlease try another type and contact AM')
    else:
        placeholder.write('Running analysis (this could take a few seconds)')

    #scale duty to suit operating speed and overwrite database
    Q_scaled, H_scaled, P_scaled, N = scale_duty(
        pump_speed,
        pd.Series(pump_db['Q']),
        pd.Series(pump_db['H']),
        pd.Series(pump_db['P']),
        pump_db['N']
    )
    pump_db['Q'] = Q_scaled
    pump_db['H'] = H_scaled
    pump_db['P'] = P_scaled
    pump_db['N'] = N
    
    idx_duty, _ = find_nearest(pump_db['Q'], Q_duty)
    duty_db = {'H': pump_db['H'][idx_duty], 'Q': pump_db['Q'][idx_duty]}
    fish_db = {'fish_type': fish_type, 'L_f': L_f, 'B_f': B_f}
    # analyse to standard
    v_strike, P_th, f_MR, P_m = NEN_analyse(pump_db, fish_db, intake, n_steps=30)
    title = f'NEN 8775 Mortality Probability\n\n Pump Type: {pump_name}\nPump Speed: {pump_speed} RPM\nFish Type: {fish_type}\nFish Length: {L_f} m'
    fig = plot_result(pump_db, P_m, duty_db, title)
    placeholder.empty()
    placeholder.write('Analysis complete!')
    # find duty mortality probability
    Q_idx, _ = find_nearest(pump_db['Q'], duty_db['Q'])
    H_idx, _ = find_nearest(pump_db['H'], duty_db['H'])
    duty_mortality = format(P_m[Q_idx, H_idx] * 100, ".2f")
    # print mortality probability and figure
    st.write(f'Mortality probability: {(duty_mortality)}%')
    st.write(fig)
    st.caption('Download figure by right clicking and selecting "Save image as..."')
