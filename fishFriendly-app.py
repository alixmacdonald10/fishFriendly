import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fishFriendly_funcs as main



def NEN_analyse(pump_db, fish_db, intake, n_steps=30):
    
    # calculate effective length of fish
    L_eff = main.length_eff(fish_db['L_f'], fish_db['B_f'], fish_db['fish_type'], intake)
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
            N_pump_scale = main.N_scale(pump_db['N'], pump_db['H'][i], H_temp)
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
                dP_th = main.collision_probability(L_eff, v_m, omega, r, pump_db['n_blade'])
                # determine angles at radial position
                r_angle_lookup = r / R_o
                r_angle_idx, _ = main.find_nearest(pd.Series(pump_db['NEN_r_array']), r_angle_lookup)
                beta_deg = pump_db['NEN_beta_array'][r_angle_idx]
                beta = np.radians(beta_deg)
                delta_deg = pump_db['NEN_delta_array'][r_angle_idx]
                delta = np.radians(delta_deg)
                # strike velocity
                v_strike = main.strike_velocity(v_m, omega, r, beta, delta)
                if v_strike > v_strike_max:
                    v_strike_max = v_strike
                # determine thickness at radial position
                r_thk_idx, _ = main.find_nearest(pd.Series(pump_db['r_imp_thk']), r * 1e3)
                d = pump_db['imp_thk'][r_thk_idx]
                # mortality factor
                df_MR = main.mortality_factor(L_eff, d, v_strike, fish_type)
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



# set title and subtitle
st.title('Fish Freindly Analysis')
st.markdown("""
            The following tool determines the probability of fish mortality as it passes
            through a pump at a specific duty according to Dutch Standard NEN 8775.\n
            
            Calculation inputs are located on the left hand toolbar.
            
            Notes on implementation:
            - Pump speed is not required. If blank then nominal speed is used.
            - Fish sizes are not required. If blank then fish sizes from NEN 8775 are used.
            - Duty points are not required. If blank then BEP will be plotted instead.
            """
)
st.caption('Missing data for AF 100 and BF 70 type pumps')

# sidebar inputs
st.sidebar.markdown("**Inputs**")
st.sidebar.caption("*Not required")
# pump name
left_column, right_column = st.sidebar.columns(2)
pump_des = left_column.selectbox(
    'Pump designation',
     ['DAF', 'DBF', 'SAF', 'SBF', 'CAF', 'CBF']
)
pump_size = right_column.selectbox(
    'Pump size',
     ['45', '60', '70', '90', '100', '120', '140']
)
if 'BF' in pump_des:
    pump_name = f'{pump_des} {pump_size}_12'
else:
    pump_name = f'{pump_des} {pump_size}_05'
# convert name to database naming convention
database_name = pump_name.lower()[1:]
# pump speed
pump_speed = st.sidebar.number_input(label="Pump speed (RPM) *", value=0)
pump_speed = float(pump_speed)
# fish type and size
fish_type = st.sidebar.selectbox('Fish type', ['eel', 'fish'])  # fish or eel
if st.sidebar.checkbox('Non-standard fish size? *'):
    L_f = st.sidebar.number_input(label="Fish length (m)", min_value=0.0, max_value=1.5, step=1.,format="%.2f")
    L_f = float(L_f)
    B_f = st.sidebar.number_input(label="Fish height (m) *", min_value=0.0, max_value=1.5, step=1.,format="%.2f")
    B_f = float(B_f)
else:
    if fish_type == 'eel':
        L_f = 0.8
        B_f = 0
    else:
        L_f = 0.23
        B_f = 0
fish_db = {'fish_type': fish_type, 'L_f': L_f, 'B_f': B_f}
# intake used?
if st.sidebar.checkbox('Type 10/BPL Intake?'):
    intake = 'Y'
else:
    intake = 'N'
# duty points
if st.sidebar.checkbox('Show Duty? *'):  
    H_duty = st.sidebar.number_input(label="Duty head (m) *", min_value=0.0, max_value=18.0, step=1.,format="%.1f")
    H_duty = float(H_duty)
    Q_duty = st.sidebar.number_input(label="Duty flow (m\N{SUPERSCRIPT THREE}/s) *", min_value=0.0, max_value=10.0, step=1.,format="%.3f")
    Q_duty = float(Q_duty)
else:
    H_duty = None
    Q_duty = None

# database info
database_path = 'database.json'

# seperator
st.markdown("""---""")
# run 
pressed = st.sidebar.button('Run Analysis')
if pressed:    
    # load in the pump data    
    pump_db = main.load_pump(database_name, database_path)
    with st.spinner(text='Running analysis (this could take up to a minute)...'):
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
        v_strike, P_th, f_MR, P_m = NEN_analyse(pump_db, fish_db, intake, n_steps=30)
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
        st.write(f'Mortality probability: {(duty_mortality)}%')

        # plot
        title = f'NEN 8775 Mortality Probability\n\n Pump Type: {pump_name}\nPump Speed: {pump_speed} RPM\nFish Type: {fish_type}\nFish Length: {L_f} m'
        fig = main.plot_result(pump_db, P_m, duty_db, title)
        st.success('Analysis complete!')

    # print mortality probability and figure
    st.write(fig)
    st.caption('Download figure by right clicking and selecting "Save image as..."')

