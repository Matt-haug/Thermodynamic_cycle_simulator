# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 09:41:03 2023

@author: Mattéo Hauglustaine
"""

import streamlit as st
from Rankine_Hirn_cycle import RH_cycle

st.set_page_config(page_title='Thermodynamic cycles simulator',layout='centered',page_icon=':fire:')
st.set_option('deprecation.showPyplotGlobalUse', False)

st.title('Thermodynamic cycles simulator :factory:')

st.divider()

cycle_choice = st.selectbox('Select the cycle you want to compute',('Rankine Hirn','More cycles upcoming ...'))
    
if cycle_choice == 'Rankine Hirn':
    
    col1,col2,col3 = st.columns(3)
    
    T_max = col1.number_input('$T_{max} \ [K]$')
    T_min = col2.number_input('$T_{min} \ [K]$')
    p_max = col3.number_input('$p_{max} \ [Pa]$')
    
    advanced_param_choice = st.multiselect('Choose wich advanced parameters you want to add',['p_min : Minimum steam pressure (condenser pressure)',
                                                                                              'x_min : Minimum vapor quality admissible for the turbine',
                                                                                              'P_out : Desired output power at the turbine shaft',
                                                                                              'qm : Desired total mass flow',
                                                                                              'eta_isent_T : Isentropic efficiency of the turbine',
                                                                                              'eta_isent_P : Isentropic efficiency of the pump',
                                                                                              'eta_mec : Mechanical efficiency (bearings, coupling, ...)',
                                                                                              'eta_VG : Vapor generator efficiency',
                                                                                              'p_reheating : Reheating(s) pressure(s)',
                                                                                              'p_bleeding : Bleeding(s) pressure(s)'])
    
    if 'p_min : Minimum steam pressure (condenser pressure)' in advanced_param_choice:
        
        p_min = st.number_input('$p_{min} \ [Pa]$')
        
    if 'x_min : Minimum vapor quality admissible for the turbine' in advanced_param_choice:
        
        x_min = st.number_input('$x_{min} \ [/]$')
        
    if 'P_out : Desired output power at the turbine shaft' in advanced_param_choice:
        
        P_out = st.number_input('$P_{out} \ [W]$')
        
    if 'qm : Desired total mass flow' in advanced_param_choice:
        
        qm = st.number_input('$q_{m} \ [kg/s]$')
    
    if 'eta_isent_T : Isentropic efficiency of the turbine' in advanced_param_choice:
        
        eta_isent_T = st.number_input('$\eta_{T_{is}} \ [/]$')
        
    if 'eta_isent_P : Isentropic efficiency of the pump' in advanced_param_choice:
        
        eta_isent_P = st.number_input('$\eta_{P_{is}} \ [/]$')
        
    if 'eta_mec : Mechanical efficiency (bearings, coupling, ...)' in advanced_param_choice:
        
        eta_mec = st.number_input('$\eta_{mec} \ [/]$')
        
    if 'eta_VG : Vapor generator efficiency' in advanced_param_choice:
        
        eta_VG = st.number_input('$\eta_{VG} \ [/]$')
    
    if 'p_reheating : Reheating(s) pressure(s)' in advanced_param_choice:
                        
        p_reheating_input = st.text_input('$p_{reheating} \ [Pa]$ (if many, write ";" between each pressure)',placeholder='Example : 2000000 ; 1.6e5')
        
        if p_reheating_input == '' or p_reheating_input == '0':
            p_reheating_input = None
        
        if p_reheating_input != None:
            p_reheating = p_reheating_input.split(';')
            p_reheating = [float(i) for i in p_reheating]
            
    if 'p_bleeding : Bleeding(s) pressure(s)' in advanced_param_choice:
        
        st.markdown('*Bleedings are yet to be implemented !*')
    
    p_min=None
    x_min=0
    P_out=None
    qm=None
    eta_isent_T=1
    eta_isent_P=1
    eta_mec=1
    eta_VG=1
    p_reheating=None
    p_bleeding=None
    
    @st.cache_data
    def compute(T_max,T_min,p_max,p_min,x_min,P_out,qm,eta_isent_T,eta_isent_P,eta_mec,eta_VG,p_reheating,p_bleeding):
            
            cycle=RH_cycle(T_max, T_min, p_max,p_min,x_min,P_out,qm,eta_isent_T,eta_isent_P,eta_mec,eta_VG,p_reheating,p_bleeding)
            cycle.evaluate()
            ts_plot,hs_plot,power_fig,states,metrics = cycle.return_data()
            
            st.pyplot(ts_plot)
            st.divider()
            st.pyplot(hs_plot)
            st.divider()
            
            if power_fig != None:    
                st.pyplot(power_fig)
            
            st.divider()
            
            if states['2'][-1] < x_min and states['2'][-1] != -1:
                st.warning('Caution, this cycle does not match the minimal vapor quality required for the turbine (x = {:.2f})'.format(x_min),icon="⚠️")
            
            col1,col2,col3,col4,col5 = st.columns(5)
            states_keys_list = list(states.keys())
            
            for i in range(len(states)):
                with col1:
                    st.subheader(':red[State {}]'.format(states_keys_list[i]))
                    st.markdown('**_entropy_**')
                    st.markdown('${:.2f}$'.format(states[states_keys_list[i]][0]*1e-3))
                    st.markdown('$kJ/kg*K$')
                with col2:
                    st.subheader('------------')
                    st.markdown('**_temperature_**')
                    st.markdown('${:.2f}$'.format(states[states_keys_list[i]][1]-273.15))
                    st.markdown('$°C$')
                with col3:
                    st.subheader('------------')
                    st.markdown('**_enthalpy_**')
                    st.markdown('${:.2f}$'.format(states[states_keys_list[i]][2]*1e-3))
                    st.markdown('$°kJ/kg$')
                with col4:
                    st.subheader('------------')
                    st.markdown('**_pressure_**')
                    st.markdown('${:.2f}$'.format(states[states_keys_list[i]][3]*1e-5))
                    st.markdown('$bar$')   
                with col5:
                    st.subheader('------------')
                    st.markdown('**_vapor quality_**')
                    if states[states_keys_list[i]][4] == -1:
                        st.markdown('$/$')
                    else:
                        st.markdown('${:.2f}$'.format(states[states_keys_list[i]][4]))
                    st.markdown('$[/]$')
            
            st.divider()
            
            col1,col2 = st.columns(2)
            with col1:
                if metrics[0] == None and metrics[1] == None:
                    st.markdown('$P_{{out}} = / \ [MW]$ ')
                    st.markdown('$q_{{m}} = / \ [kg/s]$ ')
                else:
                    st.markdown('$P_{{out}} = {:.2f} \ [MW]$ '.format(metrics[0]*1e-6))
                    st.markdown('$q_{{m}} = {:.2f} \ [kg/s]$ '.format(metrics[1]))
                st.markdown('$\eta_{{th}} = {:.2f} \ [\%]$ '.format(metrics[2]*100))
                st.markdown('$\eta_{{Carnot}} = {:.2f} \ [\%]$ '.format(metrics[3]*100))
            with col2:
                st.markdown('$W_{{turbine}} = {:.2f} \ [kW]$ '.format(metrics[4]*1e-3))
                st.markdown('$W_{{pump}} = {:.2f} \ [kW]$ '.format(metrics[5]*1e-3))
                st.markdown('$Q_{{VG}} = {:.2f} \ [kW]$ '.format(metrics[6]*1e-3))
                st.markdown('$Q_{{condenser}} = {:.2f} \ [kW]$ '.format(metrics[7]*1e-3))
                
    compute_cycle = st.button('Compute cycle',type='primary')
    
    if compute_cycle:
        
        compute(T_max, T_min, p_max,p_min,x_min,P_out,qm,eta_isent_T,eta_isent_P,eta_mec,eta_VG,p_reheating,p_bleeding)
        
st.divider()

st.caption('Made by **Mattéo Hauglustaine** - 2023 - *Click [here](https://github.com/Matt-haug/Thermodynamic_cycle_simulator) to get acces to the source code on my GitHub*')
