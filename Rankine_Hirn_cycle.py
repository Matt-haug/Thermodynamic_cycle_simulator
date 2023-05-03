# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 15:15:04 2023

@author: Mattéo Hauglustaine

Simulation of a Rankine Hirn cycle.

The purpose of this class is to compute the states of a Rankine Hirn cycle
and the major output values such as work, efficiency or power for instance.
The cycle may contain reheatings and/or bleedings. The script is based on
CoolProp, Numpy and Matplotlib libraries.

"""

### Libraries importation ###

from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt

#======================================#

class RH_cycle():
    
    """
    The basic inputs are mandatory parameters :
        - T_max : [K] Maximum steam temperature
        - T_min : [K] Minimum steam temperature (condenser temperature)
        - p_max : [Pa] Maximum steam pressure
        
    The advanced inputs are optional parameters:
        - p_min : [Pa] Minimum steam pressure (condenser pressure)
        - x_min : Minimum vapor quality admissible for the turbine
        - P_out : [W] Desired output power at the turbine shaft
        - qm : [kg/s] Desired total mass flow
        - eta_isent_T : Isentropic efficiency of the turbine
        - eta_isent_P : Isentropic efficiency of the pump
        - eta_mec : Mechanical efficiency (bearings, coupling, ...)
        - eta_VG : Vapor generator efficiency
        - p_reheating : [Pa] Reheating(s) pressure(s) (list type)
        - p_bleeding : [Pa] Bleeding(s) pressure(s) (list type)
        - makeplot : Make the plots if True
    
    The method 'evaluate' computes the cycle.
    
    The method 'return_data' allows to store the outputs into variables.
    
    """
    def __init__(self,T_max,T_min,p_max,p_min=None,x_min=0,P_out=None,qm=None,eta_isent_T=1,eta_isent_P=1,eta_mec=1,eta_VG=1,p_reheating=None,p_bleeding=None,makeplot=True):
        
        self.T_max              = T_max         # [K] Maximum steam temperature
        self.T_min              = T_min         # [K] Minimum steam temperature (condenser temperature)
        self.p_max              = p_max         # [Pa] Maximum steam pressure
        self.p_min              = p_min         # [Pa] Minimum steam pressure (condenser pressure)
        self.x_min              = x_min         # Minimum vapor quality admissible for the turbine
        self.P_out              = P_out         # [W] Desired output power of the turbine
        self.qm                 = qm            # [kg/s] Desired total mass flow
        self.eta_isent_T        = eta_isent_T   # Isentropic efficiency of the turbine
        self.eta_isent_P        = eta_isent_P   # Isentropic efficiency of the pump
        self.eta_mec            = eta_mec       # Mechanical efficiency (bearings, coupling, ...)
        self.eta_VG             = eta_VG        # Vapor generator efficiency
        self.p_reheating        = p_reheating   # [Pa] Reheating(s) pressure(s)
        self.p_bleeding         = p_bleeding    # [Pa] Bleeding(s) pressure(s)
        self.makeplot           = makeplot      # Make the plots if True
        
    def evaluate(self):
                        
        if self.P_out and self.qm != None:
            raise ValueError('Please specify either a mass flow or an output power but not both desired values')
        
        if self.T_min and self.p_min != None:
            raise ValueError('Please specify either a condenser temperature or a condenser pressure but not both desired values. Initiate the T_min parameter at "None" if you want to input a condenser pressure.')
        
        if self.p_reheating != None and type(self.p_reheating) != list:
            raise TypeError('Please make sure the "p_reheating" argument is given as a list')
        
        if self.p_bleeding != None and type(self.p_bleeding) != list:
            raise TypeError('Please make sure the "p_bleeding" argument is given as a list')
                   
        ### States ###
            
        self.states = {} # dictionnary of states in the cycle (used for plotting)
        
        # State 1
        self.p1 = self.p_max
        self.t1 = self.T_max
        self.s1 = PropsSI('S','T',self.t1,'P',self.p1,'water')
        self.h1 = PropsSI('H','T',self.t1,'P',self.p1,'water')
        self.x1 = PropsSI('Q','T',self.t1,'P',self.p1,'water')
        self.states['1'] = [self.s1,self.t1,self.h1,self.p1,self.x1]
        
        print('State 1 : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(self.p1*1e-5,self.t1-273.15,self.h1*1e-3,self.s1*1e-3,self.x1))

        if self.p_reheating == None:
            
            # State 2
            if self.p_min == None:
                self.t2 = self.T_min
                self.h2s = PropsSI('H', 'S',self.s1,'T',self.t2,'water')
                self.h2 = self.h1 - self.eta_isent_T * (self.h1 - self.h2s)
                self.p2 = PropsSI('P','S',self.s1,'H',self.h2s,'water')
                self.s2 = PropsSI('S','P',self.p2,'H',self.h2,'water')
                self.x2 = PropsSI('Q','P',self.p2,'H',self.h2,'water')
            
            else:
                self.p2 = self.p_min
                self.h2s = PropsSI('H', 'S',self.s1,'P',self.p2,'water')
                self.h2 = self.h1 - self.eta_isent_T * (self.h1 - self.h2s)
                self.t2 = PropsSI('T','S',self.s1,'H',self.h2s,'water')
                self.s2 = PropsSI('S','P',self.p2,'H',self.h2,'water')
                self.x2 = PropsSI('Q','P',self.p2,'H',self.h2,'water')
            
            self.states['2'] = [self.s2,self.t2,self.h2,self.p2,self.x2]
                            
            print('State  2 : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(self.p2*1e-5,self.t2-273.15,self.h2*1e-3,self.s2*1e-3,self.x2))
            
        else:
            
            for i in range(len(self.p_reheating)):
                
                # State 6,8,etc
                self.t6 = self.T_max 
                self.p6 = self.p_reheating[i]
                self.h6 = PropsSI('H','T',self.t6,'P',self.p6,'water')
                self.s6 = PropsSI('S','T',self.t6,'P',self.p6,'water')
                self.x6 = PropsSI('Q','T',self.t6,'P',self.p6,'water')
                self.states[str(2*i+6)] = [self.s6,self.t6,self.h6,self.p6,self.x6]
                
                if i == 0:
                    # State 5
                    self.p5 = self.p_reheating[i]
                    self.h5s = PropsSI('H', 'S',self.s1,'P',self.p5,'water')
                    self.h5 = self.h1 - self.eta_isent_T * (self.h1 - self.h5s)
                    self.t5 = PropsSI('T','P',self.p5,'H',self.h5,'water')
                    self.s5 = PropsSI('S','P',self.p5,'H',self.h5,'water')
                    self.x5 = PropsSI('Q','P',self.p5,'H',self.h5,'water')
                    self.states[str(2*i+5)] = [self.s5,self.t5,self.h5,self.p5,self.x5]
                
                else:
                    # State 7,9,etc
                    self.p7 = self.p_reheating[i]
                    self.h7s = PropsSI('H', 'S',self.states[str((2*i+5)-1)][0],'P',self.p7,'water')
                    self.h7 = self.states[str((2*i+5)-1)][2] - self.eta_isent_T * (self.states[str((2*i+5)-1)][2] - self.h7s)
                    self.t7 = PropsSI('T','P',self.p7,'H',self.h7,'water')
                    self.s7 = PropsSI('S','P',self.p7,'H',self.h7,'water')
                    self.x7 = PropsSI('Q','P',self.p7,'H',self.h7,'water')
                    self.states[str(2*i+5)] = [self.s7,self.t7,self.h7,self.p7,self.x7]
                                
                print('State {} : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(2*i+5,self.p5*1e-5,self.t5-273.15,self.h5*1e-3,self.s5*1e-3,self.x5))
 
                print('State {} : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(2*i+6,self.p6*1e-5,self.t6-273.15,self.h6*1e-3,self.s6*1e-3,self.x6))
                
            # State 2
            if self.p_min == None:
                self.t2 = self.T_min
                self.h2s = PropsSI('H', 'S',self.s6,'T',self.t2,'water')
                self.h2 = self.h6 - self.eta_isent_T * (self.h6 - self.h2s)
                self.p2 = PropsSI('P','S',self.s6,'H',self.h2s,'water')
                self.s2 = PropsSI('S','P',self.p2,'H',self.h2,'water')
                self.x2 = PropsSI('Q','P',self.p2,'H',self.h2,'water')
            
            else:
                self.p2 = self.p_min
                self.h2s = PropsSI('H', 'S',self.s6,'P',self.p2,'water')
                self.h2 = self.h6 - self.eta_isent_T * (self.h6 - self.h2s)
                self.t2 = PropsSI('T','S',self.s6,'H',self.h2s,'water')
                self.s2 = PropsSI('S','P',self.p2,'H',self.h2,'water')
                self.x2 = PropsSI('Q','P',self.p2,'H',self.h2,'water')
            
            self.states['2'] = [self.s2,self.t2,self.h2,self.p2,self.x2]
                            
            print('State 2 : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(self.p2*1e-5,self.t2-273.15,self.h2*1e-3,self.s2*1e-3,self.x2))
    
        # State 3
        self.x3 = 0
        self.t3 = self.t2
        self.p3 = PropsSI('P', 'Q',self.x3,'T',self.t3,'water')
        self.h3 = PropsSI('H', 'Q',self.x3,'T',self.t3,'water')
        self.s3 = PropsSI('S', 'Q',self.x3,'T',self.t3,'water')
        self.states['3'] = [self.s3,self.t3,self.h3,self.p3,self.x3]
        
        print('State 3 : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(self.p3*1e-5,self.t3-273.15,self.h3*1e-3,self.s3*1e-3,self.x3))
        
        # State 4
        self.p4 = self.p_max
        self.h4s = PropsSI('H', 'S',self.s3,'P',self.p4,'water')
        self.h4 = self.h3 + self.eta_isent_P * (self.h4s - self.h3)
        self.s4 = PropsSI('S', 'H',self.h4,'P',self.p4,'water')
        self.t4 = PropsSI('T', 'H',self.h4,'P',self.p4,'water')
        self.x4 = PropsSI('Q', 'H',self.h4,'P',self.p4,'water')
        self.states['4'] = [self.s4,self.t4,self.h4,self.p4,self.x4]
        
        print('State 4 : p = {:.2f} bar, T = {:.2f} °C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg*K, x = {:.2f} \n'.format(self.p4*1e-5,self.t4-273.15,self.h4*1e-3,self.s4*1e-3,self.x4))

        if self.p_bleeding != None:
            
            
            print('écrire la cas avec sous tirages')

        if self.x2 < self.x_min:
            
            print('/!\ Caution, this cycle does not match the minimal vapor quality required for the turbine (x = {:.2f}) /!\ \n'.format(self.x_min))
    
        ### Outputs ###
        
        if self.p_reheating == None:
            
            self.Wt = (self.h1 - self.h2)
            self.Wp = (self.h4 - self.h3)
            self.Qvg = (self.h1 - self.h4)
            self.eta_th = abs((self.Wt-self.Wp)/self.Qvg)
            self.Qc = self.h2 - self.h3
            self.eta_carnot = 1 - (self.t2/self.t1)
        
        else:
            
            u=0
            for i in range(len(self.p_reheating)):
                u += self.states[str(2*i+6)][2] - self.states[str(2*i+5)][2]
                            
            self.Wt = ((self.h1 - self.h2) + u)
            self.Qvg = ((self.h1 - self.h4)+ u)
            self.Wp = (self.h4 - self.h3)
            self.eta_th = abs((self.Wt-self.Wp)/self.Qvg)
            self.Qc = self.h2 - self.h3
            self.eta_carnot = 1 - (self.t2/self.t1)
        
        print('Turbine work = {:.2f} kJ/kg \n'.format(self.Wt*1e-3))
        print('Pump work = {:.2f} kJ/kg \n'.format(self.Wp*1e-3))
        print('Generator heat = {:.2f} kJ/kg \n'.format(self.Qvg*1e-3))
        print('Condenser heat = {:.2f} kJ/kg \n'.format(self.Qc*1e-3))
        print('Thermal efficiency = {:.2f} % \n'.format(self.eta_th*100))
        print('Carnot efficiency = {:.2f} % \n'.format(self.eta_carnot*100))
        
        if self.P_out != None:
            self.qm = self.P_out/(self.Wt * self.eta_mec)
            print('Output power = {:.2f} MW \n'.format(self.P_out*1e-6))
            print('Mass flow = {:.2f} kg/s \n'.format(self.qm))
            
        if self.qm != None:
            self.P_out = self.qm * self.Wt * self.eta_mec
            print('Output power = {:.2f} MW \n'.format(self.P_out*1e-6))
            print('Mass flow = {:.2f} kg/s \n'.format(self.qm))
        
        if self.makeplot == True:
            
            ### Plotting ###
            
            plt.ioff()
            
            ## T,s diagram ##
            
            self.ts_fig = plt.figure(1,figsize=(8,6))
            plt.minorticks_on()
            plt.tick_params(which='both',direction='in',top=True,right=True)
            plt.grid(which='major',linewidth=0.5)
            plt.grid(which='minor',linewidth=0.2)
            plt.ylabel('Temperature [°C]',fontsize=14,family='serif')
            plt.xlabel('Entropy [kJ/(kg*K)]',fontsize=14,family='serif')
            plt.xticks(fontsize=14,family='serif')
            plt.yticks(fontsize=14,family='serif')
            plt.ylim((-10,self.T_max-273.15+30))
            plt.xlim((-0.5,10))
            plt.title('T-s diagram',fontsize=17,family='serif',weight='bold')
            
            # Bell curve #
            
            t_bell = np.linspace(273.15,647.09,100)
            sx0 = np.empty(len(t_bell))
            sx1 = np.empty(len(t_bell))
            h_bellx0 = np.empty(len(t_bell))
            h_bellx1 = np.empty(len(t_bell))
    
            for i in range(len(t_bell)):
                sx0[i] = PropsSI('S','T',t_bell[i],'Q',0,'water')
                sx1[i] = PropsSI('S','T',t_bell[i],'Q',1,'water')
                h_bellx0[i] = PropsSI('H','T',t_bell[i],'Q',0,'water')
                h_bellx1[i] = PropsSI('H','T',t_bell[i],'Q',1,'water')
                
            plt.plot(sx0*1e-3,t_bell-273.15,sx1*1e-3,t_bell-273.15, linestyle = '-.',color = 'black',linewidth = 0.8)
            
            # Cycle curve #
            
            if self.p_reheating == None:
                
                T = np.empty((4,100))
                S = np.empty((4,100))
                            
                T[0,:] = np.linspace(self.t1,self.t2,len(T[0]))
                T[1,:] = np.ones(len(T[0]))*self.t2
                S[1,:] = np.linspace(self.s2,self.s3,len(T[0]))
                T[2,:] = np.linspace(self.t3,self.t4,len(T[0]))
                T[3,:] = np.linspace(self.t4,self.t1,len(T[0]))
                            
                for i in range(len(T[0])):
                    
                    h2s = PropsSI('H', 'S',self.s1,'T',T[0,i],'water')
                    h2 = self.h1 - self.eta_isent_T * (self.h1 - h2s)
                    p2 = PropsSI('P','S',self.s1,'H',h2s,'water')
                    S[0,i] = PropsSI('S','P',p2,'H',h2,'water')
                                                   
                    # h4s = PropsSI('H', 'S',self.s3,'T',T[2,i],'water')
                    # h4 = self.h3 + self.eta_isent_P * (h4s - self.h3)     ==> State 3 to 4 plot, unnecessary and it makes the calculation time way bigger
                    # p4 = PropsSI('P','S',self.s3,'H',h4s,'water')
                    # S[2,i] = PropsSI('S','P',p4,'H',h4,'water')
                    
                    S[3,i] = PropsSI('S','P',self.p1,'T',T[3,i],'water')
                
                S = S*1e-3
                T = T-273.15
                
                for i in range(len(T[:,0])):
                    plt.plot(S[i,:],T[i,:],color='green',linewidth=2)
                    
            else:
                
                T = np.zeros((4+2*len(self.p_reheating),100))
                S = np.zeros((4+2*len(self.p_reheating),100))
                
                for i in range(2*len(self.p_reheating)-1):
                    T[i,:] = np.linspace(self.states[str(i+5)][1],self.states[str(i+6)][1],len(T[0]))
                    
                T[-1,:] = np.linspace(self.t1,self.t5,len(T[0]))
                T[-2,:] = np.linspace(self.t6,self.t2,len(T[0]))
                T[-3,:] = np.linspace(self.t4,self.t1,len(T[0]))
                T[-4,:] = np.ones(len(T[0]))*self.t2
                S[-4,:] = np.linspace(self.s2,self.s3,len(T[0]))
                            
                for i in range(len(T[0])):
                    
                    h5s = PropsSI('H', 'S',self.s1,'T',T[-1,i],'water')
                    h5 = self.h1 - self.eta_isent_T * (self.h1 - h5s)
                    p5 = PropsSI('P','S',self.s1,'H',h5s,'water')
                    S[-1,i] = PropsSI('S','P',p5,'H',h5,'water')
                    
                    for j in range(len(self.p_reheating)):
                        S[2*j,i] = PropsSI('S','P',self.p_reheating[j],'T',T[2*j,i],'water')
                        
                    for j in range(len(self.p_reheating)-1):
                        h7s = PropsSI('H', 'S',self.states[str(2*j+6)][0],'T',T[2*j+1,i],'water')
                        h7 = self.states[str(2*j+6)][2] - self.eta_isent_T * (self.states[str(2*j+6)][2] - h7s)
                        p7 = PropsSI('P','S',self.states[str(2*j+6)][0],'H',h7s,'water')
                        S[2*j+1,i] = PropsSI('S','P',p7,'H',h7,'water')
                                                
                    h2s = PropsSI('H', 'S',self.s6,'T',T[-2,i],'water')
                    h2 = self.h6 - self.eta_isent_T * (self.h6 - h2s)
                    p2 = PropsSI('P','S',self.s6,'H',h2s,'water')
                    S[-2,i] = PropsSI('S','P',p2,'H',h2,'water')
                                   
                    # h4s = PropsSI('H', 'S',self.s3,'T',T[2,i],'water')
                    # h4 = self.h3 + self.eta_isent_P * (h4s - self.h3)     ==> State 3 to 4 plot, unnecessary and it makes the calculation time way bigger
                    # p4 = PropsSI('P','S',self.s3,'H',h4s,'water')
                    # S[2,i] = PropsSI('S','P',p4,'H',h4,'water')
                    
                    S[-3,i] = PropsSI('S','P',self.p1,'T',T[-3,i],'water')
                
                S = S*1e-3
                T = T-273.15

                for i in range(len(T[:,0])):
                    plt.plot(S[i,:],T[i,:],color='green',linewidth=2)
                        
            states_keys_list = list(self.states.keys())
            
            for i in range(len(self.states)):
                
                plt.plot(self.states[states_keys_list[i]][0]*1e-3,self.states[states_keys_list[i]][1]-273.15,'o',color='black',linewidth=4)
                plt.text(self.states[states_keys_list[i]][0]*1e-3+0.2,self.states[states_keys_list[i]][1]-273.15,states_keys_list[i],fontsize=14,family='serif')
            
            ## h,s diagram ##
            
            self.hs_fig = plt.figure(2,figsize=(8,6))
            plt.minorticks_on()
            plt.tick_params(which='both',direction='in',top=True,right=True)
            plt.grid(which='major',linewidth=0.5)
            plt.grid(which='minor',linewidth=0.2)
            plt.ylabel('Enthalpy [kJ/kg]',fontsize=14,family='serif')
            plt.xlabel('Entropy [kJ/(kg*K)]',fontsize=14,family='serif')
            plt.xticks(fontsize=14,family='serif')
            plt.yticks(fontsize=14,family='serif')
            plt.ylim((-100,self.h1*1e-3+350))
            plt.xlim((-0.5,10))
            plt.title('h-s diagram',fontsize=17,family='serif',weight='bold')
            
            # Bell curve #
                  
            plt.plot(sx0*1e-3,h_bellx0*1e-3,sx1*1e-3,h_bellx1*1e-3, linestyle = '-.',color = 'black',linewidth = 0.8)    
            
            # Cycle curve #
            
            if self.p_reheating == None:
                
                H = np.empty((4,100))
                S = np.empty((4,100))
                            
                H[0,:] = np.linspace(self.h1,self.h2,100)
                H[1,:] = np.linspace(self.h2,self.h3,100)
                H[3,:] = np.linspace(self.h4,self.h1,100)
                            
                for i in range(len(H[0])):
                    
                    h2 = self.h1 - self.eta_isent_T * (self.h1 - H[0,i])
                    p2 = PropsSI('P','S',self.s1,'H',H[0,i],'water')
                    S[0,i] = PropsSI('S','P',p2,'H',h2,'water')
                                                   
                    S[1,i] = PropsSI('S','P',self.p2,'H',H[1,i],'water')
                    
                    S[3,i] = PropsSI('S','P',self.p1,'H',H[3,i],'water')
                
                S = S*1e-3
                H = H*1e-3
                
                for i in range(len(H[:,0])):
                    plt.plot(S[i,:],H[i,:],color='green',linewidth=2)
                    
            else:
                
                H = np.zeros((4+2*len(self.p_reheating),100))
                S = np.zeros((4+2*len(self.p_reheating),100))
                                
                for i in range(2*len(self.p_reheating)-1):
                    H[i,:] = np.linspace(self.states[str(i+5)][2],self.states[str(i+6)][2],len(H[0]))
                
                H[-4,:] = np.linspace(self.h1,self.h5,len(H[0]))
                H[-3,:] = np.linspace(self.h6,self.h2,len(H[0]))
                H[-2,:] = np.linspace(self.h2,self.h3,len(H[0]))
                H[-1,:] = np.linspace(self.h4,self.h1,len(H[0]))
                            
                for i in range(len(H[0])):
                    
                    h5 = self.h1 - self.eta_isent_T * (self.h1 - H[-4,i])
                    p5 = PropsSI('P','S',self.s1,'H',H[-4,i],'water')
                    S[-4,i] = PropsSI('S','P',p5,'H',h5,'water')
                    
                    for j in range(len(self.p_reheating)):
                        S[2*j,i] = PropsSI('S','P',self.p_reheating[j],'H',H[2*j,i],'water')
                        
                    for j in range(len(self.p_reheating)-1):
                        h7 = self.states[str(2*j+6)][2] - self.eta_isent_T * (self.states[str(2*j+6)][2] - H[2*j+1,i])
                        p7 = PropsSI('P','S',self.states[str(2*j+6)][0],'H',H[2*j+1,i],'water')
                        S[2*j+1,i] = PropsSI('S','P',p7,'H',h7,'water')
                    
                    h2 = self.h6 - self.eta_isent_T * (self.h6 - H[-3,i])
                    p2 = PropsSI('P','S',self.s6,'H',H[-3,i],'water')
                    S[-3,i] = PropsSI('S','P',p2,'H',h2,'water')
                                                   
                    S[-2,i] = PropsSI('S','P',self.p2,'H',H[-2,i],'water')
                    
                    S[-1,i] = PropsSI('S','P',self.p1,'H',H[-1,i],'water')
                
                S = S*1e-3
                H = H*1e-3
                
                for i in range(len(H[:,0])):
                    plt.plot(S[i,:],H[i,:],color='green',linewidth=2)
            
            for i in range(len(self.states)):
                
                plt.plot(self.states[states_keys_list[i]][0]*1e-3,self.states[states_keys_list[i]][2]*1e-3,'o',color='black',linewidth=4)
                plt.text(self.states[states_keys_list[i]][0]*1e-3+0.2,self.states[states_keys_list[i]][2]*1e-3,states_keys_list[i],fontsize=14,family='serif')
                
            ## Power plot ##
            
            self.power_fig = None
            
            if self.qm != None:
                
                losses = np.array([abs((1-self.eta_mec)*self.Wt - self.Wp),self.Qc,(1/self.eta_VG-1)*self.Qvg]) * self.qm
                power_sum = np.concatenate((np.array([self.P_out]),losses))
                power_sum = power_sum*1e-6
                self.power_fig = plt.figure(3,figsize=(8,6))
                labels = ['Net power output\n{:.2f} MW'.format(power_sum[0]),'Mechanical losses\n{:.2f} MW'.format(power_sum[1]),'Condenser losses\n{:.2f} MW'.format(power_sum[2]),'Vapor generator losses\n{:.2f} MW'.format(power_sum[3])]
                plt.pie(power_sum,labels=labels,autopct='%1.1f%%',textprops={'fontsize':14,'family':'serif'})
                plt.title('Primary Energy flux {0:.2f} MW'.format(np.sum(power_sum)),fontsize=17,family='serif',weight='bold')
            
    def return_data(self):
        self.metrics = [self.P_out,self.qm,self.eta_th,self.eta_carnot,self.Wt,self.Wp,self.Qvg,self.Qc]
        return self.ts_fig , self.hs_fig , self.power_fig , self.states , self.metrics
        
#======================================#

#%% Test %%#

# # Example data = 'P_out' :         200e6
# #                'qm' :            None
# #                'T_max' :         540+273.15
# #                'T_min' :         30+273.15
# #                'p_max' :         140e5
# #                'x_min' :         0.89
# #                'eta_isent_T' :   0.91
# #                'eta_isent_P' :   0.85
# #                'eta_mec' :       0.98
# #                'eta_VG':         0.92
# #                'p_reheat' :      35e5
# #                'nbr_bleedings' : 5

# cycle=RH_cycle(900, 300, 180e5,P_out=500e6)
# cycle.evaluate()

# t,h,j=cycle.return_plots()

# plt.plot(t)
# plt.plot(h)



# #======================================#
