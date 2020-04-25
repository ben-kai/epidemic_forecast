# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 20:46:39 2020

@author: benjamin.kaiser@kaiserengineering.de
"""

import numpy as np
import matplotlib.pyplot as plt
##############################################################################
# country variables
# Germany 
p_count = 82000000
# number of diseased at time 0
inf_count = 100
# number of days to recover from disease
day_of_diseas = 12


# time array for the transmission parameter. How the transmission parameter changes
transm_count=np.array([[0.25   ,   0.25   ,  0.13,    0.13,   0.10,    0.07,     0.05,     0.04] , #transmission rate in % per day
                       [0.00  ,  16.0  ,     17.0,    24.0,   25.01,    31.0,   50.01,   100.0]])  #time        in days

recov_count=np.array([[0.00   ,            0.00] ,  #recovery rate in % per day
                       [0.00  ,  day_of_diseas-1]]) #time        in days
##############################################################################
# initialize SRI Model in 100% of population -> N=1
N=1
Istart=inf_count/p_count
Sstart=N-Istart
Rstart=0

# set SRI start parameters
S=Sstart
I= Istart
R=Rstart

# set time values in days
t_end=400 # Days
dt    = 0.1
t     = 0

# initial History Variables
time     = [t]
Susc     = [S*p_count]
Infe     = [I*p_count]
Reco     = [R*p_count]
Conf     = [inf_count]
dConf    = [1]
dR       = [0]

# day counter
day=1;
j=day_of_diseas
while t < t_end:
    t=t+dt
    if t > day_of_diseas:
        recov_count=np.append(recov_count,[[( dConf[-day_of_diseas])/(I*p_count)],[j]], axis=1)  
        j+=1
    transm   = np.interp(t,transm_count[1,:],transm_count[0,:])     
    recov   = np.interp(t,recov_count[1,:],recov_count[0,:])   
    
    #print(t , recov)
    S_p=-transm*S*I
    I_p=transm*S*I-recov*I
    R_p= recov*I
    S=S+S_p*dt
    I=I+I_p*dt
    R=R+R_p*dt
    #print(R*p_count,I*p_count)
    if t > day:

        # store variables
        time.append(day)
        Susc.append(S*p_count)
        Infe.append(I*p_count)
        Reco.append(R*p_count)
        Conf.append(I*p_count+R*p_count)
        dConf.append((I*p_count+R*p_count)-Conf[-2])
        dR.append(R*p_count-Reco[-2])
        #print(day, R*p_count,'hello')
        day+=1

#plt.plot(time, Susc, label='Simulation Susceptible')
plt.plot(time, Infe, label='Simulation Diseased')
plt.plot(time, Reco, label='Simulation Removed')
plt.plot(time, Conf, label='Simulation Confirmed')

plt.legend(loc='upper left')
#    plt.plot(time, force[2])
#    plt.plot(time, force[3])
plt.grid(True)
plt.show()
#    
