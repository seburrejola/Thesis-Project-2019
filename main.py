import numpy as np
import random
import time
import copy 
from datetime import datetime

import lifeparticle_force as lf_force
import classes_and_functions as caf

'''
The imported code: 'lifeparticle_force.py' has functions to calculate and to 
set the forces of interaction in the system.
The imported code: 'classes_and_functions.py' has the classes and some
functions used in the code.

The explanation of the each function will be in each code.


%%%%%%% SIMULATIONS PARAMETERS %%%%%%%
'''
dt=0.01             #Parameter for the movement of the particles
pl=[]               #List for the manage of the particles of the system
N=100               #Initial Number of particles
N_max=100           #MaxNumber of particles
steps=1000          #Number of steps for the simulation
frame=50            #Parameter for the manage of frames thats the animation will have.
mass=1              #The mas for each particle
colors=7            #Number of types of praticles, max=7
initial_colors=3    #Number of initial colors for the simlation

delta_t=50          #parameter for the meassurement of the local clustering for the mutations of the particles

theta_div= steps + 100   #Division lifecounter threshold
theta_dif= steps + 150   #Differentiation lifecounter threshold

#nu=3 el nu=3 de antes da dt*6*nu*np.pi=0.56, con nu=5 me daba 0.94
nu=0.9          #Viscosity coefficient 
boundary=20     #Size of the system
noise=0#4       #Parameter for the random motion 
d_select=1      #Parameter to select de % of forces that will be attractive


writer1=open('test.txt','w') #The txt document where the code is going to save the simulation info.

#The next 6 lines are for the meassurement of the coding time, also to aproximate the time it will take.
start=time.time()
now=datetime.now()
dt_string=now.strftime('%H:%M:%S')
print(dt_string)
print('Maximum estimated computation time in minutes:',(steps/1000)*((N_max/100)**2))
print('Maximum estimated computation time in hours:',((steps/1000)*((N_max/100)**2))/60)


'''
%%%%%%% INITIAL CONDITIONS %%%%%%%
'''

lpos=[[],[],[],[]]      #list to save the information of the particles
boundary_init=[int(boundary/3),int(boundary*(2/3))] #square where the particles are placed 

caf.initial_conditions(N,initial_colors,boundary_init,mass,pl)
caf.set_lpos(lpos,pl)

'''
%%%%%%% FORCE PARAMETERS %%%%%%%
'''

force_parameters=[]         #Lis for the parameters of each force of interaction
action_radio=1.5            #maximum radious of action
lf_force.random_parameters(force_parameters,colors,action_radio,d_select)      #Function to randomly fill force_parameters
#force_parameters=lf_force.set_parameters1_force()                             #Function to set force_parameters

force_parameters[0][0]=[0.35,0.5, -1200]
#force_parameters[1][1]=[0.35,0.5, -1200]

'''
%%%%%%% MUTATIONS FUNCTIONS %%%%%%%
'''

molecular_rad=0.4                         #Clustering check radius
mutant_colors=np.array([0,1,2,3,4,5,6])   #Types of particles that can mutate
mutate_rules=caf.set_mutation_rules()     #Randomly set the mutation rules
#mutate_rules[1][0]=random.randint(1,6)
#mutate_rules=lf_force.set_parameters1_mutate()  

def step_function(step,boundary,nu,frame,dt,N_max,delta_t,steps,theta_div,theta_dif):
        
    '''
    %%%%%%% MUTATION FUNCTIONS %%%%%%%
    In the next two 'if's we update the lists for the local clustering meassurements
    '''
    
    if len(pl)<N_max and step>100:
        nu=0.95
        e1=random.randint(0,len(pl)-1)
        if pl[e1].life_cont>=theta_div and random.uniform(0,1)<0.1:
            caf.testing_div_func(e1,pl,lpos,mass)
            if len(pl)==N_max:
                print('len(pl) =',len(pl),',| Step =',step)
                print('Remaining computation time (in minutes):',((steps-step)/1000)*((N_max/100)**2))
    
    if step%delta_t==0:#step==2:#
        for i in range(len(pl)):
            pl[i].cluster_check1_color(pl,i,molecular_rad)
    
    if (step+6)%delta_t==0:
        for i in range(len(pl)):
            pl[i].cluster_check2_color(pl,i,molecular_rad)
            if len(pl)>N_max*0.5 and pl[i].life_cont>=theta_dif+step//delta_t and random.uniform(0,1)<0.05:#and pl[i].color in mutant_colors
                
                pl[i].mutate(i,lpos,pl,colors,molecular_rad,mutate_rules)
        
    '''
    %%%%%%% MOVEMENT FUNCTIONS %%%%%%%
    In the next 'double for' the code calculates all the change in velocity due
    to the force interactions between particles and the random motion.
    
    Then in the second 'for' all the changes of the positions of the particles
    are calculated.
    '''
    
    for i in range(len(pl)):
        for j in range(len(pl)):
            if i!=j :
                radioo=np.sqrt((pl[i].x - pl[j].x)**2 +(pl[i].y - pl[j].y)**2 )
                c=-force_parameters[pl[i].color][pl[j].color][0] + 2 * force_parameters[pl[i].color][pl[j].color][1]
                if radioo < c:
                    pl[i].force_interaction(pl[j].x,pl[j].y,action_radio,pl[j].color,force_parameters,dt,radioo)
                    
        pl[i].random_vel(noise) # Random motion
        
    for l in range(len(pl)):
        pl[l].life_cont+=1
        pl[l].move_dt(dt,boundary,nu)
    
    '''
    %%%%%%% SAVING THE SIMULATON INFO %%%%%%%
    
    First we update the list with the information of the particles and then 
    we write it in the txt we opened at the begining of the code
    '''
    
    if step%frame==0:
        caf.update_scatt(lpos,pl)
        N=copy.copy(len(pl))
        writer1.write(str(N))
        writer1.write(',')
        for i in range(len(lpos)):
            writer1.write('|,')
            for j in range(len(lpos[i])):    
                writer1.write(str(lpos[i][j]))
                writer1.write(',')
        writer1.write('\n')
        
'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

# In the next loop the simulation is done.
for step in range(steps):
    step_function(step+1,boundary,nu,frame,dt,N_max,delta_t,steps,theta_div,theta_dif)



end=time.time() 
print('Computation time in minutes:',(end-start)/60) # Here we print how long the simulation took
writer1.close() # We close the document
print('len(pl) =',len(pl))
