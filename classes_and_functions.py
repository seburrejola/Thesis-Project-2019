import numpy as np
import random 
import copy
import itertools as it

import lifeparticle_force2 as lf_force

'''
Particle class and it functions.
'''
class particle:
    def __init__(self,x,y,vx,vy,mass,color):
        self.x=x
        self.y=y
        self.pos=(self.x,self.y)
        self.vx=vx
        self.vy=vy
        self.m=mass
        self.velocity=np.sqrt((vx**2)+(vy**2))
        self.R=1  #Particle radio.
        self.color=color
        self.new_is_color=[[],[]]
        self.old_is_color=[[],[]]
        self.life_cont=0

                
    def move_dt(self,dt,boundary,nu):
        '''
        This function calculates the change in the position of the particle.
        '''
        self.friction_force(nu,dt)
        if self.x + self.vx*dt > boundary*2 or self.x + self.vx*dt < 0:
            self.vx*=-1
        if self.y + self.vy*dt > boundary or self.y + self.vy*dt < 0:
            self.vy*=-1
            
        self.x += self.vx*dt
        self.y += self.vy*dt
        self.pos=(self.x,self.y)
    
    def random_vel(self,noise):
        '''
        This function calculates the change in velocitiy of the particle due to
        random motion / heat, or noise, of the system.
        '''
        velocity = random.uniform(0,1)*noise
        angle = random.uniform(0,1) 
        self.vx += velocity * np.cos( 2* np.pi * angle)
        self.vy += velocity * np.sin( 2* np.pi * angle)

          
    def force_interaction(self,x,y,action_radio,color_i,force_parameters,dt,radio):
        '''
        This function calculates the change in velocitiy of the particle due to
        interactions with an other particle of the type "color_i" and with a
        distance between them "radio"
        '''
        
        theta=np.arctan2((y - self.y),(x - self.x ))   
        a=force_parameters[self.color][color_i][0]
        b=force_parameters[self.color][color_i][1]
        d=force_parameters[self.color][color_i][2]
        force=lf_force.force_calculation(a,b,d,radio) 
        self.vx += (dt/self.m)*force*np.cos(theta)
        self.vy += (dt/self.m)*force*np.sin(theta)
                

    
    def friction_force(self,nu,dt):
        '''
        This function calculates the change in velocitiy of the particle due to
        the viscous drag.
        '''
        #Fdrag = -6 * np.pi * self.R * nu
        #dt*6*nu*np.pi=0.94 con nu=5,R=1 y mass=1
        
        #self.vx += (dt/self.m)*Fdrag*self.vx
        #self.vy += (dt/self.m)*Fdrag*self.vy
        
        self.vx += -nu*self.vx
        self.vy += -nu*self.vy
              
    def cluster_check1_color(self,pl,agent_j,molecular_rad):
        '''
        This function is the first fucntion to check the particle neighborhood 
        for the local measurment of a cluster.
        '''
        self.old_is_color=[[],[]]
        #self.old_is_color[1]=[]
        #for index in self.old_is_color[0]:
         #   self.old_is_color[1].append(pl[index].color)
                
        for i in range(len(pl)):
            if i!=agent_j:
                radio=(np.sqrt(((self.x - pl[i].x)**2)+(self.y - pl[i].y)**2))
                if radio<molecular_rad and i not in self.old_is_color[0]:
                    self.old_is_color[0].append(i)
                    self.old_is_color[1].append(pl[i].color)

    def cluster_check2_color(self,pl,agent_j,molecular_rad):
        '''
        This function is the second fucntion to check the particle neighborhood 
        for the local measurment of a cluster.
        '''
        self.new_is_color=[[],[]]
        for i in range(len(pl)):
            if i!=agent_j:
                radio=(np.sqrt(((self.x - pl[i].x)**2)+(self.y - pl[i].y)**2))
                if radio<molecular_rad and i not in self.new_is_color[0]:
                    self.new_is_color[0].append(i)
                    self.new_is_color[1].append(pl[i].color)
                    
        
    def mutate(self,agent_j,lpos,pl,colors,molecular_rad,mutate_rules):
        '''
        This function compares the two measurements of local clusters, if the 
        particle have measured the same type of praticle in both measurments, 
        then the particle will be able to mutate to a new type depending on the
        type of this particle and one of the measured particles.
        '''
        self.life_cont=0
        colors_list=['b','r','m','g','k','c','y'] 
        #self.cluster_check2_color(pl,agent_j,molecular_rad) #esto lo cambio cuando estaba trabajando en el informe
        aux2=list(set(self.new_is_color[1])&set(self.old_is_color[1]))
        if len(aux2)!=0:
            
            aux=(self.color,random.choice(aux2))
            #aux=(self.color,pl[list(set(self.new_is)&set(self.old_is))[0]].color)
            if aux in mutate_rules[0]:
                #self.color=pl[list(set(self.new_is)&set(self.old_is))[0]].color
                self.color=mutate_rules[1][mutate_rules[0].index(aux)]
           
            lpos[2][agent_j]=colors_list[self.color]

        self.old_is_color[0]=copy.copy(list(set(self.new_is_color[0])&set(self.old_is_color[0])))
        for index in self.old_is_color[0]:
            #pl[index].life_cont=0
            self.old_is_color[1].append(pl[index].color)

class Cluster_class:
    def __init__(self,particles,aparition_step):
        self.particles = particles
        self.frame_cont = 1
        self.aparition_step = aparition_step
        self.lost_agents=[]
        self.new_agents=[]
        self.old_particles=[]
        self.status='alive'
        self.updated=0
        
    def update_particles(self,cluster,step):
        dif=cluster^self.particles
        self.old_particles.append((copy.copy(self.particles),step))
        self.lost_agents.append((self.particles & dif,step))
        self.new_agents.append((cluster & dif,step))
        self.particles=cluster
        self.frame_cont+=1
        self.updated=1
        
        
def save_cluster(cl,ccl,step):
    aux_ccl=[[],[]]

    updated_clusters=[]
    for i in range(len(ccl)):#loop para generar los testings, basicamente armas cl[].step-1 en aux_ccl
        if ccl[i].status=='alive':
            aux_ccl[0].append(ccl[i].particles)
            aux_ccl[1].append(i)
            ccl[i].status='testing'
       
    
    if cl_check(aux_ccl[0])==0 and len(aux_ccl[0])>0:
        print('fuck me',step, aux_ccl[0])
    
        
    
    for j in range(len(cl)):
        not_sharing_cont=0
        if len(aux_ccl[0])==0:#por si no tengo nada alive
            break
        
        if cl[j] in aux_ccl[0]:#gate para aumentar el frame_cont si esque sigue igual que antes
            ccl[aux_ccl[1][aux_ccl[0].index(cl[j])]].frame_cont+=1
            ccl[aux_ccl[1][aux_ccl[0].index(cl[j])]].status='alive'
            updated_clusters.append(cl[j])
            
        if cl[j] not in aux_ccl[0]:
            for l in range(len(aux_ccl[0])):#loop para agregar los cluster que comparter particulas y los que no comparten ninguna
                if cl[j]&aux_ccl[0][l]!=set() and ccl[aux_ccl[1][l]].status=='testing': #compartir -> ak & ak2 != set()
                    if cl[j] not in updated_clusters: 
                        sharing_agents(cl[j],ccl,step,aux_ccl[1][l],updated_clusters)
                if cl[j]&aux_ccl[0][l]==set() and ccl[aux_ccl[1][l]].status=='testing': #no compartir -> ak & ak2 == set()
                    #ccl.append(Cluster_class(cl[j],69))
                    not_sharing_cont+=1
        
        #if not_sharing_cont==len(aux_ccl[0]):
        #    ccl.append(Cluster_class(cl[j],step+1000))
         #   updated_clusters.append(cl[j])
            #not_sharing_cont=0
            
    for m in aux_ccl[1]:#loop de muerte a quien sigue en testing
        if ccl[m].status=='testing' and ccl[m].updated==0:
            ccl[m].status='died in the step ' + str(step)
        if ccl[m].status=='testing' and ccl[m].updated==1:
            ccl[m].status='alive'
            ccl[m].updated=0
        
        
    for cluster_set in cl:
        if cluster_set not in updated_clusters:
            ccl.append(Cluster_class(cluster_set,step+1000))
        
        
    if len(ccl)==0:# or len(aux_ccl[0])==0:# and len(cl)!=0:#loop para llenar ccl
        if len(cl)!=0:
            for n in range(len(cl)):
                ccl.append(Cluster_class(cl[n],'pene'+str(step)))
                    
                    
    
    

def sharing_agents(cluster,ccl,step,j,updated_clusters):
    if len(cluster&ccl[j].particles)>=len(ccl[j].particles)*0.5:
        #print(j,step)
        ccl[j].update_particles(cluster,step)
        updated_clusters.append(cluster)
        
    if len(cluster&ccl[j].particles)<=len(ccl[j].particles)*0.5:
        #print('ctm')
        ccl.append(Cluster_class(cluster,420))
        ccl[len(ccl)-1].status='alive'#'shared agents in the step '+str(step)+' with cluster '+str(j)
        ccl[len(ccl)-1].old_particles.append((ccl[j].particles,step))
        updated_clusters.append(cluster)

   


def initial_conditions(N,colors,boundary_init,mass,pl):
    '''
    This function initialize the system. 
    
    Here we set the particle list (pl), with N amount of particles with mass
    = mass, and placed inside the rectangle with coordinates obtained with
    boundary_init. Also this particles have a type (that refers to its color)
    taking a number between 0 and the input "colors".    
    '''
    for i in range(N):
        color=random.randint(1,colors)-1
        velocity=random.uniform(0,1)
        angle=random.uniform(0,1)
        xrand=random.uniform(boundary_init[0],boundary_init[1])*2
        yrand=random.uniform(boundary_init[0],boundary_init[1])
        vx = velocity * np.cos( 2* np.pi * angle)
        vy = velocity * np.sin( 2* np.pi * angle)
        #print(vx,vy,vx**2 + vy**2 ,velocity**2)
        #mass=random.uniform(0.1,1)
        pl.append(particle(xrand,yrand,vx,vy,mass,color))
        


'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''    
def add_agent(lpos,color,mass,pl,agent_j):
    pl[agent_j].life_cont=0
    #for index in pl[agent_j].old_is_color[0]:
     #   pl[index].life_cont=0
        
    colors_list=['b','r','m','g','k','c','y']
    color=color#random.randint(1,colors)-1
    velocity=random.uniform(0,1) #podria aca poner 0 o la mitad del momento que tenia la celula madre.
    angle=random.uniform(0,1)
    e1=random.uniform(0,1)
    xrand=pl[agent_j].x+0.35*e1
    yrand=pl[agent_j].y+0.35*(1-e1)
    vx = velocity * np.cos( 2* np.pi * angle)
    vy = velocity * np.sin( 2* np.pi * angle)
    #print(vx,vy,vx**2 + vy**2 ,velocity**2)
    #mass=random.uniform(0.1,1)
    pl.append(particle(xrand,yrand,vx,vy,mass,color))
    lpos[0].append(xrand)
    lpos[1].append(yrand)
    lpos[2].append(colors_list[color])
    lpos[3].append(15)

def testing_div_func(agent_j,pl,lpos,mass):
    #if len(set(pl[agent_j].new_is_color[0])&set(pl[agent_j].old_is_color[0]))>0:
        #print('hola')
    colorr=pl[agent_j].color
    if random.uniform(0,1)<0.01:
        colorr=random.randint(0,6)
    add_agent(lpos,colorr,mass,pl,agent_j)
        
'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
def set_lpos(lpos,pl):
    '''
    This function set the list 'lpos', where lpos=[[],[],[],[]], this list is
    going to have the information of the particls per step and it will be used
    to work and save the data of the simulation.
    
    In lpos[0] we are going to have a list with all the position X of all the 
    particles, while in lpos[1] we are going to have all the positions Y.
    In lpos[2] we will have the color of each particle and in lpos[3] the size
    of the particle. 
    
    The size of the paritlces does not depend on he boundarys of the system, 
    this number just represent the size of thefigure we are going to plot with
    and scatter plot.
    '''
    for i in range(len(pl)):
        if pl[i].color==0:
            lpos[2].append('b')
            
        if pl[i].color==1:
            lpos[2].append('r')
            
        if pl[i].color==2:
            lpos[2].append('m') 
            
        if pl[i].color==3:
            lpos[2].append('g')
            
        if pl[i].color==4:
            lpos[2].append('k')
            
        if pl[i].color==5:
            lpos[2].append('c')
            
        if pl[i].color==6:
            lpos[2].append('y')
           
        lpos[0].append(pl[i].x)    
        lpos[1].append(pl[i].y)
        lpos[3].append(15)
        

        
def update_scatt(lpos,pl):
    '''
    This function updates the positions of the particles in the list lpos.
    '''
    for i in range(len(pl)):
        lpos[0][i]=pl[i].x
        lpos[1][i]=pl[i].y
    
    
    
def set_mutation_rules():
    '''
    This function randomly set the mutation rules.
    
    In this case we set a list of tuples with all the possible combinations
    of pairs pf typs of partticles. Then we set an other list of the same size, 
    that will have the resulting type for the mutating paritcle. At the end,
    this function will return a list with both lists previously explained.
    '''
    mutate_rules=[]
    mutate_rules.append(list(it.product((0,1,2,3,4,5,6),repeat=2)))
    mutate_rules.append([])
    
    cont=0
    while len(mutate_rules[1])<len(mutate_rules[0]):
        e=random.randint(1,6)
        #if e!=mutate_rules[0][cont][0] and e!=mutate_rules[0][cont][1]:
        mutate_rules[1].append(e)
        cont+=1
    return mutate_rules
    
   


def clustering_radios(cluster_list,pl):
    '''
    This function prints the average radio of the clusters in the input 
    'cluster_list'
    '''
    for i in range(len(cluster_list)):
        if len(cluster_list[i])!=1:
            radios=[]
            radio_function(radios,list(cluster_list[i])[0],cluster_list[i],pl)
            print(i,np.mean(radios))



def step_function(step,boundary,nu,frame,dt,N_max,delta_t,steps,theta_div,theta_dif,pl,lpos,mass,noise,molecular_rad,action_radio,mutate_rules,colors,force_parameters,writer1):
        
    '''
    %%%%%%% MUTATION FUNCTIONS %%%%%%%
    In the next two 'if's we update the lists for the local clustering meassurements
    '''
    
    if len(pl)<N_max and step>100:
        
        e1=random.randint(0,len(pl)-1)
        if pl[e1].life_cont>=theta_div and random.uniform(0,1)<0.1:
            testing_div_func(e1,pl,lpos,mass)
            if len(pl)==N_max:
                print('len(pl) =',len(pl),',| Step =',step)
                print('Remaining computation time (in minutes):',((steps-step)/1000)*((N_max/100)**2))
    
    if step%delta_t==0:#step==2:#
        for i in range(len(pl)):
            pl[i].cluster_check1_color(pl,i,molecular_rad)
    
    if (step+6)%delta_t==0:
        for i in range(len(pl)):
            pl[i].cluster_check2_color(pl,i,molecular_rad)
            if len(pl)>N_max*0.5 and pl[i].life_cont>=theta_dif and random.uniform(0,1)<0.05:#and pl[i].color in mutant_colors
                
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
        update_scatt(lpos,pl)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The next functions are for the code read_simulation_data.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''



def set_old_new_is(N,old_is,new_is,tipo):
    '''
    This function sets two lists for the global clustering measurment.
    '''
    if tipo=='old':
        for i in range(N):
            old_is.append([[],[]])
    if tipo=='new':
        for i in range(N):
            new_is.append([])
            



def clustering_check1(N,lpos,old_is,molecular_rad):
    '''
    This function is an adaptation of the 'particle' class fucntion: 
    particle.cluster_check1_color(...), it works in the same way but now
    this functions just need the list lpos and works outside of the 'particle'
    class object.
    '''
    
    for i in range(N):
        old_is[i][1]=[]
    for j in range(N):
        for index in old_is[j][0]:
            old_is[j][1].append(lpos[2][index])
    for l in range(N):
        for m in range(N):
            if l!=m:
                radio=np.sqrt(((lpos[0][l] - lpos[0][m])**2)+((lpos[1][l] - (lpos[1][m]))**2))
                
                if radio<molecular_rad and m not in old_is[l][0]:
                    old_is[l][0].append(m)
                    old_is[l][1].append(lpos[2][m])
                    
                    
def clustering_check2(N,lpos,new_is,molecular_rad):
    '''
    This function is an adaptation of the 'particle' class fucntion: 
    particle.cluster_check2_color(...), it works in the same way but now
    this functions just need the list lpos and works outside of the 'particle'
    class object.
    '''
    
    for i in range(N):
        new_is[i]=[[],[]]
    for l in range(N):
        for m in range(N):
            if l!=m:
                radio=np.sqrt(((lpos[0][l] - lpos[0][m])**2)+((lpos[1][l] - (lpos[1][m]))**2))
                if radio<molecular_rad and m not in new_is[l][0]:
                    new_is[l][0].append(m)
                    new_is[l][1].append(lpos[2][m])
                    
            
def clusters_means(cluster_list,lpos):
    '''
    In this functions enters a list a list where the elements represent the 
    clusters measured in the system and list that has the information of the 
    particles per step.
    
    With both lists, this functions makes a new list with the average position
    in X and Y of the clusters, also with the size of the clusters (equal to 
    the maximum radius between particles in each cluster), and it sets a color 
    for the cluster. Finally this function returns a list with the same format  
    of lpos with the purpose of ploting the clustering information.
    '''
    if len(cluster_list)==0:
        return []
    cl_means=[[],[],[],[]]
    for i in range(len(cluster_list)):
        
        if len(cluster_list[i])!=1:
            radios=[]
            xlist=[]
            ylist=[]
            for j in cluster_list[i]:
                xlist.append(lpos[0][j])
                ylist.append(lpos[1][j])
                radio_function(radios,j,cluster_list[i],lpos)
            cl_means[0].append(np.mean(xlist))
            cl_means[1].append(np.mean(ylist))
            cl_means[2].append('k')
            cte=900/(np.pi*0.25*(1.64)**2)
            cl_means[3].append(cte*(0.2+np.max(radios))**2)
        
    return cl_means
            
            
        
    
    
    
    

def radio_function(radios,agent_j,cluster,lpos):
    '''
    This function raturns a set with the radius between the particle 'agent_j'
    and all the other particles in the cluster.
    '''
    for i in cluster:
        if i!=agent_j and np.sqrt(((lpos[0][agent_j] - lpos[0][i])**2)+(lpos[1][agent_j] - lpos[1][i])**2) not in radios:
            radios.append(np.sqrt(((lpos[0][agent_j] - lpos[0][i])**2)+(lpos[1][agent_j] - lpos[1][i])**2))
    radios=set(radios) 


            
              
        
           
   
    
def clusters(N,lpos,molecular_rad,cl,old_is,new_is):
    '''
    (0).-This function builds the list 'cl' into a list of sets. Each set 
    represents a cluster, so in each set we are going to have a group of 
    numbers that are the identification number for each particle.
    
    (1).-In the first part of this function, a list 'cl_aux' is filled with the 
    neighborhood of each particle. 
    
    (2).-Then in 'cl_aux2' in each 'j' position of the list 'cl_aux' (i.e, for 
    each 'j' particle) we put the neighborhood of the particles of the 
    neighborhood of the particle 'j'.
    
    (3).-Then we start cleaning the 'cl_aux2' list. To do this we put the 
    elements of this list in a new list 'cl_aux3', but now when we are filling 
    'cl_aux3' we do not put elements that are already in the list 'cl_aux3'.
    
    (4).-Then we check if the cl_aux3 list has clusters that shares particles, 
    in other words, we check if all the clusters have different particles, to  
    do this we loop the steps (2) and (3) in orther to obtain just sets of 
    particles with each neighbor and all the neighbors of the neighbors.
    
    (5).-When the step (4) is done, 'cl_aux3' will be a list of sets with the 
    identification number of each particle, also none of these sets shares any
    particle.   
    '''
    cl_aux=[] #vecinos
    clustering_check2(N,lpos,new_is,molecular_rad)
    for i in range(N): 
        aux0=[]
        
        
        for index in new_is[i][0]:
            if index<(N-1) and lpos[2][index] in set(new_is[i][1])&set(old_is[i][1]) and index not in aux0:
                aux0.append(index)
        
        cl_aux.append(set(aux0))
    
    cl_aux2=[] #vecinos de vecinos
    for j in range(N): 
        aux_set=set()
        if cl_aux[j]!=set():
            aux_set=copy.copy(cl_aux[j])
            for index in cl_aux[j]:
                aux_set|=cl_aux[index]
            cl_aux2.append(aux_set)
            
    #if len(cl_aux2)==0:
     #   print('cl_aux len = 0')   
    
    cl_aux3=[]
    for k in range(len(cl_aux2)):
        if cl_aux2[k] not in cl_aux3:
            cl_aux3.append(cl_aux2[k])


    while cl_check(cl_aux3)==0: 
        if len(cl_aux2)==0:
            break
        
        
        for j in range(len(cl_aux2)): 
            
            aux_set=copy.copy(cl_aux2[j])
            for index in cl_aux2[j]:
                for l in range(len(cl_aux2)):
                    if index in cl_aux2[l] and l!=j:
                        aux_set|=cl_aux2[l]
            cl_aux2[j]|=aux_set        
         
        cl_aux3=[]          #clusters + fragmented clusters   
        for k in range(len(cl_aux2)):
            if cl_aux2[k] not in cl_aux3:
                cl_aux3.append(cl_aux2[k])


        
    for r in cl_aux3:       #The clusterd well measured are added to cl
        if r not in cl and len(r)>1:
            cl.append(r)
            
            

def cl_check(cl):
    '''
    The input 'cl' (cluster list) needs to have the format of a cluster list,
    in other words, a list of sets, where in this sets we are going to have
    the identification number for each particle.
    
    This function check if the elements of the list 'cl' share elements between
    them.This function will return 0 if the elements of 'cl' share elements
    between them meaning that the list 'cl' has an error, because at least one
    particle is shared between clusters.
    '''
    aux_list=[]
    if len(cl)==0:
        return 0
    if len(cl)!=0:
        for i in range(len(cl)):
            for j in range(len(cl)):
                #if j!=i and cl[i]&cl[j]==set():
                  #  return 1
                if j!=i and cl[i]&cl[j]!=set():
                    aux_list.append(i)
                    
    
    if len(aux_list)!=0:
        return 0
       



    
def save_clusters_plot(array_x,cl,step):
    '''
    This function updates the clusters per step array. 
    
    This array will be used to show how many cluster we had each time the 
    system took a measurment.
    '''
    array_x[step,0]=len(cl)

def save_freq_plot1(lista,cl):
    '''
    This function set a list with the sizes of the cluster measured.
    '''
    aux_array=np.zeros((len(cl),1))
    for i in range(len(cl)):
        aux_array[i,0]=len(cl[i])
        
    lista.append(copy.copy(aux_array))

def save_freq_plot2(lista):
    '''
    This function sets and array that is going to be used to obtain the plot
    of 'frequency of clusters' vs 'size of clusters'.
    '''
    cont=0
    cont2=0
    for i in range(len(lista)):
        cont+=len(lista[i])
        
    array_x=np.zeros((cont,1))
    for j in range(len(lista)):
        for l in range(len(lista[j])):
            array_x[cont2,0]=lista[j][l]
            cont2+=1
            
    return array_x
    
def save_freq_plot3(ccl):
    array_x=np.zeros((len(ccl),1))
    for i in range(len(ccl)):
        array_x[i,0]=copy.copy(len(ccl[i].particles))
    
    return array_x
        
        

def reading_lpos(lpos,N,lineas,colors_list,old_is,new_is,step):
    linea=lineas[step].split(',')
    linea.remove('\n')
    N=int(linea[0])
    del linea[0]
    
    
    
    #lpos2=[[],[],[],[]]
    
    set_old_new_is(N,old_is,new_is,'new')
    cont1=0
    cont2=0
    for element in linea:
        if element!='|':
            if element not in colors_list:
                lpos[cont1].append(float(element))
                cont2+=1
            if element in colors_list:
                cont2+=1
                lpos[cont1].append(element)
        if cont2==N:
            cont1+=1
            cont2=0


def average_ls_per_sim(big_ccl,life_spams,average_ls):
    cont=0
    for ccl in big_ccl:
        life_spam=[]
        for cluster in ccl:
            life_spam.append(cluster.frame_cont)
        average_ls.append([np.mean(life_spam),np.std(life_spam)])
        print('Sim',cont,'had a average and an std of', average_ls[cont],'and a max of',max(life_spam))
        cont+=1
        life_spams.append(life_spam)
        
    ak=big_ccl[0]+big_ccl[1]+big_ccl[2]+big_ccl[3]+big_ccl[4]
    ak2=[]
    for cluster in ak:
        ak2.append(cluster.frame_cont)
        
    print('Scenario average life spam',np.mean(ak2),', with an std of',np.std(ak2))
        
        
        



            