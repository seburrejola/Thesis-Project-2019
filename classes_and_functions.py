import numpy as np
import random 
import copy
import itertools as it

import lifeparticle_force as lf_force

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
        Fdrag = -6 * np.pi * self.R * nu
        self.vx += (dt/self.m)*Fdrag*self.vx
        self.vy += (dt/self.m)*Fdrag*self.vy
              
    def cluster_check1_color(self,pl,agent_j,molecular_rad):
        '''
        This function is the first fucntion to check the particle neighborhood 
        for the local measurment of a cluster.
        '''
        self.old_is_color[1]=[]
        for index in self.old_is_color[0]:
            self.old_is_color[1].append(pl[index].color)
                
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
        colors_list=['b','r','m','g','k','c','y'] 
        self.cluster_check2_color(pl,agent_j,molecular_rad)
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
            self.old_is_color[1].append(pl[index].color)


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
        


def add_agent(lpos,colors,mass,pl):
    '''
    This function adds a new particle to the system. For now it just add the 
    particle in the middle of the space, and it will have color = 6 or = 'y'.
    '''
    colors_list=['b','r','m','g','k','c','y']
    color=6#random.randint(1,colors)-1
    velocity=random.uniform(0,1)
    angle=random.uniform(0,1)
    xrand=20
    yrand=10
    vx = velocity * np.cos( 2* np.pi * angle)
    vy = velocity * np.sin( 2* np.pi * angle)
    #print(vx,vy,vx**2 + vy**2 ,velocity**2)
    #mass=random.uniform(0.1,1)
    pl.append(particle(xrand,yrand,vx,vy,mass,color))
    lpos[0].append(xrand)
    lpos[1].append(yrand)
    lpos[2].append(colors_list[color])
    lpos[3].append(15)

    

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




'''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
The next functions are for the code read_write_txt.py
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''



def set_old_new_is(N,old_is,new_is):
    '''
    This function sets two lists for the global clustering measurment.
    '''
    for i in range(N):
        old_is.append([[],[]])
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
    (0).-This function sets the list 'cl' into a list of sets. Each set 
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
            if lpos[2][index] in set(new_is[i][1])&set(old_is[i][1]) and index not in aux0:
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
        if r not in cl:
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

def save_freq_plot2(lista,cl):
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
    
                