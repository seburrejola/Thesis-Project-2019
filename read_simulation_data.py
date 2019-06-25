import matplotlib.animation as animation
import matplotlib.pyplot as plt 
import numpy as np
import copy as copy

import classes_and_functions as caf

'''
This code will open the data previously saved of a simulation and it will
animate it and also obtain the clustering data.
'''



#Here we have some parameters of the simulation.
boundary=20     
colors_list=['b','r','m','g','k','c','y'] 



#In the next linea we make the figure where we are going to do the scatter plot.
fig2 = plt.figure(figsize=(10, 5))
ax2 = fig2.add_axes([0, 0, 1, 1], frameon=False)
ax2.set_xlim(0, boundary*2), ax2.set_xticks([])#EXTREMOS DEL PLOT
ax2.set_ylim(0, boundary), ax2.set_yticks([])#EXTREMOS DEL PLOT





# We open and then we read the data
file=open('simulation1.txt','r')
lineas=file.readlines()
ims1=[]
file.close()



cl=[]                                         #Cluster list
cl_means=[]                                   #Clusters average values
clusters_per_step=np.zeros((len(lineas) ,1))  #Array to plot the amount of clusters per step.
frequency_list=[]                             #List to plot the frequency vs size of the clusters.

old_is=[]    #List1 for the local neighborhood of the particles
new_is=[]    #List2 for the local neighborhood of the particles


'''
READING AND PLOTING THE FILE
'''

molecular_rad=1.0       #Minimum radius for a cluster measurement.
lpos2=[[],[],[],[]]     #List that will have the information of the particles per step
for step in range(len(lineas)):
    
    '''
    In this section we manage to set the list lpos from the data in the
    document.
    '''
    lpos2=[[],[],[],[]]
    linea=lineas[step].split(',')
    linea.remove('\n')
    N=int(linea[0])
    del linea[0]
    
    old_is=[]
    new_is=[]
    caf.set_old_new_is(N,old_is,new_is)
    
    cont=0
    cont2=0
    for element in linea:
        if element!='|':
            if element not in colors_list:
                lpos2[cont].append(float(element))
                cont2+=1
            if element in colors_list:
                cont2+=1
                lpos2[cont].append(element)
        if cont2==N:
            cont+=1
            cont2=0
    
    '''
    Now that we have lpos setted, we proceed to do the local measurement of the 
    particle neighborhood, and then we set the cluster list 'cl'.
    
    With 'cl' setted we then obtain the average values of the cluster so we
    can plot them in the next section.
    '''
    caf.clustering_check1(N,lpos2,old_is,molecular_rad)
    caf.clustering_check2(N,lpos2,new_is,molecular_rad)
    cl=[]
    caf.clusters(N,lpos2,molecular_rad,cl,old_is,new_is)
    cl_means=copy.copy(caf.clusters_means(cl,lpos2))
    
    
    '''
    All the rest of the code is for ploting purposes.
    
    Fist we make the sactter plots, with or without the clusters measurment
    (if the cluster measrument isn't done, the len(cl_means)==0) and then we
    save the data of the sixes of the clusters and the amount fo the clusters 
    in the step.
    '''
    if len(cl_means)!=0:
        Final_plot=[[],[],[],[]]
        Final_plot[0]=lpos2[0]+cl_means[0]
        Final_plot[1]=lpos2[1]+cl_means[1]
        Final_plot[2]=lpos2[2]+cl_means[2]
        Final_plot[3]=lpos2[3]+cl_means[3]
        scat = ax2.scatter(Final_plot[0],Final_plot[1],c=Final_plot[2],s=Final_plot[3],alpha=0.30)

    if len(cl_means)==0:
        Final_plot=[[],[],[],[]]
        Final_plot[0]=lpos2[0]
        Final_plot[1]=lpos2[1]
        Final_plot[2]=lpos2[2]
        Final_plot[3]=lpos2[3]
        
        scat = ax2.scatter(Final_plot[0],Final_plot[1],c=Final_plot[2],s=Final_plot[3],alpha=0.30)    
    
    ims1.append([scat]) #We save thye scatterplot for the animation, like a frame
    caf.save_clusters_plot(clusters_per_step,cl,step) #We save the amount of clusters in this step
    caf.save_freq_plot1(frequency_list,cl) #We save the sizes of the clusters in this step.
        

    
ani1 = animation.ArtistAnimation(fig2, ims1, interval=100, blit=True,repeat_delay=0)#interval ~ tiempo entre frames
   


'''
The next few lines are to make the plots of the numbers of clusters per step
(in blue) and an histogram to obtain a plot that show the frequncy of each 
size of cluster in eh whole simulation (in green).
'''           
frequency_plot=caf.save_freq_plot2(frequency_list,cl)
    
plt.figure()
plt.subplot(2,1,1)
plt.title('Frequency Plot')
plt.hist(frequency_plot,bins=50,facecolor='green')#'bs'
plt.grid(True)

plt.subplot(2,1,2)
plt.title('Clusters per Step ')
plt.axis([0,len(lineas),0,int(max(clusters_per_step))+1])
plt.plot(clusters_per_step,c='b')
plt.grid(True)

    
            
     
    
        
