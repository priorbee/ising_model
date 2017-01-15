#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 14:30:14 2017

@author: aislingprior
"""

#!/usr/bin/env python
import numpy as np 
import matplotlib.pyplot as plt
import random as rd

#creating two lattices with 'nxn' lattice sites that have spin -1 or 1



#first = lattice w/ random spins
def rdlat(x):
        rdlat1=[]
        pos=1 #positive spins
        neg=-1 #negative

#for loop to create random spin arrays
        for i in range(x):
            rows=[]
#random number between 0 and 1, if greater than 0.5, positive spin 
            for j in range(x): 
                    if rd.random()<0.5: 
                        rows.append(neg)
                    else:
                        rows.append(pos)

            rdlat1.append(rows)

        return np.array(rdlat1)
        
    
print rdlat(5)

#second = lattice with collinear spins
def collat(x):
        collat1=[]
        spinup=np.ones(x) #array of positive spins (+1)
        spindown=-(spinup) #array of negative spins (-1)

#creating a loop that will produce x arrays that have same spin
        i=0
        while i<x:
            collat1.append(spinup)
            collat1.append(spindown)
            i=i+2
        return np.array(collat1)

print collat(5) 



# defining function that sweeps through lattice sites 
def sweeper(temp,site): 
        s=np.size(site[0]) #size of first row of lattice
        q=np.size(site)

#looping over lattice so it sweeps through each lattice point
        for i in range(q): 
        
#defining sites in the lattice w/ coordinates (a,b)
            a=rd.randint(0,s-1) #random integer between 0 and max value 
            b=rd.randint(0,s-1)
            point=site[a][b] #random lattice point
    
#finding energy difference between each lattice point and its adcacent sites    
#function for energy at point = -2*sum of interactions with neighbours
            if 0<a<s-1 and 0<b<s-1:
                    dE=-2*point*(site[a+1][b]+ site[a-1][b] +site[a][b+1]+site[a][b-1]) 
        #DeltaEnergy= -2*lattice point*(point to the right + point left + pt. up + pt. down)

#creating boundary conditions    
            elif 0<a<s-1 and b==0:
                     dE=-2*point*(site[a+1][b]+ site[a-1][b] +site[a][b+1]+ site[a][s-1]) 
                                       #s instead of b otherwise not within boundary^
            elif a==0 and 0<b<s-1:
                   dE=-2*point*(site[a+1][b]+ site[s-1][b] +site[a][b+1]+ site[a][b-1]) 
    
            elif a==0 and b==0:
                    dE=-2*point*(site[a+1][b]+ site[s-1][b] +site[a][b+1]+ site[a][b-1]) 
    
            elif a==s-1 and 0<b<s-1:
                    dE=-2*point*(site[0][b]+site[a-1][b] +site[a][b+1]+site[a][b-1]) 
     #if we increase a we are outside ^  boundary

            elif 0<a<s-1 and b==s-1:
                    dE=-2*point*(site[a+1][b]+site[a-1][b] +site[a][0]+ site[a][b-1]) 
    
            elif a==s-1 and b==s-1:
                    dE=-2*point*(site[0][b]+ site[a-1][b] +site[a][0]+ site[a][b-1])

            elif a==s-1 and b==0:
                    dE=-2*point*(site[0][b]+ site[a-1][b] +site[a][b+1]+site[a][s-1])

            elif  a==0 and b==s-1:
                    dE=-2*point*(site[a+1][b]+ site[s-1][b] +site[a][0]+ site[a][b-1]) 


# dE must be positive for a probability of the spin to flip

#calculating probability of the spin flipping
#following are correct equations but results are too small for computation due to kB being tiny
        #kB=1.38*10^23 # =boltzman's constant    
        #beta=1/(kB*Temp)
        #Pflip= np.exp-(beta*dE)

#revised method
#find probabilty in terms of boltzman's constant 
        beta1= 1/(temp)
#creating function that calculates probability of flipping, must be GT random no. to flip        
        Pflip1 = np.exp(-beta1*dE) #in terms of kB
        r=rd.random() #random number generator between 0 and 1
        
#if change in energy is negative then spin will flip
        if dE<0:
                    site[a][b]=-site[a][b]
#otherwise insert to probability of flipping equation                    
        elif Pflip1>r:
                    site[a][b]=-site[a][b]

#return site

#couldn't get method above to create working graph for some reason

#alternate method


L=10 # size of lattice (10x10)

#new sweeping function
def multisweeper(lattice, T):  #lattice size and new temperature
  
    for z in range(40000):     #for z in range(no. of sweeps), to reach equilibrium
        x = rd.randint(0, L-1) #same as previous method 
        y = rd.randint(0, L-1)       
        lattice1 = lattice[x,y] #creating lattice coordinates

#new  method to calculate boundary conditions, calculates all conditions so no sites outside lattice are chosen 
        Boundary=lattice[(x+1)%L,y]+lattice[x,(y+1)%L]+lattice[(x-1)%L,y] + lattice[x, (y-1)%L]
#new energy difference equation 
        dE2=2*lattice1*Boundary
        Pflip2 = np.exp(-dE2/T) 
# statements that flip spin when needed, similar to previous method       
        if dE2 < 0:         
            lattice1 = -lattice[x,y]       
        elif rd.random() < Pflip2: 
            lattice1 = -lattice[x,y] 
        lattice[x,y]=lattice1 
    return lattice
    

T=0.1 #starting temperature
# creating loop that calculates the magnetisation while temp is between 0.1 and 10 
while (T < 10.0): 
    
    #changing the spin values of the lattice to a list so they can be summed up
    spins=np.array(multisweeper(collat(L), T)) 
    #summing up spins to find total magnetisation of the lattice for each temperature
    magnetisation=np.sum((spins)/100.0)
    T = T + 0.1 # temp increasing in steps of 0.1
    magnetisation1=abs(magnetisation) #just need absolute value for magnetisation, the sign is arbitrary
    plt.plot(T, magnetisation1,  'o') #plotting magnetisation versus temp
    plt.xlabel('Temperature')
    plt.ylabel('Magnetisation')
    plt.title('Magnetisation vs Temperature')
   
plt.show()
