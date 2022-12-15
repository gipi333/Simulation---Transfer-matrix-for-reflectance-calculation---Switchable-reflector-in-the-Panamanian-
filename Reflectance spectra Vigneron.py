import cmath
import math

import numpy
from numpy.linalg import inv
import matplotlib.pyplot as plt


#-----------------------------------------
# Création des constantes
#-----------------------------------------

theta = 0
epsilon_sub = 2.0
epsilon_inc = 1
c = 3*(10**8)

nbr_couches = 64



#-----------------------------------------
# Création de la couche 
#-----------------------------------------

epsilon_j_non_ima = numpy.zeros(nbr_couches,dtype=complex)
epsilon_j_ima = numpy.zeros(nbr_couches,dtype=complex)
d_j = numpy.zeros(nbr_couches)


for k in range(0,nbr_couches):
    if (k%2 == 0):
        epsilon_j_non_ima[k] = (1.35)**2
        epsilon_j_ima[k] = complex((1.35)**2, 0.05)
    else: 
        epsilon_j_non_ima[k] = (1.55)**2
        epsilon_j_ima[k] = complex((1.55)**2, 0.05)
    
d_j[0:23:2] = 0.25*186*10**(-9)  
d_j[1:24:2] = 0.75*186*10**(-9)
d_j[24:49:2] = 0.25*220*10**(-9)
d_j[25:50:2] = 0.75*220*10**(-9)
d_j[48:64:2] = 0.25*270*10**(-9)
d_j[49:64:2] = 0.75*270*10**(-9)       

            



Transmitance_TE = numpy.zeros(601)
Transmitance_TM = numpy.zeros(601)
Reflectance_TM = numpy.zeros(601)
Reflectance_TE = numpy.zeros(601)


for lambdo in range(300,900):    
    
    w = 2* math.pi * (c/lambdo) * 10**(9)
    
    
    epsilon_j = numpy.zeros(nbr_couches,dtype=complex)
        
    if (lambdo < 500):   
        epsilon_j =  epsilon_j_ima    
    else:
        epsilon_j =  epsilon_j_non_ima 
    
      
    #-----------------------------------------
    # Définition de Q_inc pour le cas TE et TM
    #-----------------------------------------
    
    k_z_inc = numpy.sqrt( (w/c)**2 * epsilon_inc - (w/c)**2 * epsilon_inc**2 * math.sin( theta * (180/math.pi))**2     );
    
    
    Qinc_TE =  numpy.zeros((2,2),dtype=complex)
    
    Qinc_TE[0,0] = Qinc_TE[0,1] = 1
    Qinc_TE[1,0] = k_z_inc
    Qinc_TE[1,1] = -k_z_inc
    
    Qinc_TE_inv = inv(Qinc_TE)
    
    
    Qinc_TM = numpy.zeros((2,2),dtype=complex)
    
    Qinc_TM[0,0] = Qinc_TM[0,1] = 1
    Qinc_TM[1,0] = k_z_inc/epsilon_inc
    Qinc_TM[1,1] = -k_z_inc/epsilon_inc
    
    Qinc_TM_inv = inv(Qinc_TM)
    
    
    #-----------------------------------------
    # Définition de R_sub pour le cas TE et TM
    #-----------------------------------------
    
    k_z_sub = math.sqrt( (w/c)**2 * epsilon_sub - (w/c)**2 * epsilon_inc**2 * math.sin( theta * (180/math.pi) )**2)
    
    
    R_sub_TE = numpy.linspace(1, 4, 4).reshape((2,2))
    R_sub_TE[0,0] = R_sub_TE[0,1] = 1
    R_sub_TE[1,0] = k_z_sub
    R_sub_TE[1,1] = -k_z_sub
    
    
    R_sub_TM= numpy.linspace(1, 4, 4).reshape((2,2))
    
    R_sub_TM[0,0] = R_sub_TM[0,1] = 1
    R_sub_TM[1,0] = k_z_sub/epsilon_sub
    R_sub_TM[1,1] = -k_z_sub/epsilon_sub
    
      
    #-----------------------------------------
    # Définition de la matrice de transfert
    #-----------------------------------------
    
    
    #-----------------------------------------
    # Cas TE
    #-----------------------------------------
    
    T_init_TE = numpy.identity(2,dtype=complex)
    k_z_j = numpy.zeros(nbr_couches,dtype=complex)
    compl=complex(0,1)
    
    
    for k in range(0,nbr_couches):
        
        k_z_j[k] = cmath.sqrt( (w/c)**2 * epsilon_j[k] - (w/c)**2 * epsilon_inc**2 * math.sin(theta * (180/math.pi))**2 )
        
        Q_TE_j =  numpy.zeros((2,2),dtype=complex)
        Q_TE_j[0,0] = cmath.exp( compl * k_z_j[k] * d_j[k] )
        Q_TE_j[0,1] = cmath.exp( -compl * k_z_j[k] * d_j[k] )
        Q_TE_j[1,0] = k_z_j[k] * cmath.exp( compl * k_z_j[k] * d_j[k] )
        Q_TE_j[1,1] = -k_z_j[k] * cmath.exp( -compl * k_z_j[k] * d_j[k] )
        
        
        R_TE_j = numpy.zeros((2,2),dtype=complex)
        R_TE_j[0,0] =  R_TE_j[0,1] = 1
        R_TE_j[1,0] =  k_z_j[k]
        R_TE_j[1,1] =  - k_z_j[k]
        
        T_j_TE = numpy.dot(R_TE_j,inv(Q_TE_j))
        
        T_init_TE = numpy.dot(T_init_TE,T_j_TE)
    
    
    
    T_final_TE = Qinc_TE_inv.dot(T_init_TE).dot(R_sub_TE);
    
    Transmitance_TE[lambdo-299] = (abs(1/T_final_TE[0,0]))**2 * ((k_z_sub)/ k_z_inc) 
    Reflectance_TE[lambdo-299]  = (abs(T_final_TE[1,0]/T_final_TE[0,0]))**2
    
    
    
    
    
    T_init_TM = numpy.identity(2,dtype=complex)
    
    
    for k in range(0,nbr_couches):
        
        k_z_j[k] = cmath.sqrt( (w/c)**2 * epsilon_j[k] - (w/c)**2 * epsilon_inc**2 * math.sin(theta * (180/cmath.pi))**2 )
        
        Q_TM_j =  numpy.zeros((2,2),dtype=complex)
        Q_TM_j[0,0] = cmath.exp( compl * k_z_j[k] * d_j[k] )
        Q_TM_j[0,1] = cmath.exp( -compl * k_z_j[k] * d_j[k] )
        Q_TM_j[1,0] = (k_z_j[k] * cmath.exp( compl * k_z_j[k] * d_j[k] ))/epsilon_j[k]
        Q_TM_j[1,1] = (-k_z_j[k] * cmath.exp( -compl * k_z_j[k] * d_j[k] ))/epsilon_j[k]
        
        
        R_TM_j = numpy.zeros((2,2),dtype=complex)
        R_TM_j[0,0] =  R_TM_j[0,1] = 1
        R_TM_j[1,0] =  k_z_j[k]/epsilon_j[k]
        R_TM_j[1,1] =  - k_z_j[k]/epsilon_j[k]
        
        T_j_TM = numpy.dot(R_TM_j,inv(Q_TM_j))
        
        T_init_TM = numpy.dot(T_init_TM,T_j_TM)
    
    
    
    T_final_TM = numpy.dot(Qinc_TM_inv, numpy.dot(T_init_TM, R_sub_TM))
    
    
    Transmitance_TM[lambdo-299] =  (epsilon_inc/epsilon_sub)   * (abs(1/T_final_TM[0,0]))**2 * ((k_z_sub)/ k_z_inc) 
    Reflectance_TM[lambdo-299]  = (abs(T_final_TM[1,0]/T_final_TM[0,0]))**2
    




lambd = numpy.zeros(601)
Reflectance = numpy.zeros(601)
Transmitance = numpy.zeros(601)

for i in range(300,900):
    lambd[i-299]=i
    Reflectance[i-299] = 100*(Reflectance_TM[i-299] + Reflectance_TE[i-299])/2
    Transmitance[i-299] = 100*(Transmitance_TM[i-299] + Transmitance_TE[i-299])/2

#plt.plot(lambd,Reflectance,lambd,Transmitance)
plt.plot(lambd,Reflectance)
plt.xlim(300, 900)
plt.title('Réflectance en fonction de la longueur d onde')
plt.xlabel('longueur d onde [nm] ')
plt.ylabel('Réflectance [%]')







