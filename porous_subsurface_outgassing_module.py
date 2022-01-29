#gas_thermo_module.py

#written by George McDonald
#written 11/11/20

#From Jia et al. 2017S, outgassing velocity for sublimation of a gas under
#a porous subsurface layer. Developed in context for comet 67P.
#Equations S14, S45-S51. Note, S52 is not encoded, because that seems
#to assume Comet 67P parameters, e.g. F is missing. #Rather equation S49,
#here energy_func, is solved directly.

#Copyright 2022 George McDonald

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <https://www.gnu.org/licenses/>.

#Module for reformatting benefit report to highlight changes, remove
#miscellaneous elements etc, thus preparing for sending to companies.
#Elements borrowed from: https://stackoverflow.com/questions/40014989/python-xlrd-and-xlwt

###############################################################################

import numpy as np
import math

def Vthi_func(T_i,M):
    kB=1.380649e-23     #Boltzmann Constant [J/K]
    N_A=6.02214076e23   #Avogadro constant [mol^-1]
    m=M/N_A             #Molar mass -> molecular mass
    V_th=np.sqrt((8*kB*T_i)/(np.pi*m))
    return V_th

def fy_func(y):
    term1=np.exp(-4*y**2./np.pi)
    term2=2*y*(1-math.erf(2*y/np.sqrt(np.pi)))
    fy=term1-term2
    return fy

def rhoo_func(F,pc,yo,rhos):
    num=np.pi*(fy_func(yo)*(pc-2.)+2.)+16*yo**2.
    denom=np.pi*(fy_func(yo)*pc+4*yo)**2.
    rhoo=F*pc*(num/denom)*rhos
    return rhoo

def Vtho_func(pc,yo,T_i,M):
    num=np.pi*(fy_func(yo)*pc+4.*yo)
    denom=np.pi*(fy_func(yo)*(pc-2)+2)+16*yo**2.
    Vtho=(num/denom)*Vthi_func(T_i,M)
    return Vtho

def q(F,pc,yo,rhos,T_i,M):
    q=(1/4.)*fy_func(yo)*rhoo_func(F,pc,yo,rhos)*Vtho_func(pc,yo,T_i,M)
    return q

def energy_func(yo,F,pc,rhos,T_i,M):
    term1=(1/2.)*yo*(yo**2.+np.pi)*rhoo_func(F,pc,yo,rhos)*\
                Vtho_func(pc,yo,T_i,M)**3.
    term2=(7*np.pi/16.)*pc*((1/4.)*F*rhos*Vthi_func(T_i,M)**3.-\
                            q(F,pc,yo,rhos,T_i,M)*Vtho_func(pc,yo,T_i,M)**2.)
    zero=term1/term2-1.     #The eqn rearranged to be equal to 0.
    return zero