#saltation_module.py

#written by George McDonald
#written 6/11/20

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
from scipy.optimize import fsolve

#Velocities of splashing grains. Kok et al. 2009a Equations 17 and 24.
def calc_v_ej_avg(v_imp,g,d):
    #Average ejection velocity due to splashing.     
    a_ej_avg=0.15   #Kok 2009a Table 1.
    a=0.02          #Kok 2009a Table 1.
    v_ej_avg=(a_ej_avg*np.sqrt(g*d)/a)*(1-np.exp(v_imp/(-40*np.sqrt(g*d))))
    return v_ej_avg

def calc_P_v_ej(v_ej,v_imp,g,d):
    #Probability distribution function value of ejection velocities at
    #a particular ejection velocity.
    P_v_ej=np.exp(-1*v_ej/calc_v_ej_avg(v_imp,g,d))/calc_v_ej_avg(v_imp,g,d)
    return P_v_ej

def calc_v_ej_less(v_imp,g,d,dist_frxn):
    #Velocity below which dist_frxn of  particles are ejected
    v_ej_less=-1*calc_v_ej_avg(v_imp,g,d)*np.log(1-dist_frxn)
    return v_ej_less


#Saltation Initiation Functions. This is from Shao and Lu 2000.
#u_ft terminology is from Telfer et al. 2018.

def calc_s(rho_p,rho_a):
    return rho_p/rho_a

def calc_u_ft(rho_p,rho_a,d,g):    
    an=0.0123
    zeta=5e-4      #These are from Shao & Lu
    print_str='WARNING: To replicate Pluto version we have to move replace\n'+\
                'Shao and Lu Constant with different constant OUTSIDE square\n'+\
                'square root, NOTE: order of operations'
    print(print_str)
    
    ustar=np.sqrt(an*(d*g*calc_s(rho_p,rho_a) + zeta/(rho_a*d)))    
    return ustar


#Saltation Continuation Functions.
#These are written as per Paehtz et al. 2012 equations 75a-75d.
#u_t terminology is from Telfer et al. 2018.

def calc_u_t(mu,rho_p,rho_a,d,g,nu):
    kappa=0.4   #von Karman constant
    eta=0.1     #Efficiency of wind accel. of grains (Paehtz et al. 2012)
    
    ut=kappa*(calc_Vr(mu,rho_p,rho_a,d,g)+calc_Vo(rho_p,rho_a,d,g))/\
            ((1-eta)*np.log(calc_z_mt(mu,rho_p,rho_a,d,g)/\
              calc_z_o_grain_low_Re(mu,rho_p,rho_a,d,g)))
    return ut

#Smaller inputs into u_t
def calc_gt(rho_p,rho_a,g):
    #This is g_tilda.
    gt=g*(calc_s(rho_p,rho_a)-1)/calc_s(rho_p,rho_a)
    return gt

def calc_gf(rho_p,rho_a,d,g):
    zeta=5e-4
    gf=calc_gt(rho_p,rho_a,g)+6.*zeta/(np.pi*rho_p*d**2.)
    return gf

def calc_Vo(rho_p,rho_a,d,g):
    Vo=16.2*np.sqrt(calc_gf(rho_p,rho_a,d,g)*d)
    return Vo

#Larger terms
def calc_Vt(mu,rho_p,rho_a,d,g):
    eta=0.1
    Vt=(calc_Vo(rho_p,rho_a,d,g)+eta*calc_Vr(mu,rho_p,rho_a,d,g))/(1-eta)
    return Vt
    
def calc_z_mt(mu,rho_p,rho_a,d,g):
    alpha=1.02
    beta_t=0.095
    gamma=0.17
    
    z_mt=(alpha*beta_t*gamma*calc_Vr(mu,rho_p,rho_a,d,g)**(1/2.)*\
          calc_Vt(mu,rho_p,rho_a,d,g)**(3/2.))/calc_gt(rho_p,rho_a,g)
    return z_mt

def calc_z_o_grain_low_Re(mu,rho_p,rho_a,d,g):
    #This is the relation for particle Reynolds # < 3, equation E.3 in
    #Paehtz et al. 2012.
    z_o_grain=mu/(9.*rho_a*calc_u_ft(rho_p,rho_a,d,g))
    return z_o_grain

#For calculating the nonlinear Vr term
def Vr_func(Vr, *constants):
    mu,rho_p,rho_a,d,g=constants
    alpha=1.02
    Vr_func=-1.*(4*calc_s(rho_p,rho_a)*calc_gt(rho_p,rho_a,g)*d/(3*alpha))+\
            (Vr**2.)*(1+(32.*mu/(Vr*rho_a*d))**(2/3.))**(3/2.)
    return Vr_func

def calc_Vr(mu,rho_p,rho_a,d,g):
    Vr=np.ones(d.shape[0])*np.nan
    for di, d_val in enumerate(d):
        constants=(mu,rho_p,rho_a,d_val,g)
        Vr[di]=fsolve(Vr_func,10,args=constants)
    return Vr