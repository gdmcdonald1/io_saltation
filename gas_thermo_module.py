#gas_thermo_module.py

#written by George McDonald
#written 11/11/20

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

def clausius_clapeyron(T_val,L_s,M,P_ref,T_ref):
    R=8.314462         #Ideal gas constant[J/K/mol]
    c1=P_ref/np.exp(-L_s/(R*T_ref))
    
    P_val=c1*np.exp(-L_s/(R*T_val))
    rho_val=(P_val*M)/(R*T_val) #Gas density assuming ideal gas
    return P_val, rho_val

def sutherland_formula_rankine(T,C,T_o,mu_o):
    '''
    Calculates the dynamic viscosity for a gas as a function of
    temperature. Input constants are all specific to a gas, and
    can be obtained at: http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/visgas.html
    Note, formula is a fit to empirical data, so the correct units MUST
    be used.    
    Inputs:
        T: Desired temperature of output viscosity [Rankine]
        C: Sutherland's constant
        T_o: Reference temp [Rankine]
        mu_o: Reference viscosity at T_o [centiPoise]
    Outputs:
        mu: Viscosity at specified desired temperature T [centiPoise]   
    '''
    a=0.555*T_o+C
    b=0.55*T+C
    
    mu_cP=mu_o*(a/b)*(T/T_o)**(3/2.)
    return mu_cP