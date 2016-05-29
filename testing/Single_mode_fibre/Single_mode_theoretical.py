
# coding: utf-8

# # Single mode fibre
# This uses the theoretical result of the single mode optical fibre in order to form a point of comparison to the Fenics solver. 

# In[5]:

from __future__ import division,print_function
import numpy as np
import matplotlib.pylab as plt
from scipy.constants import c,pi
from scipy.special import jv,kv
from scipy.optimize import fsolve
from math import atan
from scipy.integrate import dblquad



def system(vec,V,Delta):


    u,w = vec
    return u**2+w**2 - V**2, jv(0,u)/u/jv(1,u) - (1-Delta)*kv(0,w)/w/kv(1,w)

def sigma(u,w,V,Delta):
    return -1 - u**2 * w**2 * Delta * jv(0,u)/(V**2 * u * jv(1,u))


def cart_to_cyl(x,y,z):
    """
    Transforms cartesian coordinates to cylyndrical coordinates
    input x, y, z
    returns tuple:  r, theta, z 
    """
    try:
        return (x**2 + y**2)**0.5, atan(y/x),z
    except ZeroDivisionError:
        
        return (x**2 + y**2)**0.5, atan(np.inf),z
def electric_field(x,y,u,w,beta,a,V,psi,n,Delta):
    n+=1
    r,theta,z = cart_to_cyl(x,y,0)
    s = sigma(u,w,V,Delta)
    if r <=a:
        er = (-1j* beta* a/u)*np.cos(n*theta +psi)*(0.5*(1-s)*jv(n-1,u*r/a) - 0.5*(1+s)*jv(n+1,u*r/a))
        etheta = (1j*beta*a/u)*np.sin(n*theta+psi)*(0.5*(1-s)*jv(n-1,u*r/a)+0.5*(1+s)*jv(n+1,u*r/a))
        ez = jv(n,u*r/a)*np.cos(n*theta+psi)
    else:
        er = -1j*beta *a *jv(n,u)/(w*kv(n,w))*(0.5*(1-s)*kv(n-1,w*r/a)+0.5*(1+s)*kv(n+1,w*r/a))*np.cos(n*theta + psi)
        etheta = 1j*beta *a *jv(n,u)/(w*kv(n,w))*(0.5*(1-s)*kv(n-1,w*r/a)-0.5*(1+s)*kv(n+1,w*r/a))*np.sin(n*theta + psi)
        ez = (jv(n,u)/w*kv(n,w))*np.cos(n*theta+psi) 
    ex = er*np.cos(theta) - etheta*np.sin(theta)
    ey =er*np.sin(theta) + etheta*np.cos(theta)
    return ex,ey,ez    


def plot_electric(x,y,u,w,beta,a,V,Delta):
    ele = np.zeros([3,len(x),len(y)],dtype=np.complex128)
    for i,xx in enumerate(x):
        for j, yy in enumerate(y):
            ele[:,i,j] = electric_field(xx,yy,u,w,beta,a,V,0,0,Delta)
    abss = (np.abs(ele[0,:,:])**2 + np.abs(ele[1,:,:])**2 + np.abs(ele[2,:,:])**2)**0.5
    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    plt.contourf(X,Y,abss)
    plt.show()
    return 0

def I(y,x,u,w,beta,a,V,psi,n,po,Delta):
    ex,ey,ez = electric_field(x,y,u,w,beta,a,V,psi,n,Delta)
    return np.abs(np.abs(ex)**2 +np.abs(ex)**2 + np.abs(ex)**2)**po


def effective_area(a,u,w,beta,V,n,Delta):
    top = dblquad(I,-2*a,2*a,lambda x : -2*a,lambda x: 2*a,args = (u,w,beta,a,V,0,n,1,Delta))
    bottom = dblquad(I,-2*a,2*a,lambda x : -2*a,lambda x: 2*a,args = (u,w,beta,a,V,0,n,2,Delta))
    A_eff = top[0]**2/bottom[0]
    return A_eff    




def main_test(n1,n0,lamda,a,guess = [1,1],plot = False):
    V  = a*2*pi/lamda*(n1**2 - n0**2)**0.5
    Delta = (n1**2 - n0**2)/(2*n1**2)
    print(V)
    print('Doing calculation for: ',n1, n0,lamda,a)

    u,w = fsolve(system,guess,args=(V,Delta))
    neff = (((n1/u)**2 + (n0/w)**2)/((1/u)**2 + (1/w)**2))**0.5
    #beta = neff*2*pi/lamda
    #x = np.linspace(-2*a,2*a,512)
    #y = np.linspace(-2*a,2*a,512)
    #if plot ==True:
    #    plot_electric(x,y,u,w,beta,a,V,Delta)

    Aeff = None #effective_area(a,u,w,beta,V,0,Delta)
    return neff,Aeff




if __name__ == '__main__':
    A = np.loadtxt('../../parameters.txt')

    n1,n0,lamda,a = 1.445-1e-4j,1.444,1.55e-6,1e-5#A
    print('core refractive index: ', n1)
    print('clading refractive index: ', n0)
    print('wavelength: ', lamda)
    print('core radius: ', a)
    neff,Aeff = main_test(n1,n0,lamda,a)
    print('resulting in:')
    print('effective index: ', neff)

    print('effective area: ', Aeff)




