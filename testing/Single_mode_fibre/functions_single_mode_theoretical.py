from scipy.special import jv,kv
from math import atan
from scipy.integrate import dblquad
from scipy.constants import c,pi
import matplotlib.pylab as plt
import numpy as np
def system(vec,V,Delta):
    u,w = vec
    return u**2+w**2 - V**2, jv(0,u)/(u*jv(1,u)) - (1-Delta)*kv(0,w)/(w*kv(1,w))


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

def electric_field(x,y,u,w,beta,a,V,Delta,psi,n):
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



def plot_electric(x,y,u,w,beta,a,V,Delta,psi,mode):
    ele = np.zeros([3,len(x),len(y)],dtype=np.complex128)
    for i,xx in enumerate(x):
        for j, yy in enumerate(y):
            ele[:,i,j] = electric_field(xx,yy,u,w,beta,a,V,Delta,psi,mode)
    abss = (np.abs(ele[0,:,:])**2 + np.abs(ele[1,:,:])**2 + np.abs(ele[2,:,:])**2)**0.5
    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    plt.contourf(X,Y,abss)
    plt.show
    return 0


def I(y,x,u,w,beta,a,V,Delta,psi,mode,po):
    ex,ey,ez = electric_field(x,y,u,w,beta,a,V,Delta,psi,mode)
    return np.abs(np.abs(ex)**2 +np.abs(ex)**2 + np.abs(ex)**2)**po


def effective_area(a,u,w,beta,V,Delta,mode):
    top = dblquad(I,-2*a,2*a,lambda x : -2*a,lambda x: 2*a,args = (u,w,beta,a,V,Delta,psi,mode,1))
    bottom = dblquad(I,-2*a,2*a,lambda x : -2*a,lambda x: 2*a,args = (u,w,beta,a,V,Delta,psi,mode,2))
    A_eff = top[0]**2/bottom[0]
    return A_eff


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


def sigma(u,w,V,Delta):
    return -1 - u**2 * w**2 * Delta * jv(0,u)/(V**2 * u * jv(1,u))


def system(vec,V,Delta):
    u,w = vec
    return u**2+w**2 - V**2, jv(0,u)/(u*jv(1,u)) - (1-Delta)*kv(0,w)/(w*kv(1,w))