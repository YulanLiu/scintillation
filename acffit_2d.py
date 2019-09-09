# Here is robust code to fit a 2D gaussian. It calculates the moments of the data to guess the initial parameters for an optimization routine. For a more complete gaussian, one with an optional additive constant and rotation, see http://code.google.com/p/agpy/source/browse/trunk/agpy/gaussfitter.py. It also allows the specification of a known error.


from numpy import *
from scipy import optimize
from scipy.optimize import curve_fit
def gaussian(a,b,c,height):
    
    return lambda f,t: height*exp(a*(f)**2+b*f*t+c*(t)**2)
    #return lambda f,t: ravel(height*exp(-a*(f)**2+b*f*t-c*(t)**2))



def fitgaussian(data,coor,opt_f,opt_t):
    acf_ini=[-opt_f,(opt_f+opt_t)/100,-opt_t,1.0]
    errorfunction = lambda p: ravel(gaussian(*p)(*coor) - data)
    p, success = optimize.leastsq(errorfunction, acf_ini)
    #p, success = curve_fit(gaussian, (coor_fit[0],coor_fit[1]), data, p0=acf_ini)
    
    return p

