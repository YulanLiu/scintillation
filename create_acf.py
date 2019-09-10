import numpy as np
import matplotlib.pyplot as plt
#import create_dysp
import scipy.signal
from scipy.optimize import curve_fit
import acffit_2d
def main(intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq,args,archive):
    intensity_masked=np.ma.masked_equal(intensity, 0)
    action = intensity_masked-np.ma.mean(intensity_masked)
    rfi_frac = np.count_nonzero(action) / float(np.size(action))
    action = np.ma.filled(action, 0)    
    
    # Correcting the sampling
    sampling_test = np.ones(np.shape(action))
    np.putmask(sampling_test, action==0, 0)   
    sampling = scipy.signal.fftconvolve(sampling_test, np.flipud(np.fliplr(sampling_test)), mode = 'full')   
    sampling[sampling == 0] = 1 #Avoids dividing by zero
    
    #Setting up array
    acf = scipy.signal.fftconvolve(action, np.flipud(np.fliplr(action)), mode = 'full')
    acf = np.divide(acf, sampling)
    #plt.imshow(acf_norm,aspect='auto')
    #plt.colorbar(use_gridspec=True)
    #plt.show()
    middle_f = acf.shape[0]/2
    middle_t = acf.shape[1]/2
    
    mhzperbin = diff / (intensity.shape[0])
    minperbin = length /60/ (intensity.shape[1])

    
    # remove spike
    acf[:,middle_t] = (acf[:,middle_t-1] + acf[:,middle_t+1]) / 2
   
    #normalizes
    midACF_freq = acf[:,middle_t] / acf[middle_f,middle_t]
    midACF_time = acf[middle_f,:] / acf[middle_f,middle_t]
    acf_norm = acf/np.amax(acf[int(middle_f):int(middle_f+middle_f*0.5),int(middle_t-middle_t*0.5):int(middle_t+middle_t*0.5)])
    
    try:
        filtered = scipy.signal.wiener(acf_norm[int(middle_f):int(middle_f+50),int(middle_t-50):int(middle_t+50)])
        residual = filtered - acf_norm[int(middle_f):int(middle_f+50),int(middle_t-50):int(middle_t+50)]
    except:
        filtered = scipy.signal.wiener(acf_norm)
        residual = filtered - acf_norm
    acfnoise = np.std(residual[int(0.25*residual.shape[0]):int(0.75*residual.shape[0]),int(0.25*residual.shape[1]):int(0.75*residual.shape[1])])
    
    return acf, middle_f, middle_t, sampling, rfi_frac, mhzperbin, minperbin, midACF_freq, midACF_time, acf_norm, acfnoise


def findend_lars(middle_f, middle_t, midACF_freq, midACF_time, mhzperbin, minperbin, sampling, acfnoise, args):
    i=5
    fend_err=[]
    while i < middle_f/2:
        if midACF_freq[middle_f+i] > 0:
            opt_f,cov_f = fit_f(i, middle_f, middle_t, midACF_freq, mhzperbin, sampling, acfnoise)
            fend_err.append([i,cov_f])
        else:
            break
        i=i+1
    fend_err=np.array(fend_err)
    try:
        best_fend=int(fend_err[np.nanargmin(fend_err[:,1]),0])
    except:
        best_fend=5
    
    j=5
    tend_err=[]
    while j< middle_t/2:
        if midACF_time[middle_t+j] > 0:
            opt_t,cov_t = fit_t(j, middle_f, middle_t, midACF_time, minperbin, sampling, acfnoise)
            tend_err.append([j,cov_t])
        else:
            break
        j=j+1
    tend_err=np.array(tend_err)
    try:
        best_tend=int(tend_err[np.nanargmin(tend_err[:,1]),0])
    except:
        best_tend=3
    
    return best_fend, best_tend 
    



    
def func_acf(x,b):
    return np.exp(-b*(x)**2)

def fit_f(f_end, middle_f, middle_t, midACF_freq, mhzperbin, sampling, acfnoise): 
    fitx_f = np.linspace(-f_end*mhzperbin,(f_end)*mhzperbin,f_end*2+1)
    #sampling_f=sampling[(middle_f-f_end):int(middle_f+f_end+1),middle_t]
    #acferr_f = np.divide(np.max(sampling_f),sampling_f)
    #acferr_f = (acferr_f**(0.5)) * acfnoise 
    vals_f = np.copy(midACF_freq)[int(middle_f-f_end):int(middle_f+f_end+1)]
    try:
        opt_f, cov_f = curve_fit(func_acf,fitx_f,vals_f,p0=0.0004,sigma=acferr_f)
    except:
        try:
            opt_f, cov_f = curve_fit(func_acf,fitx_f,vals_f,p0=0.0004)
        except:
            opt_f, cov_f = np.nan, np.nan
    
    return opt_f,cov_f

def fit_t(t_end, middle_f, middle_t, midACF_time, minperbin, sampling, acfnoise):
    fitx_t = np.linspace(-t_end*minperbin,t_end*minperbin,t_end*2+1)
    #sampling_t=sampling[middle_f,int(middle_t-t_end):int(middle_t+t_end+1)]
    #acferr_t = np.divide(np.max(sampling_t),sampling_t)
    #acferr_t = (acferr_t**0.5) * acfnoise
    vals_t = np.copy(midACF_time)[int(middle_t-t_end):int(middle_t+t_end+1)]
    try:
        opt_t, cov_t = curve_fit(func_acf,fitx_t,vals_t,p0=0.02,sigma=acferr_t)
    except:
        try:
            opt_t, cov_t = curve_fit(func_acf,fitx_t,vals_t,p0=0.02)
        except:
            opt_t, cov_t = np.nan, np.nan
    return opt_t, cov_t


def findend_slo(middle_f, middle_t,midACF_freq, midACF_time, args,intensity):
    def func_endf(x,b,c):
        return b*x+c
    try:
        endf=[]
        span_f = int (middle_f/20)
        for i in np.arange(middle_f-span_f-1):
            xcoor_f=np.linspace(i,span_f+i-1,span_f)
            ycoor_f=np.copy(midACF_freq)[int(middle_f+i):int(middle_f+span_f+i)]
            opt_f, cov_f = curve_fit(func_endf,xcoor_f,ycoor_f)
            endf.append(opt_f[0])
                #if opt_f[0]>0.15*endf[0] or midACF_freq[int(middle_f+i)] <= 0:
            if i > 6 and (opt_f[0]>0.2*endf[0] or midACF_freq[int(middle_f+i)] <=0): 
                f_end=i
                break
            else:
                f_end=int(middle_f/10)
    except:
        f_end=int(middle_f/15)

    if middle_f <=12:
        f_end = middle_f
    #print(f_end)
   

    try:
        span_t = int(middle_t/10)
        endt=[]
        for i in np.arange(middle_t-span_t-1):
            xcoor_t=np.linspace(i,span_t+i-1,span_t)
            ycoor_t=np.copy(midACF_time)[int(middle_t+i):int(middle_t+span_t+i)]
            opt_t, cov_t = curve_fit(func_endf,xcoor_t,ycoor_t)
            endt.append(opt_t[0])
            if i>3 and (opt_t[0]>0.2*endt[0] or midACF_time[int(middle_t+i)] <=0):
                t_end=i
                break
            else:
                t_end=int(middle_t/2)
    except:
        t_end=int(middle_t/3)
    if middle_t <=6:
        t_end = middle_t
    return f_end, t_end
def acf1Dfit(middle_f,middle_t,sampling,acfnoise, mhzperbin,minperbin,midACF_freq, midACF_time,f_end, t_end):
    opt_f,cov_f = fit_f(f_end, middle_f, middle_t, midACF_freq, mhzperbin, sampling, acfnoise)
    extrapolx_f = np.linspace(-f_end*mhzperbin,(f_end)*mhzperbin,f_end*2+1)
    extrapoly_f = func_acf(extrapolx_f,*opt_f)
    fitxplot_f = np.array(range(len(midACF_freq)))
    fitxplot_f = (fitxplot_f-middle_f)*mhzperbin
    
    opt_t,cov_t = fit_t(t_end, middle_f, middle_t, midACF_time, minperbin, sampling, acfnoise)
    extrapolx_t = np.linspace(-t_end*mhzperbin,(t_end)*mhzperbin,t_end*2+1)
    extrapoly_t = func_acf(extrapolx_t,*opt_t)
    fitxplot_t = np.array(range(len(midACF_time)))
    fitxplot_t = (fitxplot_t-middle_t)*minperbin
    
    cov_f=cov_f[0][0]
    cov_t=cov_t[0][0]
    opt_f=opt_f[0]
    opt_t=opt_t[0]
    return opt_f, cov_f, opt_t, cov_t, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t

def err1Dfit(diff,length,rfi_frac,freq_1D, time_1D, opt_f, cov_f, opt_t, cov_t):
    fill_factor=0.4
    stat_error = (fill_factor * ( (diff * length/60 * rfi_frac) / (freq_1D * time_1D)))**(-0.5)
   
 
    freqfit_error = 0.5*((np.log(2))**0.5)*(opt_f**(-1.5))*np.sqrt(cov_f)
    freq_error = np.sqrt(freqfit_error ** 2 + ( stat_error) ** 2)
   
    timefit_error = 0.5 * (opt_t**(-1.5)) * np.sqrt(cov_t)
    time_error = np.sqrt(timefit_error ** 2 + (stat_error) ** 2)
    
    return freq_error,time_error
def acf2Dfit(acf_mid,acf_norm,middle_f,middle_t,intensity,opt_f,opt_t,mhzperbin,minperbin,f_end,t_end):
    
    acf_center=np.copy(acf_norm[middle_f-f_end:middle_f+f_end,middle_t-t_end:middle_t+t_end])
    
    # create coordinate
    coor=np.indices(acf_center.shape)
    coor=np.asarray([(coor[0]-f_end)*mhzperbin,(coor[1]-t_end)*minperbin])
    params = acffit_2d.fitgaussian(acf_center,coor,opt_f,opt_t)

    coor_plot=np.indices(acf_mid.shape)
    coor_plot=np.asarray([(coor_plot[0]-intensity.shape[0]/2)*mhzperbin,(coor_plot[1]-intensity.shape[1]/2)*minperbin])
    fit = acffit_2d.gaussian(*params)
    acffit_2D=fit(*coor_plot)
    #print('freq_2D = %s, time_2D = %s' % (freq2D,time2D))
    #acffit.plotacf(params,acf_mid,coor)
    #plt.imshow(acf_mid,aspect='auto')
    #plt.contour(acffit_2D,aspect='auto')
    #plt.show()
    return params, acffit_2D
