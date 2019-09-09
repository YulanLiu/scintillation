import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def main(args,archive,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time,acf_mid,opt_f,opt_t,secondary, acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility,sec_axes):
    mjd = archive.rsplit('_')[1]
    pulsar = 'PSR B'+archive.rsplit('_')[2]
    diff=maxFreq-minFreq
    print('%s %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.6f %.6f'  % (mjd,freq1D,freq_error,time1D,time_error,length/60,diff,minFreq,maxFreq,Drift_rate,slope_visibility))
    
    if args.plot == "all":
        plotall(args,archive,pulsar,mjd,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time,acf_mid,opt_f,opt_t,secondary, acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility,sec_axes)
    elif args.plot == "dysp":
        plotdysp(args,archive,pulsar,mjd,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time)
    elif args.plot == "acf":
        plotacf(args,archive,pulsar,mjd,minFreq, maxFreq,length,extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t,midACF_freq, midACF_time,acf_mid,opt_f,opt_t,acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility)
    elif args.plot == "sec":
        plotsec(args,archive,pulsar,mjd,secondary,sec_axes)
    


        
def plotall(args,archive,pulsar,mjd,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time,acf_mid,opt_f,opt_t,secondary, acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility,sec_axes):
    nullfmt = NullFormatter()         # no labels
    # definitions for the axes
    left_edge=0.06
    bottom_edge=0.13
    width = 0.24
    height = 0.75
    space =0.005
    
    rect_dysp = [left_edge, bottom_edge, width, height-0.21]
    rect_acf = [left_edge+width+space,bottom_edge, width, height-0.21]
    rect_time = [left_edge+width+space, height-0.07, width, 0.26]
    rect_freq = [left_edge+2*width+space, bottom_edge , 0.12, height-0.21]
    rect_sec =[left_edge+2*width+space+space+0.18, bottom_edge , 0.25, height-0.15]
    
    plt.figure(1, figsize=(14, 6))
    axdysp = plt.axes(rect_dysp)
    axacf = plt.axes(rect_acf)
    axtime = plt.axes(rect_time)
    axfreq = plt.axes(rect_freq)
    axsec = plt.axes(rect_sec)
    
    
    axdysp.imshow(intensity,aspect='auto',extent=[0,length/60,minFreq,maxFreq],vmax=1, vmin=-0.1,cmap='jet', interpolation='None')
    axdysp.set_xlabel("time(mins)")
    axdysp.set_ylabel("frequency(Mhz)")
    axdysp.set_title("%s \n MJD %s \nTime:%.4f +- %.4f mins \n freq: %.4f +- %.4f Mhz \n observation time: %.4f mins\n \n"  % (pulsar,mjd,time1D,time_error,freq1D, freq_error,length/60), fontsize=12)

    axacf.yaxis.set_major_formatter(nullfmt)
    axacf.imshow(acf_mid,aspect='auto',extent=[-(length/60)/2,(length/60)/2,-(maxFreq-minFreq)/2,(maxFreq-minFreq)/2])
    axacf.contour(acf_mid,aspect='auto',extent=[-(length/60)/2,(length/60)/2,(maxFreq-minFreq)/2,-(maxFreq-minFreq)/2])
    axacf.contour(acffit_2D,aspect='auto',extent=[-(length/60)/2,(length/60)/2,(maxFreq-minFreq)/2,-(maxFreq-minFreq)/2])
    axacf.set_xlabel("Drift_rate=%.6f \n slope_visibility=%.6f" % (Drift_rate,slope_visibility))
    
    t_y=np.linspace(np.min(midACF_time)-0.1,1.1*np.max(midACF_time),20)
    t_x=[]
    for i in range(20):
        t_x.append(np.sqrt(1/opt_t))
    
    axtime.xaxis.set_major_formatter(nullfmt)
    axtime.scatter(fitxplot_t, midACF_time)
    axtime.set_xlim(-(length/60)/2,(length/60)/2)
    timelen=len(midACF_time)
    axtime.set_ylim(np.min(midACF_time[int(0.25*timelen):int(0.75*timelen)])-0.1,np.max(midACF_time[int(0.25*timelen):int(0.75*timelen)])+0.1)
    axtime.plot(t_x,t_y,color='red',linewidth=4.0)
    axtime.plot(extrapolx_t,extrapoly_t,linewidth=4.0)
    axtime.locator_params(axis='y',tight=True, nbins=5)
    
    f_x=np.linspace(np.min(midACF_freq)-0.1,1.2*np.max(midACF_freq),20)
    f_y=[]
    for i in range(20):
        f_y.append(-np.sqrt(np.log(2)/opt_f))
    
    axfreq.yaxis.set_major_formatter(nullfmt)
    axfreq.scatter(midACF_freq,fitxplot_f)
    freqlen=len(midACF_freq)  
    axfreq.set_xlim(np.min(midACF_freq[int(0.25*freqlen):int(0.75*freqlen)])-0.1,1.2*np.max(midACF_freq[int(0.25*freqlen):int(0.75*freqlen)]))
    axfreq.set_ylim(-(maxFreq-minFreq)/2,(maxFreq-minFreq)/2)
    #axfreq.set_xlabel("freq: \n%.4f(%.4f) Mhz" % (freq1D, freq_error))
    #axfreq.set_title("observation time: %.4f \n time: %.4f +- %.4f Min \n \n \n" % (length/60, time1D, time_error),fontsize=12)
    axfreq.plot(f_x,f_y,color='red',linewidth=4.0)
    axfreq.plot(extrapoly_f,extrapolx_f,linewidth=4.0)
    axfreq.locator_params(axis='x',tight=True, nbins=5)

    sec_mean=np.mean(secondary)
    axsec.imshow(secondary,aspect='auto',origin='lower', cmap ='binary', extent=[sec_axes[2],sec_axes[3],sec_axes[0],sec_axes[1]],interpolation='None',vmax=sec_mean+8,vmin=sec_mean+2)
    axsec.set_xlabel(r'Fringe Frequency ($10^{-3}$Hz)')
    axsec.set_ylabel(r'Delay ($\mu$s)')
    axsec.set_title("Secondary spectra")

    
    filename=archive.rsplit('.')[0]
    
   
    if args.savefigure:
        plt.savefig('%s_all.eps' % filename)
        plt.clf()
    else:
        plt.show()

def plotdysp(args,archive,pulsar,mjd,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time):
    plt.imshow(intensity,aspect='auto',extent=[0,length/60,minFreq,maxFreq],vmax=1, vmin=-0.1,cmap='jet', interpolation='None')
    plt.xlabel("time(mins)")
    plt.ylabel("frequency(Mhz)")
    plt.title("%s  MJD %s \nTime:%.2f +- %.2f mins  freq: %.2f +- %.2f Mhz  observation time: %.2f mins"  % (pulsar,mjd,time1D,time_error,freq1D, freq_error,length/60), fontsize=12)
    filename=archive.rsplit('.')[0]
    
   
    if args.savefigure:
        plt.savefig('%s_all.eps' % filename)
        plt.clf()
    else:
        plt.show()

def plotacf(args,archive,pulsar,mjd,minFreq, maxFreq,length,extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t,midACF_freq, midACF_time,acf_mid,opt_f,opt_t,acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility):
    
    nullfmt = NullFormatter()         # no labels
    # definitions for the axes
    left_edge=0.15
    bottom_edge=0.13
    width = 0.5
    height = 0.75
    space =0.005
    
    
    rect_acf = [left_edge,bottom_edge, width, height-0.21]
    rect_time = [left_edge, height-0.07, width, 0.26]
    rect_freq = [left_edge+width+space+0.01, bottom_edge , 0.26, height-0.21]
    
    
    plt.figure(1, figsize=(7, 6))
    
    axacf = plt.axes(rect_acf)
    axtime = plt.axes(rect_time)
    axfreq = plt.axes(rect_freq)
    
    
    axacf.set_ylabel("frequency(Mhz)")
    axacf.imshow(acf_mid,aspect='auto',extent=[-(length/60)/2,(length/60)/2,-(maxFreq-minFreq)/2,(maxFreq-minFreq)/2])
    axacf.contour(acf_mid,aspect='auto',extent=[-(length/60)/2,(length/60)/2,(maxFreq-minFreq)/2,-(maxFreq-minFreq)/2])
    axacf.contour(acffit_2D,aspect='auto',extent=[-(length/60)/2,(length/60)/2,(maxFreq-minFreq)/2,-(maxFreq-minFreq)/2])
    axacf.set_xlabel("time(mins) \n Drift_rate=%.6f slope_visibility=%.6f" % (Drift_rate,slope_visibility))
    
    t_y=np.linspace(np.min(midACF_time)-0.1,1.1*np.max(midACF_time),20)
    t_x=[]
    for i in range(20):
        t_x.append(np.sqrt(1/opt_t))

    axtime.set_title("%s MJD %s "  % (pulsar,mjd), fontsize=12)
    axtime.xaxis.set_major_formatter(nullfmt)
    axtime.scatter(fitxplot_t, midACF_time)
    axtime.set_xlim(-(length/60)/2,(length/60)/2)
    timelen=len(midACF_time)
    axtime.set_ylim(np.min(midACF_time[int(0.25*timelen):int(0.75*timelen)])-0.1,np.max(midACF_time[int(0.25*timelen):int(0.75*timelen)])+0.1)
    axtime.plot(t_x,t_y,color='red',linewidth=4.0)
    axtime.plot(extrapolx_t,extrapoly_t,linewidth=4.0)
    axtime.locator_params(axis='y',tight=True, nbins=5)
    
    f_x=np.linspace(np.min(midACF_freq)-0.1,1.2*np.max(midACF_freq),20)
    f_y=[]
    for i in range(20):
        f_y.append(-np.sqrt(np.log(2)/opt_f))
    
    axfreq.yaxis.set_major_formatter(nullfmt)
    axfreq.scatter(midACF_freq,fitxplot_f)
    freqlen=len(midACF_freq)  
    axfreq.set_xlim(np.min(midACF_freq[int(0.25*freqlen):int(0.75*freqlen)])-0.1,1.2*np.max(midACF_freq[int(0.25*freqlen):int(0.75*freqlen)]))
    axfreq.set_ylim(-(maxFreq-minFreq)/2,(maxFreq-minFreq)/2)
    #axfreq.set_xlabel("freq: \n%.4f(%.4f) Mhz" % (freq1D, freq_error))
    #axfreq.set_title("observation time: %.4f \n time: %.4f +- %.4f Min \n \n \n" % (length/60, time1D, time_error),fontsize=12)
    axfreq.plot(f_x,f_y,color='red',linewidth=4.0)
    axfreq.plot(extrapoly_f,extrapolx_f,linewidth=4.0)
    axfreq.locator_params(axis='x',tight=True, nbins=5)
   
    filename=archive.rsplit('.')[0]
    
   
    if args.savefigure:
        plt.savefig('%s_all.eps' % filename)
        plt.clf()
    else:
        plt.show()

def plotsec(args,archive,pulsar,mjd,secondary,sec_axes):
    sec_mean=np.mean(secondary)
    plt.imshow(secondary,aspect='auto',origin='lower', cmap ='binary', extent=[sec_axes[2],sec_axes[3],sec_axes[0],sec_axes[1]],interpolation='None',vmax=sec_mean+8,vmin=sec_mean+2)
    plt.xlabel(r'Fringe Frequency ($10^{-3}$Hz)')
    plt.ylabel(r'Delay ($\mu$s)')
    plt.title("Secondary spectra \n %s  MJD %s " % (pulsar, mjd))
  
    filename=archive.rsplit('.')[0]
    
    if args.savefigure:
        plt.savefig('%s_all.eps' % filename)
        plt.clf()
    else:
        plt.show()
