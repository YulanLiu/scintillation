import argparse
import create_dysp
import matplotlib.pyplot as plt
import create_acf
import create_sec
import output
import warnings
import numpy as np


parser=argparse.ArgumentParser(description="Select the archive.")
parser.add_argument('archives', nargs='+', help='The chosen archives')

parser.add_argument('-fill', '--fillempty', action='store_true', help='fill the empty place where the RFI had been removed')

parser.add_argument('-end', '--endmethod', type= int, default= 0, help='There have two methods to find the best end, one of the them is from lars(0), the other is useing the slope rate of ACF(1)')

parser.add_argument('-plot', '--plot', type=str, help='plot the dynamic spectrum, ACF and secendary. Options: all, dysp, acf, sec')
parser.add_argument('-s', '--savefigure', action='store_true', help='save the dysp')
parser.add_argument('-q', '--quiet', action='store_true', help='Do not print text information.')

args = parser.parse_args()
archive_list = args.archives

warnings.simplefilter('ignore', RuntimeWarning)
warnings.simplefilter('ignore', FutureWarning)

print("#mjd,freq1D,freq_error,time1D,time_error,length/60,diff,minFreq,maxFreq,Drift_rate,slope_visibility")
for archive in archive_list:
    # get the dynamic spectrum
    intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = create_dysp.main(args,archive)

    # get the ACF and ACF's noise
    acf, middle_f, middle_t, sampling, rfi_frac, mhzperbin, minperbin, midACF_freq, midACF_time, acf_norm, acfnoise= create_acf.main(intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq,args,archive)
    acf_mid=np.copy(acf_norm[middle_f-intensity.shape[0]/2:middle_f+intensity.shape[0]/2,middle_t-intensity.shape[1]/2:middle_t+intensity.shape[1]/2])
    
    # find the best end of the ACF fit
    if args.endmethod == 0:
        f_end, t_end = create_acf.findend_lars(middle_f, middle_t, midACF_freq, midACF_time, mhzperbin, minperbin, sampling, acfnoise, args)
    else:
        f_end, t_end = create_acf.findend_slo(middle_f, middle_t,midACF_freq, midACF_time, args,intensity)
    # fit the ACF and get scintillation bandwidth and timescale. There are two ways, 1D fit and 2D fit.
    opt_f, cov_f, opt_t, cov_t, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t = create_acf.acf1Dfit(middle_f,middle_t,sampling,acfnoise,mhzperbin,minperbin,midACF_freq, midACF_time,f_end, t_end)
    freq1D=np.sqrt(np.log(2)/opt_f)
    time1D=np.sqrt(1./opt_t)
    freq_error,time_error = create_acf.err1Dfit(diff,length,rfi_frac,freq1D, time1D,opt_f, cov_f, opt_t, cov_t)
    params, acffit_2D = create_acf.acf2Dfit(acf_mid,acf_norm,middle_f,middle_t,intensity,opt_f,opt_t,mhzperbin,minperbin,f_end,t_end)
    C_1,C_2,C_3,C_0=params
    Drift_rate=-(C_2/(2*C_3))
    slope_visibility=C_2/np.sqrt(np.abs(4*C_1*C_3))
    #print('%.f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.6f %.6f'  % (mjd_start,freq1D,freq_error,time1D,time_error,diff,length/60,minFreq,maxFreq,Drift_rate,slope_visibility))
    secondary,sec_axes = create_sec.main(intensity, args, minFreq, maxFreq, diff, archive, name, length, mhzperbin, minperbin)
    output.main(args,archive,intensity,minFreq, maxFreq,length, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time,acf_mid,opt_f,opt_t,secondary, acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility,sec_axes)
    #if args.plot:
        #output.main(intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq,args,archive, extrapolx_f,extrapoly_f, fitxplot_f, extrapolx_t, extrapoly_t, fitxplot_t, freq1D, time1D, freq_error,time_error, midACF_freq, midACF_time,middle_f, middle_t,acf_mid,opt_f,opt_t,secondary_log, acffit_2D,f_end,t_end,mhzperbin,minperbin,Drift_rate,slope_visibility,sec_axes)


    
    #plt.imshow(intensity,aspect='auto',extent=[0,length/60,minFreq,maxFreq],vmax=1, vmin=-0.1,cmap='jet', interpolation='None')
    #plt.imshow(acf_mid,aspect='auto')
    #plt.xlabel("time(mins)")
    #plt.ylabel("frequency(Mhz)")
    #plt.colorbar(use_gridspec=True)
    #plt.savefig('%s_.png' % archive )
    #plt.show()
    
    

    
        

