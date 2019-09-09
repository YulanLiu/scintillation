import numpy as np
from scipy.signal import medfilt
from scipy.optimize import curve_fit
def main(intensity, args, minFreq, maxFreq, diff, archive, name, length, mhzperbin, minperbin):
    intensity_masked=np.ma.masked_equal(intensity, 0)
    action = intensity_masked-np.ma.mean(intensity_masked)
    action = np.ma.filled(action, np.ma.median(action))
   
    sec_initial = np.fft.fft2(action)
    sec_i_abs = (np.abs(sec_initial))**2 
    sec_shifted = np.fft.fftshift(sec_i_abs) 
    
    shape_f = np.shape(sec_shifted)[0]
    shape_t = np.shape(sec_shifted)[1]
    t_end = int(1/7.0*shape_t)
    outer_medians = medfilt(np.median(sec_shifted[:,0:t_end], axis = 1), kernel_size=41)
    x_vals = np.linspace(-10, 10, shape_f)
    try:
        opt, cov = curve_fit(linear_sym, x_vals, outer_medians)
        med_correction = - linear_sym(x_vals, *opt) + opt[1] 
        for j in range(np.shape(sec_shifted)[0]):
            sec_shifted[j,:] = sec_shifted[j,:] + med_correction[j]
    except:
        print("not substract noise in secondary")
    sec_log = 10. * np.log10(sec_shifted / np.max(sec_shifted))
    sec_cut = sec_log[sec_log.shape[0]/2:sec_log.shape[0],:]
    sec_axes=secondary_axes(sec_cut, mhzperbin, minperbin)
    return sec_cut,sec_axes
def linear_sym(x,a,b):
    return a * np.abs(x) + b

def secondary_axes(secondary, mhzperbin, minperbin):
    conj_freq = np.fft.fftfreq(np.shape(secondary)[0], mhzperbin)
    max_conj_freq = np.max(conj_freq)
    min_conj_freq = np.min(conj_freq)
    conj_time = np.fft.fftfreq(np.shape(secondary)[1], minperbin*60) * 1000
    max_conj_time = np.max(conj_time)
    min_conj_time = np.min(conj_time)
    return min_conj_freq, max_conj_freq, min_conj_time, max_conj_time


