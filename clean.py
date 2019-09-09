import numpy as np
import datetime
import matplotlib.pyplot as plt
import scipy.optimize
import argparse
import psrchive
from scipy.optimize import leastsq
from pandas import Series
from astropy.io import fits

parser=argparse.ArgumentParser(description="Select the archive.")
parser.add_argument('archives', nargs='+', help='The chosen archives')
parser.add_argument('-rl','--removerfi_level', type=int, default=5, help='input the level of remove RFI.The smaller of input value, the more RFIs is removed.default=7,range=[1-10]')
parser.add_argument('-plot', '--plot', action='store_true', help='plot the dynamic spectrum, ACF and secendary')

args = parser.parse_args()
archive_list = args.archives

def main(args,archive_list):
    for archive in archive_list:
        ar = psrchive.Archive_load(archive)
        clean(ar, args, archive)
        name = str(archive) + '_cleaned'
        ar.unload(name)
        print ("Cleaned archive: %s" % name)


def clean(ar, args, archive):
    #Read data from psrfits files
    f=fits.open(archive)
    fref = f['HISTORY'].data['CTR_FREQ'][1]
    DM = f['SUBINT'].header['DM']
    freq = f['SUBINT'].data['DAT_FREQ'][0]
    tbin = f['HISTORY'].data['Tbin'][0]
    data=f['SUBINT'].data['DATA']
    old_weights = f['SUBINT'].data['DAT_WTS']
    # Loop through each frequency, incoherently shift in pulse phase
    I = np.zeros(data.shape, dtype=np.int)
    for i in np.arange(len(freq)):
        dt = (1 / 2.41e-4) * DM * (1. / fref**2. - 1. / freq[i]**2.)
        tshift = np.int(np.rint(dt / tbin))
        I[:, :, i, :] = np.roll(data[:, :, i, :], tshift, axis=2)
    # Subtract median intensity from each channel each bin
    I = I.squeeze()
    I = I - np.median(I, axis=-1, keepdims=True)
    
    # Create time and freq averaged profile
    template = I.mean(0).mean(0)
    template /= np.max(template)

    # get noise from area without pulse
    midp = len(template) // 2
    template_med = np.sort(template)[midp - midp / 2:midp + midp / 2]
    mean = np.mean(template_med)
    std = np.std(template_med)
    SN = (template - mean) / std

    # find the RFI
    residuals = get_residual(I,template)
    fftvals = get_fftval(residuals)
    autocorrvals = get_autocorr(I,residuals,old_weights)
    test_weights = get_weights(I,residuals,fftvals,old_weights,autocorrvals,args)
    test_weights = check(I,template,SN,test_weights,fftvals,autocorrvals)
    
    if args.plot:
        plot(I,template,test_weights,SN,ar)

    set_weights_archive(ar, test_weights)
    return ar

def get_residual(I,temp):
    nsub,nchan,nbin=I.shape
    residuals=np.empty_like(I)   
    for isub in np.arange(nsub):
        for ichan in np.arange(nchan):
            subprof=I[isub,ichan]
            err = lambda amp: amp * temp - subprof
            params,status = leastsq(err, [1.0])
            if status not in (1,2,3,4):
                print( "Bad status for least squares fit when removing profile.")
                residuals[isub,ichan]=0
            else:
                residual_sub=np.asarray(err(params))
                residuals[isub,ichan]=residual_sub
    return residuals

def get_fftval(residuals):
    residuals_masked=np.ma.masked_equal(residuals, 0)
    function = lambda resudials_masked, axis: np.max(np.abs(np.fft.rfft(resudials_masked - np.expand_dims(resudials_masked.mean(axis=axis), axis=axis),axis=axis)), axis=axis)
    diag=(function(residuals_masked, axis=-1))
    chan_scaled = np.abs(channel_scaler(diag)) / 10
    subint_scaled = np.abs(subint_scaler(diag)) / 10
    fftvals = np.max((chan_scaled, subint_scaled), axis=0)
    fftvals /= np.max(fftvals)
    return fftvals
def channel_scaler(array2d):
    """For each channel scale it.
    """
    scaled = np.empty_like(array2d)
    nchans = array2d.shape[1]
    for ichan in np.arange(nchans):
        with np.errstate(invalid='ignore', divide='ignore'):
            channel = array2d[:, ichan]
            median = np.ma.median(channel)
            channel_rescaled = channel - median
            mad = np.ma.median(np.abs(channel_rescaled))
            scaled[:, ichan] = (channel_rescaled) / mad
    return scaled

def subint_scaler(array2d):
    """For each sub-int scale it.
    """
    scaled = np.empty_like(array2d)
    nsubs = array2d.shape[0]
    for isub in np.arange(nsubs):
        with np.errstate(invalid='ignore', divide='ignore'):
            subint = array2d[isub, :]
            median = np.ma.median(subint)
            subint_rescaled = subint - median
            mad = np.ma.median(np.abs(subint_rescaled))
            scaled[isub, :] = (subint_rescaled) / mad
    return scaled
def get_autocorr(I,residuals,old_weights):
    nsub,nchan,nbin=I.shape
    autocorrvals = np.empty_like(old_weights)
    
    for isub in np.arange(nsub):
        for ichan in np.arange(nchan):
            lbin=residuals[isub, ichan]
            series=Series(lbin)
            if series.autocorr()==series.autocorr():
                autocorrvals[isub, ichan] = series.autocorr()
            else:
                autocorrvals[isub, ichan]=0
    autocorrvals /= np.max(autocorrvals)
    return autocorrvals


def get_weights(I,residuals,fftvals,old_weights,autocorrvals,args):
    level=args.removerfi_level  
    new_weights=np.ones(np.shape(old_weights))
    for (isub, ichan) in np.argwhere(fftvals >= level*np.mean(fftvals)):
        new_weights[isub, ichan]=0
    for (isub, ichan) in np.argwhere(autocorrvals >= level*np.mean(autocorrvals)):
        new_weights[isub, ichan]=0
    return new_weights


def check(I,temp,SN,test_weights,fftvals,autocorrvals):
    nsub,nchan,nbin=I.shape
    temp[SN < 7.] = 0
    #I_new = I * test_weights[:,:,np.newaxis]
    I_new = I * temp[np.newaxis, np.newaxis, :]
    inte = I_new.mean(-1)
    inte /= np.max(inte)
    inte_test = inte * test_weights
    inte_test /= np.max(inte_test)    
    
    for (isub, ichan) in np.argwhere(test_weights==0):
        if 0 < ichan < nchan-1 and inte[isub, ichan] >= 0.3 and autocorrvals[isub, ichan] <=0.4:
            if (inte_test[isub, ichan-1] >=0.2 and np.abs(inte_test[isub, ichan-1]-inte[isub, ichan])<=0.2) or (inte_test[isub, ichan+1] >=0.2 and np.abs(inte_test[isub, ichan+1]-inte[isub, ichan])<=0.2):
                test_weights[isub, ichan] = 1

                
    for (isub, ichan) in np.argwhere(test_weights==0):
        if 0 < isub < nsub-1 and inte[isub, ichan] >= 0.3 and autocorrvals[isub, ichan] <=0.4:
            if (inte_test[isub-1, ichan] >=0.2 and np.abs(inte_test[isub-1, ichan]-inte[isub, ichan])<=0.2) or (inte_test[isub+1, ichan] >=0.2 and np.abs(inte_test[isub+1, ichan]-inte[isub, ichan])<=0.2):
                test_weights[isub, ichan] = 1

    return test_weights


    


def plot(I,template,test_weights,SN,ar):
    minFreq = ar.get_centre_frequency() - (ar.get_bandwidth() / 2)
    maxFreq = ar.get_centre_frequency() + (ar.get_bandwidth() / 2)
    length = ar.integration_length()
    template[SN < 7.] = 0
    I = I * test_weights[:,:,np.newaxis]
    I *= template[np.newaxis, np.newaxis, :]
    intensity = I.mean(-1)
    intensity=intensity/np.max(intensity)
    intensity=np.rot90(intensity,1)
    plt.imshow(intensity,aspect='auto',extent=[0,length/60,minFreq,maxFreq],vmax=1, vmin=-0.1,cmap='jet', interpolation='None')
    #plt.imshow(intensity,aspect='auto',extent=[0,intensity.shape[1],0,intensity.shape[0]],cmap='jet', interpolation='None')
    #plt.xlabel("time(mins)")
    #plt.ylabel("frequency(Mhz)")
    plt.colorbar(use_gridspec=True)
    #plt.savefig('%s_.png' % archive )
    plt.show()



def set_weights_archive(archive, test_weights):
    #Apply the weigths to an archive according to the test results.
    
    for (isub, ichan) in np.argwhere(test_weights):
        integ = archive.get_Integration(int(isub))
        integ.set_weight(int(ichan), 1.0)
    for (isub, ichan) in np.argwhere(test_weights < 1):
        integ = archive.get_Integration(int(isub))
        integ.set_weight(int(ichan), 0.0)






if __name__ == "__main__":
    main(args,archive_list)
