import numpy as np
from astropy.io import fits
import psrchive
import matplotlib.pyplot as plt


def main(args,archive):
    if archive.endswith('_dysp.txt') or archive.endswith('_combine.txt'):
        intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = load_txt(archive)
    else:
        try:
            intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = load_txt('%s_dysp.txt' % (archive))
        except (IOError, ValueError):
            create(args,archive)
            intensity, minFreq, maxFreq, mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = load_txt('%s_dysp.txt' % (archive))

    return intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq
 
        

def load_txt(archive):
    with open(archive,'r') as myfile:
        lines=myfile.readlines()
    meta_dict=lines[-1]    
    _, _, minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = meta_dict.split(' ')
    minFreq = float(minFreq)
    maxFreq = float(maxFreq)
    mjd_start = float(mjd_start)
    mjd_end = float(mjd_end)
    length = float(length)
    diff = float(diff)
    midfreq = float(midfreq)
    intensity=np.genfromtxt('%s' % archive)
    intensity=intensity/np.max(intensity)
    return intensity, minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq 

            



def create(args,archive):
    arch = psrchive.Archive_load(archive)
    minFreq = arch.get_centre_frequency() - (arch.get_bandwidth() / 2)
    maxFreq = arch.get_centre_frequency() + (arch.get_bandwidth() / 2)
    midfreq = arch.get_centre_frequency()
    length = arch.integration_length()
    name = arch.get_source()
    mjd_start = float(arch.start_time().strtempo())
    mjd_end = float(arch.end_time().strtempo())
    diff = maxFreq - minFreq
    coord = arch.get_coordinates()
    ra, dec = coord.getHMSDMS().split(' ')
    site = arch.get_telescope()


    f=fits.open(archive)
    
    #Read data from psrfits files
    fref = f['HISTORY'].data['CTR_FREQ'][1]
    DM = f['SUBINT'].header['DM']
    freq = f['SUBINT'].data['DAT_FREQ'][0]
    tbin = f['HISTORY'].data['Tbin'][0]
    data=f['SUBINT'].data['DATA']
    weights = f['SUBINT'].data['DAT_WTS']
    
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
    template[SN < 7.] = 0
    I *= template[np.newaxis, np.newaxis, :]
    intensity = I.mean(-1)
    intensity /= np.max(intensity)
    if args.fillempty:
        intensity = fill(intensity,weights)
    else:
        intensity *= weights
    
    intensity=np.rot90(intensity,1)
    
    #plt.imshow(intensity,aspect='auto',extent=[0,length/60,minFreq,maxFreq],vmax=1, vmin=-0.1,cmap='jet', interpolation='None')
    #plt.imshow(intensity,aspect='auto',extent=[0,intensity.shape[1],0,intensity.shape[0]],cmap='jet', interpolation='None')
    #plt.xlabel("time(mins)")
    #plt.ylabel("frequency(Mhz)")
    #plt.colorbar(use_gridspec=True)
    #plt.savefig('%s_.png' % archive )
    #plt.show()
    
    
    np.savetxt('%s_dysp.txt' % archive,intensity, fmt='%1.6f')
    with open('%s_dysp.txt' % archive,"a") as myfile:
        myfile.write("\n # %f %f %f %f %s %s %s %s %f %f %f" % (minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq))
    print('Created %s_dysp.txt' % (archive))


def fill(intensity,weights):
    nsub,nchan=np.shape(intensity)
    intensity_w = intensity * weights
    for (isub, ichan) in np.argwhere(weights==0):
        if intensity[isub,ichan]>0.3 and 0 < ichan < nchan-1 and (intensity_w[isub,ichan-1]>0.3 or intensity_w[isub,ichan-1]>0.3):
            intensity_w[isub,ichan]=(intensity_w[isub,ichan-1]+intensity_w[isub,ichan+1])/2.0
    for (isub, ichan) in np.argwhere(weights==0):
        if intensity[isub,ichan]>0.3 and 0 < isub < nsub-1 and (intensity_w[isub-1,ichan]>0.3 or intensity_w[isub-1,ichan]>0.3):
            intensity_w[isub,ichan]=(intensity_w[isub-1,ichan]+intensity_w[isub+1,ichan])/2.0
    return intensity_w
    



