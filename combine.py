import argparse
import numpy as np
import create_dysp
parser=argparse.ArgumentParser(description="Select the archive.")
parser.add_argument('archives', nargs='+', help='The chosen archives')


args = parser.parse_args()
archive_list = args.archives
switch = 1
for archive in archive_list:
    print('adding the %s' % archive)
    intensity,minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length,diff,midfreq = create_dysp.main(args,archive)
    if switch:
        Ic = np.copy(intensity)
        length_c=length
        switch=0
    else:
        Ic = np.concatenate((Ic, intensity), axis=1)
        length_c =length_c + length
print(Ic.shape,length_c/60)
filename='J'+archive.rsplit('_')[2] + '_' + str(int(midfreq)) + '_combine.txt'
np.savetxt('%s' % filename, Ic, fmt='%1.6f')
with open('%s' % filename,"a") as myfile:
    myfile.write('\n # %f %f %f %f %s %s %s %s %f %f %f' % (minFreq, maxFreq,mjd_start,mjd_end,site,ra,dec,name,length_c,diff,midfreq))
