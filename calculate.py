import numpy as np
import pandas as pd
import time
import astropy.io.fits as ft
import astropy.wcs as wcs
import subprocess
import warnings


def dist(x, y, ct):
    return (ct-x)**2+(ct-y)**2


def deal(path, name, center):
    file = path+name+'_r.fits'
    print('FILE IS --> %s' % file)
    hdu = ft.open(file)
    # get header and flux values
    header = hdu[0].header
    data = hdu[0].data
    # convert ra and dec to x and y
    wx, wy = center[0], center[1]
    px, py = np.round(wcs.WCS(header).wcs_world2pix(wx, wy, 0))
    # use the software sextractor to dectect the file and get the segmentation map
    conf_sex = '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME seg.fits'
    _sp = subprocess.Popen('%s %s %s' % (sex, file, conf_sex), shell=True, executable=shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    _sp.wait()
    # if necessary, the following code will get the rms value of the background
    # output = _sp.stdout.readlines()
    # rms = float(output[output.index(b'\x1b[1M> Scanning image\n') - 1].split()[4])
    # print('RMS = %f' % rms)
    # collect the pixels belong to the galaxy
    seg = ft.open('seg.fits')[0].data
    y, x = np.where(seg == seg[py][px])
    # arg = np.argsort(data[y, x])[::-1]
    f = np.sort(data[y, x])[::-1]
    ff = np.copy(f)
    for k in range(1, len(ff)):
        ff[k] += ff[k-1]
    # TODO try to use pixels belong to the galaxy to calculate petrosion radius
    return


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    begin_time = time.clock()

    settings = pd.read_table('linux_setting.param', sep=' ', header=0, index_col='obj')
    type1_fits_directory = settings.ix['type1', 'path']
    type2_fits_directory = settings.ix['type2', 'path']
    shell = settings.ix['shell', 'path']
    ds9 = settings.ix['ds9', 'path']
    sex = settings.ix['sex', 'path']

    catalog = pd.read_csv('list.csv')
    catalog = catalog[catalog.Z1 < 0.05]
    catalog.index = range(len(catalog))

    for i in range(1):
        subprocess.call('rm *.fits catalog.txt', shell=True, executable=shell)
        ctl = catalog.ix[i]
        deal(type2_fits_directory, ctl.NAME2, [ctl.RA2, ctl.DEC2])
    end_time = time.clock()
    print('@The function takes %f seconds to complete' % (end_time - begin_time))
