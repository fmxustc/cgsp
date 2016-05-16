import numpy as np
import pandas as pd
import time
import astropy.io.fits as ft
import astropy.wcs as wcs
import subprocess
import warnings
from termcolor import cprint as log
import platform


def cal(path, name, center):

    hdu = ft.open(path + name + '_r.fits')

    # get header and flux values
    header = hdu[0].header
    data = np.copy(hdu[0].data)

    # convert ra and dec to x and y
    wx, wy = center[0], center[1]
    px, py = np.round(wcs.WCS(header).wcs_world2pix(wx, wy, 0))

    # use the software sextractor to dectect the file and get the segmentation map and background map
    conf_sex = '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME seg.fits'
    _sp = subprocess.Popen('%s %s %s' %
                           (sex, path + name + '_r.fits', conf_sex),
                           shell=True,
                           executable=shell,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    _sp.wait()
    conf_sex = '-CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME bg.fits'
    _sp = subprocess.Popen('%s %s %s' %
                           (sex, path + name + '_r.fits', conf_sex),
                           shell=True,
                           executable=shell,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
    _sp.wait()
    # if necessary, the following code will get the rms value of the background
    output = _sp.stdout.readlines()
    rms = float(output[output.index(b'\x1b[1M> Scanning image\n') - 1].split()[
        4])

    # collect the pixels belong to the galaxy
    seg = np.copy(ft.open('seg.fits')[0].data)
    bg = np.copy(ft.open('bg.fits')[0].data)
    y, x = np.where(seg == seg[py][px])

    # calculate quasi-petrosion-isophote
    sgf = data[y, x]
    arg = np.argsort(sgf)[::-1]
    f = sgf[arg]
    n = len(f)
    ff = np.copy(f)
    for k in range(1, n):
        ff[k] += ff[k - 1]
    eta = np.array([f[l] - 0.2 * ff[l] / (l + 1) for l in range(n)])
    try:
        qpr = int(np.argwhere(eta < 0)[0])
    except IndexError:
        qpr = n - 1
    f = f[:qpr]
    ff = ff[:qpr]

    # calculate distance between every pixel and the center
    df = np.array([np.sqrt((px - x[l])**2 + (py - y[l])**2)
                   for l in arg[:qpr]])

    # calculate cumulative flux, surface brightness, mean surface brightness
    radius = np.ceil(max(df)) + 1
    sb = np.zeros(radius)
    for k in range(qpr):
        sb[np.ceil(df[k])] += f[k]
    cf = np.copy(sb)
    for t in np.arange(1, radius):
        cf[t] += cf[t - 1]
    sb /= 2 * np.arange(radius) + 1
    msb = cf / np.arange(radius)**2

    # get the value of r50 depend on the cf array
    rcf = cf/cf[-1]
    r50 = int(np.argwhere(rcf > 0.5)[0])

    # calculate gini index
    gf = f[::-1]
    gini = np.sum((2 * np.arange(qpr) - qpr + 1) * gf) / (
        np.mean(gf) * qpr * (qpr - 1))

    # calculate concentration index
    r20 = int(np.argwhere(rcf > 0.2)[0])
    r80 = int(np.argwhere(rcf > 0.8)[0])
    concentration = r80/r20

    # calculate moment index
    mt = (df**2) * f
    rff = ff/ff[-1]
    moment = np.sum(mt[:np.argwhere(rff > 0.2)[0]])/np.sum(mt)

    # calculate asymmetry index:
    sda = (data[y, x] - bg[y, x])
    sx, sy = np.int_(np.round(2*px-x)), np.int_(np.round(2*py-y))
    sdb = (data[sy, sx] - bg[sy, sx])
    asymmetry = np.sum([abs(sda[l]-sdb[l]) if seg[sy[l]][sx[l]] == seg[py][px] else abs(sda[l]) for l in arg[:qpr]])/abs(np.sum(f))
    return radius, gini, moment, asymmetry, concentration


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    start_time = time.clock()

    if platform.system() == 'Linux':
        settings = pd.read_table('linux_setting.param',
                                 sep=' ',
                                 header=0,
                                 index_col='obj')
    else:
        settings = pd.read_table('mac_setting.param',
                                 sep=' ',
                                 header=0,
                                 index_col='obj')

    type1_fits_directory = settings.ix['type1', 'path']
    type2_fits_directory = settings.ix['type2', 'path']
    shell = settings.ix['shell', 'path']
    ds9 = settings.ix['ds9', 'path']
    sex = settings.ix['sex', 'path']
    blue = 'blue'
    red = 'red'
    green = 'green'

    catalog = pd.read_csv('list.csv')
    catalog = catalog[catalog.Z1 < 0.05]
    catalog.index = range(len(catalog))

    calculated_set1 = {}
    calculated_set2 = {}

    for i in range(len(catalog)):
        t0 = time.clock()
        ctl = catalog.ix[i]
        if ctl.NAME1 not in calculated_set1:
            r, g, m, a, c = cal(type1_fits_directory, ctl.NAME1, [ctl.RA1,
                                                                  ctl.DEC1])
            catalog.at[i, 'G1'] = g
            catalog.at[i, 'M1'] = m
            catalog.at[i, 'A1'] = a
            catalog.at[i, 'C1'] = c
            catalog.at[i, 'R1'] = r
            calculated_set1[ctl.NAME1] = i
        else:
            j = calculated_set1[ctl.NAME1]
            catalog.at[i, 'G1'] = catalog.at[j, 'G1']
            catalog.at[i, 'M1'] = catalog.at[j, 'M1']
            catalog.at[i, 'A1'] = catalog.at[j, 'A1']
            catalog.at[i, 'C1'] = catalog.at[j, 'C1']
            catalog.at[i, 'R1'] = catalog.at[j, 'R1']
        if ctl.NAME2 not in calculated_set2:
            r, g, m, a, c = cal(type2_fits_directory, ctl.NAME2, [ctl.RA2,
                                                                  ctl.DEC2])
            catalog.at[i, 'G2'] = g
            catalog.at[i, 'M2'] = m
            catalog.at[i, 'A2'] = a
            catalog.at[i, 'C2'] = c
            catalog.at[i, 'R2'] = r
            calculated_set2[ctl.NAME2] = i
        else:
            j = calculated_set2[ctl.NAME2]
            catalog.at[i, 'G2'] = catalog.at[j, 'G2']
            catalog.at[i, 'M2'] = catalog.at[j, 'M2']
            catalog.at[i, 'A2'] = catalog.at[j, 'A2']
            catalog.at[i, 'C2'] = catalog.at[j, 'C2']
            catalog.at[i, 'R2'] = catalog.at[j, 'R2']
        t1 = time.clock()
        log('INDEX==> %d' % i, 'cyan', end='   ', attrs=['bold'])
        log('OBJECT==> %s %s' % (ctl.NAME1, ctl.NAME2), 'green', end='    ')
        log('processed in %f seconds' % (t1 - t0), 'blue')

    catalog.to_csv('data2.csv',
                   columns=['NAME1', 'R1', 'G1', 'M1', 'A1', 'C1', 'NAME2',
                            'R2', 'G2', 'M2', 'A2', 'C2'],
                   index_label=['INDEX'],
                   sep=' ',
                   float_format='%e')
    end_time = time.clock()
    log('@The function takes %f seconds to complete' % (end_time - start_time),
        'grey',
        attrs=['bold'])
