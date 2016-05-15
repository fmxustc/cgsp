import numpy as np
import pandas as pd
import time
import astropy.io.fits as ft
import astropy.wcs as wcs
import subprocess
import warnings
from termcolor import cprint as log
import seaborn as sns
import matplotlib.pyplot as plt
import platform


def cal(path, name, center):

    def dist_ceil(tx, ty):
        return np.ceil(np.sqrt((px - tx) ** 2 + (py - ty) ** 2))

    file = path+name+'_r.fits'
    hdu = ft.open(file)

    # get header and flux values
    header = hdu[0].header
    data = hdu[0].data

    # convert ra and dec to x and y
    wx, wy = center[0], center[1]
    px, py = np.round(wcs.WCS(header).wcs_world2pix(wx, wy, 0))

    # use the software sextractor to dectect the file and get the segmentation map and background map
    conf_sex = '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME seg.fits'
    _sp = subprocess.Popen('%s %s %s' % (sex, file, conf_sex), shell=True, executable=shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    _sp.wait()
    conf_sex = '-CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME bg.fits'
    _sp = subprocess.Popen('%s %s %s' % (sex, file, conf_sex), shell=True, executable=shell, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    _sp.wait()

    # if necessary, the following code will get the rms value of the background
    output = _sp.stdout.readlines()
    rms = float(output[output.index(b'\x1b[1M> Scanning image\n') - 1].split()[4])

    # collect the pixels belong to the galaxy
    seg = ft.open('seg.fits')[0].data
    y, x = np.where(seg == seg[py][px])

    # calculate quasi-petrosion-isophote
    arg = np.argsort(data[y, x])[::-1]
    f = np.sort(data[y, x])[::-1]
    n = len(f)

    ff = np.copy(f)
    for k in range(1, n):
        ff[k] += ff[k-1]
    eta = f-0.2*np.array(list(map(lambda t: ff[t]/(t+1), range(n))))
    try:
        qpr = float(np.argwhere(eta < 0)[0])
    except IndexError:
        qpr = n

    # calculate gini index
    gf = f[:qpr][::-1]
    gini = sum(list(map(lambda l: (2*(l+1)-qpr-1)*gf[l], np.arange(qpr))))/(gf.mean()*qpr*(qpr-1))

    # calculate concentration index
    dst = list(map(dist_ceil, x, y))
    radius = np.ceil(max(dst) + 1)
    sb = np.zeros(radius)
    ssb = np.copy(sb)
    for k in np.arange(qpr):
        d = dst[arg[k]]
        sb[d] += f[k] / (2 * d + 1)
    ssb = sb * [2 * d + 1 for d in np.arange(len(sb))]
    for k in range(1, len(ssb)):
        ssb[k] += ssb[k-1]
    try:
        cr = np.argwhere(sb < rms * 2)[0]
        concentration = ssb[0.3 * cr] / ssb[cr]
    except IndexError:
        concentration = np.nan

    # calculate moment index
    mk = int(np.argwhere(ff > ff[qpr]*0.2)[0])
    mf = list(map(lambda l: f[l]*((y[arg[l]]-py)**2+(x[arg[l]]-px)**2), np.arange(qpr)))
    moment = np.sum(mf[:mk])/np.sum(mf)

    # calculate asymmetry index
    asymmetry = 0
    bg = ft.open('bg.fits')[0].data
    ad = data - bg
    asum = f[:qpr].sum()
    for k in np.arange(qpr):
        ix, iy = x[arg[k]], y[arg[k]]
        jx, jy = py*2-iy, px*2-ix
        if seg[jy][jx] == seg[py][px]:
            ad[jy][jx] *= 1
        else:
            ad[jy][jx] *= 0
        asymmetry += abs(ad[iy][ix]-ad[jy][jx])/asum

    radius = np.ceil(np.sqrt(radius))
    return radius, gini, moment, asymmetry, concentration


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    start_time = time.clock()

    if platform.system() == 'Linux':
        settings = pd.read_table('linux_setting.param', sep=' ', header=0, index_col='obj')
    else:
        settings = pd.read_table('mac_setting.param', sep=' ', header=0, index_col='obj')

    type1_fits_directory = settings.ix['type1', 'path']
    type2_fits_directory = settings.ix['type2', 'path']
    shell = settings.ix['shell', 'path']
    ds9 = settings.ix['ds9', 'path']
    sex = settings.ix['sex', 'path']

    catalog = pd.read_csv('list.csv')
    catalog = catalog[catalog.Z1 < 0.05]
    catalog.index = range(len(catalog))

    calculated_set1 = {}
    calculated_set2 = {}

    for i in range(len(catalog)):
        t0 = time.clock()
        ctl = catalog.ix[i]
        if ctl.NAME1 not in calculated_set1:
            r, g, m, a, c = cal(type1_fits_directory, ctl.NAME1, [ctl.RA1, ctl.DEC1])
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
            r, g, m, a, c = cal(type2_fits_directory, ctl.NAME2, [ctl.RA2, ctl.DEC2])
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
        log('processed in %f seconds' % (t1-t0), 'blue')

    catalog.to_csv('data.csv', columns=['NAME1', 'R1', 'G1', 'M1', 'A1', 'C1', 'NAME2', 'R2', 'G2', 'M2', 'A2', 'C2'],
                   index_label=['INDEX'], sep=' ', float_format='%e')
    end_time = time.clock()
    log('@The function takes %f seconds to complete' % (end_time - start_time), 'grey', attrs=['bold'])
