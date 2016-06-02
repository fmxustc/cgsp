import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
import subprocess
import platform


def normalize(a):
    mn, mx = a.min(), a.max()
    return (a - mn) / (mx - mn)


def do_nothing(a):
    return a


def draw(d, item, ax, label, func=do_nothing, model='kde'):
    x = item + '1'
    y = item + '2'
    a = np.sort(func(d[x].values))
    b = np.sort(func(d[y].values))
    mn, mx = max(a.min(), b.min()), min(a.max(), b.max())
    bins_hist = np.linspace(min(a.min(), b.min()), max(a.max(), b.max()), 10)
    bins_kde = np.linspace(mn, mx, 100)
    if model == 'hist':
        ax.set(xlabel=name[item], ylabel='Count')

        sns.distplot(a,
                     bins=bins_hist,
                     ax=ax,
                     kde=False,
                     hist=True,
                     hist_kws={"histtype": "step",
                               "linewidth": 2,
                               "alpha": 1,
                               "color": 'b',
                               "label": label + '1'})
        sns.distplot(b,
                     bins=bins_hist,
                     ax=ax,
                     kde=False,
                     hist=True,
                     hist_kws={"histtype": "step",
                               "linewidth": 2,
                               "alpha": 1,
                               "color": 'r',
                               "label": label + '2'})
    else:
        ax.set(xlabel=name[item], ylabel='Gaussian Kernel Density')
        sns.distplot(a,
                     bins=bins_kde,
                     ax=ax,
                     label=label + '1',
                     color='b',
                     kde=True,
                     hist=False)
        sns.distplot(b,
                     bins=bins_kde,
                     ax=ax,
                     label=label + '2',
                     color='r',
                     kde=True,
                     hist=False)


def pic():
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
    sample = data0[(data0.A2 > 0.05) & (data0.A2 < 0.10)]
    sample.index = range(len(sample))
    objs = ''
    for k in range(len(sample)):
        x = sample.ix[k]
        objs = objs + type1_fits_directory + x.NAME1 + '_r.fits '
        # objs = objs + type2_fits_directory + x.NAME2 + '_r.fits '

    view = '-geometry 1920x1080 -view layout vertical -view panner no -view buttons no -view info yes -view magnifier no -view colorbar no'
    stts = '-invert -cmap value 1.75 0.275 -zoom 0.5 -minmax -log'
    subprocess.Popen('%s %s %s %s' % (ds9, view, stts, objs), shell=True, executable='/bin/zsh')
    exit(0)

warnings.filterwarnings('ignore')

# get data
data0 = pd.read_csv('data0.csv', sep=' ')
data1 = pd.read_csv('data1.csv', sep=' ')
data2 = pd.read_csv('data2.csv', sep=' ')
data3 = pd.read_csv('data3.csv', sep=' ')
# picture settings
fig, ar = plt.subplots(2, 2, figsize=(16, 8), sharex=False, sharey=False)
fig.suptitle('R Band Differences')

name = {
    'G': 'Gini Index',
    'M': 'Moment_20 Index',
    'A': 'Asymmetry Index',
    'C': 'Concentration Index',
}
# pic(), func=np.log10
# print(len(data[data.A1 >0.6]),len(data[data.A2 >0.6]))
draw(data0, 'M', ar[0, 0], 'type', model='kde', func=np.log10)
draw(data1, 'M', ar[0, 1], 'type', model='kde', func=np.log10)
draw(data2, 'M', ar[1, 0], 'type', model='kde', func=np.log10)
draw(data3, 'M', ar[1, 1], 'type', model='kde', func=np.log10)
# draw(data, 'M', ar[0, 1], 'type', model='hist', func=np.log10)
# draw(data, 'A', ar[1, 0], 'type', model='kde')
# draw(data, 'C', ar[1, 1], 'type', model='kde', func=np.log10)

plt.show()
