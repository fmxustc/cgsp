import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings


def normalize(a):
    mn, mx = a.min(), a.max()
    return (a-mn)/(mx-mn)


def do_nothing(a):
    return a


def draw(d, item, ax, label, func=do_nothing, model='kde'):
    x = item+'1'
    y = item+'2'
    a = np.sort(func(d[x].values))
    b = np.sort(func(d[y].values))
    mn, mx = max(a.min(), b.min()), min(a.max(), b.max())
    bins_hist = np.linspace(mn, mx, 20)
    bins_kde = np.linspace(mn, mx, 100)
    if model == 'hist':
        ax.set(xlabel=name[item], ylabel='Gaussian Kernel Density')

        sns.distplot(a, bins=bins_hist, ax=ax, kde=False, hist=True,
                     hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1, "color": 'b', "label": label+'1'})
        sns.distplot(b, bins=bins_hist, ax=ax, kde=False, hist=True,
                     hist_kws={"histtype": "step", "linewidth": 2, "alpha": 1, "color": 'r', "label": label+'2'})
    else:
        ax.set(xlabel=name[item], ylabel='Count')
        sns.distplot(a, bins=bins_kde, ax=ax, label=label+'1', color='b', kde=True, hist=False)
        sns.distplot(b, bins=bins_kde, ax=ax, label=label+'2', color='r', kde=True, hist=False)


# R Band figure
warnings.filterwarnings('ignore')

# get data
data = pd.read_csv('data2.csv', sep=' ')
# picture settings
fig, ar = plt.subplots(2, 2, figsize=(16, 8), sharex=False, sharey=False)
fig.suptitle('R Band Differences')

name = {
    'G': 'Gini Index',
    'M': 'log(Moment_20) Index',
    'A': 'log(Asymmetry) Index',
    'C': 'Concentration Index',
}
draw(data, 'G', ar[0, 0], 'type', model='kde')
draw(data, 'M', ar[0, 1], 'type', model='kde', func=np.log10)
draw(data, 'A', ar[1, 0], 'type', model='kde')
draw(data, 'C', ar[1, 1], 'type', model='kde')
plt.show()
