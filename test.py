import subprocess
import pandas as pd
import platform
import warnings
import yapf


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

warnings.filterwarnings('ignore')
data = pd.read_csv('data2.csv', sep=' ', index_col='INDEX')
sample = data[data.A2 > 0.5]
print(len(sample))
# sample = sample.drop_duplicates('NAME1')
sample.index = range(len(sample))
t = 25
# print(sample.ix['J143318.48+344404.4'].NAME1)
ll = ['J143318.48+344404.4', 'J144545.12+513450.9', 'J130515.26+045234.3', 'J074535.90+243443.7', 'J131517.27+442425.6']
# print()
# for k in range(1):
#
#     objs = ''
#     for i in range(len(ll)):
#         x = sample.ix[ll[i]]
#         if str(x.index[0]) != 'NAME1':
#             x = x.ix[0]
#             objs = objs+'/Users/franky/Desktop/type2cut/'+x.NAME2+'_r.fits '
#         else:
#             objs = objs+'/Users/franky/Desktop/type2cut/'+x.NAME2+'_r.fits '
#     cmd = '~/bin/ds9'
#     view = '-geometry 1920x1080 -view layout vertical -view panner no -view buttons no -view info yes -view magnifier no -view colorbar no'
#     settings = '-invert -cmap value 1.75 0.275 -zoom 0.5 -minmax -log'
#     subprocess.Popen('%s %s %s %s' % (cmd, view, settings, objs), shell=True, executable='/bin/zsh')
# # for i in range(len(ll)):
# #     x = sample.ix[ll[i]]
# #     name = x.index[0] if str(x.index[0]) != 'NAME1' else x.NAME2
# #     conf_sex = '-CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME %s.fits' % name
# #
# #     _sp = subprocess.Popen('%s %s %s' %
# #                            (sex, type2_fits_directory + name + '_r.fits', conf_sex),
# #                            shell=True,
# #                            executable=shell,
# #                            stdout=subprocess.PIPE,
# #                            stderr=subprocess.STDOUT)
# #     _sp.wait()
objs = ''
for i in range(len(sample)):
    x = sample.ix[i]
    objs = objs+'/Users/franky/Desktop/type2cut/'+x.NAME2+'_r.fits '
    cmd = '~/bin/ds9'
    view = '-geometry 1920x1080 -view layout vertical -view panner no -view buttons no -view info yes -view magnifier no -view colorbar no'
    settings = '-invert -cmap value 1.75 0.275 -zoom 0.5 -minmax -log'
    subprocess.Popen('%s %s %s %s' % (cmd, view, settings, objs), shell=True, executable='/bin/zsh')