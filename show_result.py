import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# R Band figure

# get data
data = pd.read_csv('data.csv', sep=' ')
# picture settings
fr, ar = plt.subplots(2, 2, figsize=(16, 8), sharex=False, sharey=False)
fr.suptitle('R Band Differences')

# Gini Coefficient
g1 = sns.distplot(data['G1'].values, bins=25, ax=ar[0, 0], hist=True, kde=False,
                  hist_kws={"histtype": "step", "linewidth": 2.5, "alpha": 1, 'color': 'r'})
g1.set(xlabel='Gini Coefficient', ylabel='Count', xlim=[0.4, 0.65])
# ar[0, 0].legend('ssss')
g2 = sns.distplot(data['G2'].values, bins=25, ax=ar[0, 0], hist=True, kde=False,
                  hist_kws={"histtype": "step", "linewidth": 2.5, "alpha": 1, 'color': 'b'})


#  Concentration Index
print(len(data[data.C1 > 0.2]))
rc = sns.distplot(np.log10(data[data.C1 > 0.2]['C1'].values), bins=100, color='r', ax=ar[0, 1], hist=False, kde=True, label='type1')
rc.set(xlabel='Concentration Index')
sns.distplot(np.log10(data[data.C1 > 0.2]['C2'].values), bins=100, color='b', ax=ar[0, 1], hist=False, kde=True, label='type2')

#  Moment Index
rm = sns.distplot(np.log10(data[data.M2 < 0.03]['M1'].values), bins=100, color='r', ax=ar[1, 0], hist=False, kde=True, label='type1')
rm.set(xlabel='Moment Index', ylabel='Gaussian Kernel Density')
sns.distplot(np.log10(data[data.M2 < 0.03]['M2'].values), bins=100, color='c', ax=ar[1, 0], hist=False, kde=True, label='type2')
# Asymmetry Index
rm2 = sns.distplot(np.log10(data['A1'].values), bins=100, color='b', ax=ar[1, 1], hist=False, kde=True, label='type1')
rm2.set(xlabel='Asymmetry Index', ylabel='Gaussian Kernel Density')
sns.distplot(np.log10(data['A2'].values), bins=100, color='c', ax=ar[1, 1], hist=False, kde=True, label='type2')
# rm = sns.distplot(np.log10(data['M1'].values), bins=100, color='g', ax=ar[1, 0], hist=False, kde=True, label='type1')
# rm.set(xlabel='Moment Index', ylabel='Gaussian Kernel Density')
# sns.distplot(np.log10(data['M2'].values), bins=100, color='c', ax=ar[1, 0], hist=False, kde=True, label='type2')
# rm2 = sns.distplot(np.log10(data2['M1'].values), bins=100, color='g', ax=ar[1, 1], hist=False, kde=True, label='type1')
# rm2.set(xlabel='Moment Index (agn removed)', ylabel='Gaussian Kernel Density')
# sns.distplot(np.log10(data2['M2'].values), bins=100, color='c', ax=ar[1, 1], hist=False, kde=True, label='type2')

# #  Asymmetry Index
# rm = sns.distplot(data[data.A1 < 0.6]['A1'].values, bins=100, color='brown', ax=ar[1, 1], hist=False, kde=True, label='type1')
# rm.set(xlabel='Asymmetry Index')
# sns.distplot(data[data.A2 < 0.6]['A2'].values, bins=100, color='hotpink', ax=ar[1, 1], hist=False, kde=True, label='type2')

plt.show()
