# %%

# Plot allele frequency for TP13 in STP

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# %%
palette = sns.light_palette("seagreen")
sns.set(font='Arial', font_scale=0.80)
sns.set_palette(palette)
sns.set_style('ticks')
# %%
# read in data for male and female mosquitoes
df_female = pd.read_csv("~/SoftwarePaper3/STP/ANALYZED/E_15_00500_00450000000_000100000000_0017500_0011700_0000000_0098700_0098700/F_Mean_0001.csv")
df_male = pd.read_csv("~/SoftwarePaper3/STP/ANALYZED/E_15_00500_00450000000_000100000000_0017500_0011700_0000000_0098700_0098700/M_Mean_0001.csv")
# %%
# rescale so everything is in terms of frequency
time = df_female['Time']
df_female = df_female.drop('Time', axis=1)
df_male = df_male.drop('Time', axis=1)




# %%
df_female = df_female.div(df_female.sum(axis=1), axis=0)
df_male = df_male.div(df_male.sum(axis=1), axis=0)
# %%
# Convert to long format
df_female['Time'] = time
df_male['Time'] = time
df_male_melt = df_male.melt('Time', var_name='cols', value_name='vals')
df_male_melt.rename(columns={'cols': 'Genotype'}, inplace=True)
df_female_melt = df_female.melt('Time', var_name='cols', value_name='vals')
df_female_melt.rename(columns={'cols': 'Genotype'}, inplace=True)

# %%
# Plot
fig, axes = plt.subplots(1, 2, figsize=(20,8))
fig.suptitle('Entomological Dynamics - Allele Frequency')
rel_times = [730, 737, 744, 751, 758, 765, 772, 779]

g = sns.lineplot(ax=axes[0],x="Time", y="vals", hue='Genotype', data=df_male_melt)
g.set(xlim=(400, 1000))
axes[0].set_title('Male Mosquitoes')
axes[0].set(xlabel='Time (day)', ylabel='Allele frequency')
for i in rel_times:
    x = axes[0].axvline(i) 
    x.set_zorder(0)

g = sns.lineplot(ax=axes[1],x="Time", y="vals", hue='Genotype', data=df_female_melt)
g.set(xlim=(400, 1000))
axes[1].set_title('Female Mosquitoes')
axes[1].set(xlabel='Time (day)', ylabel='Allele frequency')
for i in rel_times:
    x = axes[1].axvline(i) 
    x.set_zorder(0)

# save
plt.savefig('allele_freq.tiff', dpi=300)    


# %%
