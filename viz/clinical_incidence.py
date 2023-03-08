# %%
# Plot clinical incidence (cases) for TP13 releases
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %%
palette = sns.light_palette("seagreen")
sns.set(font='Arial', font_scale=0.80)
# sns.set_palette(palette)
sns.set_style('ticks')
# %%
# read in human states data
df = pd.read_csv(
    '~/SoftwarePaper3/STP/ANALYZED/E_15_00500_00450000000_000100000000_0017500_0011700_0000000_0098700_0098700/H_Mean_0001.csv')
# df = pd.read_csv('/Users/agastyamondal/H_Mean_0001.csv')
# extract labels for prevalence states
labels = "^clin_inc"
time = df['Time']
df = df.filter(regex=labels)

# calculate all-ages prevalence
age_idx = ["00_01", "01_02", "02_03", "03_04", "04_05", "05_06"]
age_idx = {
    '00_01': '0-5',
    '01_02': '5-17',
    '02_03': '17-40',
    '03_04': '40-60',
    '04_05': '60-99',
    '05_06': '99+'
}
STP_pop = 223000
clininc_df = pd.DataFrame()
for idx in age_idx.keys():
    clininc_df[age_idx[idx]] = df.filter(
        regex=idx).sum(axis=1) * (STP_pop/1000)

clininc_df['All-ages'] = clininc_df.sum(axis=1)
clininc_df['Time'] = time

# scale by population/1000


# %%
rel_times = [730, 737, 744, 751, 758, 765, 772, 779]
prev_df_melt = clininc_df.melt('Time', var_name='cols', value_name='vals')
prev_df_melt.rename(columns={'cols': 'Age Group (yrs)'}, inplace=True)

g = sns.lineplot(x="Time", y="vals", hue='Age Group (yrs)', data=prev_df_melt)
g.set_title('Epidemiological Dynamics - Clinical Incidence of Malaria')
g.set(xlabel='Time (day)', ylabel='Cases/1000')
g.set(xlim=(400))
for i in rel_times:
    x = g.axvline(i, linewidth=0.075)
    x.set_zorder(0)


# %%
g.figure.savefig('clin_inc.tiff', dpi=300)


# %%
