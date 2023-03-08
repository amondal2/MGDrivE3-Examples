# %%
# Plot spatial distribution of nodes on STP

import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt


points_df = pd.read_csv('../data/STP_grid_Sites.csv')
# %%

fig = px.scatter_mapbox(points_df,
                        lat="y", 
                        lon="x", 
                        color="TrapsType",
                        color_continuous_scale=px.colors.diverging.Tealrose,
                        zoom=9.5,
                        height=800,
                        width=800
)

fig.update_layout(mapbox_style="basic")
fig.update_layout(margin={"r":0,"t":0,"l":0,"b":0})
fig.update_layout(showlegend=False)
fig.write_image('stp_spatial.svg')


# %%
