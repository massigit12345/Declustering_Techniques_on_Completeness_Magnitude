import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import contextily as ctx
import cartopy.mpl.gridliner as gridliner

earthquakes = pd.read_csv('Catalogo.csv')

colors = ['lightgreen','darkgreen','orange', 'red','darkred']  # Green to Red
cmap = mcolors.LinearSegmentedColormap.from_list('depth_cmap', colors)

min_lon = earthquakes['longitude'].min()
max_lon = earthquakes['longitude'].max()
min_lat = earthquakes['latitude'].min()
max_lat = earthquakes['latitude'].max()

source_3 = ctx.providers.Esri.WorldPhysical
source_8 = ctx.providers.USGS.USImagery
source_10 = ctx.providers.USGS.USTopo

sources = [source_3, source_8, source_10]

for srcs in sources:
    fig, ax = plt.subplots(figsize=(20, 18))
    scatter = ax.scatter(
        earthquakes['longitude'],
        earthquakes['latitude'],
        c=earthquakes['depth'],
        s=earthquakes['magnitude'] ** 3.4,  
        cmap=cmap,
        edgecolor='black',  
        alpha=0.7  
    )

    cbar = plt.colorbar(scatter, shrink=0.64, anchor=(0.0, 0.18))
    cbar.set_label('Depth')


    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Base Catalogue')

    ctx.add_basemap(ax,crs=ccrs.PlateCarree(), source=srcs, alpha=1)

    magnitudes = [1, 2, 3, 4, 5]

    max_magnitude = earthquakes['magnitude'].max()
    sizes = [size ** 3.4 for size in magnitudes]


    legend_handles = [plt.scatter([], [], s=size, marker='o', edgecolor='black', facecolor='white') for size in sizes]
    legend_labels = [f'Magnitude {magnitude}.0' for magnitude in magnitudes]

    legend = plt.legend(legend_handles, legend_labels, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                        borderaxespad=0.0, borderpad=2.5, labelspacing=2.0, frameon=True)
    legend.get_frame().set_facecolor('white')  # Set legend background color
    legend.get_frame().set_edgecolor('black')
    legend.get_frame().set_linewidth(2)
    ax.set_xlim(min_lon, max_lon)
    ax.set_ylim(min_lat, max_lat)

    ax.grid(linewidth=0.5, color='black', alpha=1, linestyle='--')

    plt.show()
