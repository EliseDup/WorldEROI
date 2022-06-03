import matplotlib.pyplot as plt
import model_methods
import geopandas as gpd
from shapely.geometry import Point, Polygon


def plot_e_out_eroi_wind(df):
    plt.figure()
    df_onshore = model_methods.df_cum_eout_eroi(df, "wind_onshore_e", "wind_onshore_eroi")
    df_offshore = model_methods.df_cum_eout_eroi(df, "wind_offshore_e", "wind_offshore_eroi")
    df_total = model_methods.df_cum_eout_eroi(df, "wind_e", "wind_eroi")
    plt.plot(df_total['e_cum'], df_total['eroi'], label='Total')
    plt.plot(df_onshore['e_cum'], df_onshore['eroi'], label='Onshore')
    plt.plot(df_offshore['e_cum'], df_offshore['eroi'], label='Offshore')
    plt.grid(True, color="#93a1a1", alpha=0.3)
    plt.legend(loc='upper right', fancybox=True, shadow=True)
    plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
    plt.ylabel(r'EROI', labelpad=6, fontsize=12)
    plt.show()


# e_label = pv_e or csp_e
# eroi_label = pv_eroi or csp_eroi
def plot_e_out_eroi_solar(df, e_label, eroi_label):
    df_solar = model_methods.df_cum_eout_eroi(df, e_label, eroi_label)
    plt.figure()
    plt.plot(df_solar['e_cum'], df_solar['eroi'])
    plt.grid(True, color="#93a1a1", alpha=0.3)
    plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
    plt.ylabel(r'EROI', labelpad=6, fontsize=12)
    plt.show()


# Geoplots
def geo_plot(df, column, label):
    map = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    crs = {'init': 'EPSG:4326'}
    geometry = [Point(xy) for xy in zip(df['Lon'], df['Lat'])]
    geo_df = gpd.GeoDataFrame(df, crs=crs, geometry=geometry)

    fig, ax = plt.subplots(figsize=(10, 10))
    map.boundary.to_crs(epsg=4326).plot(ax=ax, color='black')
    geo_df.plot(column=column, ax=ax, cmap='rainbow',
                legend=True, legend_kwds={'shrink': 0.3},
                markersize=10)
    ax.set_title(label)

    plt.show()
