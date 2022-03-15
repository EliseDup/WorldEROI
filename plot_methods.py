import matplotlib.pyplot as plt
from pandas import DataFrame
from pandas import concat
import model_methods


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
    plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)
    plt.show()


def plot_e_out_eroi_solar(df):
    df_solar = model_methods.df_cum_eout_eroi(df, 'solar_e', 'solar_eroi')
    #df_plot_solar_rooftop = DataFrame(index=rooftop_PV.index)
    #df_plot_solar_rooftop[['E_out', 'EROI']] = rooftop_PV.loc[:, ['Total E_out [EJ/year]', 'EROI']]
    #df_plot_solar = concat([df_plot_solar, df_plot_solar_rooftop])
    #df_plot_solar.drop(df_plot_solar[df_plot_solar['EROI'] <= 0].index, inplace=True)
    #df_plot_solar = df_plot_solar.sort_values(by=['EROI'], ascending=False)
    #df_plot_solar['cum_E_out'] = df_plot_solar['E_out'].cumsum()

    plt.figure()
    plt.plot(df_solar['e_cum'], df_solar['eroi'])
    plt.grid(True, color="#93a1a1", alpha=0.3)
    plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
    plt.ylabel(r'EROI [/]', labelpad=6, fontsize=12)
    plt.show()