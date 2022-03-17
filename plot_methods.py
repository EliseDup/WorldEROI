import matplotlib.pyplot as plt
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
    plt.ylabel(r'EROI', labelpad=6, fontsize=12)
    plt.show()


def plot_e_out_eroi_pv(df):
    df_solar = model_methods.df_cum_eout_eroi(df, 'pv_e', 'pv_eroi')
    plt.figure()
    plt.plot(df_solar['e_cum'], df_solar['eroi'])
    plt.grid(True, color="#93a1a1", alpha=0.3)
    plt.xlabel(r'Cumulated Production [EJ/year]', labelpad=8, fontsize=12)
    plt.ylabel(r'EROI', labelpad=6, fontsize=12)
    plt.show()