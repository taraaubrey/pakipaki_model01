import pandas as pd
import matplotlib.pyplot as plt

def open_sw(filename, date_col, value_col):
    df = pd.read_csv(
    filename,
    index_col=0,
    parse_dates=True
    )
    df['Date'] = pd.to_datetime(df[date_col], format='%Y-%m-%d %H:%M:%S')
    df['Level'] = pd.to_numeric(df[value_col], errors='coerce') /1000
    df.set_index('Date', inplace=True)
    df = df['Level'].squeeze()
    df = df.resample('D').mean()

    # normalize to zero mean
    df = df - df.mean()

    return df

def open_hbrc(filename, col):
    # Read the CSV file
    df = pd.read_csv(filename)
    
    # Parse the DateTime column to proper datetime format
    df['DateTime'] = pd.to_datetime(df['DateTime'], format='%d/%m/%Y %H:%M')
    
    # Set DateTime as index and convert to Series
    series = df.set_index('DateTime')[col].squeeze()
    series.index = series.index.normalize()  # Keep as datetime
    
    # Fill missing dates
    full_date_range = pd.date_range(start=series.index.min(), 
                                   end=series.index.max(), 
                                   freq='D')
    series = series.reindex(full_date_range)
    
    return series



def run_model(head, rain, evap, params, days=None):
    import pastas as ps

    ml = ps.Model(head)
    ml.add_noisemodel(ps.ArNoiseModel())

    # pumping start date
    d = rain-evap
    d[d>0] = 0
    # d[d<0] = 1
    d = d.abs()
    # Set May-October to 0 (months 5-10)
    d[d.index.month.isin([5, 6, 7, 8, 9, 10, 11])] = 0

    rain1d = rain.shift(1).copy()
    # add recharge model
    if params['with_evap']:
        rch = ps.rch.FlexModel()
        rm = ps.RechargeModel(rain1d, evap, recharge=rch, rfunc=ps.Gamma(), name="fast")
        ml.add_stressmodel(rm)
    else:
        rf = ps.StressModel(rain1d, rfunc=ps.Gamma(), name="fast", up=True, settings='prec')
        ml.add_stressmodel(rf)

    # add recharge model (slow recharge)
    if params['with_slow']:
        srain = rain.copy()
        c = params['rainfall_dmax']
        srain[srain>c] = c
        if params['slow_func']=='sum':
            inrf = srain.shift(days).rolling(f'{days}D').sum()
        elif params['slow_func']=='mean':
            inrf = srain.shift(days).rolling(f'{days}D').mean()
        else:
            raise ValueError("Invalid slow_func. Use 'sum' or 'mean'.")
        rm = ps.StressModel(inrf, rfunc=ps.Gamma(), name="slow", up=True, settings='prec')
        ml.add_stressmodel(rm)

    # add pumping model
    stm = ps.StressModel(d, rfunc=ps.Hantush(), name='pumping', up=False, settings='prec')
    ml.add_stressmodel(stm)

    # ml.solve(tmax='2023-12-06')
    ml.solve(
        tmin=rain.index[0].strftime('%Y-%m-%d'),
        )
    
    results = ml.plots.results()
    diagnostics = ml.plots.diagnostics()

    # save figures
    results[0].get_figure().savefig(f"pastas/output/{params['name']}_{days}D_results.png", dpi=300)
    diagnostics[0].get_figure().savefig(f"pastas/output/{params['name']}_{days}D_diagnostics.png", dpi=300)

    plot_reconstruction(ml, params, days)
    plt.close('all')
    return ml

    
def plot_reconstruction(ml, params, days):
    mlparams = ml.parameters

    # Create subplots
    fig, axes = plt.subplots(len(ml.stressmodels)+1, 1, figsize=(12, 8), sharex=True)

    i = 0
    dat = []
    for sm_name, sm in ml.stressmodels.items():
        p = []
        for param in sm.parameters.initial.index:  # Changed from rf.parameters to rm.parameters
            p_value = mlparams.loc[param, 'optimal']
            p.append(p_value)
        sm_dat = sm.simulate(p=p, tmin=sm.get_stress().index[0].strftime('%Y-%m-%d'))
        axes[i].plot(sm_dat.index, sm_dat.values)
        axes[i].set_title(f'{sm_name}')
        axes[i].set_ylabel('value')
        axes[i].grid(True)
        i += 1
        dat.append(sm_dat)

    # sum all the datasets
    head = ml.oseries.series
    initial_val = head.values[0]

    tot_dat = pd.concat(dat, axis=1).sum(axis=1)
    nat_dat = pd.concat(dat[:-1], axis=1).sum(axis=1)

    dat_i = tot_dat[head.index[0]]
    nat_i = nat_dat[head.index[0]]
    tot_dat = tot_dat + (initial_val - dat_i)
    nat_dat = nat_dat + (initial_val - nat_i)

    axes[i].plot(tot_dat.index, tot_dat.values)
    axes[i].set_title(f'Reconstructed water level')
    axes[i].set_ylabel('value')
    axes[i].grid(True)

    axes[i].plot(nat_dat.index, nat_dat.values, c='black')
    axes[i].set_title(f'Reconstructed water level')
    axes[i].set_ylabel('value')
    axes[i].grid(True)

    axes[i].plot(head.index, head.values, c='green')
    axes[i].set_title(f'Reconstructed water level')
    axes[i].set_ylabel('value')
    axes[i].grid(True)

    plt.tight_layout()
    # Optional: Save the figure
    plt.savefig(f"pastas/output/{params['name']}_{days}D_reconstruction.png", dpi=300)

    # save series
    dfout = pd.concat(dat, axis=1)
    dfout['level_est'] = tot_dat
    dfout['level_nat'] = nat_dat
    dfout['pct_rfrch'] = (dat[0] / (dat[0] + dat[1]))
    dfout.to_csv(f"pastas/output/{params['name']}_{days}D_series.csv")

    plt.close()