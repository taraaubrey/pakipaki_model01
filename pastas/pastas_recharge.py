from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import pastas as ps
from helpers import *

rain = pd.read_csv(
    r'data\40256rf_mm_D_clean.csv', 
    index_col=0, 
    parse_dates=True
    ).squeeze()

rain1h = pd.read_csv(
    r'data\40256rf_mm_clean.csv', 
    index_col=0, 
    parse_dates=True
    ).squeeze()

rfaw = open_hbrc(
    r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\a_source\clm_climate_data\nz_hawkes_bay\HBRC_rainfall\2009_2024_daily-rainfall-for-AwanuiStreamFlume.csv",
    col='Rainfall_mm'
    )

evap = open_hbrc(
    r"C:\Users\tfo46\OneDrive - University of Canterbury\Tara_PhD\c_PhD\c_Data\a_source\clm_climate_data\nz_hawkes_bay\HBRC_rainfall\BridgePa_PET.csv",
    col='PET (y)'
    )

# lag rain by 1 day
rain1d = rain.shift(1)
# normalize rainfall from 0 to 1
rain1dn = (rain1d - rain1d.min()) / (rain1d.max() - rain1d.min())

# evap = pd.read_csv(
#     r'data\40256et_mm_D_clean.csv', 
#     index_col=0, 
#     parse_dates=True
#     ).squeeze()

head = pd.read_csv(
    r'data\pk4_LEVEL_D_clean.csv', 
    index_col=0, 
    parse_dates=True
    ).squeeze()

pw = open_sw(
    r'data\river-level-for-poukawa.csv',
    'Poukawa Stream at Douglas Road (x)',
    'Poukawa Stream at Douglas Road (y)')

aw = open_sw(
    r"data/river-level-for-awanui-s.csv", 
    'Awanui Stream at Flume (x)', 
    'Awanui Stream at Flume (y)')

hconf = open_hbrc(r"data/1501_Harrisons.csv", col='LEVEL')


# MODEL B ############################
"""
Same as model A but with recharge model.
"""
# run3: slow component
# run4: no slow component
# run6: sum of rainfall over days and low rainfall_dmax (5 mm)
params = {
    'name': 'run8',
    'rfdays': list(np.arange(3, 21)),
    'rainfall_dmax': 5,
    'with_evap': False,
    'with_slow': True,
    'slow_func': 'mean',  # 'sum' or 'mean'
    'input_head': 'aw',
    }
params_df = pd.DataFrame(list(params.items()), columns=['Parameter', 'Value'])
params_df.to_csv(f"pastas/output/{params['name']}_setup_params.csv", index=False)

def dictdf_to_df(dictdf, colname='optimal'):
    # concatonate dictionary into dataframe
    optimal_params = {}
    for days, params_df in dictdf.items():
        optimal_params[days] = params_df[colname]
    # Convert to DataFrame with days as columns and parameters as rows
    return pd.DataFrame(optimal_params).T

def getmlstats(ml):
    istats = {
        'rmse': ml.stats.rmse(),
        'r2': ml.stats.rsq(),
        'mae': ml.stats.mae(),
        'kge': ml.stats.kge()
    }
    return pd.DataFrame.from_dict(istats, orient='index', columns=['stats'])

if params['input_head'] == 'aw':
    # mask = (aw-aw.shift(1)) < 0.15 # remove really large spikes
    in_df = aw['2015-01-01':].copy()
else:
    in_df = head.copy()

if params['with_slow']:
    if isinstance(params['rfdays'], list):
        ml_params_dict = {}
        ml_stats_dict = {}
        for days in params['rfdays']:
            ml = run_model(in_df, rfaw, evap, params, days=days)
            ml_params_dict[days] = ml.parameters.copy()
            ml_stats_dict[days] = getmlstats(ml)
        ml_params = dictdf_to_df(ml_params_dict, colname='optimal')
        ml_stats = dictdf_to_df(ml_stats_dict, colname='stats')

    else:
        ml = run_model(in_df, rfaw, evap, params, days=params['rfdays'])
        ml_params = ml.parameters.copy()
        ml_stats = getmlstats(ml)
else:
    ml = run_model(in_df, rfaw, evap, params)
    ml_params = ml.parameters.copy()
    ml_stats = getmlstats(ml)

ml_params.to_csv(f"pastas/output/{params['name']}_params.csv")
ml_stats.to_csv(f"pastas/output/{params['name']}_stats.csv")

# MODEL C ############################
"""
Same as model B: no linear trend.
"""
tmin = '2024-01-15'
tmax = head.index[-1].strftime('%Y-%m-%d')
ml = ps.Model(head)
ml.add_noisemodel(ps.ArNoiseModel())

# add recharge model
rch = ps.rch.FlexModel()
rm = ps.RechargeModel(rain1d, evap, recharge=rch, rfunc=ps.Gamma(), name="recharge")
ml.add_stressmodel(rm)

# add streamflow model
swm = ps.StressModel(aw, rfunc=ps.Polder(), name="awanui", settings="waterlevel")
ml.add_stressmodel(swm)

# ml.solve(tmax='2023-12-06')
ml.solve(
    # tmin=tmin
    )
results = ml.plots.results()

params = ml.parameters.copy()
# get output series
outputs = ml.get_output_series()

# save figure
fig = results[0].get_figure()
fig.savefig(f"pastas/output/ModelC.png", dpi=300)


# MODEL Da (RECHARGE) ############################
"""
1. Run model Da to get initial recharg parameters (only fitting to late time series)
2. Run model Db to get pumping parameters using fixed recharge parameters
3. Run optimization to find RMSE when pumping is likley occuring.
"""

well = rain.copy()
well.name = 'well'
well[:] = 0
start = pd.Timestamp('2023-11-16')
end = start + pd.Timedelta(days=30)
well[start:end] = 1

tmin = '2024-01-15'
tmax = head.index[-1].strftime('%Y-%m-%d')
ml = ps.Model(head)
ml.add_noisemodel(ps.ArNoiseModel())

# add recharge model
evap = rain1d.copy()
evap.name = 'evap'
evap[:] = 0
rch = ps.rch.FlexModel()
rm = ps.RechargeModel(rain1d, evap, recharge=rch, rfunc=ps.Gamma(), name="recharge")
ml.add_stressmodel(rm)

# add pumping model
stm = ps.StressModel(well, rfunc=ps.Hantush(), name='pumping', up=False, settings='prec')
ml.add_stressmodel(stm)

# ml.solve(tmax='2023-12-06')
ml.solve(
    tmin=tmin
)
results = ml.plots.results()

# get recharge
rch_sim = ml.get_stress("recharge")


# get output series
outputs = ml.get_output_series()

# save figure
fig = results[0].get_figure()
fig.savefig(f"pastas/output/ModelDa.png", dpi=300)

# check against whole series
params = ml.parameters.copy()

# rainfall
recharge = ml.get_stress('recharge')
# recharge = ml.get_contributions()[0]

# recharge.sum()/rainfall.sum() # 2.2% recharge
# recharge.sum()/137 # average 0.035 mm/d over entire period
#########################################################
"""
2. Optimize pumping days
"""
#########################################################
import numpy as np
start_dates = pd.date_range('2023-11-05', '2023-12-01')
durations = np.arange(30, 91, 5)  # 30 to 90 days, step 5

best_rmse = np.inf
best_start = None
best_duration = None

for start in start_dates:
    prev_rmse = np.inf
    for duration in durations:
        end = start + pd.Timedelta(days=duration)
        well[:] = 0
        well[start:end] = 1

        ml = ps.Model(head)
        ml.add_stressmodel(rm)
        stm = ps.StressModel(well, rfunc=ps.Hantush(), name='pumping', up=False, settings='prec')
        ml.add_stressmodel(stm)
        for p, val in params['optimal'].items():
            if p in ['recharge_A', 'recharge_n', 'recharge_a', 'recharge_f']:
                ml.set_parameter(p, val, vary=False)
        ml.solve()
        res = ml.residuals()
        rmse = np.sqrt(np.mean(res**2))
        print(f"Start: {start.date()}, Duration: {duration} days, RMSE: {rmse:.3f}")

        if rmse < best_rmse:
            best_rmse = rmse
            best_start = start
            best_duration = duration

        # Early stopping if RMSE increases
        if rmse > prev_rmse:
            break
        prev_rmse = rmse

print(f"Best start: {best_start}, duration: {best_duration} days, RMSE: {best_rmse:.3f}")

#########################################################
"""
3. Now run with optimal parameters as starting values
"""
#########################################################

well[:] = 0
start = pd.Timestamp('2023-11-08')
end = start + pd.Timedelta(days=35)
well[start:end] = 1

ml = ps.Model(head)
ml.add_stressmodel(rm)

stm = ps.StressModel(well, rfunc=ps.Hantush(), name='pumping', up=False, settings='prec')
ml.add_stressmodel(stm)

# # fix recharge parameters
# for p, val in params['optimal'].items():
#     if p in ['recharge_A', 'recharge_n', 'recharge_a', 'recharge_srmax', 'recharge_lp', 'recharge_ks', 'recharge_gamma', 'recharge_kv', 'recharge_simax']:
#         ml.set_parameter(p, val, vary=False)

ml.solve()
ml.plots.results()

print('here')


#########################################################
"""
Model E
"""
#########################################################

well[:] = 0
start = pd.Timestamp('2023-11-08')
end = start + pd.Timedelta(days=35)
well[start:end] = 1

ml = ps.Model(head)
rch = ps.rch.FlexModel()
rm = ps.RechargeModel(rain1d, evap, recharge=rch, rfunc=ps.Gamma(), name="recharge")
ml.add_stressmodel(rm)

stm = ps.StressModel(well, rfunc=ps.Hantush(), name='pumping', up=False, settings='prec')
ml.add_stressmodel(stm)
ml.solve()
ml.plots.results()

recharge = ml.get_stress('recharge')
print(recharge.sum()) # precipitation over period is 286.4 mm
# Berendrecht: 252.64 mm
# Peterson: 270.0 mm
# Flex Model: 291.6 mm

