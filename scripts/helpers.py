import pandas as pd
import numpy as np
from pathlib import Path
from setup import *


def extract_budget(model_name=None):
    import flopy
    import numpy as np
    
    def get_gwf(model_name, gwf=None):
        if gwf is None:
            import flopy

            sim = flopy.mf6.MFSimulation.load(sim_ws='.')
            gwf = sim.get_model(model_name)
        return sim, gwf
    
    # _, gwf = get_gwf(model_name)
    # gwf.ghb[0].name

    lst_path = f"{model_name}.lst"
    lst = flopy.utils.Mf6ListBudget(lst_path)
    incremental, cumulative = lst.get_dataframes(diff=True, start_datetime=None)
    # inc.columns = inc.columns.map(lambda x: x.lower().replace("_","-"))
    col_names = {
        'ghb': 'awanui',
        'ghb2': 'confined',
        'ghb3': 'poukawa',
        'ghb4': 'spring',
        'rcha': 'recharge',
        'sto-ss': 'stoss',
        'total': 'total',
        'wel': 'mbr',
        'wel2': 'inflow',
        'wel3': 'outflow',
        'in-out': 'inout'
    }
    # replace any _ in column names
    incremental.columns = incremental.columns.map(lambda x: x.lower().replace("_",""))
    incremental.rename(columns=col_names, inplace=True)
    incremental.index = np.arange(1, len(incremental) + 1)
    incremental.index.name = "kper"

    incremental['sw'] = incremental['awanui'] + incremental['poukawa'] + incremental['spring']
    incremental['inflow'] = incremental['recharge'] + incremental['confined']

    # in present day: awanui and spring seperate
    # incremental['awanui-spring'] = incremental['awanui'] + incremental['spring']
    # inc.to_csv(f"{MODEL_DIR}/inc.csv")
    incremental.to_csv(f"output.budget.csv")
    return incremental


# def extract_spring_obs(
#         gwf=None, model_name=None, 
#         samples_path=None, 
#         conf_h_path=None,
#         aw_h_path=None,
#         pw_h_path=None,
#         spring_h_path=None, 
#         pw_q_path=None
#         ):
#     import geopandas as gpd
#     import pandas as pd

#     def extract_sample_heads(
#             sample_locations, 
#             gwf=None, 
#             kstpkper=None, 
#             conf_ts=None, 
#             aw_ts=None, 
#             pw_ts=None, 
#             spring_ts=None,
#             ):
#         import pandas as pd

#         xy_spring = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['spring']['x'],
#                 y=sample_locations.loc['spring']['y'],
#                 local=False,
#                 forgive=True)
#         xy_aw = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['aw']['x'],
#                 y=sample_locations.loc['aw']['y'],
#                 local=False,
#                 forgive=True)
#         xy_pw = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['pw']['x'],
#                 y=sample_locations.loc['pw']['y'],
#                 local=False,
#                 forgive=True)
#         xyz_pk4 = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['pk4']['x'],
#                 y=sample_locations.loc['pk4']['y'],
#                 z=sample_locations.loc['pk4']['z'],
#                 local=False,
#                 forgive=True)
#         times = gwf.output.head().get_times()

#         pdata = {}
#         for i in range(len(kstpkper)):
#             n = kstpkper[i]
#             kper = n[1]
#             kstp = n[0]
#             heads = gwf.output.head().get_data(kstpkper=n)
#             time = times[i]
#             pk4_h = heads[xyz_pk4[0], xyz_pk4[1], xyz_pk4[2]]
#             aw_gw_h = heads[0, xy_aw[0], xy_aw[1]]
#             pw_gw_h = heads[0, xy_pw[0], xy_pw[1]]

            
#             # ghb[0] = awanui
#             # ghb[1] = conf
#             # ghb[2] = poukawa
#             # ghb[3] = spring
#             if kper == 0:
#                 # index [kper][1st item][index, val, ...]
#                 conf_h = gwf.ghb[1].stress_period_data.get_data()[kper][0][1]
#                 aw_sw_h = gwf.ghb[0].stress_period_data.get_data()[kper][0][1]
#                 pw_sw_h = gwf.ghb[2].stress_period_data.get_data()[kper][0][1]
                
#                 spring_df = pd.DataFrame(gwf.ghb[3].stress_period_data.get_data()[kper])
#                 mask = spring_df['cellid'] == (0, xy_spring[0], xy_spring[1])
#                 spring_h = spring_df[mask]['bhead'].values[0]

#             elif kper == 1:
#                 conf_h = conf_ts.loc[int(time)].values[0]
#                 spring_h = spring_ts.loc[int(time)][-1]
#                 aw_sw_h = aw_ts.loc[int(time)][-1]
#                 pw_sw_h = pw_ts.loc[int(time)][-1]

#             elif kper in [2, 3]:
#                 # index [kper][1st item][index, val, ...]
#                 conf_h = gwf.ghb[1].stress_period_data.get_data()[kper][0][1]
#                 aw_sw_h = gwf.ghb[0].stress_period_data.get_data()[kper][0][1]
#                 pw_sw_h = gwf.ghb[2].stress_period_data.get_data()[kper][0][1]
#                 # spring head from awanui ghb
#                 spring_df = pd.DataFrame(gwf.ghb[0].stress_period_data.get_data()[kper])
#                 spring_h = spring_df[spring_df['cellid'] == (0, xy_spring[0], xy_spring[1])]['bhead'].values[0]
        
#             conf_diff = pk4_h - conf_h
#             spring_diff = pk4_h - spring_h
#             aw_gw_diff = pk4_h - aw_gw_h
#             pw_gw_diff = pk4_h - pw_gw_h
#             aw_sw_diff = pk4_h - aw_sw_h
#             pw_sw_diff = pk4_h - pw_sw_h
#             aw_swgw_diff = aw_sw_h - aw_gw_h
#             pw_swgw_diff = pw_sw_h - pw_gw_h

#             pdata[n] = [
#                 time, 
#                 pk4_h, 
#                 conf_h, 
#                 aw_sw_h, pw_sw_h, 
#                 spring_h, 
#                 aw_gw_h, pw_gw_h, 
#                 conf_diff, 
#                 spring_diff,  
#                 aw_sw_diff, pw_sw_diff, 
#                 aw_gw_diff, pw_gw_diff,
#                 aw_swgw_diff, pw_swgw_diff
#                 ]
#         # create a DataFrame with the results
#         cols = [
#             'time', 
#             'pk4head', 
#             'confhead', 
#             'awswhead', 'pwswhead', 
#             'springhead', 
#             'awgwhead', 'pwgwhead',
#             'pk4confdiff', 
#             'pk4springdiff', 
#             'pk4awswdiff', 'pk4pwswdiff',
#             'pk4awgwdiff', 'pk4pwgwdiff',
#             'awswgwdiff', 'pwswgwdiff'
#         ]
#         head_results = pd.DataFrame.from_dict(pdata, orient='index', columns=cols)
#         head_results['kstp'] = [i[0]+1 for i in head_results.index]
#         head_results['kper'] = [i[1]+1 for i in head_results.index]

#         head_results.reset_index(drop=True, inplace=True)

#         time_start1 = int(head_results.loc[(head_results['kstp'] == 1) & (head_results['kper'] == 2)]['time'].values[0]-1)
#         time_stop1 = int(head_results.loc[(head_results['kstp'] == head_results['kstp'].max()) & (head_results['kper'] == 2)]['time'].values[0]+1)
#         time_start2 = int(head_results.loc[(head_results['kstp'] == 1) & (head_results['kper'] == 4)]['time'].values[0]-1)

#         head_results['firstpostday'] = -1
#         negs = head_results.loc[time_start1:time_stop1, 'pk4springdiff'].lt(0)
#         if negs.any():
#             head_results.loc[time_start1:time_stop1, 'firstpostday'] = head_results.loc[time_start1:time_stop1, 'pk4springdiff'].lt(0).idxmax()
#         else:
#             head_results.loc[time_start1:time_stop1, 'firstpostday'] = -1
        
#         negs = (time_start2) - head_results.loc[time_start2:, 'pk4springdiff'].lt(0)
#         if negs.any():
#             head_results.loc[time_start2:, 'firstpostday'] = (time_start2) - head_results.loc[time_start2:, 'pk4springdiff'].lt(0).idxmax()
#         else:
#             head_results.loc[time_start2:, 'firstpostday'] = -1

#         return head_results[['time', 'kper', 'kstp'] + cols + ['firstpostday']].reset_index(drop=True)
        
#     def extract_spring_flux(sample_locations, gwf=None, kstpkper=None):
#         import geopandas as gpd
#         import pandas as pd

#         xy_spring = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['spring']['x'],
#                 y=sample_locations.loc['spring']['y'],
#                 local=False,
#                 forgive=True)
        
#         # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
#         cbb = gwf.output.budget()

#         dfs = []
#         for idx in range(cbb.get_nrecords()):
#             rec = cbb.recordarray[idx]
#             package_type = rec['text'].strip().decode()
#             package_name = rec['paknam2'].strip().decode()
#             kper = rec['kper']
#             kstp = rec['kstp']
#             data = cbb.get_record(idx)
#             if package_name in ['GHB_SP', 'GHB_AW']:
#                 df = pd.DataFrame(data)
#                 kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
#                 df['index'] = kij
#                 df['kper'] = kper
#                 df['kstp'] = kstp
#                 df['package'] = package_name
#                 df.rename(columns={'q': 'springcum'}, inplace=True)
#                 if (kper in [1, 2]) and package_name == 'GHB_SP':
#                     dfs.append(df[['package', 'kper', 'kstp', 'index', 'springcum']])
#                 elif (kper in [3, 4]) and package_name == 'GHB_AW':
#                     dfs.append(df[['package', 'kper', 'kstp', 'index', 'springcum']])
        
#         # combine all dataframes
#         ghb_df = pd.concat(dfs, ignore_index=True)
#         mask = ghb_df['package'] == 'GHB_SP'
        
#         spring_cum = ghb_df[mask][['kper', 'kstp', 'springcum']].groupby(['kper', 'kstp']).sum().reset_index()
#         spring_q = ghb_df[ghb_df['index'] == (0, xy_spring[0], xy_spring[1])]
#         spring_q.rename(columns={'springcum': 'springq'}, inplace=True)
        
#         # merge spring flux and spring cumulative
#         merged = pd.merge(spring_q, spring_cum, on=['kper', 'kstp'], how='left').fillna(0)
#         return merged[['kper', 'kstp', 'springcum', 'springq']]
        
    
    
#     def extract_pw_flux(sample_locations, gwf=None, kstpkper=None, pw_ts=None):
#         import numpy as np
#         import pandas as pd

#         xy_pw = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['pw']['x'],
#                 y=sample_locations.loc['pw']['y'],
#                 local=False,
#                 forgive=True)
        
#         # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
#         cbb = gwf.output.budget()
#         times = cbb.get_times()

#         dfs = []
#         for idx in range(cbb.get_nrecords()):
#             rec = cbb.recordarray[idx]
#             package_type = rec['text'].strip().decode()
#             package_name = rec['paknam2'].strip().decode()
#             kper = rec['kper']
#             kstp = rec['kstp']
#             data = cbb.get_record(idx)
#             if package_name in ['GHB_PW']:
#                 df = pd.DataFrame(data)
#                 kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
#                 df['index'] = kij
#                 df['kper'] = kper
#                 df['kstp'] = kstp
#                 df['package'] = package_name
#                 df.rename(columns={'q': 'pwcum'}, inplace=True)
#                 dfs.append(df[['package', 'kper', 'kstp', 'index', 'pwcum']])
        
#         pw_kper1 = []
#         for t in times:
#             if (t >= 2) & (t < 53):
#                 pw_q = pw_ts.loc[t,:]['sm_m3d']
#                 pw_kper1.append(pw_q)
#         pw_kper1 = pd.Series(pw_kper1) * -1
#         pw_kper1.index = np.arange(1, len(pw_kper1) + 1)
#         pw_kper1.name = 'pwinflow'

#         # combine all dataframes
#         ghb_df = pd.concat(dfs, ignore_index=True)

#         pw_cum = ghb_df[['kper', 'kstp', 'pwcum']].groupby(['kper', 'kstp']).sum().reset_index()
#         pw_cum = pd.merge(pw_cum, pw_kper1, left_index=True, right_index=True, how='left')
#         # get first value from series and backfill
#         pw_cum.fillna(pw_ts.iloc[0]['sm_m3d'] * -1, inplace=True)

#         # real pw (inflow + actual local flow)
#         pw_cum['pwtot'] = pw_cum['pwcum'] + pw_cum['pwinflow']

#         pw_q = ghb_df[ghb_df['index'] == (0, xy_pw[0], xy_pw[1])]
#         pw_q.rename(columns={'pwcum': 'pwq'}, inplace=True)
#         merged = pd.merge(pw_cum, pw_q, on=['kper', 'kstp'], how='left').fillna(0)
#         return merged[['kper', 'kstp', 'pwtot', 'pwcum','pwq']]
    
#     def extract_aw_flux(sample_locations, gwf=None, kstpkper=None):
#         import geopandas as gpd
#         import pandas as pd

#         xy_aw = gwf.modelgrid.intersect(
#                 x=sample_locations.loc['aw']['x'],
#                 y=sample_locations.loc['aw']['y'],
#                 local=False,
#                 forgive=True)
        
#         # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
#         cbb = gwf.output.budget()

#         dfs = []
#         for idx in range(cbb.get_nrecords()):
#             rec = cbb.recordarray[idx]
#             package_type = rec['text'].strip().decode()
#             package_name = rec['paknam2'].strip().decode()
#             kper = rec['kper']
#             kstp = rec['kstp']
#             data = cbb.get_record(idx)
#             if package_name in ['GHB_AW']:
#                 df = pd.DataFrame(data)
#                 kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
#                 df['index'] = kij
#                 df['kper'] = kper
#                 df['kstp'] = kstp
#                 df['package'] = package_name
#                 df.rename(columns={'q': 'awcum'}, inplace=True)
#                 dfs.append(df[['package', 'kper', 'kstp', 'index', 'awcum']])
        
#         # combine all dataframes
#         ghb_df = pd.concat(dfs, ignore_index=True)
#         aw_cum = ghb_df[['kper', 'kstp', 'awcum']].groupby(['kper', 'kstp']).sum().reset_index()
#         aw_q = ghb_df[ghb_df['index'] == (0, xy_aw[0], xy_aw[1])]
#         aw_q.rename(columns={'awcum': 'awq'}, inplace=True)
#         merged = pd.merge(aw_cum, aw_q, on=['kper', 'kstp'], how='left').fillna(0)
#         return merged[['kper', 'kstp', 'awcum', 'awq']]
    
#     if gwf is None:
#         import flopy

#         sim = flopy.mf6.MFSimulation.load(sim_ws='.')
#         gwf = sim.get_model(model_name)

#     sample_locations = gpd.read_file(samples_path)
#     sample_locations.set_index('obsnme', inplace=True)

#     # open ts files
#     conf_ts = pd.read_csv(conf_h_path, index_col=0)
#     awh_ts = pd.read_csv(aw_h_path, index_col=0)
#     pwh_ts = pd.read_csv(pw_h_path, index_col=0)
#     spring_ts = pd.read_csv(spring_h_path, index_col=0)
#     pwq_ts = pd.read_csv(pw_q_path, index_col=0)

#     kstpkper = gwf.output.head().get_kstpkper()

#     obs_results = extract_sample_heads(sample_locations, gwf, kstpkper, conf_ts, awh_ts, pwh_ts, spring_ts)
#     spring_flux = extract_spring_flux(sample_locations, gwf, kstpkper)
#     pw_flux = extract_pw_flux(sample_locations, gwf, kstpkper, pwq_ts)
#     aw_flux = extract_aw_flux(sample_locations, gwf, kstpkper)

#     merged = pd.merge(obs_results, spring_flux, on=['kper', 'kstp'])
#     merged = pd.merge(merged, pw_flux, on=['kper', 'kstp'], how='left')
#     merged = pd.merge(merged, aw_flux, on=['kper', 'kstp'], how='left')
#     merged['awtot'] = merged['awcum'] + merged['pwtot'] + merged['springcum']
#     merged.fillna(0, inplace=True)
#     # save the results to a CSV file
#     merged.to_csv(f"obs_results.csv", index=True)
#     return


def extract_ghb_fluxes(model_name, gwf=None, save=True, sample_path=None):
    from functools import reduce

    def get_gwf(gwf, model_name):
        if gwf is None:
            import flopy

            sim = flopy.mf6.MFSimulation.load(sim_ws='.')
            gwf = sim.get_model(model_name)
        return sim, gwf

    def get_sample_xy(gwf, samples_path):
        import geopandas as gpd

        sample_locations = gpd.read_file(samples_path)
        sample_locations.set_index('obsnme', inplace=True)
        
        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        xy_aw = gwf.modelgrid.intersect(
                x=sample_locations.loc['aw']['x'],
                y=sample_locations.loc['aw']['y'],
                local=False,
                forgive=True)
        xy_pw = gwf.modelgrid.intersect(
                x=sample_locations.loc['pw']['x'],
                y=sample_locations.loc['pw']['y'],
                local=False,
                forgive=True)
        xyz_pk4 = gwf.modelgrid.intersect(
                x=sample_locations.loc['pk4']['x'],
                y=sample_locations.loc['pk4']['y'],
                z=sample_locations.loc['pk4']['z'],
                local=False,
                forgive=True)

        samples = {
            'spring': (0, int(xy_spring[0]), int(xy_spring[1])),
            'awanui': (0, int(xy_aw[0]), int(xy_aw[1])),
            'poukawa': (0, int(xy_pw[0]), int(xy_pw[1])),
            'pk4': (int(xyz_pk4[0]), int(xyz_pk4[1]), int(xyz_pk4[2]))
        }

        return samples

    def df_to_array(df, value_col, idomain, fill_value=np.nan):
        """
        Convert dataframe with index tuples to 3D array
        
        Parameters:
        - df: dataframe with index column containing (layer, row, col) tuples
        - value_col: column name to use as values (e.g., 'AWq', 'CONFq')
        - fill_value: value to use for missing locations
        """
        # Extract layer, row, col from index tuples
        # coords = df['index'].apply(lambda x: eval(x) if isinstance(x, str) else x)
        layers = list(df['k'])
        rows = list(df['i'])
        cols = list(df['j'])
        
        # Create empty array
        arr = np.full(idomain.shape, fill_value)
        
        # Fill array with values
        for i, (lay, row, col) in enumerate(zip(layers, rows, cols)):
            if not pd.isna(df[value_col].iloc[i]):
                arr[lay, row, col] = df[value_col].iloc[i]
        
        return arr

    _, gwf = get_gwf(gwf, model_name)
    cbb = gwf.output.budget()
    idomain = gwf.dis.idomain.get_data()
    idomain = np.where(idomain > 0, 1, 0)

    dfs = {
        'GHB_AW': [],
        'GHB_CONF': [],
        'GHB_PW': [],
        'GHB_SP': []
    }
    arrs = {
        'GHB_AW': {},
        'GHB_CONF': {},
        'GHB_PW': {},
        'GHB_SP': {},
    }
    samples_ts = {
        'GHB_SP': []
    }

    if sample_path:
        all_samples = get_sample_xy(gwf, sample_path)
        samples = {name: loc for name, loc in all_samples.items() if name =='spring'}
    else:
        samples = None

    for idx in range(cbb.get_nrecords()):
        rec = cbb.recordarray[idx]
        package_type = rec['text'].strip().decode()
        package_name = rec['paknam2'].strip().decode()
        kper = rec['kper']
        kstp = rec['kstp']
        time = rec['totim']
        data = cbb.get_record(idx)

        # if GHB
        if package_type == 'GHB':
            df = pd.DataFrame(data)
            kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
            # df['index'] = kij
            df['kper'] = kper
            df['kstp'] = kstp
            df['time'] = time
            df['k'] = kij.apply(lambda x: x[0])
            df['i'] = kij.apply(lambda x: x[1])
            df['j'] = kij.apply(lambda x: x[2])
            df['package'] = package_name

            # reformat column name
            suffix = package_name.split('_')[-1]
            df.rename(columns={'q': f'{suffix}q'}, inplace=True)
            dfs[package_name].append(df[['kper', 'kstp', 'time', 'k', 'i', 'j', f'{suffix}q']])

            if isinstance(samples, dict) & (package_name == 'GHB_SP'):
                arr = df_to_array(df, f'{suffix}q', idomain)
                sample_fluxes = {}
                for name, (lay, row, col) in samples.items():
                    sample_fluxes[name] = arr[lay, row, col]
                sample_fluxes['kper'] = kper
                sample_fluxes['kstp'] = kstp

                samples_ts[package_name].append(pd.DataFrame(sample_fluxes, index=[time]))
    
    
    if samples:
        sample_dfs = []
        for pkg, df_list in samples_ts.items():
            pkg_df = pd.concat(df_list, ignore_index=False)
            pkg_df.reset_index(inplace=True)
            pkg_df.rename(columns={'index': 'time'}, inplace=True)
            pkg_df = pkg_df.melt(
                id_vars=['time', 'kper', 'kstp'],
                value_vars=list(samples.keys()), 
                var_name='sample', 
                value_name='flux')
            pkg_df['package'] = pkg
            sample_dfs.append(pkg_df.dropna())
        # combine all
        sample_results = pd.concat(sample_dfs)

        # fluxes = sample_results.groupby('time').sum()
        # save
        if save:
            sample_results.to_csv(f"output.spring_fluxes.csv")

    # combine all dataframes
    combined_dfs = {}
    for pkg, df_list in dfs.items():
        combined_dfs[pkg] = pd.concat(df_list, ignore_index=True)

        if save:
            combined_dfs[pkg].to_csv(f"output.{pkg}_fluxes.csv", index=False)
    
    return combined_dfs


def to_structured(data_dict):
    # Get the shape from the first array (assuming all arrays have same shape)
    first_arr = next(iter(data_dict.values()))
    shape = first_arr.shape
    
    # Create dtype with field names as keys
    dtype = [(key, arr.dtype) for key, arr in data_dict.items()]
    
    # Create structured array with correct shape
    structured_arr = np.zeros(shape, dtype=dtype)
    
    # Fill each field
    for key, arr in data_dict.items():
        structured_arr[key] = arr
    
    return structured_arr


def get_sample_xy(gwf, samples_path):
    import geopandas as gpd

    sample_locations = gpd.read_file(samples_path)
    sample_locations.set_index('obsnme', inplace=True)
    
    xy_spring = gwf.modelgrid.intersect(
            x=sample_locations.loc['spring']['x'],
            y=sample_locations.loc['spring']['y'],
            local=False,
            forgive=True)
    xy_aw = gwf.modelgrid.intersect(
            x=sample_locations.loc['aw']['x'],
            y=sample_locations.loc['aw']['y'],
            local=False,
            forgive=True)
    xy_pw = gwf.modelgrid.intersect(
            x=sample_locations.loc['pw']['x'],
            y=sample_locations.loc['pw']['y'],
            local=False,
            forgive=True)
    xyz_pk4 = gwf.modelgrid.intersect(
            x=sample_locations.loc['pk4']['x'],
            y=sample_locations.loc['pk4']['y'],
            z=sample_locations.loc['pk4']['z'],
            local=False,
            forgive=True)

    samples = {
        'spring': (0, int(xy_spring[0]), int(xy_spring[1])),
        'awanui': (0, int(xy_aw[0]), int(xy_aw[1])),
        'poukawa': (0, int(xy_pw[0]), int(xy_pw[1])),
        'pk4': (int(xyz_pk4[0]), int(xyz_pk4[1]), int(xyz_pk4[2]))
    }

    return samples


def df_to_array(df, value_col, idomain, fill_value=np.nan):
    """
    Convert dataframe with index tuples to 3D array
    
    Parameters:
    - df: dataframe with index column containing (layer, row, col) tuples
    - value_col: column name to use as values (e.g., 'AWq', 'CONFq')
    - fill_value: value to use for missing locations
    """
    # Extract layer, row, col from index tuples
    # coords = df['index'].apply(lambda x: eval(x) if isinstance(x, str) else x)
    layers = list(df['k'])
    rows = list(df['i'])
    cols = list(df['j'])
    
    # Create empty array
    arr = np.full(idomain.shape, fill_value)
    
    # Fill array with values
    for i, (lay, row, col) in enumerate(zip(layers, rows, cols)):
        if not pd.isna(df[value_col].iloc[i]):
            arr[lay, row, col] = df[value_col].iloc[i]
    
    return arr

def get_gwf(gwf, model_name):
    if gwf is None:
        import flopy

        sim = flopy.mf6.MFSimulation.load(sim_ws='.')
        gwf = sim.get_model(model_name)
    return sim, gwf

def extract_model_heads(model_name, gwf=None, sample_path=None):
    import pandas as pd

    def get_gwf(gwf, model_name):
        if gwf is None:
            import flopy

            sim = flopy.mf6.MFSimulation.load(sim_ws='.')
            gwf = sim.get_model(model_name)
        return sim, gwf    

    def get_sample_xy(gwf, samples_path):
        import geopandas as gpd

        sample_locations = gpd.read_file(samples_path)
        sample_locations.set_index('obsnme', inplace=True)
        
        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        xy_aw = gwf.modelgrid.intersect(
                x=sample_locations.loc['aw']['x'],
                y=sample_locations.loc['aw']['y'],
                local=False,
                forgive=True)
        xy_pw = gwf.modelgrid.intersect(
                x=sample_locations.loc['pw']['x'],
                y=sample_locations.loc['pw']['y'],
                local=False,
                forgive=True)
        xyz_pk4 = gwf.modelgrid.intersect(
                x=sample_locations.loc['pk4']['x'],
                y=sample_locations.loc['pk4']['y'],
                z=sample_locations.loc['pk4']['z'],
                local=False,
                forgive=True)

        samples = {
            'spring': (0, int(xy_spring[0]), int(xy_spring[1])),
            'awanui': (0, int(xy_aw[0]), int(xy_aw[1])),
            'poukawa': (0, int(xy_pw[0]), int(xy_pw[1])),
            'pk4': (int(xyz_pk4[0]), int(xyz_pk4[1]), int(xyz_pk4[2]))
        }

        return samples

    def compare_ghb_heads(gwf, i, ghb_col, all_samples):
        ghb_h = pd.DataFrame(gwf.ghb[i].ts.timeseries.get_data())
        ghb_h.set_index('ts_time', inplace=True)

        df = pd.merge(all_samples, ghb_h, how='left', right_index=True, left_on='time')
        return df[ghb_col] - df['pk4']

    _, gwf = get_gwf(gwf, model_name)
    kstpkper = gwf.output.head().get_kstpkper()

    if sample_path:
        samples = get_sample_xy(gwf, sample_path)
    else:
        samples = None
    
    samples_ts = {}
    arrs = {}
    for i in range(len(kstpkper)):
        n = kstpkper[i]
        kper = n[1]
        kstp = n[0]
        heads = gwf.output.head().get_data(kstpkper=n)
        time = gwf.output.head().get_times()[i]
        
        arrs[f'{kper}.{kstp}'] = heads
        for k in range(heads.shape[0]):
            # save each layer separately
            np.savetxt(f"output.heads_kper{kper+1}_kstp{kstp+1}_time{time}.dat", heads[k])
        
        if samples:
            sample_heads = {}
            for name, (lay, row, col) in samples.items():
                sample_heads[name] = heads[lay, row, col]
            sample_heads['kper'] = kper + 1  # make 1-based
            sample_heads['kstp'] = kstp + 1  # make 1-based

            samples_ts[f'{kper}.{kstp}'] = pd.DataFrame(sample_heads, index=[time])
    
    if samples:
        # combine dfs
        all_samples = pd.concat(samples_ts.values(), ignore_index=False)
        all_samples = all_samples.reset_index().rename(columns={'index': 'time'})

        # add derived columns
        all_samples['pk4-spr-diff'] = compare_ghb_heads(gwf, 3, 'ts_array_11', all_samples)


        all_samples['pk4-conf-diff'] = compare_ghb_heads(gwf, 1, 'ts_array_0', all_samples)
        # all_samples['d-losing'] = days_losing_regime(sim, all_samples)

        all_samples.fillna(-999, inplace=True)

        recession_rates = all_samples.groupby(['kper']).last() - all_samples.groupby(['kper']).first()
        
        # get days to regime switch
        recession_rates['fliptime'] = -1
        for kper in recession_rates.index:
            mask = (all_samples['kper'] == kper) & (all_samples['pk4-spr-diff'] > 0)
            i_days = all_samples[mask].index.to_list()
            if len(i_days) > 0:
                day = i_days[0]
            else:
                day = -1
            recession_rates.loc[kper, 'fliptime'] = day
        
        recession_rates.to_csv(f"output.sample_recession_rates.csv")

        all_samples.to_csv(f"output.sample_heads.csv", index=False)
        return all_samples


def compare_ghb_heads(gwf, i, ghb_col, all_samples):
    ghb_h = pd.DataFrame(gwf.ghb[i].ts.timeseries.get_data())
    ghb_h.set_index('ts_time', inplace=True)

    df = pd.merge(all_samples, ghb_h, how='left', right_index=True, left_index=True)
    return df[ghb_col] - df['pk4']


# def days_losing_regime(sim, all_samples):
#     df = all_samples.copy()
#     tdis = pd.DataFrame(sim.tdis.perioddata.get_data())

#     end_i = 0
#     for i in range(len(tdis)):
#         perlen = tdis.loc[i, 'perlen']
#         start_i = end_i + 1
#         end_i = start_i + perlen - 1
        
#         if i in [0, 2]:
#             value = [0]
#         else:
#             mask = ~all_samples.loc[start_i:end_i, 'pk4-spr-diff'].lt(0)
#             if mask.all():
#                 value = [0] * perlen

#             min_index = mask.idxmin()




# def extract_model_heads_truth(truth_dir, model_name, gwf=None):
#     import pandas as pd
#     from pathlib import Path

#     if gwf is None:
#         import flopy

#         sim = flopy.mf6.MFSimulation.load(sim_ws='.')
#         gwf = sim.get_model(model_name)

#     kstpkper = gwf.output.head().get_kstpkper()
    
#     head_data = {}
#     for i in range(len(kstpkper)):
#         n = kstpkper[i]
#         kper = n[1]
#         kstp = n[0]
#         heads = gwf.dis.top.get_data()

#         # save all heads
#         heads_df = pd.DataFrame(heads).stack().reset_index()
#         heads_df.columns = ['i', 'j', 'head']
#         heads_df['kper'] = kper + 1 # make kper 1-based
#         heads_df['kstp'] = kstp + 1 # make kstp 1-based
#         head_data[n] =  heads_df
    
#     # merge all dataframes
#     all_heads = pd.concat(head_data.values(), ignore_index=True)

#     # save to csv
#     all_heads.to_csv(Path(truth_dir,f"{model_name}_heads_truth.csv"), index=False)
#     return

def create_awghb_truth(ghb_dfs, dir):
    """
    Replace all values in the kper, kstp dataframes with the truth values.
    Only valid for kper 1 and kstp 1.
    """

    df = ghb_dfs['GHB_AW'].copy()

    kstps = df['kstp'].unique()
    kstp_include = np.arange(min(kstps), max(kstps), 10)  # include every 10th kstp

    # replace all values with GHB_Q
    df['weight'] = 0
    df['std'] = GHB_Qstd
    df['AWq'] = GHB_Q
    # set std only on kper 1 and kstp 1
    # df.loc[(df['kper'] == 1) & (df['kstp'] == 1), 'weight'] = 1/GHB_Qstd
    df['weight'] = np.where(
        (df['kper'] < 3) & (df['kstp'].isin(kstp_include)), 
        1/GHB_Qstd, 0)

    df.to_csv(Path(dir, f"output.GHB_AW_fluxes.truth.csv"), index=False)
    return

def create_sprghb_truth(ghb_dfs, dir):
    """
    Replace all values in the kper, kstp dataframes with the truth values.
    Only valid for kper 1 and kstp 1.
    """

    df = ghb_dfs['GHB_SP'].copy()

    kstps = df['kstp'].unique()
    kstp_include = np.arange(min(kstps), max(kstps), 10)  # include every 10th kstp

    # replace all values with GHB_Q
    df['weight'] = 0
    df['std'] = GHB_Qstd
    df['SPq'] = GHB_SPRING_Q
    # set std only on kper 1 and kstp 1
    # df.loc[(df['kper'] == 1) & (df['kstp'] == 1), 'weight'] = 1/GHB_Qstd
    df['weight'] = np.where(
        (df['kper'] < 2) & (df['kstp'].isin(kstp_include)), 
        1/GHB_Qstd, 0)

    df.to_csv(Path(dir, f"output.GHB_SP_fluxes.truth.csv"), index=False)
    return

def create_confghb_truth(ghb_dfs, dir):
    """
    Replace all values in the kper, kstp dataframes with the truth values.
    Only valid for kper 1 and kstp 1.
    """
    df = ghb_dfs['GHB_CONF'].copy()
    
    kstps = df['kstp'].unique()
    kstp_include = np.arange(min(kstps), max(kstps), 10)  # include every 10th kstp
    
    # replace all values with GHB_Q
    df['std'] = GHB_Qstd
    df['CONFq'] = 0
    df['weight'] = np.where(
        df['kstp'].isin(kstp_include), 
        1/GHB_Qstd, 0)
    # set std only on kper 1 and kstp 1

    df.to_csv(Path(dir, f"output.GHB_CONF_fluxes.truth.csv"), index=False)
    return


def create_head_truth(arr, TRUTHREL_DIR):
    np.savetxt(Path(TRUTHREL_DIR, f"output.heads.truth.dat"), arr + HEAD_offset)

    # std
    std_arr = np.ones_like(arr) * HEAD_std
    np.savetxt(Path(TRUTHREL_DIR, f"output.heads.std.dat"), std_arr)

    # weight
    weight_arr = np.ones_like(arr) * 1/HEAD_std
    np.savetxt(Path(TRUTHREL_DIR, f"output.heads.weight.dat"), weight_arr)

    return

def samples_truth(gwf, TRUTHREL_DIR):
    import pandas as pd

    recession_fn = Path(f"output.sample_recession_rates.csv")
    pk4_path = Path(f'{MODEL_NAME}.pk4_level.csv')
    spr_path = Path(f'{MODEL_NAME}.ts_spring_raw.csv')
    conf_path = Path(f'{MODEL_NAME}.ts_confined_raw.csv')
    
    df = pd.read_csv(pk4_path, index_col=0)
    sp_df = pd.read_csv(spr_path, index_col=0)
    conf_df = pd.read_csv(conf_path, index_col=0)
    rec_df = pd.read_csv(recession_fn, index_col='kper')

    times = list(df.index)
    times_include = list(np.arange(1, max(times), 10)) + [2]  # include every 10th kstp
    
    df.rename(columns={'sm_level_mRL': 'pk4'}, inplace=True)
    df['spring'] = sp_df['sm_level_mRL'].values
    df['conf'] = conf_df['level_mRL'].values

    df['std'] = np.where(~df['pk4'].isna(), PK4_std, 0)
    df['weight'] = np.where(
        (df.index.isin(times_include)) & (~df['pk4'].isna()), 1/PK4_std, 0)

    df['pk4-spr-diff'] = df['spring'] - df['pk4']
    df['pk4-conf-diff'] = df['conf'] - df['pk4']

    df.fillna(0, inplace=True)

    # truth values
    rec_df.loc[2, 'pk4'] = df['pk4'].iloc[-1] - df['pk4'].iloc[1]
    rec_df.loc[2, 'pk4-spr-diff'] = df['pk4-spr-diff'].iloc[-1] - df['pk4-spr-diff'].iloc[1]
    rec_df.loc[2, 'pk4-conf-diff'] = df['pk4-conf-diff'].iloc[-1] - df['pk4-conf-diff'].iloc[1]
    
    # default values
    rec_df['std'] = 0.
    rec_df['weight'] = 0.
    rec_df['fliptime'] = -1.

    # get fliptime
    mask = df['pk4-spr-diff'] > 0
    i_days = df[mask].index.to_list()
    if len(i_days) > 0:
        day = int(np.min(i_days))
    else:
        day = -1
    
    rec_df.loc[2, 'fliptime'] = day
    rec_df.loc[2, 'std'] = PK4_std
    rec_df.loc[2, 'weight'] = 1/PK4_std

    rec_df[['time','kstp', 'pk4', 'pk4-spr-diff', 'pk4-conf-diff', 'fliptime', 'std', 'weight']].to_csv(Path(TRUTHREL_DIR, f"output.sample_recession_rates.truth.csv"))
    
    df = df[['pk4', 'pk4-spr-diff', 'pk4-conf-diff', 'std', 'weight']].reset_index().rename(columns={'index': 'time'})
    df.to_csv(Path(TRUTHREL_DIR, f"output.sample_heads.truth.csv"), index=False)
    
    return

def create_budget_truth(budget, TRUTHREL_DIR):
    budget['confined'] = 0. # must be greater than
    budget['sw'] = budget['recharge'] * -1 # must be less than

    # default values
    budget['conf-std'] = 0.
    budget['conf-weight'] = 0.
    budget['inflow-std'] = 0.
    budget['inflow-weight'] = 0.

    budget['confined'] = 0. # must be greater than
    budget.loc[:, 'conf-std'] = 3000
    budget.loc[:, 'conf-weight'] = (1/budget.loc[:, 'conf-std']).values

    # sum of ghb conf and recharge can't be greater than awanui outflow
    budget.loc[1, 'inflow'] = 73000 # must be less
    budget.loc[2, 'inflow'] = 63000 # must be less than
    budget.loc[1, 'inflow-std'] = budget.loc[1, 'inflow'] * 0.05
    budget.loc[2, 'inflow-std'] = budget.loc[2, 'inflow'] * 0.05
    budget.loc[:, 'inflow-weight'] = (1/budget.loc[:, 'inflow-std']).values
    
    budget.replace(np.inf, 0, inplace=True)
    budget.to_csv(Path(TRUTHREL_DIR, f"output.budget.truth.csv"))