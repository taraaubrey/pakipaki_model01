from setup import *


def extract_heads_and_budget(model_name=None):
    import flopy
    import numpy as np
    
    hds_path = f"{model_name}.hds"
    hds = flopy.utils.HeadFile(hds_path)
    d = hds.get_data()  # get the head data from the file
    for k, dlay in enumerate(d):
        tmp_hds_fn = f"{model_name}_hdslay{k+1}.txt"
        np.savetxt(tmp_hds_fn, dlay, fmt="%15.6E")

    lst_path = f"{model_name}.lst"
    lst = flopy.utils.Mf6ListBudget(lst_path)
    incremental, cumulative = lst.get_dataframes(diff=True, start_datetime=None)
    # inc.columns = inc.columns.map(lambda x: x.lower().replace("_","-"))
    col_names = {
        'ghb': 'bghbaw',
        'ghb2': 'bghbconf',
        'ghb3': 'bghbpw',
        'ghb4': 'bghbspring',
        'rcha': 'brcha',
        'sto-ss': 'bstoss',
        'total': 'btotal',
        'wel': 'bwelmbr',
        'wel2': 'bwelin',
        'wel3': 'bwelout',
    }
    # replace any _ in column names
    incremental.columns = incremental.columns.map(lambda x: x.lower().replace("_",""))
    incremental.rename(columns=col_names, inplace=True)
    incremental.index = np.arange(1, len(incremental) + 1)
    incremental.index.name = "kper"

    # in present day: awanui and spring seperate
    incremental['bawtotal'] = incremental['bghbaw'] + incremental['bghbspring']
    # inc.to_csv(f"{MODEL_DIR}/inc.csv")
    incremental.to_csv(f"incremental_budget.csv")
    return


def extract_spring_obs(
        gwf=None, model_name=None, 
        samples_path=None, 
        conf_h_path=None,
        aw_h_path=None,
        pw_h_path=None,
        spring_h_path=None, 
        pw_q_path=None
        ):
    import geopandas as gpd
    import pandas as pd

    def extract_sample_heads(
            sample_locations, 
            gwf=None, 
            kstpkper=None, 
            conf_ts=None, 
            aw_ts=None, 
            pw_ts=None, 
            spring_ts=None,
            ):
        import pandas as pd

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
        times = gwf.output.head().get_times()

        pdata = {}
        for i in range(len(kstpkper)):
            n = kstpkper[i]
            kper = n[1]
            kstp = n[0]
            heads = gwf.output.head().get_data(kstpkper=n)
            time = times[i]
            pk4_h = heads[xyz_pk4[0], xyz_pk4[1], xyz_pk4[2]]
            aw_gw_h = heads[0, xy_aw[0], xy_aw[1]]
            pw_gw_h = heads[0, xy_pw[0], xy_pw[1]]

            
            # ghb[0] = awanui
            # ghb[1] = conf
            # ghb[2] = poukawa
            # ghb[3] = spring
            if kper == 0:
                # index [kper][1st item][index, val, ...]
                conf_h = gwf.ghb[1].stress_period_data.get_data()[kper][0][1]
                aw_sw_h = gwf.ghb[0].stress_period_data.get_data()[kper][0][1]
                pw_sw_h = gwf.ghb[2].stress_period_data.get_data()[kper][0][1]
                
                spring_df = pd.DataFrame(gwf.ghb[3].stress_period_data.get_data()[kper])
                mask = spring_df['cellid'] == (0, xy_spring[0], xy_spring[1])
                spring_h = spring_df[mask]['bhead'].values[0]

            elif kper == 1:
                conf_h = conf_ts.loc[int(time)].values[0]
                spring_h = spring_ts.loc[int(time)][-1]
                aw_sw_h = aw_ts.loc[int(time)][-1]
                pw_sw_h = pw_ts.loc[int(time)][-1]

            elif kper in [2, 3]:
                # index [kper][1st item][index, val, ...]
                conf_h = gwf.ghb[1].stress_period_data.get_data()[kper][0][1]
                aw_sw_h = gwf.ghb[0].stress_period_data.get_data()[kper][0][1]
                pw_sw_h = gwf.ghb[2].stress_period_data.get_data()[kper][0][1]
                # spring head from awanui ghb
                spring_df = pd.DataFrame(gwf.ghb[0].stress_period_data.get_data()[kper])
                spring_h = spring_df[spring_df['cellid'] == (0, xy_spring[0], xy_spring[1])]['bhead'].values[0]
        
            conf_diff = pk4_h - conf_h
            spring_diff = pk4_h - spring_h
            aw_gw_diff = pk4_h - aw_gw_h
            pw_gw_diff = pk4_h - pw_gw_h
            aw_sw_diff = pk4_h - aw_sw_h
            pw_sw_diff = pk4_h - pw_sw_h
            aw_swgw_diff = aw_sw_h - aw_gw_h
            pw_swgw_diff = pw_sw_h - pw_gw_h

            pdata[n] = [
                time, 
                pk4_h, 
                conf_h, 
                aw_sw_h, pw_sw_h, 
                spring_h, 
                aw_gw_h, pw_gw_h, 
                conf_diff, 
                spring_diff,  
                aw_sw_diff, pw_sw_diff, 
                aw_gw_diff, pw_gw_diff,
                aw_swgw_diff, pw_swgw_diff
                ]
        # create a DataFrame with the results
        cols = [
            'time', 
            'pk4head', 
            'confhead', 
            'awswhead', 'pwswhead', 
            'springhead', 
            'awgwhead', 'pwgwhead',
            'pk4confdiff', 
            'pk4springdiff', 
            'pk4awswdiff', 'pk4pwswdiff',
            'pk4awgwdiff', 'pk4pwgwdiff',
            'awswgwdiff', 'pwswgwdiff'
        ]
        head_results = pd.DataFrame.from_dict(pdata, orient='index', columns=cols)
        head_results['kstp'] = [i[0]+1 for i in head_results.index]
        head_results['kper'] = [i[1]+1 for i in head_results.index]

        head_results.reset_index(drop=True, inplace=True)

        time_start1 = int(head_results.loc[(head_results['kstp'] == 1) & (head_results['kper'] == 2)]['time'].values[0]-1)
        time_stop1 = int(head_results.loc[(head_results['kstp'] == head_results['kstp'].max()) & (head_results['kper'] == 2)]['time'].values[0]+1)
        time_start2 = int(head_results.loc[(head_results['kstp'] == 1) & (head_results['kper'] == 4)]['time'].values[0]-1)

        head_results['firstnegday'] = -1
        negs = head_results.loc[time_start1:time_stop1, 'pk4springdiff'].lt(0)
        if negs.any():
            head_results.loc[time_start1:time_stop1, 'firstnegday'] = head_results.loc[time_start1:time_stop1, 'pk4springdiff'].lt(0).idxmax()
        else:
            head_results.loc[time_start1:time_stop1, 'firstnegday'] = -1
        
        negs = (time_start2) - head_results.loc[time_start2:, 'pk4springdiff'].lt(0)
        if negs.any():
            head_results.loc[time_start2:, 'firstnegday'] = (time_start2) - head_results.loc[time_start2:, 'pk4springdiff'].lt(0).idxmax()
        else:
            head_results.loc[time_start2:, 'firstnegday'] = -1

        return head_results[['time', 'kper', 'kstp'] + cols + ['firstnegday']].reset_index(drop=True)
        
    def extract_spring_flux(sample_locations, gwf=None, kstpkper=None):
        import geopandas as gpd
        import pandas as pd

        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        
        # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
        cbb = gwf.output.budget()

        dfs = []
        for idx in range(cbb.get_nrecords()):
            rec = cbb.recordarray[idx]
            package_type = rec['text'].strip().decode()
            package_name = rec['paknam2'].strip().decode()
            kper = rec['kper']
            kstp = rec['kstp']
            data = cbb.get_record(idx)
            if package_name in ['GHB_SP', 'GHB_AW']:
                df = pd.DataFrame(data)
                kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
                df['index'] = kij
                df['kper'] = kper
                df['kstp'] = kstp
                df['package'] = package_name
                df.rename(columns={'q': 'springcum'}, inplace=True)
                if (kper in [1, 2]) and package_name == 'GHB_SP':
                    dfs.append(df[['package', 'kper', 'kstp', 'index', 'springcum']])
                elif (kper in [3, 4]) and package_name == 'GHB_AW':
                    dfs.append(df[['package', 'kper', 'kstp', 'index', 'springcum']])
        
        # combine all dataframes
        ghb_df = pd.concat(dfs, ignore_index=True)
        mask = ghb_df['package'] == 'GHB_SP'
        spring_cum = ghb_df[mask][['kper', 'kstp', 'springcum']].groupby(['kper', 'kstp']).sum().reset_index()
        spring_q = ghb_df[ghb_df['index'] == (0, xy_spring[0], xy_spring[1])]
        spring_q.rename(columns={'springcum': 'springq'}, inplace=True)
        merged = pd.merge(spring_q, spring_cum, on=['kper', 'kstp'], how='left').fillna(0)
        return merged[['kper', 'kstp', 'springcum', 'springq']]
        
    def extract_pw_flux(sample_locations, gwf=None, kstpkper=None, pw_ts=None):
        import numpy as np
        import pandas as pd

        xy_pw = gwf.modelgrid.intersect(
                x=sample_locations.loc['pw']['x'],
                y=sample_locations.loc['pw']['y'],
                local=False,
                forgive=True)
        
        # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
        cbb = gwf.output.budget()
        times = cbb.get_times()

        dfs = []
        for idx in range(cbb.get_nrecords()):
            rec = cbb.recordarray[idx]
            package_type = rec['text'].strip().decode()
            package_name = rec['paknam2'].strip().decode()
            kper = rec['kper']
            kstp = rec['kstp']
            data = cbb.get_record(idx)
            if package_name in ['GHB_PW']:
                df = pd.DataFrame(data)
                kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
                df['index'] = kij
                df['kper'] = kper
                df['kstp'] = kstp
                df['package'] = package_name
                df.rename(columns={'q': 'pwcum'}, inplace=True)
                dfs.append(df[['package', 'kper', 'kstp', 'index', 'pwcum']])
        
        pw_kper1 = []
        for t in times:
            if (t >= 2) & (t < 53):
                pw_q = pw_ts.loc[t,:]['sm_m3d']
                pw_kper1.append(pw_q)
        pw_kper1 = pd.Series(pw_kper1) * -1
        pw_kper1.index = np.arange(1, len(pw_kper1) + 1)
        pw_kper1.name = 'pwinflow'

        # combine all dataframes
        ghb_df = pd.concat(dfs, ignore_index=True)

        pw_cum = ghb_df[['kper', 'kstp', 'pwcum']].groupby(['kper', 'kstp']).sum().reset_index()
        pw_cum = pd.merge(pw_cum, pw_kper1, left_index=True, right_index=True, how='left')
        # get first value from series and backfill
        pw_cum.fillna(pw_ts.iloc[0]['sm_m3d'] * -1, inplace=True)

        # real pw (inflow + actual local flow)
        pw_cum['pwtot'] = pw_cum['pwcum'] + pw_cum['pwinflow']

        pw_q = ghb_df[ghb_df['index'] == (0, xy_pw[0], xy_pw[1])]
        pw_q.rename(columns={'pwcum': 'pwq'}, inplace=True)
        merged = pd.merge(pw_cum, pw_q, on=['kper', 'kstp'], how='left').fillna(0)
        return merged[['kper', 'kstp', 'pwtot', 'pwcum','pwq']]
    
    def extract_aw_flux(sample_locations, gwf=None, kstpkper=None):
        import geopandas as gpd
        import pandas as pd

        xy_aw = gwf.modelgrid.intersect(
                x=sample_locations.loc['aw']['x'],
                y=sample_locations.loc['aw']['y'],
                local=False,
                forgive=True)
        
        # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
        cbb = gwf.output.budget()

        dfs = []
        for idx in range(cbb.get_nrecords()):
            rec = cbb.recordarray[idx]
            package_type = rec['text'].strip().decode()
            package_name = rec['paknam2'].strip().decode()
            kper = rec['kper']
            kstp = rec['kstp']
            data = cbb.get_record(idx)
            if package_name in ['GHB_AW']:
                df = pd.DataFrame(data)
                kij = df['node'].apply(lambda x: gwf.modelgrid.get_lrc(x-1)[0])
                df['index'] = kij
                df['kper'] = kper
                df['kstp'] = kstp
                df['package'] = package_name
                df.rename(columns={'q': 'awcum'}, inplace=True)
                dfs.append(df[['package', 'kper', 'kstp', 'index', 'awcum']])
        
        # combine all dataframes
        ghb_df = pd.concat(dfs, ignore_index=True)
        aw_cum = ghb_df[['kper', 'kstp', 'awcum']].groupby(['kper', 'kstp']).sum().reset_index()
        aw_q = ghb_df[ghb_df['index'] == (0, xy_aw[0], xy_aw[1])]
        aw_q.rename(columns={'awcum': 'awq'}, inplace=True)
        merged = pd.merge(aw_cum, aw_q, on=['kper', 'kstp'], how='left').fillna(0)
        return merged[['kper', 'kstp', 'awcum', 'awq']]
    
    if gwf is None:
        import flopy

        sim = flopy.mf6.MFSimulation.load(sim_ws='.')
        gwf = sim.get_model(model_name)

    sample_locations = gpd.read_file(samples_path)
    sample_locations.set_index('obsnme', inplace=True)

    # open ts files
    conf_ts = pd.read_csv(conf_h_path, index_col=0)
    awh_ts = pd.read_csv(aw_h_path, index_col=0)
    pwh_ts = pd.read_csv(pw_h_path, index_col=0)
    spring_ts = pd.read_csv(spring_h_path, index_col=0)
    pwq_ts = pd.read_csv(pw_q_path, index_col=0)

    kstpkper = gwf.output.head().get_kstpkper()

    obs_results = extract_sample_heads(sample_locations, gwf, kstpkper, conf_ts, awh_ts, pwh_ts, spring_ts)
    spring_flux = extract_spring_flux(sample_locations, gwf, kstpkper)
    pw_flux = extract_pw_flux(sample_locations, gwf, kstpkper, pwq_ts)
    aw_flux = extract_aw_flux(sample_locations, gwf, kstpkper)

    merged = pd.merge(obs_results, spring_flux, on=['kper', 'kstp'])
    merged = pd.merge(merged, pw_flux, on=['kper', 'kstp'], how='left')
    merged = pd.merge(merged, aw_flux, on=['kper', 'kstp'], how='left')
    merged['awtot'] = merged['awcum'] + merged['pwtot'] + merged['springcum']
    merged.fillna(0, inplace=True)
    # save the results to a CSV file
    merged.to_csv(f"obs_results.csv", index=True)
    return

def extract_model_heads(model_name, gwf=None):
    import pandas as pd

    if gwf is None:
        import flopy

        sim = flopy.mf6.MFSimulation.load(sim_ws='.')
        gwf = sim.get_model(model_name)

    kstpkper = gwf.output.head().get_kstpkper()
    
    head_data = {}
    for i in range(len(kstpkper)):
        n = kstpkper[i]
        kper = n[1]
        kstp = n[0]
        heads = gwf.output.head().get_data(kstpkper=n)

        # save all heads
        heads_df = pd.DataFrame(heads[0]).stack().reset_index()
        heads_df.columns = ['i', 'j', 'head']
        heads_df['kper'] = kper + 1 # make kper 1-based
        heads_df['kstp'] = kstp + 1 # make kstp 1-based
        head_data[n] =  heads_df
    
    # merge all dataframes
    all_heads = pd.concat(head_data.values(), ignore_index=True)

    # save to csv
    all_heads.to_csv(f"{model_name}_heads.csv", index=False)
    return

