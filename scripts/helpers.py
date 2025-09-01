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
        'ghb': 'ghb_aw',
        'ghb2': 'ghb_conf',
        'ghb3': 'ghb_pw',
        'ghb4': 'ghb_spring',
        'rcha': 'rcha',
        'sto-ss': 'sto_ss',
        'total': 'total',
        'wel': 'wel_mbr',
        'wel2': 'wel_in',
        'wel3': 'wel_out',
    }
    cumulative.rename(columns=col_names, inplace=True)
    # inc.index.name = "totim"
    cumulative.index.name = "time"
    # inc.to_csv(f"{MODEL_DIR}/inc.csv")
    cumulative.to_csv(f"cum.csv")
    return


def extract_spring_obs(gwf=None, model_name=None, samples_path=None, conf_h_path=None, spring_h_path=None):
    import geopandas as gpd
    import pandas as pd

    def extract_sample_heads(sample_locations, model_name, gwf=None, kstpkper=None, conf_ts=None, spring_ts=None):
        import pandas as pd

        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        # xyz_pk2 = gwf.modelgrid.intersect(
        #         x=sample_locations.loc['pk2']['x'],
        #         y=sample_locations.loc['pk2']['y'],
        #         z=sample_locations.loc['pk2']['z'],
        #         local=False,
        #         forgive=True)
        xyz_pk4 = gwf.modelgrid.intersect(
                x=sample_locations.loc['pk4']['x'],
                y=sample_locations.loc['pk4']['y'],
                z=sample_locations.loc['pk4']['z'],
                local=False,
                forgive=True)

        pdata = {}
        for n in kstpkper:
            kper = n[1]
            kstp = n[0]
            heads = gwf.output.head().get_data(kstpkper=n)
            
            pk4_head = heads[xyz_pk4[0], xyz_pk4[1], xyz_pk4[2]]
            if kper == 1: #period
                conf_h = conf_ts.loc[kstp+1].values[0]
                spring_h = spring_ts.loc[kstp+1][-1]
            else:
                conf_h = gwf.ghb[1].stress_period_data.get_data()[kper][kstp][1]
                spring_df = pd.DataFrame(gwf.ghb[3].stress_period_data.get_data()[0])
                spring_h = spring_df[spring_df['cellid'] == (0, xy_spring[0], xy_spring[1])]['bhead'].values[0]
            conf_diff = pk4_head - conf_h
            spring_diff = pk4_head - spring_h

            pdata[n] = [pk4_head, conf_h, spring_h, conf_diff, spring_diff]
        # create a DataFrame with the results
        cols = ['pk4_head', 'conf_head', 'spring_head', 'pk4conf_diff', 'pk4spring_diff']
        head_results = pd.DataFrame.from_dict(pdata, orient='index', columns=cols)
        head_results['kstp'] = [i[0]+1 for i in head_results.index]
        head_results['kper'] = [i[1]+1 for i in head_results.index]

        return head_results[['kper', 'kstp'] + cols].reset_index(drop=True)

    def extract_sample_fluxes(sample_locations, model_name, gwf=None, kstpkper=None):
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
                df.rename(columns={'q': 'spring_q'}, inplace=True)
                dfs.append(df[['package', 'kper', 'kstp', 'index', 'spring_q']])
        
        # combine all dataframes
        ghb_df = pd.concat(dfs, ignore_index=True)
        return ghb_df[ghb_df['index'] == (0, xy_spring[0], xy_spring[1])]
    
    if gwf is None:
        import flopy

        sim = flopy.mf6.MFSimulation.load(sim_ws='.')
        gwf = sim.get_model(model_name)

    sample_locations = gpd.read_file(samples_path)
    sample_locations.set_index('obsnme', inplace=True)

    # open ts files
    conf_ts = pd.read_csv(conf_h_path, index_col=0)
    spring_ts = pd.read_csv(spring_h_path, index_col=0)

    kstpkper = gwf.output.head().get_kstpkper()

    obs_results = extract_sample_heads(sample_locations, model_name, gwf, kstpkper, conf_ts, spring_ts)
    spring_flux = extract_sample_fluxes(sample_locations, model_name, gwf, kstpkper)

    merged = pd.merge(obs_results, spring_flux[['kper', 'kstp', 'spring_q']], on=['kper', 'kstp'])

    # save the results to a CSV file
    merged.to_csv(f"obs_results.csv", index=True)
    return

