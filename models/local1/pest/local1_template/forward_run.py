import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import pyemu

# function added thru PstFrom.add_py_function()
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
    inc, cum = lst.get_dataframes(diff=True, start_datetime=None)
    # inc.columns = inc.columns.map(lambda x: x.lower().replace("_","-"))
    cum.columns = cum.columns.map(lambda x: x.lower().replace("_", "-"))
    # inc.index.name = "totim"
    cum.index.name = "totim"
    # inc.to_csv(f"{MODEL_DIR}/inc.csv")
    cum.to_csv(f"cum.csv")
    return




# function added thru PstFrom.add_py_function()
def extract_spring_obs(gwf=None, model_name=None, samples_path=None):
    import geopandas as gpd
    import pandas as pd

    def extract_sample_heads(sample_locations, gwf=None):
        import pandas as pd

        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        xyz_pk2 = gwf.modelgrid.intersect(
                x=sample_locations.loc['pk2']['x'],
                y=sample_locations.loc['pk2']['y'],
                z=sample_locations.loc['pk2']['z'],
                local=False,
                forgive=True)
        xyz_pk4 = gwf.modelgrid.intersect(
                x=sample_locations.loc['pk4']['x'],
                y=sample_locations.loc['pk4']['y'],
                z=sample_locations.loc['pk4']['z'],
                local=False,
                forgive=True)
        
        heads = gwf.output.head().get_data()
        spring_head = heads[0, xy_spring[0], xy_spring[1]]
        pk2_head = heads[xyz_pk2[0], xyz_pk2[1], xyz_pk2[2]]
        pk4_head = heads[xyz_pk4[0], xyz_pk4[1], xyz_pk4[2]]
        # create a DataFrame with the results

        head_results = pd.DataFrame({
            'spring_head': [spring_head],
            'pk2_head': [pk2_head],
            'pk4_head': [pk4_head],
            'spring_pk4_diff': [spring_head - pk4_head],
            'spring_pk2_diff': [spring_head - pk2_head],
            'pk4_pk2_diff': [pk4_head - pk2_head]
        }, index=[0])

        return head_results

    def extract_sample_fluxes(sample_locations, gwf=None):
        import geopandas as gpd
        import pandas as pd

        xy_spring = gwf.modelgrid.intersect(
                x=sample_locations.loc['spring']['x'],
                y=sample_locations.loc['spring']['y'],
                local=False,
                forgive=True)
        
        # spring_node = gwf.modelgrid.get_node((0, xy_spring[0], xy_spring[1]))
        drn_q = gwf.output.budget().get_data(text='DRN')[0]
        kij = [gwf.modelgrid.get_lrc(node-1) for node in drn_q['node'].tolist()] # nodes are not zero-indexed need to minus 1
        lays = [i[0][0] for i in kij]
        rows = [j[0][1] for j in kij]
        cols = [i[0][2] for i in kij]
        drn_df = pd.DataFrame(drn_q)
        drn_df['k'] = lays
        drn_df['i'] = rows
        drn_df['j'] = cols
        # to dataframe
        drn_df.set_index(['k', 'i', 'j'], inplace=True)

        try:
            spring_q = drn_df.loc[(0, xy_spring[0], xy_spring[1]), 'q']
        except:
            spring_q = -9999

        return spring_q
    
    if gwf is None:
        import flopy

        sim = flopy.mf6.MFSimulation.load(sim_ws='.')
        gwf = sim.get_model(model_name)

    sample_locations = gpd.read_file(samples_path)
    sample_locations.set_index('obsnme', inplace=True)

    obs_results = extract_sample_heads(sample_locations, gwf)
    spring_flux = extract_sample_fluxes(sample_locations, gwf)

    obs_results['spring_flux'] = spring_flux
    obs_results.index.name = 'kper'

    # save the results to a CSV file
    obs_results.to_csv(f"obs_results.csv", index=True)
    
    return


def main():

    try:
       os.remove(r'cum.csv')
    except Exception as e:
       print(r'error removing tmp file:cum.csv')
    try:
       os.remove(r'obs_results.csv')
    except Exception as e:
       print(r'error removing tmp file:obs_results.csv')
    pyemu.helpers.apply_list_and_array_pars(arr_par_file='mult2model_info.csv',chunk_len=50)
    pyemu.os_utils.run(r'mf6')

    extract_heads_and_budget(model_name='local1')
    extract_spring_obs(gwf=None, model_name='local1', samples_path=r'..\..\..\..\data\sample_locations.shp')

if __name__ == '__main__':
    mp.freeze_support()
    main()

