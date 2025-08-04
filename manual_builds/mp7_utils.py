# modified from Wes Kitlasten

import os
import numpy as np
# import flopy
import pandas as pd
# import pyemu
#import geopandas as gpd
# from nzmf6.build.src import utils


def samples_to_pars(gwf, obs_data=None, npartob=100, radius=20, nan_height=2, gclass=None,
                    input_factor=1.0, par_type=None, site_group='none'):
    '''tracer: dict with key = tracer name
                obs_file: csv with tracer observations (obsnme,date,nztme,tr,site_id,nztmn,tritium_code,sigtr,boredepth)
                tracer_input_file: filename of csv with tracer concentration with time
                input_factor: latitude correction factor for tracer'''

    par = utils.samples_to_obs(gwf, obs_data=obs_data, gclass = gclass,
                               obs_type=par_type, site_group=site_group)
    # some bore characteristics... maybe better in clean_chem util?
    if 'radius' in obs_data.keys() and 'radius' not in par.columns:
        par.loc[:, 'radius'] = obs_data['radius']
    else:
        par.loc[:, 'radius'] = radius
    if 'npartob' in obs_data.keys() and 'npartob' not in par.columns:
        par.loc[:, 'npartob'] = obs_data['npartob']
    else:
        par.loc[:, 'npartob'] = npartob
    if 'height' in obs_data.keys() and 'height' not in par.columns:
        par.loc[:, 'height'] = obs_data['height']
    else:
        par.loc[:, 'height'] = nan_height
    par['height'].fillna(nan_height, inplace=True)
    par['z'].clip(nan_height, None)
    # absolutely necessary to sort if going between different grids with same pst obs
    par.sort_index(inplace=True)
    if 'input_factor' not in par.columns:
        par.loc[:, 'input_factor'] = input_factor
    par.sort_index(inplace=True)
    return par


def build_kdtree(gwf, mask=None):
    # needs to be inside points_to_partloc for forward_run.py
    # kdtree approx method, but

    from scipy.spatial import KDTree
    print('building KDTree')
    nlist = []
    sidxlist = []
    if mask is None:
        idom = gwf.dis.idomain.array
    else:
        idom = np.where(mask>0,1,0)
    x, y, z = gwf.modelgrid.xyzcellcenters    
    xr = (x*idom).ravel()
    yr = (y*idom).ravel()
    zr = (z*idom).ravel()
    nodes = np.stack([xr, yr, zr], axis=0).T
        
    sj = np.concatenate([np.arange(x.shape[1]) for _ in range(0, x.shape[0])])
    si = np.concatenate([np.ones(x.shape[1]) * _ for _ in range(0, x.shape[0])])
    for i in range(0,idom.shape[0]):
        ir = idom[i].ravel()
        # need aligned array for indexes
        sjr = sj * ir
        sir = si * ir
        skr = i * ir
        sidx = np.stack([sjr, sir, skr], axis=0).T
        sidxlist.append(sidx)
    sidxs = np.concatenate(sidxlist)
    tree = KDTree(nodes, leafsize=50, compact_nodes=True, balanced_tree=False)
    return tree, sidxs

def points_to_rloc(df, gwf):
    bottom = gwf.dis.botm.array.copy()
    top = gwf.dis.top.array.copy()
    idom = gwf.dis.idomain.array.copy()
    delr = gwf.dis.delr.data.copy()
    delc = gwf.dis.delc.data.copy()
    xorigin = gwf.dis.xorigin.data
    yorigin = gwf.dis.yorigin.data
    # thickness as 3D array, no ravel
    t = [(top - bottom[0])]
    for i, b in enumerate(bottom[1:]):
        t.append(bottom[i] - b)
    thick = np.stack(t, axis=0) * idom
    df['i'] = df['i'].astype(int)
    df['j'] = df['j'].astype(int)
    df['k'] = df['k'].astype(int)
    # origin is lower left, positive y is up
    # nztmn increases to the north
    df['fj'] = df['j'].apply(lambda x: xorigin + np.sum(delr[:x]))
    df['fi'] = df['i'].apply(lambda y: yorigin + np.sum(delc) - np.sum(delc[:y+1]))
    df['fk'] = df.loc[:,['k', 'i', 'j']].apply(lambda x: bottom[x[0], x[1], x[2]], axis=1)

    # now local cell position, rloc is in pos x, y, z directions
    df['rj'] = (df['x'] - df['fj']) / delr[df['j']]
    df['ri'] = (df['y'] - df['fi']) / delc[df['i']]
    df['rk'] = (df['elev'] - df['fk']).div(thick[df['k'], df['i'], df['j']])

    # clip it for subtle rounding errors
    df['rj'] = df['rj'].clip(1e-6,1-1e-6)
    df['ri'] = df['ri'].clip(1e-6,1-1e-6)
    df['rk'] = df['rk'].clip(1e-6,1-1e-6)
    return df

# function added thru PstFrom.add_py_function()
def cylinder_partloc(gwf, indf=None, sim_ws='.', mpsim_name='', wrkrs=1,
                     write_partcoord=False,
                     radius=20, height=2, npartob=100,
                     gclass=None, gclass_keys=None):
    '''assume loc_df has id, gk, gi, gj (grid coordinates
    as floats (i.e. fractional lay, row, col)
    Find nearest cell in 3D
    
    indf: dataframe with columns [radius, height, npartob, nztme, nztmn, boredepth] for all locations
    '''
    import numpy as np

    idom = gwf.dis.idomain.array.copy()
    top = gwf.dis.top.array.copy()
    if gclass is None or gclass_keys is None:
        gclass = idom
    else:
        if len(gclass_keys) > 0:
            temp = np.zeros(idom.shape)
            for gk in gclass_keys:
                temp = np.where(gclass==gk, gclass, temp)
            gclass = temp * idom
    gclass = gclass.astype(int)


    # get parameters for particle locations
    if 'radius' not in indf.columns:
        indf['radius'] = radius
    if 'height' not in indf.columns:
        indf['height'] = height
    if 'npartob' not in indf.columns:
        indf['npartob'] = npartob
    # indf['z'] is par, not 'boreelev' but need elev for intersection
    indf['ij'] = indf[['x', 'y']].apply(lambda x: gwf.modelgrid.intersect(x[0], x[1], z=None,
                                                                                    local=False,
                                                                                    forgive=True), axis=1)
    indf['i'] = indf['ij'].apply(lambda x: x[0])
    indf['j'] = indf['ij'].apply(lambda x: x[1])
    indf['top'] = top[tuple(indf[['i', 'j']].values.T)]

    tree = False
    # df of part locations for each site
    print(f'creating particles for {len(indf)} sites in {mpsim_name}')
    
    all_particles = []
    for site in indf.index:
        radius = indf.loc[site, 'radius']
        height = indf.loc[site, 'height']
        npartob = indf.loc[site, 'npartob']
        
        parts = uniform_parts(radius=radius, height=height, npartob=npartob)
        parts['x'] = parts['dx'] + indf.loc[site,'x']
        parts['y'] = parts['dy'] + indf.loc[site, 'y']
        parts['site'] = site

        # original boreelev is bottom, add dz to move up along cyl
        if np.isnan(indf.loc[site, 'z']):
            parts['elev'] = indf.loc[site, 'top'] - 0.5 + parts['dz'] # 0.5 is an offset from the top
        else:
            parts['elev'] = indf.loc[site, 'top'] - indf.loc[site, 'z'] + parts['dz']
        parts['kij'] = parts[['x', 'y', 'elev']].apply(lambda x: \
                                                           gwf.modelgrid.intersect(x[0], x[1], z=x[2], local=False,
                                                                                   forgive=True), axis=1)
        parts['k'] = parts['kij'].apply(lambda x: x[0])
        parts['i'] = parts['kij'].apply(lambda x: x[1])
        parts['j'] = parts['kij'].apply(lambda x: x[2])



        # take care of those sneaking out of domain, returning nan
        notin = [_ for _ in parts.index if np.isnan(parts.loc[_,'i']) or
                 np.isnan(parts.loc[_,'j']) or
                 np.isnan(parts.loc[_,'k'])]
        isin = [_ for _ in parts.index if _ not in notin]
        if len(isin)>0:
            parts.loc[isin, ['k','i','j','rk','ri','rj']] = points_to_rloc(parts.loc[isin,:], gwf).copy()

        if len(notin)>0:
            print('calculating nearest {} points for site {} not in bounds or gclass_keys'.format(len(notin), site))
            if not tree:
                tree, sidxs = build_kdtree(gwf, mask=gclass)
            # KDtree will find nearest cell center
            # even if current point in active cell already
            # so only look for "notin" subset
            outpoints = np.array([parts.loc[notin, 'x'],
                               parts.loc[notin, 'y'],
                               parts.loc[notin, 'elev']]).transpose()
            dist, indices = tree.query(outpoints, k=1, workers=wrkrs)
            parts.loc[notin, 'j'] = sidxs[indices, 0]
            parts.loc[notin, 'i'] = sidxs[indices, 1]
            parts.loc[notin, 'k'] = sidxs[indices, 2]
            parts.loc[notin, 'rj'] = 0.5
            parts.loc[notin, 'ri'] = 0.5
            parts.loc[notin, 'rk'] = 0.5
            parts.loc[notin,['k','i','j','rk','ri','rj']] = parts.loc[notin, ['k','i','j','rk','ri','rj']].copy()

        # time offset and drape, dum for now
        parts['dum'] = 0
        # for zero based
        parts[['i', 'j', 'k']] = parts[['i', 'j', 'k']] + 1
        
        
        fname = os.path.join(sim_ws, f'part_loc.{site}.csv')
        with open(fname, 'w+') as f:
            f.write('# input style (1); locationstyle (1); count, id_option (1); k,i,j (Cellnumber), LocalX, LocalY, LocalZ, TimeOffset, Drape\n')
            f.write('1\n')
            f.write('1\n')
            f.write('{} {}\n'.format(len(parts), 1))
        parts.loc[:, ['k', 'i', 'j', 'rj', 'ri', 'rk', 'dum', 'dum']].to_csv(fname, mode='a', sep=' ', header=False)
        if write_partcoord:
            fname = f'part_coord.{site}.csv'
            parts.to_csv(os.path.join(sim_ws, fname))
        all_particles.append(parts[['site', 'k', 'i', 'j', 'rj', 'ri', 'rk']].copy())
    
    # merge dfs
    if len(all_particles) > 0:
        all_particles = pd.concat(all_particles)
        # reset index
        all_particles.reset_index(inplace=True, drop=True)
        

    return all_particles


def rando_parts(radius, height, npartob):
    ang = np.random.uniform(0, 360, int(npartob))
    # relative cell location
    x = radius * np.cos(ang)
    y = radius * np.sin(ang)
    z = np.random.normal(0, height, int(npartob))
    df = pd.DataFrame([x, y, z], index=['dx', 'dy', 'dz']).T
    df.index.name = 'pnum'
    df.to_csv(os.path.join(sim_ws, 'random_parts.csv'))
    return ()


def uniform_parts(radius=2, height=2, npartob=100):
    import numpy as np
    import pandas as pd
    num_ang = 5
    num_z = np.ceil(npartob / num_ang)
    ang = np.linspace(0, 2 * np.pi, int(num_ang) + 1)[:-1]
    # relative cell location
    z = np.linspace(0, height, int(num_z))
    df_list = []
    np.random.seed(666)
    for h in z:
        da = np.random.uniform()*2*np.pi
        x = radius * np.cos([_ + da for _ in ang])
        y = radius * np.sin([_ + da for _ in ang])
        df = pd.DataFrame([x, y], index=['dx', 'dy']).T
        df['dz'] = h
        df_list.append(df)
    df = pd.concat(df_list)
    df.reset_index(inplace=True, drop=True)
    df.index.name = 'pnum'
    return df
    

def samples_to_mp7(gwf, par, mpsim_name, model_ws='.',
                   geoclass=None, gclass_keys=[]):
    # import numpy as np
    # import pandas as pd
    # import os
    # import flopy

    # if type(gwf) == str:
    #     # get model stuff
    #     sim = flopy.mf6.MFSimulation.load(sim_ws=model_ws, load_only=['dis'])
    #     gwf = sim.get_model(gwf)
    # if geoclass is None:
    #     gclist = [_ for _ in os.listdir(model_ws) if 'geoclass' in _]
    #     if len(gclist) == 1:
    #         gclass = np.stack([np.read_txt(os.path.join(model_ws, gclist[0]))] * gwf.dis.nlay.data)
    #     else:
    #         gclass = np.zeros(gwf.dis.botm.array.shape)
    #         for f in gclist:
    #             if 'geoclass_layer' in f:
    #                 k = int(f.split('geoclass_layer')[-1].split('.')[0])-1
    #                 gclass[k] = np.loadtxt(os.path.join(model_ws, f))
    # if len(gclass_keys)==0:
    #     gclass_keys = [int(_) for _ in np.unique(gclass)]
    # if type(par) == str:
    #     par = pd.read_csv(os.path.join(model_ws, par), index_col=0)
    # if 'obsnme' in par.columns:
    #     par.index = par.obsnme
    cylinder_partloc(gwf, indf=par, mpsim_name=mpsim_name, sim_ws=model_ws,
                     gclass=gclass, gclass_keys=gclass_keys)
    # write mp7 files
    with open(os.path.join(model_ws, f'{mpsim_name}.mpnam'), 'w+') as f:
        f.write('# Name file for MODPATH 7, generated by Wes\n')
        f.write(f'MPBAS      {gwf.name}.mpbas\n')
        f.write('GRBDIS     {}.dis.grb\n'.format(gwf.name))
        f.write('TDIS       {}\n'.format(gwf.simulation.tdis.filename))
        f.write('HEAD       {}.hds\n'.format(gwf.name))
        f.write('BUDGET     {}.cbc\n'.format(gwf.name))

    with open(os.path.join(model_ws, f'{gwf.name}.mpbas'), 'w+') as f:
        f.write('# MPSIM package for MODPATH 7, generated by Wes\n')
        # see Abrams et al 2012, sfr ANY is my opinion
        f.write('0\n')
        #f.write('RCH\n6\n')
        # f.write('SFR\n5\n')
        for k in range(0, gwf.dis.nlay.data):
            f.write("OPEN/CLOSE '{}.porosity_layer{}.arr' 1.0 (FREE)\n".format(gwf.name, k + 1))

    # endtime=365*1000000 #1million max, stop time option = 3
    with open(os.path.join(model_ws, f'{mpsim_name}.mpsim'), 'w+') as f:
        f.write('# MPSIM package for MODPATH 7, generated by Wes\n')
        f.write(f'{mpsim_name}.mpnam\n')
        f.write(f'{mpsim_name}.mplst\n')
        f.write(
            '1 2 2 2 1 0\n')  # 3 and 4 are pass through(1) or not (2) weak sink and source resp. last needs to be 0 or cells need to be listed below, no error issued\n')
        f.write(f'{mpsim_name}.mpend\n')
        f.write('0\n')
        f.write('1 #reference time option\n')
        f.write('0.0\n')
        f.write('2 # stop time option, time below\n')
        # f.write('{}\n'.format(endtime))
        f.write('0 #zone data option\n')
        f.write('0\n')
        f.write('{}\n'.format(len(par.index)))
        for i in par.index:
            f.write('{}\n'.format(i.split('-')[-1]))
            f.write('1\n')
            f.write('0.0\n')
            f.write("EXTERNAL 'part_loc.{}.csv'\n".format(i))
    return


def zones_to_mp7(model_name, porosity_files, zone_files, model_ws='.',
                 skip_lays=[], zone=1, ndiv=10, mpsim_name=False, tdis_name=False):
    '''issue with MODPATH, more particles than expected causing access violation'''
    assert len(porosity_files) == len(zone_files), f'pososity files not equal to zone files'
    nrow, ncol = np.loadtxt(os.path.join(model_ws,zone_files[0])).shape
    if not mpsim_name:
        mpsim_name = 'cellage_{}.backward'.format(zone)
    if not tdis_name:
        tdis_name = '{}.tdis'.format(model_name)
    # write mp7 files
    with open(os.path.join(model_ws, f'{mpsim_name}.mpnam'), 'w+') as f:
        f.write('# Name file for MODPATH 7, generated by Wes\n')
        f.write(f'MPBAS      {model_name}.mpbas\n')
        f.write('GRBDIS     {}.dis.grb\n'.format(model_name))
        f.write('TDIS       {}\n'.format(tdis_name))
        f.write('HEAD       {}.hds\n'.format(model_name))
        f.write('BUDGET     {}.cbc\n'.format(model_name))

    with open(os.path.join(model_ws, f'{model_name}.mpbas'), 'w+') as f:
        f.write('# MPBAS package for MODPATH 7, generated by Wes\n')
        # see Abrams et al 2012, sfr ANY is my opinion since my streams are often near the bottom of the cell
        f.write('0\n')
        # f.write('RCH\n6\n')
        # f.write('SFR\n5\n')
        for pf in porosity_files:
            f.write(f"OPEN/CLOSE '{os.path.basename(pf)}' 1.0 (FREE)\n")

    # endtime=365*1000000 #1million max, stop time option = 3
    with open(os.path.join(model_ws, f'{mpsim_name}.mpsim'), 'w+') as f:
        f.write('# MPSIM package for MODPATH 7, generated by Wes\n')
        f.write(f'{mpsim_name}.mpnam\n')
        f.write(f'{mpsim_name}.mplst\n')
        f.write(
            '1 2 2 2 1 0\n')  # 3 and 4 are pass through(1) or not (2) weak sink and source, resp. last needs to be 0 or cells need to be listed below, no error issued\n')
        f.write(f'{mpsim_name}.mpend\n')
        f.write('0\n')
        f.write('1 #reference time option\n')
        f.write('0.0\n')
        f.write('2 # stop time option, time below\n')
        # f.write('{}\n'.format(endtime))
        f.write('0 #zone data option\n')
        f.write('0\n')
        f.write('1\n')
        f.write('zone_{}\n'.format(zone))
        f.write('1\n')
        f.write('0.0\n')
        f.write('INTERNAL\n')
        f.write('4\n')
        f.write('2 0\n')
        f.write('1\n')
        f.write('{}\n'.format(zone))
        for k in range(0, len(zone_files)): # for each layer
            if k in skip_lays:
                blank = os.path.join(model_ws, 'blank.arr')
                if not os.path.exists(blank):
                    np.savetxt(blank, np.zeros((nrow,ncol)), fmt='%i')
                f.write("OPEN/CLOSE '{}' 1 (FREE) 0\n".format('blank.arr'))
            else:
                f.write(f"OPEN/CLOSE '{zone_files[k]}' 1 (FREE) 0\n")
        f.write(f'{ndiv[0]} {ndiv[1]} {ndiv[2]}\n')
    return


def cells_to_mp7(gwf, arr, lays=[], pg='loc', zone=1, ndiv=10, model_ws='.', tdis_filename=False, mpsim_name=False):
    '''lays: zero-based'''
    if not tdis_filename:
        tdis_filename = f'{gwf.name}.tdis'
    idx = np.where(arr == zone)
    if len(lays) > 0:
        rows, cols = idx
        lays = lays * len(rows)
    else:
        if len(idx) == 3:
            lays, rows, cols = idx
        else:
            rows, cols = idx
            # lays needs to be same len as rows and cols
            lays = list(range(0, gwf.dis.nlay.data)) * len(rows)

    if not mpsim_name:
        mpsim_name = 'cellage_{}.backward'.format(zone)
    # tdis in sim not gwf, passing sim is stupid
    for f in os.listdir(model_ws):
        if f.endswith('.tdis'):
            tdis_name = f
    # write mp7 files
    with open(os.path.join(model_ws, f'{mpsim_name}.mpnam'), 'w+') as f:
        f.write('# Name file for MODPATH 7, generated by Wes\n')
        f.write(f'MPBAS      {gwf.name}.mpbas\n')
        f.write('GRBDIS     {}.dis.grb\n'.format(gwf.name))
        f.write('TDIS       {}\n'.format(tdis_filename))
        f.write('HEAD       {}.hds\n'.format(gwf.name))
        f.write('BUDGET     {}.cbc\n'.format(gwf.name))

    with open(os.path.join(model_ws, f'{gwf.name}.mpbas'), 'w+') as f:
        f.write('# MPSIM package for MODPATH 7, generated by Wes\n')
        # see Abrams et al 2012, sfr ANY is my opinion
        f.write('0\n')
        # f.write('RCH\n6\n')
        #f.write('SFR\n5\n')
        for k in range(0, gwf.dis.nlay.data):
            f.write("OPEN/CLOSE '{}.porosity_layer{}.arr' 1.0 (FREE)\n".format(gwf.name, k + 1))

    # endtime=365*1000000 #1million max, stop time option = 3
    with open(os.path.join(model_ws, f'{mpsim_name}.mpsim'), 'w+') as f:
        f.write('# MPSIM package for MODPATH 7, generated by Wes\n')
        f.write(f'{mpsim_name}.mpnam\n')
        f.write(f'{mpsim_name}.mplst\n')
        f.write(
            '1 2 2 2 1 0\n')  # 3 and 4 are pass through(1) or not (2) weak sink and source, resp. last needs to be 0 or cells need to be listed below, no error issued\n')
        f.write(f'{mpsim_name}.mpend\n')
        f.write('0\n')
        f.write('1 #reference time option\n')
        f.write('0.0\n')
        f.write('2 # stop time option, time below\n')
        # f.write('{}\n'.format(endtime))
        f.write('0 #zone data option\n')
        f.write('0\n')
        f.write('1\n')
        f.write('{}_{}\n'.format(pg, zone))
        f.write('1\n')
        f.write('0.0\n')
        f.write('INTERNAL\n')
        f.write('2\n')
        f.write('1 {}\n'.format(len(rows)))
        f.write('1 {} 0\n'.format(len(rows)))
        # TemplateSubdivisionType=1 (above) requires vert and hor ndiv for all faces,
        # e.g. 0,0,0,0,0,0,1,1,0,0 for 1 division on bottom face
        f.write('0 0 0 0 0 0 0 0 {} {} 0 0\n'.format(ndiv, ndiv))
        for k, i, j in zip(lays, rows, cols):
            f.write('{0} {1} {2} {0} {1} {2}\n'.format(k + 1, i + 1, j + 1))
    return


def df_to_ins(df, ins_pth, columns=[], poi=[]):
    df = df.copy()
    df['l1'] = 'l1'
    for c in columns:
        if c in poi:
            df['t_' + c] = df[columns[0]].apply(lambda x: '!{}_poi_{}!'.format(x, c))
        else:
            df['t_' + c] = '~,~'
    nc = ['l1'] + ['t_' + i for i in columns]

    with open(ins_pth, 'w+') as f:
        f.write('{}\n{}\n'.format('pif ~', 'l1'))
    df.loc[:, nc] \
        .to_csv(ins_pth,
                index=False, header=False, sep=' ', mode='a')
    return ()


def mp7_random_samples(gwf, gclass_file=False, gsgp=[191], model_ws='.', poly_path=None,
                       sample_size=0.5, height=2, min_depth=0, max_depth=80, tracer=None):
    # setting random state should ensure each PstFrom build generates same
    # different states for depth and date
    depth_rand_state = np.random.RandomState(1234567890)
    date_rand_state = np.random.RandomState(666666666)
    # gwf stuff
    if type(gwf) == str:
        sim = flopy.mf6.MFSimulation.load(sim_ws=model_ws, load_only=['dis', 'ghb', 'wel'], verify_data=False)
        gwf = sim.get_model()
    idom = gwf.dis.idomain
    depth = (gwf.dis.top.array - gwf.dis.botm.array) * idom
    x, y, z = gwf.modelgrid.xyzcellcenters
    # geologic class or zone array
    if not gclass_file:
        geoclass = np.loadtxt(os.path.join(model_ws, f'{gwf.name}.geoclass.arr'))
    else:
        geoclass = np.loadtxt(os.path.join(model_ws, gclass_file))
    geoclass = geoclass * idom

    # new points and simulation for each zone
    mpsim_names = []
    for gs in gsgp:
        mpsim_name = f'{tracer}.rand{gs}'
        mpsim_names.append(mpsim_name)
        # 3d array of idom exclude
        # TODO: exclude wel and ghb boudnary cells
        idx = [_ for _ in np.argwhere((geoclass == gs) & (depth < max_depth))]
        # all points of gs
        df = pd.DataFrame(idx, columns=['k', 'i', 'j'])
        df['x'] = x[df.i, df.j]
        df['y'] = y[df.i, df.j]
        df['height'] = height
        df['depth'] = depth[df.k, df.i, df.j]
        df['uni_depth'] = depth_rand_state.uniform(0, 1, len(idx))
        df['uni_date'] = date_rand_state.uniform(0, 1, len(idx))
        df['min_depth'] = min_depth
        df['max_depth'] = max_depth
        df.loc[:, 'z'] = df.loc[:, 'min_depth'] + \
                                  df.loc[:, 'uni_depth'] * (df.loc[:, ['max_depth', 'depth']].min(axis=1))
        df['nom_depth'] = df.apply(lambda x: int(round(x['z'])), axis=1)
        df['node'] = df.apply(lambda x: int(gwf.modelgrid.get_node(x[['k', 'i', 'j']])[0]), axis=1)

        # sorted by e, n , depth
        # necessary to sort AFTER DROPPING OOB if going between different grids with same pst obs
        df.sort_values(['x', 'y', 'z'], inplace=True)
        df.reset_index(inplace=True, drop=True)
        df['part_gp'] = df.index
        df['obsnme'] = df.apply(lambda x: f"pg_{int(x['part_gp'])}_{x['nom_depth']}", axis=1)
        df.index = df.obsnme
        df.sort_index(inplace=True)
        if 'site_id' not in df.columns and 'obsnme' in df.columns:
            df['site_id'] = df['obsnme']

        # add group if poly_path included
        if poly_path is not None:
            df = utils.poly_groups(df, poly_path)
            df_list = []
            for rc in df.gnme.unique():
                # cell center sampling screws this up
                if sample_size <= 1:
                    npoints = int(df[df.gnme == rc].shape[0] * sample_size)
                else:
                    npoints = int(sample_size)
                ridx = depth_rand_state.choice(df[df.gnme == rc].index, size=npoints, replace=False)
                df_list.append(df.loc[ridx, :])
            df = pd.concat(df_list)
        else:
            if sample_size <= 1:
                npoints = int(df.shape[0] * sample_size)
            else:
                npoints = int(sample_size)
            ridx = depth_rand_state.choice(df.index, size=npoints, replace=False)
            df = df.loc[ridx, :]
        # need some date to meet requirements of df, related to tr
        df['min_date'] = pd.to_datetime('01/01/1940', format='%d/%m/%Y')
        df['max_date'] = pd.to_datetime('01/01/2020', format='%d/%m/%Y')
        df[tracer] = 0
        df['sigtr'] = 0
        if 'date' not in df.columns:
            df['datetime'] = (df['min_date'] + df['uni_date'] * (df['max_date'] - df['min_date']))
            df['date'] = df['datetime'].dt.date
        df[['site_id', 'date', 'x', 'y', 'z',
            'height', tracer, 'sigtr']] \
            .to_csv(os.path.join(model_ws, f'{mpsim_name}.obs.csv'), index='obsnme')
    return mpsim_names


# def mp7_previous_samples(obs_data, model_ws='.', mpsim_name='', model_name=''):
#     # get obs locations from previous run for consistency
#     if type(obs_data) == str:
#         obs_data = pd.read_csv(obs_data)
#     # need bottom of bore for samples_to_mp7, not mean
#     height = 2
#     obs_data['z'] = obs_data['mean_depth'] + height / 2
#     df = samples_to_mp7(model_ws, model_name, x_field='mean_nztme', y_field='mean_nztmn',
#                         z_field='z',
#                         obs_data=obs_data.loc[:, ['mean_nztme', 'mean_nztmn', 'z']],
#                         npartob=100, radius=10, height=height, nan_depth=10, cellidx=True,
#                         min_depth=0.2, mpsim_name=mpsim_name, tracer='age')
#     cylinder_partloc(indf=df, mpsim_name=mpsim_name, sim_ws=model_ws)
#     return


def mpsim_to_obsnme(mpsim, mpsim_name, sim_ws='.'):
    # requires mpsim_name to trigger how obsnme is determined
    with open(os.path.join(sim_ws, f'{mpsim_name}.mpsim'), 'r') as f:
        data = f.read()
    # if 'boreelev' not in mpsim.columns:
        # mpsim['boreelev'] = mpsim['z0']
    if 'external' in data.lower():
        # particle groups defined by external part_loc file order in mpsim
        pg = 0
        with open(os.path.join(sim_ws, f'{mpsim_name}.mpsim'), 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'external' in [_.lower() for _ in line.strip().split()]:
                mpsim.loc[mpsim['particlegroup'] == pg, 'obsnme'] = \
                    line.strip().split()[-1].split('.',1)[-1].rsplit('.',1)[0]
                pg = pg + 1
    else:
        # all particles same mp group, define obsnme by particle
        # derive particle obsnme from order in mpend file
        # doesn't guarantee the same name due to rounding... what to do?
        if 'z' not in mpsim.columns:
            mpsim['obsnme'] = mpsim[['particleid', 'k0']]. \
                apply(lambda x: f"{mpsim_name.replace('.', '_')}_{int(x[0])}_lay_{int(x[1])}", axis=1)
            mpsim['z'] = mpsim['model_top'] - mpsim['z0'] + mpsim['zloc']
        else:
            mpsim['obsnme'] = mpsim[['particleid', 'z']].\
                apply(lambda x: f"{mpsim_name.replace('.','_')}_{int(x[0])}_{int(x[1])}", axis=1)
        mpsim.sort_values(by='obsnme', inplace=True)
    return mpsim

def mp_to_table(mpsim_name, model_name=None,
              file_type='obs', obs_type=None, sim_ws='.'):
    '''df from mpsim results, depth from simulation'''
    if model_name is None:
        with open(os.path.join(sim_ws,'mfsim.nam'), 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'gwf6' in line.lower():
                model_name = line.split()[1].replace('.nam','')
    if obs_type is None:
        obs_type = mpsim_name.split('.')[0]
    outfile = f'{mpsim_name}.{file_type}.csv'
    # get sim stuff
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, load_only=['dis'])
    gwf = sim.get_model(model_name)
    top = gwf.dis.top.array
    mg = gwf.modelgrid
    fpth = os.path.join(sim_ws, f'{mpsim_name}.mpend')
    eobj = flopy.utils.EndpointFile(fpth)
    e = eobj.get_alldata()
    mpsim = pd.DataFrame.from_records(e)
    mpsim.fillna(99999, inplace=True)
    mpsim['kij'] = mg.get_lrc(list(mpsim.node0.values))

    # 0 based for array compatibility
    mpsim.loc[:, 'k'] = mpsim.loc[:, 'kij'].apply(lambda x: x[0])
    mpsim.loc[:, 'i'] = mpsim.loc[:, 'kij'].apply(lambda x: x[1])
    mpsim.loc[:, 'j'] = mpsim.loc[:, 'kij'].apply(lambda x: x[2])

    mpsim['model_top'] = top[tuple(mpsim[['i', 'j']].values.T)]
    mpsim['x'] = mpsim['x0'] + mg.xoffset
    mpsim['y'] = mpsim['y0'] + mg.yoffset
    mpsim['boreelev'] = mpsim['z0']
    #mpsim['z'] = mpsim['model_top'] - mpsim['boreelev']

    mpsim = mpsim_to_obsnme(mpsim, mpsim_name, sim_ws=sim_ws)
    mpsim.set_index('obsnme', inplace=True, drop=True)
    if 'date' not in mpsim.columns:
        # arbitrary sample date
        mpsim['date'] = pd.to_datetime('20/02/2020', format='%d/%m/%Y')
    mpsim[obs_type] = 1e-9
    # 1 based for compatibility with model input files
    mpsim.loc[:, 'k'] = mpsim.loc[:, 'k'] + 1
    mpsim.loc[:, 'i'] = mpsim.loc[:, 'i'] + 1
    mpsim.loc[:, 'j'] = mpsim.loc[:, 'j'] + 1
    mpsim.sort_index(inplace=True)
    mpsim[[obs_type, 'x', 'y', 'model_top', 'boreelev', 'z', 'date', 'i', 'j', 'k']]. \
        to_csv(os.path.join(sim_ws, outfile), date_format='%d/%m/%Y')
    return outfile

def mp_to_age(mpsim_name, model_name=None, sim_ws='.', gclass_dict={},
              obs_file=None, outname=None, make_figs=False,
              quart_mult=1.5, make_obs=False):
    # assumes simulation results are in days
    if model_name is None:
        with open(os.path.join(sim_ws, 'mfsim.nam'), 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'gwf6' in line.lower():
                model_name = line.split()[1].replace('.nam', '')
    if obs_file is None:
        obs_file = f'{mpsim_name}.obs.csv'
    obs = pd.read_csv(os.path.join(sim_ws, obs_file))
    # try:
    #     obs['date'] = pd.to_datetime(obs['date']).dt.date
    # except:
    #     obs['date'] = pd.to_datetime(obs['date'], format='%d/%m/%Y').dt.date
    if outname is None:
        outname = obs_file.replace('.obs', '.mrt.sim')
    # get sim stuff
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, load_only=['dis'])
    gwf = sim.get_model(model_name)
    idom = gwf.dis.idomain.array
    botm = gwf.dis.botm.array
    top = gwf.dis.top.array
    # center of layers
    depth = [(top - botm[0]) / 2]
    for k in range(1, botm.shape[0]):
        thick = (botm[k - 1] - botm[k]) / 2
        depth.append(top - botm[k] + thick)
    depth = np.stack(depth, axis=0) * idom
    mg = gwf.modelgrid
    fpth = os.path.join(sim_ws, f'{mpsim_name}.mpend')
    eobj = flopy.utils.EndpointFile(fpth)
    e = eobj.get_alldata()
    mpsim = pd.DataFrame.from_records(e)
    mpsim = mpsim.fillna(99999)
    # drop outliers if enough particles in each unique group
    if mpsim.groupby('particlegroup').count().particleid.gt(10).all() and \
            mpsim.groupby('particlegroup').count().particleid.ne(len(mpsim)).all():
        # drop outlier particles within particle group
        # drop outlier particles within particle group
        med = mpsim.groupby('particlegroup').median()['time']
        uq = mpsim.groupby('particlegroup').quantile(0.75)['time']
        lq = mpsim.groupby('particlegroup').quantile(0.25)['time']
        iqr = uq - lq
        for mpg in med.index:
            group = [_ for _ in mpsim.index if mpsim.loc[_, 'particlegroup'] == mpg]
            sel = [_ for _ in group if
                   (mpsim.loc[_, 'time'] > (med[mpg] + iqr[mpg] * quart_mult) or
                   mpsim.loc[_, 'time'] < (med[mpg] - iqr[mpg] * quart_mult))]
            if len(sel) >= 1 and len(sel) != len(group):
                mpsim.drop(sel,inplace=True)
            print(f'dropping {len(sel)} from {len(group)} as outliers for particle group {mpg}')
        
    # initial cell
    mpsim['kij0'] = mg.get_lrc(list(mpsim.node0.values))
    mpsim.loc[:, 'k0'] = mpsim.loc[:, 'kij0'].apply(lambda x: x[0])
    mpsim.loc[:, 'i0'] = mpsim.loc[:, 'kij0'].apply(lambda x: x[1])
    mpsim.loc[:, 'j0'] = mpsim.loc[:, 'kij0'].apply(lambda x: x[2])
    # model top at initial location
    top = gwf.dis.top.array
    mpsim['model_top'] = top[tuple(mpsim[['i0', 'j0']].values.T)]
    # exit cell
    mpsim['kij'] = mg.get_lrc(list(mpsim.node.values))
    mpsim.loc[:, 'k'] = mpsim.loc[:, 'kij'].apply(lambda x: x[0])
    mpsim.loc[:, 'i'] = mpsim.loc[:, 'kij'].apply(lambda x: x[1])
    mpsim.loc[:, 'j'] = mpsim.loc[:, 'kij'].apply(lambda x: x[2])

    # get obsnme from mpend header
    mpsim = mpsim_to_obsnme(mpsim, mpsim_name, sim_ws=sim_ws)
    # merge sim with sample to get date, time needs to be in time delta days
    mpsim = mpsim.merge(obs, on='obsnme', suffixes=[None, 'obs'])
    obs.set_index('obsnme', inplace=True, drop=True)
    # age for all sfr, rch cells is near 0 at surface
    # 1 based nodal ages tables with k, i, j, node
    bcage_list = []
    for bc in ['rch', 'sfr', 'mbr']:
        if os.path.exists(os.path.join(sim_ws, f"bcage.{bc}.csv")):
            bcage_list.append(
                pd.read_csv(os.path.join(sim_ws, f"bcage.{bc}.csv")))
    assert len(bcage_list) > 0
    bcage = pd.concat(bcage_list, axis=0)
    bcage.dropna(axis=0, subset=['k', 'i', 'j'], inplace=True)
    # convert i,j in bcage to 0 based for merge
    bcage[['k', 'i', 'j']] -= 1
    bcage['node'] = (bcage.k * gwf.dis.ncol.data * gwf.dis.nrow.data +
                     bcage.i * gwf.dis.ncol.data +
                     bcage.j)    
    # node already 0 based derived from 0 based index above
    # bcage['node'] -= 1
    bcage['node'] = bcage['node'].astype(int)
    bcnm = bcage.groupby('node').mean(numeric_only=True)
    
    # start working in years
    bcnm['bcage'] = bcnm['age']/365.25
    mpsim['time'] = mpsim['time'] / 365.25
    
    # combine on exit node
    # flopy.utils.EndpointFile is 0 based as should be mg.get_lrc
    assert mpsim['node'].isna().sum() == 0
    mpsim['node'] = mpsim['node'].astype(int)
    mpsim.set_index('node', drop=True, inplace=True)
    mpsim = bcnm[['bcage']].merge(mpsim, on='node', how='right')
    # some particles ending at wonky places in stochastic reals, e.g. GHB
    mpsim['bcage'].fillna(0, inplace=True)    # fillna(1000*365.25, inplace=True)#tr_df['time'].max(), inplace=True)
    mpsim['bcage'].clip(1e-8, None, inplace=True)
    mpsim['mrtsim'] = mpsim['time'] + mpsim['bcage']
    assert mpsim['mrtsim'].ge(0).all()

    # get stats
    df_list = []
    stats = ['mean', 'median', 'std', 'min', 'max']
    cols = []
    # groupby obliterates non numeric fields
    for q in mpsim.select_dtypes(include=np.number).columns:
        df_list.append(mpsim[['obsnme', q]].groupby('obsnme').mean(numeric_only=True))  # .clip(1e-40, None))
        df_list.append(mpsim[['obsnme', q]].groupby('obsnme').median(numeric_only=True))  # .clip(1e-40, None))
        df_list.append(mpsim[['obsnme', q]].groupby('obsnme').std(numeric_only=True))  # .fillna(0).clip(1e-40, None))
        df_list.append(mpsim[['obsnme', q]].groupby('obsnme').min(numeric_only=True))  # .clip(1e-40, None))
        df_list.append(mpsim[['obsnme', q]].groupby('obsnme').max(numeric_only=True))  # .clip(1e-40, None))
        cols = cols + [f"{q}_{_}" for _ in stats]
    mdf = pd.concat(df_list, axis=1)
    mdf.columns = cols
    mdf['nztme_mean'] = mg.xoffset + mdf['x0_mean']
    mdf['nztmn_mean'] = mg.yoffset + mdf['y0_mean']
    # particles based on bottom of bore
    mdf[['boreelev_mean','model_top_mean']] = mpsim[['obsnme', 'z0', 'model_top']].groupby(by='obsnme').min(numeric_only=True)
    mdf['boredepth_mean'] = mdf['model_top_mean'] - mdf['boreelev_mean']
    mdf = mdf.copy()
    # the 64 bit integer limits determine the Timedelta limits, ~106751 days
    # so input_date is short, only use to look up tritium input, not decay
    mdf['tunits'] = pd.to_timedelta(mdf['mrtsim_mean'].clip(1, 100000), 'd')

    # group by gclass in gclass_dict
    # TODO: pass key instead of implicit 'mrt'
    if len(gclass_dict.keys())>0:
        if 'gclass' in obs.columns:
            obs['gclass'] = obs['gclass'].astype(int)
        sim_vals = {}
        obs_vals = {}
        for name in gclass_dict.keys():
            # mean, min, max of mean of all particles at each site
            sm = mdf.loc[obs['gclass'].isin(gclass_dict[name]),'mrtsim_mean'].mean()
            sn = mdf.loc[obs['gclass'].isin(gclass_dict[name]),'mrtsim_mean'].min()
            sx = mdf.loc[obs['gclass'].isin(gclass_dict[name]),'mrtsim_mean'].max()
            sim_vals[name] = {'gclass_mean': sm, 'gclass_min': sn, 'gclass_max': sx}
        gsim = pd.DataFrame.from_dict(sim_vals, orient='index')
        gsim.index.name = 'geoname'
        gsim.sort_index(inplace=True)
        gsim.to_csv(os.path.join(sim_ws, outname.replace('mrt.sim', 'mrtgclass.sim')))
        if 'mrt' in obs.columns:
            for name in gclass_dict.keys():
                om = obs.loc[obs['gclass'].isin(gclass_dict[name]), 'mrt'].mean()
                on = obs.loc[obs['gclass'].isin(gclass_dict[name]), 'mrt'].min()
                ox = obs.loc[obs['gclass'].isin(gclass_dict[name]), 'mrt'].max()
                obs_vals[name] = {'gclass_mean': om, 'gclass_min': on, 'gclass_max': ox}
            gobs = pd.DataFrame.from_dict(obs_vals, orient='index')
            gobs.index.name = 'geoname'
            gobs.sort_index(inplace=True)
            gobs.to_csv(os.path.join(sim_ws, outname.replace('mrt.sim', 'mrtgclass.obs')))

    mdf['log_mrtsim'] = np.log10(mdf['mrtsim_mean'])
    mdf['log_mrtsim'].fillna(-50, inplace=True)
    outcols = ['mrtsim_mean', 'mrtsim_median', 'mrtsim_min', 'mrtsim_max', 'mrtsim_std',
               'bcage_mean', 'bcage_median', 'bcage_min', 'bcage_max', 'bcage_std',
               'nztme_mean', 'nztmn_mean', 'boredepth_mean', 'boreelev_mean', 'model_top_mean',
               'log_mrtsim', 'time_mean', 'tunits']
    mdf.sort_index(inplace=True)
    mdf[outcols]. \
        to_csv(os.path.join(sim_ws, outname), date_format='%d/%m/%Y')
    if make_figs:
        for site in mdf.index:
            mpsim.loc[site,['mrtsim']].hist(bins=100)
            plt.title(f"age for site {site}\n,"
                      f"min = {mdf.loc[site, f'mrtsim_min']},\n"
                      f"mean = {mdf.loc[site, f'mrtsim_mean']},\n"
                      f"max = {mdf.loc[site, f'mrtsim_max']},\n"
                      f"mean bcage = {mdf.loc[site, f'bcage_mean']}")
            if not os.path.exists(os.path.join(sim_ws, '..', 'figures')):
                os.mkdir(os.path.join(sim_ws, '..', 'figures'))
            plt.savefig(os.path.join(sim_ws, '..', 'figures', f'{site}.png'))
            plt.close()
    return mdf

def mp_to_tracer(mpsim_name, model_name=None, sim_ws='.',
                 obs_file=None,
                 par_file=None,
                 tracer_input_file=None,
                 obs_type=None, outname=None,
                 make_figs=False):
    import matplotlib.pyplot as plt
    # assumes simulation results are in days
    if model_name is None:
        with open(os.path.join(sim_ws,'mfsim.nam'), 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'gwf6' in line.lower():
                model_name = line.split()[1].replace('.nam','')
    if obs_type is None:
        obs_type = mpsim_name.split('.')[0]
    if obs_file is None:
        obs_file = f'{mpsim_name}.obs.csv'
    obs = pd.read_csv(os.path.join(sim_ws, obs_file), index_col=0)
    if par_file is None:
        par_file = f'{mpsim_name}.par.csv'
    par = pd.read_csv(os.path.join(sim_ws, par_file), index_col=0)
    if outname is None:
        outname = obs_file.replace('.obs','.sim')

    # get sim stuff
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, load_only=['dis'])
    gwf = sim.get_model(model_name)
    idom = gwf.dis.idomain.array
    botm = gwf.dis.botm.array
    top = gwf.dis.top.array
    # center of layers
    depth = [(top - botm[0])/2]
    for k in range(1,botm.shape[0]):
        thick = (botm[k-1] - botm[k])/2
        depth.append(top - botm[k] + thick)
    depth = np.stack(depth, axis=0)*idom

    # TODO: just read the file
    mpsim = mp_to_age(mpsim_name, model_name=model_name, sim_ws=sim_ws)
    mpsim['date'] = pd.to_datetime(obs['date'], format='%d/%m/%Y')
    mpsim['input_date'] = mpsim['date'] - mpsim['tunits']
    if tracer_input_file is not None:
        # read tracer input data, merge to sim data with sample date (mpsim) on input date
        tracer_input = pd.read_csv(tracer_input_file)
        # concat to get all dates, drop NaT, watch format!!
        tracer_input['input_date'] = pd.to_datetime(tracer_input.input_date, format='%d/%m/%Y')
        tracer_input = pd.concat([tracer_input[['input_date',obs_type, 'input_factor']],
            mpsim[['input_date']]],axis=0)
        tracer_input.dropna(subset=['input_date'], inplace=True)
        tracer_input.drop_duplicates(subset=['input_date'], inplace=True)
        tracer_input.set_index('input_date', inplace=True, drop=True)
        tracer_input.sort_index(inplace=True)
        # interpolate where possible, then backfill the rest
        tracer_input.interpolate(method='time', inplace=True)
        tracer_input.bfill(inplace=True)
        if obs_type == 'tritium':
            # N(t) = No * e**(-lam/t)
            # e(-lam*t) is decay term with lam being ln(2)/T1/2
            # (radioactive decay of tritium with a half-life 12.32 yrs),
            c0 = mpsim.loc[:, 'input_date'].\
                apply(lambda x: tracer_input.loc[x, obs_type]*(tracer_input.loc[x, 'input_factor']))
            c0.index = mpsim.index
            mpsim[f'init_{obs_type}'] = c0
            mpsim.loc[:, f"{obs_type}sim"] = mpsim.loc[:, f'init_{obs_type}'] * \
                                              np.exp(-(np.log(2) / (12.32)) * mpsim.loc[:, 'mrtsim_mean'])
            cols = [f"{obs_type}sim",f'init_{obs_type}']
    mpsim[f"{obs_type}sim"].clip(1e-10, None, inplace=True)
    # log transform after stats on particles
    mpsim[f"log_{obs_type}sim"] = np.log10(mpsim[f"{obs_type}sim"])
    mpsim[f"log_{obs_type}sim"].fillna(-50, inplace=True)
    cols = cols + [f"log_{obs_type}sim"]

    mpsim['log_bcage_mean'] = np.log10(mpsim['bcage_mean'])
    mpsim['log_bcage_mean'].fillna(-50, inplace=True)
    cols = cols + ['log_bcage_mean']

    mpsim[f"log_{obs_type}"] = np.log10(obs[obs_type])
    mpsim[f"log_{obs_type}"].fillna(-50, inplace=True)
    mpsim[f"log_{obs_type}_mean_res"] = mpsim[f"log_{obs_type}sim"] - mpsim[f"log_{obs_type}"]
    cols = cols + [f"log_{obs_type}",f"log_{obs_type}_mean_res"]#,f"log_{obs_type}_median_res"]

    mpsim.sort_index(inplace=True)
    mpsim[cols +
        ['input_date', 'date', 'nztme_mean', 'nztmn_mean', 'boreelev_mean', 'boredepth_mean']].\
        to_csv(os.path.join(sim_ws, outname))
    return

def rep_with(og, f='site_id', rep_with={}, part=True, axis=0):
    # common markers used interchangably, including double space
    def rep(k, axis=0, part=True):
        for i in rep_with.keys():
            for j in rep_with[i]:
                k = [c.replace(j, i) for c in k if \
                     (part or str(c) == j)]
        return k

    if type(og) == pd.DataFrame:
        og.columns = og.columns.str.lower()
        if axis == 0:  # replace in index
            for i in rep_with.keys():
                for j in rep_with[i]:
                    og[f] = og[f] \
                        .apply(lambda x: str(x).strip().replace(j, i).lower() if \
                        (part or str(x) == j) else str(x))
        elif axis == 1:  # replace in column
            lst = [c.strip() for c in og.columns]
            cols = rep(lst, axis=axis)
            og.columns = cols
            if 'id' in og.columns and 'site_id' not in og.columns:
                og['site_id'] = og['id']
                og.drop('id', inplace=True, axis=1)
        else:
            AssertionError(axis in [0, 1])
    elif type(og) == list:
        og = [str(_).lower() for _ in og]
        og = rep(og)
        if 'id' in og and 'site_id' not in og:
            idx = og.index('id')
            og[idx] = 'site_id'
    else:
        AssertionError(type(og) not in [list, pd.DataFrame])
    return og



def clean_chem(fpth, sample_file, sample_req, loc_req, loc_file=False, out_file=False,
               rep_with_part={'_': [' ', '/', '-', '__'], '@': [' @', '@ ', '_@', '@_']},
               rep_with_whole={'nan': ['NaN', 'Nan', 'NA']},
               add_cols=['z', 'sigtr'],
               attr_by_keys={'z': {0: ['river', 'stream', 'creek', 'culvert',
                                                'lake', 'pond', 'pool', 'drain',
                                                'spring', 'seep', 'swamp', 'marsh']},
                             2: ['shallow']}):
    '''attempt to clean messy chem files
        fpth: path to files and where new files will be written
        sample_file: file containing fields in 'sample_req'
        loc_file: file containing fields in 'loc_req'
        sample_req: should be original file column names
        loc_req: should be original file column names
        attr_by_keys: will reset field identified by level 0 key to
                value identified by level 1 key if site_id contains
                one of the keywords in the list
        rep_with: replace all occurances of chars in cols and id with key
        TODO: allow "cols" to be original or anticipated?
        TODO: map some fields to others? (like id to site_id)
        TODO: add_cols is in "cleaned" format, inconsistent with *_req
        Example: mdf=utils.clean_chem(data_dir,sample_file,
                                      loc_file=loc_file,
                                      loc_req=['ID','NZTM_E','NZTM_N'],
                                      sample_req=sample_req,
                                      add_cols=['z','sigtr','tritium_code'])
        mdf=utils.clean_chem(r'..\\..\\data\\chemistry\\wairau','GNS-SR_2019-063 Marlborough All Results from LIMS with Coordinates.csv',
                             loc_file='MDC__AllResults_for_Conny.csv',sample_req=['Site ID','tritium','Date'],loc_req=['ID','NZTM_E','NZTM_N'],
                             add_cols=['z','sigtr'],out_file='wairau_tritium_clean.csv')
        '''
    # sample file
    sample_df = pd.read_csv(os.path.join(fpth, sample_file), usecols=sample_req)
    assert set(sample_req).issubset(set(sample_df.columns))
    sample_df = rep_with(sample_df, rep_with=rep_with_part, axis=1)
    sample_df = rep_with(sample_df, rep_with=rep_with_part, axis=0)
    sample_req = rep_with(sample_req, rep_with=rep_with_part)
    assert set(sample_req).issubset(set(sample_df.columns))
    sample_df = sample_df.replace('nan', np.nan)

    # location file
    loc_df = pd.read_csv(os.path.join(fpth, loc_file))
    assert set(loc_req).issubset(set(loc_df.columns))
    loc_df = rep_with(loc_df, rep_with=rep_with_part, axis=1)
    loc_df = rep_with(loc_df, rep_with=rep_with_part, axis=0)
    loc_req = rep_with(loc_req, rep_with=rep_with_part)
    assert set(loc_req).issubset(set(loc_df.columns))
    loc_df = loc_df.replace('nan', np.nan)
    add_cols = rep_with(add_cols, rep_with=rep_with_part)
    out_cols = list(set(loc_req + sample_req + add_cols))

    for ac in add_cols:
        if ac not in sample_df.columns:
            sample_df[ac] = np.nan
        if ac not in loc_df.columns:
            loc_df[ac] = np.nan

    # some rows are missing coords, but same site_id has them further down (e.g. P28W/4402)
    # outer merge plus dropna leaves all rows with valid loc_req from both
    # gives all site_id : location combinations
    gloc = sample_df[loc_req].merge(loc_df[loc_req], how='outer').dropna().drop_duplicates()
    # merge loc_df cols with gloc
    loc_df = gloc.merge(loc_df, how='left', on=loc_req)

    # need to update nan locations in sample_df with same site_id from loc_df
    nansams = [sample_df.loc[_,'site_id'] for _ in sample_df.index if np.isnan(sample_df.loc[_, 'x'])]
    nansams = nansams + [sample_df.loc[_,'site_id'] for _ in sample_df.index if np.isnan(sample_df.loc[_, 'y'])]
    for nansam in nansams:
        if nansam in loc_df.site_id.values:
            sample_df.loc[sample_df['site_id'] == nansam, ['x', 'y']] = loc_df.loc[loc_df['site_id'] == nansam, ['x', 'y']]
    nansams = [sample_df.loc[_, 'site_id'] for _ in sample_df.index if np.isnan(sample_df.loc[_, 'z'])]
    for nansam in nansams:
        if nansam in loc_df.site_id.values:
            depth = np.mean(loc_df.loc[loc_df['site_id'] == nansam, 'z'].values)
            sample_df.loc[sample_df['site_id'] == nansam, 'z'] = depth

    # merge sample_df cols with loc_df
    # on site_id will blend all data for same site id, huge issue for wairau_river
    # this is what screwed me!!!
    # # # # # mdf = loc_df.merge(sample_df, how='left', on='site_id')

    # drop where sample_req is nan (other tracers, no date, etc)
    idx = sample_df[sample_df[sample_req].isna().any(axis=1)].index
    sample_df.drop(idx, inplace=True)
    # assign some missing depths based on key words
    for k in attr_by_keys.keys():
        if k in sample_df.columns:
            v_s = attr_by_keys[k]
            for v in v_s.keys():
                mask = []
                for key in v_s[v]:
                    mask = mask + [i for i in sample_df.index if np.isnan(sample_df.loc[i, k]) \
                                   and key in sample_df.loc[i, 'site_id']]
                sample_df.loc[mask, k] = v
    sample_df['date'] = pd.to_datetime(sample_df['date'], format='%d/%m/%Y')
    # absolutely necessary if going between different grids with same pst obs
    # sort by postion before assigning idx
    sample_df.sort_values(['x', 'y', 'z'], inplace=True)
    sample_df.reset_index(inplace=True, drop=True)
    # mp7 name restricted to 16 char
    sample_df['obsnme'] = sample_df[['site_id', 'z']].apply(lambda x: x[0] + '_nan' if np.isnan(x[1]) \
        else x[0] + '_' + str(int(x[1])), axis=1)
    # last 12 chars leaving 4 for unique and '_'
    sample_df['obsnme'] = sample_df['obsnme'].apply(lambda x: x[-12:])
    # blind truncation in previous command may cause non-unique sites (in addition to non-unq from dates)
    # make unique obsnme fields
    sample_df['obsnme'] = sample_df['obsnme'] + '_' + sample_df.groupby('obsnme').cumcount(numeric_only=True).astype(str)
    sample_df.index = sample_df.obsnme
    sample_df.sort_index(inplace=True)
    if out_file:
        sample_df[['obsnme'] + out_cols].to_csv(os.path.join(fpth, out_file), index=False)
    return sample_df[out_cols]