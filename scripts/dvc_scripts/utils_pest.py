import os
import pyemu

def define_mult_array(pf, ws,
          tag='local1.recharge',
          ib=None,
          sr=None,
          grid_gs=None,
          lb=0.2, ub=5.0,
          ulb=0.01, uub=100,
          add_coarse=True,
          lays=[0, 1, 2, 3, 4],
          pp_space=10):
    
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    
    for i, f in enumerate(files):
        if i in lays:
            if isinstance(f,str):
                base = f.split(".")[1].replace("_","")
            else:
                base = f[0].split(".")[1]
            
            # # grid (fine) scale parameters
            # pf.add_parameters(
            #     f,
            #     zone_array=ib[i],
            #     par_type="grid", #specify the type, these will be unique parameters for each cell
            #     geostruct=grid_gs, # the gestatisical structure for spatial correlation 
            #     par_name_base=base+"gr", #specify a parameter name base that allows us to easily identify the filename and parameter type. "_gr" for "grid", and so forth.
            #     pargp=base+"gr", #likewise for the parameter group name
            #     lower_bound=lb, upper_bound=ub, #parameter lower and upper bound
            #     ult_ubound=uub, ult_lbound=ulb # The ultimate bounds for multiplied model input values. Here we are stating that, after accounting for all multipliers, Kh cannot exceed these values. Very important with multipliers
            #     )
                            
            # # pilot point (medium) scale parameters
            # pargp=base
            # prefix_dict = {i:[pargp]}
            # df_pp = pyemu.pp_utils.setup_pilotpoints_grid(
            #     sr=sr,  # model spatial reference
            #     ibound=ib, # to which cells to setup ppoints
            #     prefix_dict=prefix_dict, #prefix to add to parameter names
            #     pp_dir=ws, 
            #     tpl_dir=ws, 
            #     every_n_cell=5) # pilot point spacing
            # pp_file = os.path.join(ws,f"{pargp}pp.dat")
            
            pf.add_parameters(
                f,
                zone_array=ib[i],
                par_type="pilotpoints", #specify the type, these will be unique parameters for each cell
                geostruct=grid_gs, # the gestatisical structure for spatial correlation 
                par_name_base=base+"pp", #specify a parameter name base that allows us to easily identify the filename and parameter type. "_gr" for "grid", and so forth.
                pargp=base+"pp", #likewise for the parameter group name
                lower_bound=lb, upper_bound=ub, #parameter lower and upper bound
                ult_ubound=uub, ult_lbound=ulb,
                pp_options= {
                    "pp_space": pp_space, # specify the spacing of the pilot points
                    }
                )

            # # write the pilot point file
            # pst = pf.build_pst()
            # par_pp = pst.add_parameters(pp_file+".tpl")
            # pst.parameter_data.loc[par_pp.parnme, ['parval1','parlbnd','parubnd', 'pargp']] = parval1, lb, ub, pargp
            
            if add_coarse==True:
                # constant (coarse) scale parameters
                pf.add_parameters(f,
                                    zone_array=ib[i],
                                    par_type="constant",
                                    geostruct=grid_gs,
                                    par_name_base=base+"cn",
                                    pargp=base+"cn",
                                    lower_bound=lb, upper_bound=ub,
                                    ult_ubound=uub, ult_lbound=ulb)
    return


def wel(pf, ws, name='wel', tag='local1.wel_stress_period_data', grid_gs=None,
        q_bounds=[0.1, 10], q_ultbounds=[0.01, 10]):
    name = name + '_cond'
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    for f in files:
        name = name + '_cond'
        pf.add_parameters(f,
                            par_type="grid",
                            geostruct=grid_gs,
                            par_name_base=name+"gr",
                            pargp=name+"gr",
                            index_cols=[0,1,2], #column containing lay,row,col
                            use_cols=[3], #column containing conductance values
                            lower_bound=q_bounds[0],
                            upper_bound=q_bounds[1],
                            ult_lbound=q_ultbounds[0],
                            ult_ubound=q_ultbounds[1])
        pf.add_parameters(f,
                            par_type="constant",
                            geostruct=grid_gs,
                            par_name_base=name+"cn",
                            pargp=name+"cn",
                            index_cols=[0,1,2],
                            use_cols=[3],  
                            lower_bound=q_bounds[0],
                            upper_bound=q_bounds[1],
                            ult_lbound=q_ultbounds[0],
                            ult_ubound=q_ultbounds[1])
    return

def drn(pf, ws, name='drn', tag='local1.drn_stress_period_data', grid_gs=None, cond_bounds=[0.1, 10], cond_ultbounds=[0.01, 100], head_bounds=[32.5, 42], head_ultbounds=[None, None]):
    name = name + '_cond'
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    for f in files:
        pf.add_parameters(f,
                            par_type="grid",
                            geostruct=grid_gs,
                            par_name_base=name+"gr",
                            pargp=name+"gr",
                            index_cols=[0,1,2], #column containing lay,row,col
                            use_cols=[4], #column containing conductance values
                            lower_bound=cond_bounds[0],
                            upper_bound=cond_bounds[1],
                            ult_lbound=cond_ultbounds[0],
                            ult_ubound=cond_ultbounds[1])
        pf.add_parameters(f,
                            par_type="constant",
                            geostruct=grid_gs,
                            par_name_base=name+"cn",
                            pargp=name+"cn",
                            index_cols=[0,1,2],
                            use_cols=[4],  
                            lower_bound=cond_bounds[0],
                            upper_bound=cond_bounds[1],
                            ult_lbound=cond_ultbounds[0],
                            ult_ubound=cond_ultbounds[1])

        # constant and grid scale additive head parameters
        name = name + '_head'
        pf.add_parameters(f,
                            par_type="grid",
                            geostruct=grid_gs,
                            par_name_base=name+"gr",
                            pargp=name+"gr",
                            index_cols=[0,1,2],
                            use_cols=[3],   # column containing head values
                            par_style="a", # specify additive parameter
                            transform="none", # specify not log-transform
                            lower_bound=head_bounds[0],
                            upper_bound=head_bounds[1],
                            ult_lbound=head_ultbounds[0],
                            ult_ubound=head_ultbounds[1])
        pf.add_parameters(f,
                            par_type="constant",
                            geostruct=grid_gs,
                            par_name_base=name+"cn",
                            pargp=name+"cn",
                            index_cols=[0,1,2],
                            use_cols=[3],
                            par_style="a", 
                            transform="none",
                            lower_bound=head_bounds[0],
                            upper_bound=head_bounds[1],
                            ult_lbound=head_ultbounds[0],
                            ult_ubound=head_ultbounds[1])
    return

def chd(pf, ws, name='chd',
        tag='local1.chd_stress_period_data',
        grid_gs=None,
        head_bounds=[-2, 2],
        head_ultbounds=[None, None]):
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    for f in files:
        # constant and grid scale additive head parameters
        name = name + '_head'
        pf.add_parameters(f,
                            par_type="grid",
                            geostruct=grid_gs,
                            par_name_base=name+"gr",
                            pargp=name+"gr",
                            index_cols=[0,1,2],
                            use_cols=[3],   # column containing head values
                            par_style="a", # specify additive parameter
                            transform="none", # specify not log-transform
                            lower_bound=head_bounds[0],
                            upper_bound=head_bounds[1],
                            ult_lbound=head_ultbounds[0],
                            ult_ubound=head_ultbounds[1])
        pf.add_parameters(f,
                            par_type="constant",
                            geostruct=grid_gs,
                            par_name_base=name+"cn",
                            pargp=name+"cn",
                            index_cols=[0,1,2],
                            use_cols=[3],
                            par_style="a", 
                            transform="none",
                            lower_bound=head_bounds[0],
                            upper_bound=head_bounds[1],
                            ult_lbound=head_ultbounds[0],
                            ult_ubound=head_ultbounds[1])
    return