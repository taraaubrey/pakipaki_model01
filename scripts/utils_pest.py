import os
import pyemu
import copy

def define_mult_array(pf, ws,
          tag='local1.recharge',
          ib=None,
          sr=None,
          grid_gs=None,
          lb=0.2, ub=5.0,
          add_coarse=True,
          lays=[0, 1, 2, 3, 4],
          pp_space=10,
          fine=True,
          fine_gs=None,
          ):
    
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    
    for i, f in enumerate(files):
        if i in lays:
            if isinstance(f,str):
                base = f.split(".")[1].replace("_","")
            else:
                base = f[0].split(".")[1]
            
            pf.add_parameters(
                f,
                zone_array=ib[0],
                par_type="pilotpoints", #specify the type, these will be unique parameters for each cell
                geostruct=grid_gs, # the gestatisical structure for spatial correlation 
                par_name_base=base+"-pp", #specify a parameter name base that allows us to easily identify the filename and parameter type. "_gr" for "grid", and so forth.
                pargp=base+"pp", #likewise for the parameter group name
                lower_bound=lb, upper_bound=ub, #parameter lower and upper bound
                # ult_ubound=uub, ult_lbound=ulb,
                pp_options= {
                    "pp_space": pp_space, # specify the spacing of the pilot points
                    }
                )
            
            if add_coarse==True:
                # constant (coarse) scale parameters
                pf.add_parameters(f,
                                    zone_array=ib[0],
                                    par_type="constant",
                                    geostruct=grid_gs,
                                    par_name_base=base+"-cn",
                                    pargp=base+"cn",
                                    lower_bound=lb, upper_bound=ub,
                                    # ult_ubound=uub, ult_lbound=ulb
                                    )
    return


def rcha(
        pf, ws, ib,
        name='rch', 
        tag='local1.wel_stress_period_data', 
        fine_gs=None,
        grid_gs=None,
        constant_gs=None,
        rch_bounds=[0.1, 10],
        grid=True,
        constant=True,
        pp_space=10,
        ):
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    # for f in files:
    if grid:
        if fine_gs:
            pf.add_parameters(
                files,
                zone_array=ib[0],
                par_type="grid",
                geostruct=fine_gs,
                par_name_base=name+"-fngr",
                pargp=name+"fngr",
                lower_bound=rch_bounds[0],
                upper_bound=rch_bounds[1],
                )
        pf.add_parameters(
            files,
            zone_array=ib[0],
            par_type="pilotpoints",
            geostruct=grid_gs,
            par_name_base=name+"-pp",
            pargp=name+"pp",
            lower_bound=rch_bounds[0],
            upper_bound=rch_bounds[1],
            pp_options= {
                "pp_space": pp_space, # specify the spacing of the pilot points
                }
        )
    if constant:
        pf.add_parameters(
            files,
            zone_array=ib[0],
            par_type="constant",
            geostruct=constant_gs,
            par_name_base=name+"-cn",
            pargp=name+"cn",
            lower_bound=rch_bounds[0],
            upper_bound=rch_bounds[1],
        )
    return


def wel(
        pf, ws, 
        name='wel', 
        tag='local1.wel_stress_period_data', 
        fine_gs=None,
        grid_gs=None,
        constant_gs=None,
        q_bounds=[0.1, 10],
        grid=True,
        constant=True,
        ):
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    # for f in files:
    if grid:
        if fine_gs:
            pf.add_parameters(
                files,
                par_type="grid",
                geostruct=fine_gs,
                par_name_base=name+"-fngr",
                pargp=name+"fngr",
                index_cols=[0,1,2],
                use_cols=[3],
                lower_bound=q_bounds[0],
                upper_bound=q_bounds[1],
                )
        pf.add_parameters(
            files,
            par_type="grid",
            geostruct=grid_gs,
            par_name_base=name+"-crgr",
            pargp=name+"crgr",
            index_cols=[0,1,2],
            use_cols=[3],
            lower_bound=q_bounds[0],
            upper_bound=q_bounds[1],
        )
    if constant:
        pf.add_parameters(
            files,
            par_type="constant",
            geostruct=constant_gs,
            par_name_base=name+"-cn",
            pargp=name+"cn",
            index_cols=[0,1,2],
            use_cols=[3],  
            lower_bound=q_bounds[0],
            upper_bound=q_bounds[1],
        )
    return

def ghb(
        pf, ws, 
        name='drn', 
        tag='local1.drn_stress_period_data',
        fine_gs=None, 
        grid_gs=None, 
        constant_gs=None,
        cond_bounds=[0.1, 10], 
        head_bounds=[32.5, 42], 
        grid=True,
        constant=True,
        ):
    
    if tag.split('.')[-1].startswith('ghbspr'):
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
        kper_files = kper_files[:-2]  # exclude last two stress periods (no ghb)
    else:
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]

    # conductances all the same for each stress period
    if grid:
        if fine_gs:
            iname = name + f'-cond-fngr'
            # every cell scale (fine scale)
            pf.add_parameters(
                kper_files,
                par_type="grid",
                geostruct=fine_gs,
                par_name_base=iname,
                pargp=iname,
                index_cols=[0,1,2],
                use_cols=[4],
                lower_bound=cond_bounds[0],
                upper_bound=cond_bounds[1],
                )
        # every cell (larger scale)
        iname = name + f'-cond-crgr'
        pf.add_parameters(
            kper_files,
            par_type="grid",
            geostruct=grid_gs,
            par_name_base=iname,
            pargp=iname,
            index_cols=[0,1,2],
            use_cols=[4],
            lower_bound=cond_bounds[0],
            upper_bound=cond_bounds[1],
            )
    if constant:
        # constant and grid scale additive head parameters
        pf.add_parameters(
            kper_files,
            par_type="constant",
            geostruct=constant_gs,
            par_name_base=name+"-cond-cn",
            pargp=name+"-cond-cn",
            index_cols=[0,1,2],
            use_cols=[4],  
            lower_bound=cond_bounds[0],
            upper_bound=cond_bounds[1],
            )


    # for f in kper_files:
    #     kper = f.split(".")[1][-1]
    #     if kper not in ['1']:

    #         if grid:
    #             if fine_gs:
    #                 iname = name + f'-head-fngr'
    #                 pf.add_parameters(
    #                     f,
    #                     par_type="grid",
    #                     geostruct=fine_gs,
    #                     par_name_base=iname,
    #                     pargp=iname,
    #                     index_cols=[0,1,2],
    #                     use_cols=[3],   
    #                     par_style="a", 
    #                     transform="none",
    #                     lower_bound=head_bounds[0],
    #                     upper_bound=head_bounds[1],
    #                     )
    #             iname = name + f'-head-crgr'
    #             pf.add_parameters(
    #                 f,
    #                 par_type="grid",
    #                 geostruct=grid_gs,
    #                 par_name_base=iname,
    #                 pargp=iname,
    #                 index_cols=[0,1,2],
    #                 use_cols=[3],   
    #                 par_style="a", 
    #                 transform="none",
    #                 lower_bound=head_bounds[0],
    #                 upper_bound=head_bounds[1],
    #                 )
                
    #         if constant:
    #         # constant and grid scale additive head parameters
    #             pf.add_parameters(f,
    #                                 par_type="constant",
    #                                 geostruct=constant_gs,
    #                                 par_name_base=name+"-head-cn",
    #                                 pargp=name+"-head-cn",
    #                                 index_cols=[0,1,2],
    #                                 use_cols=[3],
    #                                 par_style="a", 
    #                                 transform="none",
    #                                 lower_bound=head_bounds[0],
    #                                 upper_bound=head_bounds[1],
    #                                 )
    return


def add_ts_parameters(pf, name, TEMP_DIR, f, lb, ub):
    import numpy as np
    import pandas as pd
    
    df = pd.read_csv(os.path.join(TEMP_DIR, f))
    pf.add_parameters(
        f,
        par_type="constant",
        par_name_base=df.columns.tolist()[1:],
        pargp=name+'-cn',
        index_cols=[0],
        use_cols=np.arange(1, len(df.columns)).tolist(),
        par_style="a", 
        transform="none",
        lower_bound=lb,
        upper_bound=ub,
        )

# def chd(pf, ws, name='chd',
    #     tag='local1.chd_stress_period_data',
    #     grid_gs=None,
    #     head_bounds=[-2, 2],
    #     head_ultbounds=[None, None]):
    # files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    # for f in files:
    #     # constant and grid scale additive head parameters
    #     name = name + '_head'
    #     pf.add_parameters(f,
    #                         par_type="grid",
    #                         geostruct=grid_gs,
    #                         par_name_base=name+"gr",
    #                         pargp=name+"gr",
    #                         index_cols=[0,1,2],
    #                         use_cols=[3],   # column containing head values
    #                         par_style="a", # specify additive parameter
    #                         transform="none", # specify not log-transform
    #                         lower_bound=head_bounds[0],
    #                         upper_bound=head_bounds[1],
    #                         # ult_lbound=head_ultbounds[0],
    #                         # ult_ubound=head_ultbounds[1]
    #                         )
    #     pf.add_parameters(f,
    #                         par_type="constant",
    #                         geostruct=grid_gs,
    #                         par_name_base=name+"cn",
    #                         pargp=name+"cn",
    #                         index_cols=[0,1,2],
    #                         use_cols=[3],
    #                         par_style="a", 
    #                         transform="none",
    #                         lower_bound=head_bounds[0],
    #                         upper_bound=head_bounds[1],
    #                         # ult_lbound=head_ultbounds[0],
    #                         # ult_ubound=head_ultbounds[1]
    #                         )
    # return