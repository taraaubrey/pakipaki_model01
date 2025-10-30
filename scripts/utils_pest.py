import os
import pyemu
import copy

def define_mult_array(
        pf, ws,
        tag='local1.recharge',
        ib=None,
        grid_gs=None,
        grid_suffix='-fngr',
        lb=0.2, ub=5.0,
        lays=[0, 1, 2, 3, 4],
        pp_space=None,
        pp_gs=None,
        ):
    
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    
    for i, f in enumerate(files):
        if i in lays:
            if isinstance(f,str):
                iname = f.split(".")[1].replace("_","")
            else:
                iname = f[0].split(".")[1]
            
            if pp_space:
                pf.add_parameters(
                    f,
                    zone_array=ib[0],
                    par_type="pilotpoints",
                    par_name_base=iname+'-pp',
                    geostruct=pp_gs,
                    pargp=iname+'-pp', 
                    lower_bound=lb, upper_bound=ub,
                    transform="log",
                    par_style="multiplier",
                    pp_options= {
                        "pp_space": pp_space, # specify the spacing of the pilot points
                        }
                    )
            
            if grid_gs:
                # constant (coarse) scale parameters
                pf.add_parameters(
                    f,
                    zone_array=ib[0],
                    par_type="grid",
                    par_name_base=iname+grid_suffix,
                    geostruct=grid_gs,
                    pargp=iname+grid_suffix,
                    lower_bound=lb, upper_bound=ub,
                    transform="log",
                    par_style="multiplier",
                    )
            
            # constant scale parameters
            pf.add_parameters(
                f,
                zone_array=ib[0],
                par_type="constant",
                par_name_base=iname+'-cn',
                pargp=iname+'-cn',
                lower_bound=lb, upper_bound=ub,
                transform="log",
                par_style="multiplier",
                )
            
    return


def rcha(
        pf, ws, ib,
        name='rch', 
        tag='local1.wel_stress_period_data', 
        grid_gs=None,
        grid_suffix='-fngn',
        rch_bounds=[0.1, 10],
        pp_gs=None,
        pp_space=None,
        time_gs = None,
        datetimes=None,
        ):
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    # for f in files:
    if grid_gs:
        pf.add_parameters(
            files,
            zone_array=ib[0],
            par_type="grid",
            geostruct=grid_gs,
            par_name_base=name + grid_suffix,
            pargp=name + grid_suffix,
            lower_bound=rch_bounds[0],
            upper_bound=rch_bounds[1],
            )
    if pp_space:
        pf.add_parameters(
            files,
            zone_array=ib[0],
            par_type="pilotpoints",
            geostruct=pp_gs,
            par_name_base=name+"-pp",
            pargp=name+"-pp",
            lower_bound=rch_bounds[0],
            upper_bound=rch_bounds[1],
            pp_options= {
                "pp_space": pp_space, # specify the spacing of the pilot points
                }
        )

    pf.add_parameters(
        files,
        zone_array=ib[0],
        par_type="constant",
        par_name_base=name+"-cn",
        pargp=name+"-cn",
        lower_bound=rch_bounds[0],
        upper_bound=rch_bounds[1],
    )

    if time_gs:
        pf.add_parameters(
            files,
            zone_array=ib[0],
            par_type="constant",
            par_name_base=name+"-tcn",
            pargp=name+"-tcn",
            lower_bound=rch_bounds[0],
            upper_bound=rch_bounds[1],
            datetime=datetimes,
            geostruct=time_gs,
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
        name='ghb', 
        tag='local1.drn_stress_period_data',
        grid_gs=None, 
        head_gs=None,
        cond_bounds=None, 
        head_bounds=None,
        ):
    
    if tag.split('.')[-1].startswith('ghbspr'):
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith("0.txt")]
    else:
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("0.txt") or f.endswith("2.txt"))]

    # conductances all the same for each stress period
    if cond_bounds:
        if grid_gs:
            if isinstance(grid_gs, dict):
                for suffix, gs in grid_gs.items():
                    iname = name + f'-cond-{suffix}'
                    # every cell scale (fine scale)
                    pf.add_parameters(
                        kper_files,
                        par_type="grid",
                        geostruct=gs,
                        par_name_base=iname,
                        pargp=iname,
                        index_cols={'k':0, 'i':1, 'j':2},
                        use_cols=[4],
                        lower_bound=cond_bounds[0],
                        upper_bound=cond_bounds[1],
                        )

        pf.add_parameters(
            kper_files,
            par_type="constant",
            par_name_base=name+"-cond-cn",
            pargp=name+"-cond-cn",
            index_cols={'k':0, 'i':1, 'j':2},
            use_cols=[4],  
            lower_bound=cond_bounds[0],
            upper_bound=cond_bounds[1],
            )

    # conductances all the same for each stress period
    if head_bounds:
        if head_gs:
            if isinstance(head_gs, dict):
                for suffix, gs in head_gs.items():
                    iname = name + f'-head-{suffix}'
                    # every cell scale (fine scale)
                    pf.add_parameters(
                        kper_files,
                        par_type="grid",
                        geostruct=gs,
                        par_name_base=iname,
                        pargp=iname,
                        index_cols={'k':0, 'i':1, 'j':2},
                        use_cols=[3],   
                        par_style="a", 
                        transform="none",
                        lower_bound=head_bounds[0],
                        upper_bound=head_bounds[1],
                        )

        pf.add_parameters(
            kper_files,
            par_type="constant",
            par_name_base=name+"-head-cn",
            pargp=name+"-head-cn",
            index_cols={'k':0, 'i':1, 'j':2}, # {'k':0, 'i':1, 'j':2}
            use_cols=[3],   
            par_style="a", 
            transform="none",
            lower_bound=head_bounds[0],
            upper_bound=head_bounds[1],
            )

    return


def conf_ghb(
        pf, ws,
        name='ghb', 
        tag='local1.drn_stress_period_data',
        cond_gs=None, 
        head_gs=None,
        cond_bounds=None, 
        head_bounds=None,
        pp_space=5,
        ):
    
    if tag.split('.')[-1].startswith('ghbspr'):
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith("0.txt")]
    else:
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("0.txt") or f.endswith("2.txt"))]

    # conductances all the same for each stress period
    if cond_bounds:
        pf.add_parameters(
            kper_files,
            par_type="grid",
            geostruct=cond_gs,
            par_name_base=name+"-cond-gr",
            pargp=name+"-cond-gr",
            index_cols={'k':0, 'i':1, 'j':2},
            use_cols=[4],
            lower_bound=cond_bounds[0],
            upper_bound=cond_bounds[1],
            )

        pf.add_parameters(
            kper_files,
            par_type="constant",
            par_name_base=name+"-cond-cn",
            pargp=name+"-cond-cn",
            index_cols={'k':0, 'i':1, 'j':2},
            use_cols=[4],  
            lower_bound=cond_bounds[0],
            upper_bound=cond_bounds[1],
            )

    # conductances all the same for each stress period
    if head_bounds:
        pf.add_parameters(
            kper_files,
            par_type="grid",
            par_name_base=name+"-head-gr",
            pargp=name+"-head-gr",
            index_cols={'k':0, 'i':1, 'j':2}, # {'k':0, 'i':1, 'j':2}
            use_cols=[3],   
            par_style="a", 
            transform="none",
            lower_bound=head_bounds[0],
            upper_bound=head_bounds[1],
            )

        pf.add_parameters(
            kper_files,
            par_type="constant",
            par_name_base=name+"-head-cn",
            pargp=name+"-head-cn",
            index_cols={'k':0, 'i':1, 'j':2}, # {'k':0, 'i':1, 'j':2}
            use_cols=[3],   
            par_style="a", 
            transform="none",
            lower_bound=head_bounds[0],
            upper_bound=head_bounds[1],
            )

    return


def add_ts_parameters(
        pf, ws,
        datetimes,
        name,
        f,
        bounds,
        parnames,
        use_cols,
        time_gs=None,
        ):
    import numpy as np
    import pandas as pd

    # constant
    iname = name + '-cn'

    pf.add_parameters(
        f,
        par_type="constant",
        par_name_base=parnames,
        pargp=[iname] * len(parnames),
        index_cols=0,
        use_cols=use_cols,
        # use_rows=use_rows_indexes,
        par_style="a",
        transform="none",
        lower_bound=bounds[0],
        upper_bound=bounds[1],
        geostruct=None,  # explicitly set geostruct for constant parameters
        )
    
    return
        
    # if time_gs:
    #     if isinstance(time_gs, dict):
    #         for suffix, gs in time_gs.items():
    #             iname = name + f'-{suffix}'
    #             # every cell scale (fine scale)
    #             pf.add_parameters(
    #                 f,
    #                 par_type="grid",
    #                 geostruct=gs,
    #                 par_name_base=parnames,
    #                 pargp=iname,
    #                 index_cols=0,
    #                 use_cols=use_cols,
    #                 par_style="a", 
    #                 transform="none",
    #                 lower_bound=bounds[0],
    #                 upper_bound=bounds[1],
    #                 )
