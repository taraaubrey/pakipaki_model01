import os
import pyemu
import copy

def define_mult_array(
        pf, ws,
        tag='local1.recharge',
        ib=None,
        p_ins = [],
        lb=0.2, ub=5.0,
        ulb=0.1, uub=10.,
        lays=[0, 1, 2, 3, 4],
        ):
    
    files = [f for f in os.listdir(ws) if tag in f.lower() and f.endswith(".txt")]
    
    for i, f in enumerate(files):
        for ins in p_ins:
            par_type = ins.get('type', 'pilotpoints')
            gs = ins.get('gs', None)
            pp_space = ins.get('pp_space', None)
            suffix = ins.get('name_suffix', '')
            zone_array = ins.get('zone_array', ib[0])

            if i in lays:
                if isinstance(f,str):
                    iname = f.split(".")[1].replace("_","")
                else:
                    iname = f[0].split(".")[1]
            
                if par_type == 'pilotpoints':
                    pf.add_parameters(
                        f,
                        zone_array=zone_array,
                        par_type="pilotpoints",
                        par_name_base=iname+suffix,
                        geostruct=gs,
                        pargp=iname+suffix, 
                        lower_bound=lb, upper_bound=ub,
                        ult_lbound=ulb, ult_ubound=uub,
                        transform="log",
                        par_style="multiplier",
                        pp_options= {
                            "pp_space": pp_space, # specify the spacing of the pilot points
                            }
                        )
            
                if par_type == 'grid':
                    # constant (coarse) scale parameters
                    pf.add_parameters(
                        f,
                        zone_array=zone_array,
                        par_type="grid",
                        par_name_base=iname+suffix,
                        geostruct=gs,
                        pargp=iname+suffix,
                        lower_bound=lb, upper_bound=ub,
                        ult_lbound=ulb, ult_ubound=uub,
                        transform="log",
                        par_style="multiplier",
                        )
                
                if par_type == 'constant':
                    # constant scale parameters
                    pf.add_parameters(
                        f,
                        zone_array=zone_array,
                        par_type="constant",
                        par_name_base=iname+suffix,
                        pargp=iname+suffix,
                        lower_bound=lb, upper_bound=ub,
                        ult_lbound=ulb, ult_ubound=uub,
                        transform="log",
                        par_style="multiplier",
                        )
            
    return


def rcha(
        pf, ws,
        tag='local1.recharge',
        ib=None,
        p_ins = [],
        lb=0.2, ub=5.0,
        ulb=0.1, uub=10.,
        ):
    ss_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("0.txt") or f.endswith("2.txt"))]
    tr_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("1.txt") or f.endswith("3.txt"))]

    ins_files = {
        'rchss': ss_files,
        'rchtr': tr_files,
    }

    for name, files in ins_files.items():
        for ins in p_ins:
            par_type = ins.get('type', 'pilotpoints')
            gs = ins.get('gs', None)
            pp_space = ins.get('pp_space', None)
            suffix = ins.get('name_suffix', '')
            zone_array = ins.get('zone_array', ib[0])

            if par_type == 'pilotpoints':
                pf.add_parameters(
                    files,
                    zone_array=zone_array,
                    par_type="pilotpoints",
                    par_name_base=name+suffix,
                    geostruct=gs,
                    pargp=name+suffix, 
                    lower_bound=lb, upper_bound=ub,
                    ult_lbound=ulb, ult_ubound=uub,
                    transform="log",
                    par_style="multiplier",
                    pp_options= {
                        "pp_space": pp_space, # specify the spacing of the pilot points
                        }
                    )
        
            if par_type == 'grid':
                # constant (coarse) scale parameters
                pf.add_parameters(
                    files,
                    zone_array=zone_array,
                    par_type="grid",
                    par_name_base=name+suffix,
                    geostruct=gs,
                    pargp=name+suffix,
                    lower_bound=lb, upper_bound=ub,
                    ult_lbound=ulb, ult_ubound=uub,
                    transform="log",
                    par_style="multiplier",
                    )
            
            if par_type == 'constant':
                # constant scale parameters
                pf.add_parameters(
                    files,
                    zone_array=zone_array,
                    par_type="constant",
                    par_name_base=name+suffix,
                    pargp=name+suffix,
                    lower_bound=lb, upper_bound=ub,
                    ult_lbound=ulb, ult_ubound=uub,
                    transform="log",
                    par_style="multiplier",
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

def ghb_heads(
        pf, ws,
        name='ghb', 
        tag='local1.drn_stress_period_data',
        head_gs=None,
        head_bounds=[None, None],
        ult_bounds=[None, None],
        ):
    
    kper_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("0.txt") or f.endswith("2.txt"))]

    name = name.replace('-', '')
    # head but only in the kper 1, and 3 (ss) 
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
                        ult_lbound=ult_bounds[0],
                        ult_ubound=ult_bounds[1],
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
            ult_lbound=ult_bounds[0],
            ult_ubound=ult_bounds[1],
            )

        return

def ghb_cond(
        pf, ws,
        name='ghb', 
        tag='local1.drn_stress_period_data',
        grid_gs=None,
        cond_bounds=[None, None], 
        ult_bounds=[None, None],
        ):
    
    if tag.split('.')[-1].startswith('ghbspr'):
        kper_files = [f for f in os.listdir(ws) if tag in f.lower() and (f.endswith("0.txt") or f.endswith("1.txt"))]
    else:
        kper_files = [f for f in os.listdir(ws) if tag in f.lower()]
    name = name.replace('-', '')
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
                        ult_lbound=ult_bounds[0],
                        ult_ubound=ult_bounds[1],
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
            ult_lbound=ult_bounds[0],
            ult_ubound=ult_bounds[1],
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
