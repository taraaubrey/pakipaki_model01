import os
import multiprocessing as mp
import numpy as np
import pandas as pd
import pyemu
def main():

    pyemu.helpers.apply_list_and_array_pars(arr_par_file='mult2model_info.csv',chunk_len=50)

if __name__ == '__main__':
    mp.freeze_support()
    main()

