"""
Functions intended to read and write model outputs. Basic dataframes are saved
as csv files. Also, a copy of the model is generated.
"""

import os, sys
import numpy as np
import pandas as pd
import cantera as ct
import matplotlib.pyplot as plt
from shutil import copy2, rmtree, copytree, ignore_patterns

def loader(directory_loc):
    list_dir = os.listdir(directory_loc)
    
    df_i_file = [f for f in list_dir if 'df_i' in f]
    df_p_file = [f for f in list_dir if 'df_p' in f]
    df_f_file = [f for f in list_dir if 'df_f' in f]
    df_y_file = [f for f in list_dir if 'df_y' in f]
    df_r_file = [f for f in list_dir if 'df_r' in f]
    
    df_i_loc = directory_loc+'/'+df_i_file[0]
    df_p_loc = directory_loc+'/'+df_p_file[0]
    df_f_loc = directory_loc+'/'+df_f_file[0]
    df_y_loc = directory_loc+'/'+df_y_file[0]
    df_r_loc = directory_loc+'/'+df_r_file[0]
    
    df_i = pd.read_csv(df_i_loc)
    df_i.columns = df_i.columns.map(int)
    
    df_p = pd.read_csv(df_p_loc)
    df_f = pd.read_csv(df_f_loc)
    df_y = pd.read_csv(df_y_loc)
    df_r = pd.read_csv(df_r_loc)
    
    print('\nDataFrames loaded from load_folder')
    
    return df_i,df_p,df_f,df_y,df_r

def save_model(cwd,save_folder):
    
    try:
        os.chdir(cwd+'/Saved_Results')
    except:
        os.mkdir(cwd+'/Saved_Results')
        os.chdir(cwd+'/Saved_Results')
        
    if os.path.exists(save_folder):
        print('\nWARNING: save_folder already exists. Approve to overwrite.')
        print('\n"Enter" to continue or "Ctrl+c" to cancel.')
        
        user_in = input()   
        if user_in == KeyboardInterrupt:
            sys.exit(0)
        else:
            rmtree(save_folder)
            
    # Ignore subfolders
    ignore = ignore_patterns('__pycache__','Old_cti_version')
    
    os.mkdir(save_folder)
    copy2(cwd+'/pemfc_runner.py',save_folder)
    copytree(cwd+'/Inputs',save_folder+'/Inputs')
    copytree(cwd+'/Shared_Funcs',save_folder+'/Shared_Funcs',ignore=ignore)
        
    os.chdir(cwd)
    
    return print('\nModel saved...')

def save_dfs(cwd,save_folder,df_i,df_p,df_f,df_y,df_r):
    
    os.mkdir(cwd+'/Saved_Results/'+save_folder+'/Saved_dfs')
    os.chdir(cwd+'/Saved_Results/'+save_folder+'/Saved_dfs')    
    
    df_i.to_csv('df_i.csv',index=False)    
    df_p.to_csv('df_p.csv',index=False)
    df_f.to_csv('df_f.csv',index=False)
    df_y.to_csv('df_y.csv',index=False)
    df_r.to_csv('df_r.csv',index=False)
    
    os.chdir(cwd)
    
    return print('\nDataFrames saved...')

def save_plot(cwd,save_folder,fig_name):
    
    try: os.mkdir(cwd+'/Saved_Results/'+save_folder+'/Figures')
    except: None
    
    os.chdir(cwd+'/Saved_Results/'+save_folder+'/Figures')
    plt.savefig(fig_name)
    os.chdir(cwd)
    
    return None