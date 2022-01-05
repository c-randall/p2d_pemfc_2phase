""" Import needed modules """
"-----------------------------------------------------------------------------"
import os
from Shared_Funcs.post_pak_dfs import *
from Shared_Funcs.post_pak_plots import *
from Shared_Funcs.read_and_write import *

""" Generate or load dataframes """
"-----------------------------------------------------------------------------"    
if post_only == 2:
    directory = os.getcwd()+'/Saved_Results/'+load_folder+'/Saved_dfs'
    df_i,df_p,df_f,df_y,df_r = loader(directory)
    save = 0
    
else:
    
    if debug: 
        df_t = df_debug(sol,ca,gdl,cl)
    if polar == 'no' and i_OCV == 0:
        i_find = 0

    df_f,df_p,df_y,df_r = df_tool(df_i,df_p,ca,gdl,cl,i_find,yamlfile) 

""" Save model files and outputs """
"-----------------------------------------------------------------------------" 
if save:
    cwd = os.getcwd()
    save_model(cwd,save_folder)
    
    save_dfs(cwd,save_folder,df_i,df_p,df_f,df_y,df_r)
    
""" Produce requested figures """
"-----------------------------------------------------------------------------"
if debug:
    debug_plts(df_t,ca,gdl,cl,debug,save,save_folder)
    
if grads:
    grad_plts(df_y,ca,gdl,cl,grads,save,save_folder)
    
if radial:
    radial_plts(df_r,ca,gdl,cl,radial,save,save_folder)
    
if polar and polar != 'no':
    polar_plt(df_p,ca,gdl,cl,polar,data,save,save_folder)
    
if over_p:
    over_plt(df_p,ca,gdl,cl,save,save_folder)
    
if i_ver:
    verification(df_i,ca,gdl,cl,gdl_cl,i_find)