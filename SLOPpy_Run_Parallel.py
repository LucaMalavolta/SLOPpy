import yaml, argparse, sys, SLOPpy, os, glob
from joblib import Parallel, delayed
import numpy as np

if __name__ == '__main__':

#     SLOPpy.sloppy_run()

    sloppyLogFolder='./sloppyLogs'  

    if not os.path.exists(sloppyLogFolder): os.mkdir(sloppyLogFolder)


    parser = argparse.ArgumentParser(prog='SLOPpy_Run_Parallel', description='SLOPpy runner in parallel')
    parser.add_argument('config_file', type=str, nargs=1, help='config file')

    args = parser.parse_args()
    file_conf = args.config_file[0]

    with open(file_conf, 'r') as file:
        config = yaml.safe_load(file)
    file.close()

    orig_stdout = sys.stdout
    
    
    def launch_sloppy_thread(night,config,file_conf):#this function runs on the individual nights. it is parallelized by Parallel/delayed
        
        #building a configuration file with only one night to process
        night_config=config.copy()
        del night_config['nights']#all the nights removed from the dict
        
        night_config['nights']={}
        night_config['nights'][night]=config['nights'][night]#re-adding only one night to the dict
        
        night_config_file=night+'_'+file_conf
        with open(night_config_file, 'w') as file:
            yaml.dump(night_config, file)#writing the dict in yaml format

        f = open(sloppyLogFolder+'/'+night_config_file+'.log', 'w')#redirecting the standard output to a log file dedicated to the night
        sys.stdout = f
        
        SLOPpy.sloppy_run(file_conf=night_config_file)#launching SLOPpy on the night

        f.close()
        
        return()
    

    # Parallel(n_jobs=len(config['nights'].keys()))(delayed(launch_sloppy_thread)(night,config,file_conf) for night in config['nights'].keys())

    sys.stdout = orig_stdout
    
    #This is a workaround to make SLOPpy perform the combined analysis.
    #For each night, the parallelized analysis produces pickle files whose names are preceded by the night string: <night>+<pickle filename>
    #The following loop creates a symlink to the <night>+<pickle file name> with name <pickle file name>: this way the combined analysis is performed using the results from the parallelized night-wise analysis
    for night in config['nights'].keys():
        lista_all=np.sort(glob.glob(night+'*'+night+'*.p'))
        for l in lista_all:
            # print(l)
            # print(l[11:])
            if not os.path.islink(l[11:]):
                if ('doublet' in l[11:]) and (night in l[11:]): continue #to skip the combined analysis run on individual nights
                os.symlink(l,l[11:])
            
  
    # SLOPpy.sloppy_run(file_conf=config)#launching SLOPpy on the original config file (maybe to run the combined fit?)

    sys.stdout = orig_stdout