import yaml, argparse, sys, SLOPpy, os
from joblib import Parallel, delayed


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

    # print(config['nights'].keys())

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
    

    Parallel(n_jobs=len(config['nights'].keys()))(delayed(launch_sloppy_thread)(night,config,file_conf) for night in config['nights'].keys())
    
    SLOPpy.sloppy_run(file_conf=config)#launching SLOPpy on the original config file (maybe to run the combined fit?)

    sys.stdout = orig_stdout