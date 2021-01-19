import subprocess
import sys
from joblib import Parallel, delayed
import joblib

import numpy as np
import pandas as pd
from sklearn.model_selection import ParameterGrid


def save(filename,obj):
    with open(filename,'wb') as f:
        joblib.dump(obj,f)
        
def load(filename):
    with open(filename,'rb') as f:
        obj = joblib.load(f)
    return obj


def run_sncf(A,Z,c0=0,c1=0,w0=0):
    args = ['./sncf.x',str(A),str(Z),str(c0),str(c1),str(w0)]
    #outfile = open('out.out','w'); errfile = open('err.out','w')
    #result = subprocess.run(args,check=True,stdout=outfile,stderr=errfile)
    process = subprocess.run(args,check=False,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )
    #process.stdout.close()
    print (process)
    output = []; line = ""
    with open('temp.out') as f:
        line = f.readlines()[0]
        #line = line.strip()
    line = [ float(x) for x in line.strip().split() ]
    #print (line)
    output.append(line[2]) 
    output.append(line[3])
    #f.close(); outfile.close(); errfile.close()
    # energy, ch. radius
    return output

"""
def run_cycle(edf='nnlosat_256',c0=0,c1=0,w0=0):
    args = ['./my_script.sh',str(edf),str(c0),str(c1),str(w0)]
    result = subprocess.run(args,check=True)
    return result
"""


class EDF(object):

    def __init__(self,c0=0,c1=0,w0=0,edf='nnlosat_256'):
        self.c0=c0
        self.c1=c1
        self.w0=w0
        # todo check if file exists
        self.edf = edf
        filename = "Interactions/{ee}.in".format(ee=edf)
        args = ['cp',filename, 'interaction.in']
        result = subprocess.run(args,check=True)
        self.sncf = lambda a,z: run_sncf(a,z,self.c0,self.c1,self.w0)

    def __call__(self,az):
        az = np.array(az)
        if az.ndim == 1:
            return self.sncf(az[0],az[1])
        elif az.ndim > 1:
            n_out = len(self.sncf(az[0,0],az[0,1]))
            output = np.zeros( shape=(az.shape[0],n_out) )
            for k in range(output.shape[0]):
                output[k,:] = self.sncf(az[k,0],az[k,1])
            return output

    """
    y: (n_samples, n_features)
    Now: n_features=2, (E, rch)
    Return:
    
    """
    def _loss(self,y_obs,y_pred):
        loss = np.zeros(y_obs.shape[1])
        loss[0] = np.sum( (y_obs[:,0] - y_pred[:,0])**2 )
        # Radii when rch(exp) is known
        loss[1] = np.where(y_obs[:,1]>0, (y_obs[:,1] - y_pred[:,1])**2, 0. ).sum()
        loss[0] = np.sqrt( loss[0]/y_obs.shape[0] )
        loss[1] = np.sqrt( loss[1]/y_obs[ y_obs[:,1]>0 ].shape[0] )
        return loss


def evaluate(inputs,exp,params):
    print (params)
    edf = EDF(**params)
    y = edf(inputs)
    error = edf._loss(exp,y)
    #print ("Error: {e:.3f}\t{r:.3f}".format(e=error[0],r=error[1]))
    return error

def full_run(inputs,exp,grid,n_jobs=-1):
    return Parallel(n_jobs=n_jobs, verbose=10)(delayed(evaluate)(inputs,exp,params) for params in grid) 


def load_dataset():
    # Read data
    df = pd.read_csv('tab_exp.dat', index_col=0,delim_whitespace=True,engine='python' )
    mask = [ st.startswith('#')==False for st in df.index.tolist() ]
    df = df[mask]
    inputs = df.loc[:, ['A','Z'] ].to_numpy()
    exp = df.loc[:, 'Energy':].to_numpy()
    return df


if __name__=="__main__":
    
    # Data sets
    df_full = load_dataset()
    df_o16 = df_full.loc[['16o','40ca','48ca','90zr','132sn','208pb'],:]
    df_sami = df_o16.drop(['16o'])
    df_sami.loc['132sn','R(ch)'] = 0.
    print ("Sami: ", df_sami)
    dfs = (df_sami,df_o16,df_full)
    names = ("sami","o16","full")

    # Grid
    param_grid = {'c0': np.arange(0,-45,-5), 'c1': np.arange(0,45,5), 'w0': np.arange(30,150,10) }
    #param_grid = {'c0': (-28, -26,-25,-24,-22), 'c1': np.arange(0,30,5), 'w0': np.arange(50,140,10) }
    param_grid = ParameterGrid(param_grid)
    
    n_jobs=6
    grids = [ list(param_grid)[i:i +n_jobs] for i in range(0, len(param_grid), n_jobs)]
 
    
    # Loop over datasets
    for name, df in zip(names,dfs):
        # Input and target
        inputs = df_sami.loc[:, ['A','Z'] ].to_numpy()
        exp = df_sami.loc[:, 'Energy':].to_numpy()
        filename = 'results_{}'.format(name)

        results = []
        for grid in grids:
            res = full_run(inputs,exp,grid,n_jobs=n_jobs)
            results.extend(res)
            #pass
        print ("\n\nDone {}!\n\n".format(name))
        save('{}.pkl'.format(filename), results)

        results = load('{}.pkl'.format(filename) )
        results = np.array(results)
        print (results.shape)


        df = pd.DataFrame.from_dict(param_grid)
        df['err_e'] = results[:,0]
        df['err_r'] = results[:,1]
        print (df)
        save('{}.pkl'.format(filename), df)
        df.to_csv('grid_search_{}.csv'.format(name),index=False) 
