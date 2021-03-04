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

"""
Flag: True->custom interaction; False->Skyrme 
"""
def run_sncf(A,Z,c0=0,c1=0,w0=0,flag=True):
    args = ['./sncf.x',str(A),str(Z),str(c0),str(c1),str(w0)] if flag else ['./sncf.x',str(A),str(Z)]
    process = subprocess.run(args,check=False,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )
    
    output = []; line = ""
    with open('temp.out') as f:
        line = f.readlines()[0]
        #line = line.strip()
    line = [ float(x) for x in line.strip().split() ]
    #print (line)
    output.append(line[2])       # energy 
    output.append(line[3])       # rch
    output.append(line[2]/float(A) )   # E/A
    #f.close(); outfile.close(); errfile.close()
    # energy, ch. radius
    return output



class EDF(object):

    def __init__(self,c0=0,c1=0,w0=0,edf='nnlosat_256'):
        self.c0=c0
        self.c1=c1
        self.w0=w0
        # todo check if file exists
        self.edf = ''.join(edf.split())
        #print (self.edf)
        if len(self.edf)>0:
            filename = "Interactions/{ee}.in".format(ee=self.edf)
            args = ['cp',filename, 'interaction.in']
            result = subprocess.run(args,check=True)
            self.flag = True
        else:
            with open('interaction.in','r') as f, open('sncf.in','r') as fi:
                l2 = str( f.readlines()[1] )
                if int(l2)<=0:
                    self.edf = fi.readlines()[1].strip()
                    self.flag = False
        print ("Using ",self.edf)
        # sncf
        self.sncf = lambda a,z: run_sncf(a,z,self.c0,self.c1,self.w0,self.flag)


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
    Now: n_features=3, (E, rch, E/A)
    Return:    
    """
    def _loss(self,y_obs,y_pred):
        loss = np.zeros(y_obs.shape[1])
        # Energy
        loss[0] = np.sum( (y_obs[:,0] - y_pred[:,0])**2 )
        loss[0] = np.sqrt( loss[0]/y_obs.shape[0] )
        # Radii when rch(exp) is known
        nr = y_obs[ y_obs[:,1]>0 ].shape[0]
        loss[1] = np.where(y_obs[:,1]>0, (y_obs[:,1] - y_pred[:,1])**2, 0. ).sum()
        loss[1] = np.sqrt( loss[1]/nr )
        # E/A
        loss[2] = np.sum( (y_obs[:,2] - y_pred[:,2])**2 )
        loss[2] = np.sqrt( loss[2]/y_obs.shape[0] )
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
    df['Energy_per_nucleon'] = df['Energy']/df['A']
    return df


def load_sami():
    df_full = load_dataset()
    df_o16 = df_full.loc[['16o','40ca','48ca','90zr','132sn','208pb'],:]
    df_sami = df_o16.drop(['16o'])
    df_sami.loc['132sn','R(ch)'] = 0.
    inputs = df_sami.loc[:, ['A','Z'] ].to_numpy()
    exp = df_sami.loc[:, 'Energy':].to_numpy()
    return (inputs,exp,df_sami)



if __name__=="__main__":
    
    # Data sets
    df_full = load_dataset()
    df_o16 = df_full.loc[['16o','40ca','48ca','90zr','132sn','208pb'],:]
    df_sami = df_o16.drop(['16o'])
    df_sami.loc['132sn','R(ch)'] = 0.
    #print ("Sami\n", df_sami)
    dfs = (df_sami,)   #(df_sami,df_o16,df_full)
    names = ("sami",)  # ("sami","o16","full")

    # Grid
    # param_grid = {'c0': np.arange(0,-45,-5), 'c1': np.arange(0,45,5), 'w0': np.arange(30,150,10) }
    param_grid = {'edf':('av4p_256',), 'c0': np.arange(0,-160,-5), 'c1': np.arange(0,50,5), 'w0': np.arange(0,150,10) }
    # param_grid = {'edf':['av4p_1256',], 'c0': [0.,], 'c1':[0.,], 'w0':[0.,] }
    param_grid = ParameterGrid(param_grid)
    print ("N. models:\t", len(list(param_grid)) ) 
    n_jobs= 24
    grids = [ list(param_grid)[i:i +n_jobs] for i in range(0, len(param_grid), n_jobs)]
 
    # Loop over datasets
    for name, df in zip(names,dfs):
        # Input and target
        inputs = df_sami.loc[:, ['A','Z'] ].to_numpy()
        exp = df_sami.loc[:, 'Energy':].to_numpy()
        filename = 'results_{}'.format(name)
        print ("Filename :", filename)
 
        st_point = 0
        grids = grids[st_point:]

        results = []
        for k, grid in enumerate(grids):
            res = full_run(inputs,exp,grid,n_jobs=n_jobs)
            results.extend(res)
            # save often (perhaps too often)
            save('{}.pkl'.format(filename), results)
            print ("Saved block {}".format(k+st_point) )
            #pass
        print ("\n\nDone {}!\n\n".format(name))
        

        results = load('{}.pkl'.format(filename) )
        results = np.array(results)
        print (results.shape)


        df = pd.DataFrame.from_dict(param_grid)
        df['err_e'] = results[:,0]
        df['err_r'] = results[:,1]
        df['err_e_a'] = results[:,2]
        print (df)
        save('{}.pkl'.format(filename), df)
        df.to_csv('grid_search_{}.csv'.format(name),index=False) 
