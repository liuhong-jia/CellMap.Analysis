import scanpy as sc
import os
import pandas as pd
import numpy as np

from scipy.stats import pearsonr
from collections.abc import Iterable
from scipy.spatial.distance import jensenshannon
def jsd(v1,v2):
    return jensenshannon(v1-v2)
from sklearn.metrics import mean_squared_error
def rmse(x1,x2):
    return mean_squared_error(x1,x2,squared=False)

def ssim(im1,im2,M=1):
    im1, im2 = im1/im1.max(), im2/im2.max()
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1 = np.sqrt(((im1 - mu1) ** 2).mean())
    sigma2 = np.sqrt(((im2 - mu2) ** 2).mean())
    sigma12 = ((im1 - mu1) * (im2 - mu2)).mean()
    k1, k2, L = 0.01, 0.03, M
    C1 = (k1*L) ** 2
    C2 = (k2*L) ** 2
    C3 = C2/2
    l12 = (2*mu1*mu2 + C1)/(mu1 ** 2 + mu2 ** 2 + C1)
    c12 = (2*sigma1*sigma2 + C2)/(sigma1 ** 2 + sigma2 ** 2 + C2)
    s12 = (sigma12 + C3)/(sigma1*sigma2 + C3)
    ssim = l12 * c12 * s12
    return ssim

def str_remove_spe(string):
    s_list = []
    for s in string:
        if s.isalnum():
            s_list.append(s)
        else:
            s_list.append('_')
    new_string = ''.join(s_list)
    return new_string

def get_RCTD_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results
    
def get_tangram_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def get_seurat_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results = results.iloc[:,1:-1]
    results.columns = [str_remove_spe(c.split('prediction.score.')[1]) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results
    
def get_c2l_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c.split('q05cell_abundance_w_sf_')[1]) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    results = (results.T/results.sum(axis=1)).T
    return results

def get_stereo_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def get_spotlight_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results
    
def get_spatialdwls_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results = results.iloc[:,1:]
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def get_destvi_results(results_path,celltype=None):
    if (not os.path.exists(results_path)) and (celltype is not None):
        results = pd.DataFrame(np.zeros((1000,len(celltype))),index=range(1000),columns=celltype)
        results = results.loc[:,np.unique(results.columns)]
        print(f'not found: {results_path}')
    else: 
        results = pd.read_csv(results_path,index_col=0)
        results.columns = [str_remove_spe(c) for c in results.columns]
        if celltype is not None:
            results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
        results = results.loc[:,np.unique(results.columns)]
    return results
	
def get_stride_results(results_path,celltype=None):
    if (not os.path.exists(results_path)) and (celltype is not None):
        results = pd.DataFrame(np.zeros((1000,len(celltype))),index=range(1000),columns=celltype)
        results = results.loc[:,np.unique(results.columns)]
        print(f'not found: {results_path}')
    else:    
        results = pd.read_csv(results_path,index_col=0,sep='\t')
        results = results.fillna(0)
        results.columns = [str_remove_spe(c) for c in results.columns]
        if celltype is not None:
            results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
            results = results.loc[:,celltype]
        results = results.loc[:,np.unique(results.columns)]
    return results
	
def get_cellmap_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def get_redeconve_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results = results.T.reset_index(drop=True)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    results = (results.T/results.sum(axis=1)).T
    return results

def read_results(result_path, func, index, celltype):
    if os.path.exists(result_path):
        results = func(result_path, celltype=celltype)
    else:
        results = pd.DataFrame(0, index=index, columns=celltype)
    return results

def get_dataset_results(dataset_number):
    st = sc.read_h5ad(f'/mnt/disk1/hongjia/work/xiu/Nature_Methods/SimualtedSpatalData/dataset{dataset_number}/Spatial.h5ad')
    gd = st.uns['density']
    gd.columns = [str_remove_spe(c) for c in gd.columns]
    gd = (gd.T/gd.sum(axis=1)).T
    gd = gd.loc[:,np.unique(gd.columns)]
    index = gd.index
    celltype = gd.columns
    result_path = f'/mnt/disk1/hongjia/work/xiu/Nature_Methods/Simulated_results/dataset{dataset_number}'
    RCTD_results = read_results(os.path.join(result_path, 'RCTD_result.txt'), func=get_RCTD_results, index=index, celltype=celltype)
    c2l_results = read_results(os.path.join(result_path, 'Cell2location_result.txt'), func=get_c2l_results, index=index, celltype=celltype)
    tangram_results = read_results(os.path.join(result_path, 'Tangram_result.txt'), func=get_tangram_results, index=index, celltype=celltype)
    seurat_results = read_results(os.path.join(result_path, 'Seurat_result.txt'), func=get_seurat_results, index=index, celltype=celltype)
    stereo_results = read_results(os.path.join(result_path, 'Stereoscope_result.txt'), func=get_stereo_results, index=index, celltype=celltype)
    spotlight_results = read_results(os.path.join(result_path, 'SPOTlight_result.txt'), func=get_spotlight_results, index=index, celltype=celltype)
    spatialdwls_results = read_results(os.path.join(result_path,'SpatialDWLS_result.txt'), func=get_spatialdwls_results, index=index, celltype=celltype)
    destvi_results = read_results(os.path.join(result_path, 'DestVI_result.txt'), func=get_destvi_results, index=index, celltype=celltype)
    stride_results = read_results(os.path.join(result_path, 'STRIDE_result.txt'), func=get_stride_results, index=index, celltype=celltype)
    redeconve_results = read_results(os.path.join(result_path, 'Redeconve.result.txt'), func=get_redeconve_results, index=index, celltype=celltype)
    cellmap_results = read_results(os.path.join(result_path, 'CellMap_result.txt'), func=get_cellmap_results, index=index, celltype=celltype)
    return gd, [RCTD_results, c2l_results, tangram_results, seurat_results, stereo_results, spotlight_results, spatialdwls_results, destvi_results, stride_results, redeconve_results, cellmap_results]

def compare_results(gd,result_list,metric='pcc',columns=None,axis=1):
    if metric=='pcc':
        func = pearsonr
        r_ind = 0
    if metric=='jsd':
        func = jensenshannon
        r_ind = None
    if metric=='rmse':
        func = rmse
        r_ind = None
    if metric=='ssim':
        func = ssim
        r_ind = None
    if isinstance(result_list, pd.DataFrame):
        c_list = []
        if axis == 1:
            print('axis: ',1)
            for i,c in enumerate(gd.columns):
                r = func(gd.iloc[:,i].values, np.nan_to_num(np.clip(result_list.iloc[:,i],0,1)))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        else:
            print('axis: ',0)
            for i,c in enumerate(gd.index):
                r = func(gd.iloc[i,:].values, np.nan_to_num(np.clip(result_list.iloc[i,:],0,1)))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        df = pd.DataFrame(c_list,index=gd.columns,columns=columns)
    else:
        df_list = []
        for n,res in enumerate(result_list):
            # print('- '+str(n))
            c_list = []
            if axis == 1:
                for i,c in enumerate(gd.columns):
                    r = func(gd.iloc[:,i].values, np.nan_to_num(np.clip(res.iloc[:,i],0,1)))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list,index=gd.columns)
            else:
                for i,c in enumerate(gd.index):
                    r = func(gd.iloc[i,:].values, np.nan_to_num(np.clip(res.iloc[i,:],0,1)))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list,index=gd.index)
            df_list.append(df_tmp)
        df = pd.concat(df_list,axis=1)
        df.columns = columns
    return df

def cal_index(method='pcc'): 
    df_list = []
    for dataset_number in range(1,33):
        print(f'dataset number: {dataset_number}')
        # dataset_number = 10
        gd, results = get_dataset_results(dataset_number)
        df = compare_results(
            gd,
            results,
            columns = ['RCTD','Cell2location','Tangram','Seurat','Stereoscope','SPOTlight','SpatialDWLS','DestVI','STRIDE','Redeconve', 'CellMap'],
            axis=1,
            metric=method
        )
        df['dataset'] = dataset_number
        df_list.append(df)
    all_clusters = pd.concat(df_list,axis=0)
    return all_clusters

def get_as_score(Result):
    tools_num = Result.shape[0]
    Tools_score=[]
    methods = list(Result.index)
    score_col = []
    list_up = list(range(1,Result.shape[1]+1))
    list_down  = list(range(Result.shape[1],0,-1))
    
    for method in methods:
        if method == 'PCC' or method == 'SSIM':
            Tools_score.append(pd.Series(list_down, index=Result.loc[method,:].sort_values(ascending=False).index))
            
        if method == 'JS' or method == 'RMSE':
            Tools_score.append(pd.Series(list_up, index=Result.loc[method,:].sort_values(ascending=False).index))
        score_col.append(method)
        
    score=pd.concat([m for m in Tools_score],axis=1)
    score.columns = score_col
    score = score/Result.shape[1]
    return score

def make_score(dataset_all, Tools, path_input, path_output):
    for dataset in dataset_all:
        Tools_data=[x for x in range(len(Tools))]
        for i in range(len(Tools)):
            File = f'{path_input}/dataset{dataset}_{Tools[i]}_Metrics.txt'
            if os.path.isfile(File): 
                Tools_data[i] = pd.read_table(File,sep='\t',index_col=0, header=0)
                Tools_data[i] = Tools_data[i].mean()
                Tools_data[i]['Tool'] = Tools[i]
            else:
                print(f'files not found in dataset{dataset} for {Tools[i]}')
                Tools_data[i] = pd.DataFrame([-2, -2, 10, 10], columns = ['Genes'], index = ['PCC', 'SSIM', 'RMSE', 'JS']).T
                Tools_data[i] = Tools_data[i].mean()
                Tools_data[i]['Tool'] = Tools[i]
        Result = pd.concat([m for m in Tools_data], axis = 1)
        Result.columns = Result.loc[["Tool"],:].values.flatten()
        Result.drop('Tool', axis = 0, inplace = True)
        
        score = get_as_score(Result)
        score.to_csv(f'{path_output}/dataset{dataset}_score.txt', header = 1, index = 1)

def make_all_score(dataset_all, path_input):
    score_all = pd.DataFrame()
    for dataset in dataset_all:
        a = pd.read_csv(f'{path_input}/Accuracy_Rank/dataset{dataset}_score.txt',header=0,index_col=0)
        score_all = pd.concat([score_all,a],axis=0)
    return score_all

if __name__ == "__main__":

    ## step 1
    outpath = '/mnt/disk1/hongjia/work/xiu/Nature_Methods/Score'
    PCC = cal_index(method='pcc')
    PCC.to_csv(os.path.join(outpath, 'PCC.csv'))

    RMSE = cal_index(method='rmse')
    RMSE.to_csv(os.path.join(outpath, 'RMSE.csv'))

    SSIM = cal_index(method='ssim')
    SSIM.to_csv(os.path.join(outpath, 'SSIM.csv'))

    JS = cal_index(method='jsd')
    JS.to_csv(os.path.join(outpath, 'JS.csv'))

    ### step 2
    path = '/mnt/disk1/hongjia/work/xiu/Nature_Methods/Score'
    path_output = '/mnt/disk1/hongjia/work/xiu/Nature_Methods/Score'
    PCC = pd.read_csv(os.path.join(path, 'PCC.csv'), sep = ',', header = 0, index_col = 0)
    RMSE = pd.read_csv(os.path.join(path, 'RMSE.csv'), sep = ',', header = 0, index_col = 0)
    SSIM = pd.read_csv(os.path.join(path, 'SSIM.csv'), sep = ',', header = 0, index_col = 0)
    JS = pd.read_csv(os.path.join(path, 'JS.csv'), sep = ',', header = 0, index_col = 0)

    Datasets = list(set(PCC['dataset']))
    RScore = None
    Tools = ['RCTD','Cell2location','Tangram','Seurat','Stereoscope','SPOTlight','SpatialDWLS','DestVI','STRIDE','DSTG', 'Redeconve', 'CellMap']
    for Data in Datasets:
        PCCUse = PCC[PCC['dataset'] == Data]
        PCCUse = PCCUse.dropna(axis = 1, how = 'all')
        PCCUse = PCCUse.fillna(0)
        
        RMSEUse = RMSE[RMSE['dataset'] == Data]
        RMSEUse = RMSEUse.dropna(axis = 1, how = 'all')
        RMSEUse = RMSEUse.fillna(1)
        
        JSUse = JS[JS['dataset'] == Data]
        JSUse = JSUse.dropna(axis = 1, how = 'all')
        JSUse = JSUse.fillna(1)
        
        SSIMUse = SSIM[SSIM['dataset'] == Data]
        SSIMUse = SSIMUse.dropna(axis = 1, how = 'all')
        SSIMUse = SSIMUse.fillna(0)
        
        Tools = PCCUse.columns.tolist()
        for Tool in Tools[0:-1]:
            Metric = pd.concat([PCCUse[Tool], SSIMUse[Tool], RMSEUse[Tool], JSUse[Tool]], axis = 1)
            Metric.columns = ['PCC', 'SSIM', 'RMSE', 'JS']
            Metric.to_csv(f'{path_output}/dataset{Data}_{Tool}_Metrics.txt', sep = '\t')
            Result = pd.DataFrame(Metric.mean())
            Result.columns = ['AverageScore']
            Result['Tool'] = Tool
            Result['Data'] = Data
            if RScore is None:
                RScore = Result
            else:
                RScore = pd.concat([RScore, Result])
    RScore['Metrics'] = RScore.index

    ### step 3
    path_input = '/mnt/disk1/hongjia/work/xiu/Nature_Methods/Score'
    path_output = path_input + '/Accuracy_Rank/'
    Datasets = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,26,27,28,29,30,31,32]
    Tools = ['RCTD','Cell2location','Tangram','Seurat','Stereoscope','SPOTlight','SpatialDWLS','DestVI','STRIDE','DSTG','Redeconve', 'CellMap']
    make_score(Datasets, Tools, path_input, path_output)
    Score = make_all_score(Datasets, path_input)
    AS = pd.DataFrame(Score.sum(axis=1).values,index = Score.index.values,columns=['merge'])
    AS = AS/4
    print(Score)
    print(AS.groupby(AS.index).median())

    

import scanpy as sc
import os
import pandas as pd
import numpy as np

import argparse
from scipy.stats import pearsonr
from collections.abc import Iterable
from scipy.spatial.distance import jensenshannon
def jsd(v1,v2):
    return jensenshannon(v1-v2)
from sklearn.metrics import mean_squared_error
def rmse(x1,x2):
    return mean_squared_error(x1,x2,squared=False)

def ssim(im1,im2,M=1):
    im1, im2 = im1/im1.max(), im2/im2.max()
    mu1 = im1.mean()
    mu2 = im2.mean()
    sigma1 = np.sqrt(((im1 - mu1) ** 2).mean())
    sigma2 = np.sqrt(((im2 - mu2) ** 2).mean())
    sigma12 = ((im1 - mu1) * (im2 - mu2)).mean()
    k1, k2, L = 0.01, 0.03, M
    C1 = (k1*L) ** 2
    C2 = (k2*L) ** 2
    C3 = C2/2
    l12 = (2*mu1*mu2 + C1)/(mu1 ** 2 + mu2 ** 2 + C1)
    c12 = (2*sigma1*sigma2 + C2)/(sigma1 ** 2 + sigma2 ** 2 + C2)
    s12 = (sigma12 + C3)/(sigma1*sigma2 + C3)
    ssim = l12 * c12 * s12
    return ssim

def str_remove_spe(string):
    s_list = []
    for s in string:
        if s.isalnum():
            s_list.append(s)
        else:
            s_list.append('_')
    new_string = ''.join(s_list)
    return new_string

def get_cellmap_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def read_results(result_path, func, celltype):
    if os.path.exists(result_path):
        results = func(result_path, celltype=celltype)
    # else:
    #     results = pd.DataFrame(0, index=index, columns=celltype)
        return results

def get_dataset_results(dataset_number, path_result):
    # st = sc.read_h5ad(f'/mnt/disk1/hongjia/work/xiu/Nature_Methods/SimualtedSpatalData/dataset{dataset_number}/Spatial.h5ad')
    # gd = st.uns['density']
    gd = read_results('/mnt/disk1/hongjia/work/HER2+/shaozong/RCTD.txt', func=get_cellmap_results, celltype=None)
    gd.columns = [str_remove_spe(c) for c in gd.columns]
    gd = (gd.T/gd.sum(axis=1)).T
    gd = gd.loc[:,np.unique(gd.columns)]   
    index = gd.index
    celltype = gd.columns
    cellmap_results = read_results('/mnt/disk1/hongjia/work/HER2+/shaozong/cellmap.txt', func=get_cellmap_results, celltype=celltype)
    return gd, cellmap_results

def get_RCTD_results(results_path,celltype=None):
    results = pd.read_csv(results_path,index_col=0)
    results.columns = [str_remove_spe(c) for c in results.columns]
    if celltype is not None:
        results.loc[:,np.setdiff1d(celltype,results.columns)] = 0
    results = results.loc[:,np.unique(results.columns)]
    return results

def compare_results(gd,result_list,metric='pcc',columns=None,axis=1):
    if metric=='pcc':
        func = pearsonr
        r_ind = 0
    if metric=='jsd':
        func = jensenshannon
        r_ind = None
    if metric=='rmse':
        func = rmse
        r_ind = None
    if metric=='ssim':
        func = ssim
        r_ind = None
    if isinstance(result_list, pd.DataFrame):
        c_list = []
        if axis == 1:
            print('axis: ',1)
            for i,c in enumerate(gd.columns):
                r = func(gd.iloc[:,i].values, np.nan_to_num(np.clip(result_list.iloc[:,i],0,1)))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        else:
            print('axis: ',0)
            for i,c in enumerate(gd.index):
                r = func(gd.iloc[i,:].values, np.nan_to_num(np.clip(result_list.iloc[i,:],0,1)))
                if isinstance(result_list, Iterable):
                    if r_ind is not None:
                        r = r[r_ind]
                c_list.append(r)
        df = pd.DataFrame(c_list,index=gd.columns,columns=columns)
    else:
        df_list = []
        for n,res in enumerate(result_list):
            # print('- '+str(n))
            c_list = []
            if axis == 1:
                for i,c in enumerate(gd.columns):
                    r = func(gd.iloc[:,i].values, np.nan_to_num(np.clip(res.iloc[:,i],0,1)))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list,index=gd.columns)
            else:
                for i,c in enumerate(gd.index):
                    r = func(gd.iloc[i,:].values, np.nan_to_num(np.clip(res.iloc[i,:],0,1)))
                    if isinstance(res, Iterable):
                        if r_ind is not None:
                            r = r[r_ind]
                    c_list.append(r)
                df_tmp = pd.DataFrame(c_list,index=gd.index)
            df_list.append(df_tmp)
        df = pd.concat(df_list,axis=1)
        df.columns = columns
    return df

def cal_index(method, dataset_number, path_result): 

    gd, results = get_dataset_results(dataset_number, path_result)
    df = compare_results(
        gd,
        results,
        columns = ['CellMap'],
        axis=1,
        metric=method
    )
    return df

if __name__ == "__main__":

    path_result = '/mnt/disk1/hongjia/work/HER2+/shaozong'
    dataset_number = 1
    PCC = cal_index('pcc', dataset_number, path_result)
    PCC = PCC.dropna(axis = 1, how = 'all')
    PCC = PCC.fillna(0)

    RMSE = cal_index('rmse', dataset_number, path_result)
    RMSE = RMSE.dropna(axis = 1, how = 'all')
    RMSE = RMSE.fillna(1)

    SSIM = cal_index('ssim', dataset_number, path_result)
    SSIM = SSIM.dropna(axis = 1, how = 'all')
    SSIM = SSIM.fillna(0)

    JS = cal_index('jsd', dataset_number, path_result)
    JS = JS.dropna(axis = 1, how = 'all')
    JS = JS.fillna(1)

    Metric = pd.concat([PCC['CellMap'], SSIM['CellMap'], RMSE['CellMap'], JS['CellMap']], axis = 1)
    Metric.columns = ['PCC', 'SSIM', 'RMSE', 'JS']
    print(Metric.mean(axis=0))