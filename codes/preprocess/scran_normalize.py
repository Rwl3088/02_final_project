def scran_normalize_JZ(adata):
    import subprocess
    import numpy as np
    import pandas as pd
    #import rpy2.robjects as ro
    #from rpy2.robjects import numpy2ri
    #from rpy2.robjects.packages import importr
    #importr('scran')
    from scipy.io import mmwrite
    mmwrite('temp.mtx',adata.X.T)
    #numpy2ri.activate()
    #ro.r.assign('mat', adata.X.T)
    #ro.r('mat <- Matrix::readMM("temp.mtx")')
    #qclust_params = 'mat'
    # qclust_params = f'mat, min.size={min_size}, max.size={max_size}'
    #ro.reval(f'cl <- quickCluster({qclust_params})')
    #csf_params = f'mat, clusters=cl' 
    # csf_params = f'mat, clusters=cl, min.mean={min_mean}' 
    #sf = np.asarray(ro.reval(
    #    f'computeSumFactors({csf_params})'
    #))
    subprocess.call("Rscript /B_ALL/script/scran_normalize.R",shell=True)
    sf = pd.read_csv("sf.csv",index_col = 0).to_numpy()
    adata.obs['sf'] = sf   
    adata.layers['counts'] = adata.X.copy()
    adata.X /= adata.obs['sf'].values[:, None]
    #numpy2ri.deactivate()
    #subprocess.call('rm temp.mtx',shell=True)
    return adata
