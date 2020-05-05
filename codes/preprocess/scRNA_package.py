def print_Cs(adata, prefix='PC', n_genes=100, n_comps=20):
    import pandas as pd
    for i in range(n_comps):
        print(f'> {prefix}{i}')
        pc = pd.DataFrame(
            adata.varm[f'{prefix}s'][:,i],
            index=adata.var_names,
            columns=[prefix],
        )
        print('\n'.join(
            pc.abs().sort_values(
                by=prefix,
                ascending=False,
            ).head(n_genes).index.tolist())
        )

def plot_silhouette(
        adata, 
        n_comps,
        fname,
        rep='X_pca',
        algs=None,
        res_array=None
    ):
    # NOTE Please check silhouette score first argument
    import matplotlib.pyplot as plt
    from sklearn.metrics import silhouette_score
    fig = plt.figure(dpi=150, figsize=(10,5))
    for a in algs:
        sil_res = [
            r for r in res_array if adata.obs[f'{a}_{r}'].unique().shape[0] > 1
        ]
        sil = [
            silhouette_score(
                adata.obsm[rep][:,:n_comps],
                adata.obs[f'{a}_{r}']
            ) for r in sil_res 
        ]
        ks = [
            adata.obs[f'{a}_{r}'].astype(int).max() + 1 for r in sil_res
        ]
        plt.plot(sil_res, sil, label=f'{a}')
        for r, s, k in zip(sil_res, sil, ks):
            plt.text(r, s*1.01, k, fontdict={'fontsize': 6})
    plt.legend()
    fig.savefig(fname, bbox_inches='tight')

def ica(adata, n_components=None, inplace=True, **kwargs): 
    from sklearn.decomposition import FastICA 
    ica_transformer = FastICA(n_components=n_components, **kwargs) 
    x_ica = ica_transformer.fit_transform(adata.X) 
    if inplace:
        adata.obsm['X_ica'] = x_ica 
        adata.varm['ICs'] = ica_transformer.components_.T 
    else:
        return ica_transformer 

def nmf(adata, n_components=None, inplace=True, **kwargs): 
    from sklearn.decomposition import NMF 
    nmf_transformer = NMF(n_components=n_components, **kwargs) 
    x_nmf = nmf_transformer.fit_transform(adata.X) 
    if inplace:
        adata.obsm['X_nmf'] = x_nmf 
        adata.varm['NMFs'] = nmf_transformer.components_.T 
    else:
        return nmf_transformer 

def scran_normalize(adata):
    import numpy as np
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.packages import importr
    importr('scran')
    numpy2ri.activate()
    ro.r.assign('mat', adata.X.T)
    qclust_params = 'mat'
    # qclust_params = f'mat, min.size={min_size}, max.size={max_size}'
    ro.reval(f'cl <- quickCluster({qclust_params})')
    csf_params = f'mat, clusters=cl' 
    # csf_params = f'mat, clusters=cl, min.mean={min_mean}' 
    sf = np.asarray(ro.reval(
        f'computeSumFactors({csf_params})'
    ))
    adata.obs['sf'] = sf   
    adata.layers['counts'] = adata.X.copy()
    adata.X /= adata.obs['sf'].values[:, None]
    numpy2ri.deactivate()
    return adata

def mast_de(adata, key, perc=-1, covs=''):
    import numpy as np
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr
    ro.reval('rm(list=ls())')
    importr('MAST')
    numpy2ri.activate()
    pandas2ri.activate()
    ro.reval('options(mc.cores=2)')
    ro.r.assign(
        'mat',
        pd.DataFrame(
            adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names
        )
    )
    print('Filtering genes')
    if perc > 0:
        ro.reval('mat <- mat[,colSums(mat > 0) > length(colnames(mat)) * ' + str(perc) + ']; dim(mat);')
    else:
        ro.reval('mat <- mat[,colSums(mat > 0) > ' + str(adata.obs[key].value_counts().min() - 1) + ']; dim(mat);')
    if adata.obs.shape[1] == 0:
        adata.obs['Barcode'] = list(adata.obs.index)
    ro.r.assign('cdat', adata.obs)
    ro.reval('fdat <- as.data.frame(colnames(mat)); \
        row.names(fdat) <- fdat[,1]')
    ro.reval('raw <- FromMatrix(t(mat), cdat, fdat)')
    print('Data loaded')
    de = None
    for group in adata.obs[key].cat.categories:
        print(f'Group {group}')
        cmd = 'group <- colData(raw)$' + key + '; \
            levels(group) <- factor(c(unique(group), "-1")); \
            group[group != "' + group + '"] = "-1"; \
            colData(raw)$group <- group; \
            zlmCond <- zlm(~group + n_genes' + covs + ', raw); \
            summaryCond <- summary(zlmCond, doLRT="group' + group + '"); \
            summaryDT <- summaryCond$datatable; \
            fcHurdle <- merge(summaryDT[contrast=="group' + group + '" & component=="H",.(primerid, `Pr(>Chisq)`)], summaryDT[contrast=="group' + group + '" & component=="logFC", .(primerid, coef, ci.hi, ci.lo)], by="primerid"); \
            fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, "fdr")];'
        ro.reval(cmd)
        index = list(ro.reval('data.frame(fcHurdle)$primerid'))
        coef = list(ro.reval('data.frame(fcHurdle)$coef'))
        fdr = list(ro.reval('data.frame(fcHurdle)$fdr'))
        group_de = pd.DataFrame(
            {f'{group}_coef': coef, f'{group}_fdr': fdr},
            index=index,
        )
        de = pd.concat([de, group_de], axis=1)
    return de

def slingshot(adata, start, n_pcs=5, cl=None):
    import numpy as np
    import pandas as pd
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri, pandas2ri
    from rpy2.robjects.packages import importr
    importr('slingshot')
    numpy2ri.activate()
    pandas2ri.activate()
    ro.r.assign('pca', adata.obsm['X_pca'][:,:n_pcs]) 
    ro.r.assign('cl', adata.obs[cl])
    ro.reval('sds <- newSlingshotDataSet(pca, cl)') 
    ro.reval(f'sce <- slingshot(sds, cl, start.clus="{start}")')
    pt = pd.DataFrame(np.asarray(ro.reval('slingPseudotime(sce)')), index=adata.obs_names)
    pt.columns = [f'{cl}_lineage_{c}' for c in pt.columns]
    try:
        adata.obs = adata.obs.drop(pt.columns, axis=1)
    except KeyError:
        print('PT keys not dropped in obs dataframe: Not found.')
    adata.obs = pd.concat([adata.obs, pt], axis=1)
    adata.uns['slingshot'] = {} 
    adata.uns['slingshot']['lineages'] = {} 
    lineages = np.asarray(np.asarray(ro.reval('sce@lineages')))
    for i, l in enumerate(lineages): 
        adata.uns['slingshot']['lineages'][i] = list(np.asarray(l))
    numpy2ri.deactivate()
    pandas2ri.deactivate()
    return adata

def plot_pseudotime(adata, save, embedding='umap'):
    import copy
    from matplotlib.cm import viridis
    import scanpy as sc
    if embedding == 'umap':
        efunc = sc.pl.umap
    elif embedding == 'tsne':
        efunc = sc.pl.tsne
    cmap = copy.copy(viridis)
    cmap.set_under(color='lightgray', alpha=.5)
    efunc(
        adata,
        color=list(adata.obs.columns),
        cmap=cmap,
        save=save,
        vmin=0,
        legend_loc='on data',
    ) 

def plot_pseudotime_genes(df, cells, genes, row_colors, col_colors,
        filepath, figsize=(30,10), dpi=100):
    import seaborn as sns
    g = sns.clustermap(df.loc[cells, genes], cmap='viridis', col_cluster=False, row_cluster=False, yticklabels=False, figsize=(30,10), row_colors=row_colors, col_colors=col_colors, z_score=1)
    g.ax_heatmap.set_xlabel('')
    g.ax_heatmap.set_ylabel('')
    g.cax.set_visible(False)
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)
    g.savefig(filepath, bbox_inches='tight', dpi=dpi)

def de_genes(adata, groupby, lin):
    import pandas as pd
    import scanpy as sc
    sc.tl.rank_genes_groups(adata=adata, groupby=groupby, use_raw=False)
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    de_df = pd.DataFrame({group: result[key][group]
        for group in groups for key in ['names']}).head(10)
    de_df.index.name = 'Genes'
    de_df.columns.name = 'Groups'
    genes = de_df.T.stack().loc[lin].reindex(lin, level=0)
    return genes

def fullpath_closure(data_dir):
    import os
    def fullpath(path):
        return os.path.join(data_dir, path)
    return fullpath

def cell_entropy(adata):
    import scipy.stats as ss
    df = adata.to_df()
    df_p = df.div(df.sum(axis=1), axis=0)
    ent = df_p.apply(ss.entropy, axis=1)
    ent.name = 'entropy'
    return ent

def get_proliferation(org):
    # Compiled from Tirosh et al., Science, (2015)
    # Downloaded from https://github.com/theislab/scanpy_usage/blob/master/180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt
    proliferation = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8', 'HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    if org == 'hsapiens':
        pass
    else:
        proliferation = orthologs(proliferation, 'hsapiens', org)
    return proliferation

def get_cell_cycle(org):
    # FIXME use automatic conversion
    # Compiled from Whitfield et al, Molecular biology of the cell, (2002)
    # http://genome-www.stanford.edu/Human-CellCycle/HeLa/supplement.shtml
    # Supplemental figure 15
    # Converted to orthologs with g:Profiler
    cc_genes = {}
    if org == 'hsapiens':
        cc_genes = {
            'G1/S': ['CCNE1', 'E2F1', 'CDC6', 'PCNA'],
            'S': ['RFC4', 'DHFR', 'RRM2', 'RAD51'],
            'G2': ['CDC2', 'TOP2A', 'CCNF', 'CCNA2'],
            'G2/M': ['STK15', 'BUB1', 'CCNB1', 'PLK1'],
            'M/G1': ['PTTG1', 'RAD21', 'VEGFC', 'CDKN3']
        }
    elif org == 'mmusculus':
        cc_genes = {
            'G1/S': ['Ccne1', 'E2f1', 'Cdc6', 'Pcna', 'Pcna-ps2'],
            'S': ['Rfc4', 'Dhfr', 'Rrm2', 'Rad51'],
            'G2': ['Cdk1', 'Top2a', 'Ccnf', 'Ccna2'],
            'G2/M': ['Bub1', 'Ccnb1', 'Plk1'],
            'M/G1': ['Pttg1', 'Rad21', 'Vegfc', 'Cdkn3'] 
        }
    
    return cc_genes

def upper_quartile(df, q=0.75):
    quart = df.quantile(q, axis=1)
    factor = quart / quart.median()
    df_n = df.div(factor, axis=0)
    return df_n

#### These should be together
from numba import njit
@njit
def c(a, b):
    import numpy as np
    n, k = a.shape
    m, k = b.shape
    mu = (k - 1) / 2
    sig = ((k - 1) * (k + 1) / 12) ** .5
    out = np.empty((n, m))
    a = a - mu
    b = b - mu
    for i in range(n):
        for j in range(m):
            out[i, j] = a[i] @ b[j] / k / sig ** 2
    return out

def rank(a):
    import numpy as np
    i, j = np.meshgrid(*map(np.arange, a.shape), indexing='ij')
    s = a.argsort(1)
    out = np.empty_like(s)
    out[i, s] = j
    return out

def spearman_df(a, b):
    import pandas as pd
    return pd.DataFrame(c(rank(a.values), rank(b.values)), a.index, b.index)
####

def light_color(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def HSC_palette():
    import seaborn as sns
    pal = sns.color_palette('tab20', n_colors=20)[:20] + [(.85, .85, .85)]
    pal = [pal[i] for i in [7, -1, 0]]
    pal = [light_color(pal[i], 1.3) for i in range(len(pal))]
    lut = dict(zip(['non-HSC', 'undefined', 'HSC'], pal))
    return pal, lut

def adata_project(adata, adata_h):
    import scanpy.api as sc
    adata_p = adata.copy()
    sc.tl.pca(adata_p, random_state=42, n_comps=2)
    sc.tl.tsne(adata_p, random_state=42, n_pcs=2)
    adata_p.obsm['X_tsne'] = adata_h.obsm['X_tsne'].copy()
    return adata_p

def adata_var_filter(_adata, prefix_l=[], gene_l=[]):
    n_vars_before = _adata.n_vars
    for prefix in prefix_l:
        _adata = _adata[:, _adata.var.index[~_adata.var.index.str.startswith(prefix)]]
    if _adata.n_vars == n_vars_before:
        print('* No variables filtered using prefix.') 
    else:
        print(f'* {n_vars_before - _adata.n_vars} variables filtered using prefix.')
    n_vars_before = _adata.n_vars
    _adata = _adata[:, _adata.var.index[~_adata.var.index.isin(gene_l)]]
    if _adata.n_vars == n_vars_before:
        print('* No variables filtered using gene name.') 
    else:
        print(f'* {n_vars_before - _adata.n_vars} variables filtered using gene name.')
    return _adata

def adata_var_get(_adata, prefix_l=[], gene_l=[]):
    vars_l = set()
    for prefix in prefix_l:
        vars_l |= set(_adata.var.index[_adata.var.index.str.startswith(prefix)])
    vars_l |= set(_adata.var.index[_adata.var.index.isin(gene_l)])
    return _adata[:, list(vars_l)]

def rows_permutation(df_inner):
    import numpy as np
    return df_inner.astype(np.float64).apply(np.random.permutation, axis=0)

def get_noise_ratio(df_inner, n=100):
    import numpy as np
    from sklearn.decomposition import PCA
    noise_ratio = 0.
    for i in range(n):
        pca = PCA(n_components=1)
        pca = pca.fit(rows_permutation(df_inner))
        noise_ratio = max(noise_ratio, pca.explained_variance_ratio_[0])
        print('.', end='')
    print(n)
    return noise_ratio

def get_topev(df, n=100):
    import numpy as np
    from sklearn.decomposition import PCA
    evs = []
    for i in range(n):
        pca = PCA(n_components=1)
        pca = pca.fit(rows_permutation(df))
        evs += [pca.explained_variance_ratio_[0]]
        print('.', end='')
    print(n)
    return evs

def get_topev_parallel(df, n=1000):
    from joblib import Parallel, delayed
    import multiprocessing as mp
    import os
    import pandas as pd
    size = df.shape[0]
    # Set to the maximum number of CPUs the process can use
    if os.cpu_count() > 1:
        n_jobs = os.cpu_count() // 2
    topev = Parallel(n_jobs=n_jobs)(
        delayed(get_topev)(df, size) for _ in range(n)
    )
    return topev

def denoise(df_orig, nr):
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    mu = np.mean(df_orig, axis=0)
    pca = PCA()
    pca.fit(df_orig)
    nComp = (pca.explained_variance_ratio_ > nr).sum()
    df_hat = np.dot(pca.transform(df_orig)[:,:nComp], pca.components_[:nComp,:])
    df_hat = pd.DataFrame(df_hat, index=df_orig.index, columns=df_orig.columns)
    df_hat += mu
    return df_hat

def remove_pcs(df_orig, pcs=[]):
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    mu = np.mean(df_orig, axis=0)
    pca = PCA()
    pca.fit(df_orig)
    comp = [i for i in range(df_orig.shape[0]) if i not in pcs]
    df_hat = np.dot(pca.transform(df_orig)[:, comp], pca.components_[comp,:])
    df_hat = pd.DataFrame(df_hat, index=df_orig.index, columns=df_orig.columns)
    df_hat += mu
    return df_hat

def bootstrap(df_orig, n_samples):
    return apply_parallel(df_orig, bootstrap_func, n_samples)

def bootstrap_func(df_orig, size):
    import numpy as np
    return df_orig.apply(np.random.choice, size=size, replace=True).sum()

def apply_parallel(df_orig, func, n_executions):
    ''' Repeat func on df_orig for n_executions, in parallel '''
    from joblib import Parallel, delayed
    import multiprocessing as mp
    import os
    import pandas as pd
    size = df_orig.shape[0]
    # Set to the maximum number of CPUs the process can use
    n_jobs = os.cpu_count()
    bs = Parallel(n_jobs=n_jobs)(delayed(bootstrap_func)(df_orig, size)
        for _ in range(n_executions))
    return pd.concat(bs, axis=1)

def load_cellranger(mtx, genes, barcodes):
    ''' Load CellRanger output
        Parameters:
        mtx : filepath for mtx matrix with sparse data
        genes : filepath for genes list corresponding to the first
          dimension of the matrix
        barcodes : filepath for barcodes list corresponding to the
          second dimension of the matrix
    '''
    import pandas as pd
    from scipy.io import mmread
    # Read the matrix
    mat = mmread(mtx)
    df = pd.SparseDataFrame(mat)
    genes = pd.read_csv(genes, sep='\t', header=None)
    barcodes = pd.read_csv(barcodes, sep='\t', header=None)[0]
    df.index = genes[[0, 1]]
    df.index = pd.MultiIndex.from_tuples(df.index)
    df.columns = list(barcodes)
    print('For the same gene symbol there may be duplicates. Selecting the one with the highest median signal.')
    df['Median'] = df.median(axis=1)
    df = df.sort_values(by=['Median'], ascending=False, na_position='last')
    df = df.drop(columns=['Median'])
    df = df.fillna(0).to_dense()
    df = df.groupby(level=1).first()
    # Transpose in the form cells x genes
    df = df.T
    return df

def var_names_make_unique(df, method='median'):
    df = df.T
    if method == 'median':
        df['Median'] = df.median(axis=1)
        df = df.sort_values(by=['Median'], ascending=False, na_position='last')
        df = df.drop(columns=['Median'])
        df = df.groupby(level=0).first()
        # Transpose in the form cells x genes
        df = df.T
    else:
        raise NotImplementedError()
    return df

def subsample_category(adata, category, size=100):
    ''' Get subsample of size for each group defined by category '''
    import numpy as np
    grouped = adata.obs.loc[:, [category]].groupby(category)
    subsampled = adata[[x for a in grouped.groups.values() for x in np.random.choice(a, size=size)], :]
    subsampled.obs_names = [f's{i}' for i in range(subsampled.shape[0])]
    return subsampled

def violin_null(adata_null, color, category, n_cols=4, q_thr=1e-3, es_thr=1):
    ''' Violin plot for each group in category, null distribution included.
        Kolgomorov-Smirnov of each group against the NULL
    '''
    import matplotlib.pyplot as plt
    import scipy.stats as ss
    import seaborn as sns
    from statsmodels.stats.multitest import fdrcorrection
    with plt.style.context('seaborn-white'):
        n_rows = len(color) // n_cols
        if len(color) % n_cols != 0:
            n_rows += 1
        fig, axs = plt.subplots(n_rows, n_cols, figsize=(n_cols*4, n_rows*4))
        if n_cols != 1:
            axs = axs.ravel()
        elif len(color) == 1:
            axs = [axs]
        for i, c in enumerate(color):
            sns.violinplot(data=adata_null.obs, x=category, y=c, linewidth=.1, inner=None, scale='width', ax=axs[i])
            sns.stripplot(data=adata_null.obs, x=category, y=c, linewidth=0, size=1, color='k', alpha=.5, ax=axs[i])
            ylim = axs[i].get_ylim()
            ylim_star = ylim[1] - .08 * (ylim[1]-ylim[0])

            groups = list(adata_null.obs[category].cat.categories)[:-1]
            es = []
            pval = []
            for g in groups:
                x = adata_null.obs.loc[adata_null.obs[category] == g, c]
                y = adata_null.obs.loc[adata_null.obs[category] == 'N', c]
                es += [abs(x.median() - y.median())]
                pval += [ss.ks_2samp(x, y).pvalue]
            qval = fdrcorrection(pval, method='negcorr')[1]

            for j, g in enumerate(groups):
                if qval[j] < q_thr and es[j] > es_thr:
                    axs[i].text(j, ylim_star, '*', ha='center')
        plt.subplots_adjust(wspace=.3)
    
    return fig

def ontologies(go_terms, org):
    import requests as rq
    response = rq.post(
        url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/',
        json={
            'organism': org,
            'target': 'ENSG',
            'query': go_terms,
        },
    )
    return list(set([d['name'] for d in response.json()['result']]))

def orthologs(genes, org, target):
    import requests as rq
    response = rq.post(
        url='https://biit.cs.ut.ee/gprofiler/api/orth/orth/',
        json={
            'organism': org,
            'target': target,
            'query': genes,
        },
    )
    return [d['name'] for d in response.json()['result'] if d['name'] != 'N/A']


def hippo_signature(top=50):
    hippo_oligo = [
        'Mobp', 'Opalin', 'Tmem125',
    ]
    hippo_opcs = [
        'Pdgfra', 'Sema3d', 'Bmp4',
    ]
    hippo_astro = [
        'Aldoc', 'Aldh1l1', 'Gfap',
    ]
    # hippo_npcs = [
        # 'Hopx', 'Gfap',
    # ]
    hippo_homeo_microglia = [
        'Gpr34', 'Fcrls', 'Cx3cr1',
    ]
    hippo_act_microglia = [
        'Ccl2', 'Gpr84',
    ]
    hippo_vasc_epith = [
        'Vwf', 'Tmem252', 'Pecam1',
    ]
    hippo_macro_immune = [
        'Hmcn1',
    ]
    hippo_pyramid = [
        'Camk2a', 'Camk2b', 'Grin2a', 'Tbr1', 'Trank1', 'Adcy1', 'Grin2b',
    ]
    hippo_interneurons = [
        'Tbr1', 'Igfbpl1', 'Vip', 'Sst', 'Pvalb', 'Lamp5',
    ]
    hippo_smooth_muscle = [
        'Rbpms2', 'Dhx58os',
    ]
    hippo_pericytes = [
        'Tmem45a', 'Slc6a20a', 'Art3',
    ]
    hippo_endothelial = [
        'Abcb1a',
    ]

    loc = locals()
    signature = {
        m: loc[m][:top] for m in loc.keys() if m.startswith('hippo_')
    }

    return signature


def abcam_signature(top=50):
    abcam_neuroepith = [
        'Nes', 'Sox2', 'Notch1', 'Hes1', 'Hes3', 'Ocln', 'Cdh1', 'Sox10',
    ]
    abcam_radialglia = [
        'Vim', 'Pax6', 'Hes1', 'Hes5', 'Gfap', 'Slc1a3', 'Fabp7', 'Tnc', 'Cdh2', 'Nes', 'Sox2',
    ]
    abcam_ipc = [
        'Eomes', 'Ascl1',
    ]
    abcam_immatureneurons = [
        'Dcx', 'Neurod1', 'Tbr1', 'Tubb3', 'Stmn1',
    ]
    abcam_matureneurons = [
        'Rbfox3', 'Map2', 'Nefm', 'Nefh', 'Syp', 'Dlg4',
    ]
    abcam_glutamatergic = [
        'Slc17a7', 'Slc17a6', 'Grin1', 'Grin2b', 'Gls', 'Glul',
    ]
    abcam_gabaergic = [
        'Slc6a1', 'Gabbr1', 'Gabbr2', 'Gad2', 'Gad1',
    ]
    abcam_dopaminergic = [
        'Th', 'Slc6a3', 'Foxa2', 'Kcnj6', 'Nr4a2', 'Lmx1b',
    ]
    abcam_serotonergic = [
        'Tph1', 'Slc6a4', 'Fev',
    ]
    abcam_cholinergic = [
        'Chat', 'Slc18a3', 'Ache',
    ]
    abcam_microglia = [
        'Itgam', 'Ptprc', 'Aif1', 'Adgre1', 'Cd68', 'Cd40',
    ]
    abcam_opcs = [
        'Pdgfra', 'Cspg4',
    ]
    abcam_oligo = [
        'Olig1', 'Olig2', 'Olig3', 'Cldn11', 'Mbp', 'Mog', 'Sox10',
    ]
    abcam_astro = [
        'Gfap', 'Slc1a3', 'Slc1a2', 'Glul', 'S100b', 'Aldh1l1',
    ]
    abcam_schwannprecursors = [
        'Sox10', 'Gap43', 'Fabp7', 'Mpz', 'Dhh', 'Ngfr',
    ]
    abcam_myelinatingschwann = [
        'Sox10', 'S100b', 'S100a1', 'Egr2', 'Mbp', 'Mpz',
    ]
    abcam_nonmyelinatingschwann = [
        'Sox10', 'S100b', 'S100a1', 'Gap43', 'Ncam1', 'Ngfr', 
    ]

    loc = locals()
    signature = {
        m: loc[m][:top] for m in loc.keys() if m.startswith('abcam_')
    }

    return signature


def RS007_signature(top=50):
    # Whole liver (RS007) markers + PF
    untreated_endo_markers = ['Gpihbp1','Clec4g','Ptprb','Aqp1','S100a16','Egfl7','Plpp3','Gpr182','Dnase1l3','Kdr','Cd300lg','Fcgr2b','Sdpr','Fam167b','Cldn5','Adgrf5','Nrp1','Gng11','Bmp2','Tspan7','Ramp2','Tinagl1','Mmrn2','Adgrl4','Stab2','Meis2','Cd59a','Cyp4b1','Tm4sf1','Oit3','Tpbgl','Npr1','Sema6a','Arhgap31','Flt4','Timp3','Pecam1','Vamp5','Arhgef15','Pde2a','Stab1','Snrk','Eng','Plpp1','Tmem2','Myct1','Clec14a','Tie1','Elk3','Cmtm8']
    untreated_kuppfer_markers = ['Ctss','Lyz2','Laptm5','Spi1','Ctsh','Ccl6','Tyrobp','Mpeg1','H2-Ab1','Aif1','Rgs2','Cd74','H2-Aa','Cybb','Clec4a1','Pld4','Ctsc','Cd53','H2-Eb1','Wfdc17','Unc93b1','Lst1','Ly86','Sat1','Csf1r','Ptpn6','Cd300a','Clec4a3','Pirb','AF251705','Lcp1','Ckb','Ms4a6c','Fcer1g','Hck','Cd44','Clec12a','Cyba','Pilra','Tep1','Adgre1','Fam129a','Il6ra','C1qa','Ifi30','C1qb','Fgd2','Hk2','Axl','Nfam1']
    untreated_T_markers = ['Hcst','Lck','Trbc1','Cd3g','Ptprcap','Ltb','Gimap3','Rac2','Trbc2','Rasal3','Selplg','Skap1','Sept1','Thy1','Ms4a4b','B4galnt1','Tbc1d10c','Coro1a','Nkg7','Tmsb10','Limd2','Cd52','Ptpn22','Bcl2','Zap70','Il2rb','Ptpn18','Epsti1','Ptprc','Fam189b','Btg1','Dusp2','Rinl','Txk','Pglyrp1','Sept6','Ablim1','Lat','Evl','Klk8','Tesc','Ubash3a','Sh3bgrl3','Stat4','Tespa1','Itk','Il2rg','Cd3e','Fxyd5','Lsp1']
    untreated_hep1_markers = ['Serpina1c','Serpina1a','Apoc3','Fabp1','Serpina1b','Chchd10','Ttr','Serpina1d','Apoa2','Mup20','Mup10','Car3','Bhmt','Apoc4','Apoa1','Mgst1','Ttc36','Gnmt','Gc','Mup3','Aldob','Hpd','Trf','Sephs2','Ahsg','Alb','Phyh','Serpina3k','Mup7','Dbi','H2-Q10','Nudt7','Apoc1','Ddt','Scp2','Cat','Gstp1','Stard10','Mup11','Selenbp2','Cdo1','Uox','Wfdc21','Hp','Hpx','Sult2a8','Pcbd1','Adh1','Hrsp12','Msrb1']
    untreated_hsc_markers = ['Cxcl12','Dcn','Lrat','Sod3','Col14a1','Colec11','Gdf2','Reln','Rgs5','Prelp','Ecm1','Lhfp','Tmem56','Hand2','Rbp1','Bmp5','Plvap','Angptl2','Tgfbi','Vipr1','Ntn1','Col1a2','Ifitm1','Cygb','G0s2','Gucy1a3','Angptl6','Steap4','Ramp1','Sept4','Gdf10','Cryab','Ank3','Pth1r','Col3a1','Mustn1','Colec10','Arvcf','Abcc9','Mylk','Bgn','Raph1','Rbms3','Lrp1','Nrxn1','Ryk','Cped1','Scarf2','Spry1','Fgfr2']
    untreated_B_markers = ['Vpreb3','Iglc3','Cd79a','Iglc2','Cd79b','Serpinb1a','Ly6d','Tnfrsf13c','Igkc','Blnk','Cd19','Trim7','Ebf1','Ighm','H2-Eb1','Ighd','H2-DMb2','Ms4a1','Apobec3','H2-Ob','Cd74','Siglecg','Fcrl1','Fcmr','Ptprcap','H2-Aa','Cd37','Cyp4f18','Rhoh','Fcrla','Napsa','Fcer2a','Cd22','Gm43291','H2-Ab1','S1pr4','Lamb3','H2-Oa','Cd52','Tspan32','Dusp2','Pou2f2','Trim59','Bank1','Coro1a','Cd2','Ms4a4c','Spib','Ralgps2','Capg']
    untreated_hep2_markers = ['Hpd','Apoc4','Apoh','Car3','Serpina1d','Mup11','Ahsg','Ttc36','Gc','Bhmt','H2-Q10','Serpina3k','Serpina1c','Sephs2','Azgp1','Serpina1a','Mgst1','Ttr','Cat','Stard10','Serpina1b','Gnmt','Apoa1','Apoc3','Hp','Fgb','Aldob','Cyp2d9','Phyh','F2','Ambp','Trf','Pgrmc1','Serpinc1','Chchd10','Nudt7','Ddt','Alb','Mup20','Ces3a','Vtn','Mup7','Fgg','Mup3','Apoa2','Hpx','C3','C8g','Adh1','Kng1']
    untreated_dendrit_markers = ['Spib','Siglech','Iglc3','Cox6a2','Rpgrip1','Smim5','Irf8','Ighm','Bcl11a','Lair1','Pld4','Sell','Atp1b1','Ly6d','Pacsin1','Ly6c2','Rnase6','Gm5547','Cd7','Nucb2','Tspan13','St8sia4','Prkca','Mpeg1','Rgs10','Ctsh','Dnajc7','Amica1','Ppfia4','Cd200','Ly86','Fcrla','Blnk','Cybb','Alox5ap','BC147527','P2ry14','Ptpn6','Clec12a','Card11','Mctp2','Lrrc16a','Pafah1b3','Lefty1','Plac8','Bloc1s2','Sdc4','Cyth4','Tyrobp','Gna15']

    loc = locals()
    signature = {
        m: loc[m][:top] for m in loc.keys() if m.endswith('_markers')
    }

    return signature

def misc_mouse_signature(top=50):

    pf_markers = ['Upk3b', 'Gpm6a', 'Krt7', 'Msln']
    chol_markers = ['Epcam', 'Sox9', 'P2rx2', 'Dmbt1', 'Gm609', 'Cldn4', 'Itgb6', 'Pak3', 'Cftr', 'Krt19', 'Prom1', 'St14', 'Foxj1', 'Ddr1', 'Pak6']
    kuppfer_markers = ['Clec4f']

    loc = locals()
    signature = {
        m: loc[m][:top] for m in loc.keys() if m.endswith('_markers')
    }

    return signature

def dissociation_signature():

    dissociation_stress = 'Actg1 Btg1 Cxcl1 Dnajb4 Errfi1 H3f3b Hspb1 Irf1 Klf6 Mir22hg Nfkbia Pcf11 Pxdc1 Sdc4 Srf Tpm3 Usp2 Gadd45g Ankrd1 Btg2 Cyr61 Dusp1 Fam132b Hipk3 Hsph1 Irf8 Klf9 Mt1 Nfkbiz Pde4b Rap1b Serpine1 Srsf5 Tppp3 Wac Hspe1 Arid5a Ccnl1 Dcn Dusp8 Fos Hsp90aa1 Id3 Itpkc Litaf Mt2 Nop58 Per1 Rassf1 Skil Srsf7 Tra2a Zc3h12a Ier5 Atf3 Ccrn4l Ddx3x Egr1 Fosb Hsp90ab1 Idi1 Jun Lmna Myadm Nppc Phlda1 Rhob Slc10a6 Stat3 Tra2b Zfand5 Kcne4 Atf4 Cebpb Ddx5 Egr2 Fosl2 Hspa1a Ier2 Junb Maff Myc Nr4a1 Pnp Rhoh Slc38a2 Tagln2 Trib1 Zfp36 Bag3 Cebpd Des Eif1 Gadd45a Hspa1b Ier3 Jund Mafk Myd88 Odc1 Pnrc1 Ripk1 Slc41a1 Tiparp Tubb4b Zfp36l1 Bhlhe40 Cebpg Dnaja1 Eif5 Gcc1 Hspa5 Ifrd1 Klf2 Mcl1 Nckap5l Osgin1 Ppp1cc Sat1 Socs3 Tnfaip3 Tubb6 Zfp36l2 Brd2 Csrnp1 Dnajb1 Erf Gem Hspa8 Il6 Klf4 Midn Ncoa7 Oxnad1 Ppp1r15a Sbno2 Sqstm1 Tnfaip6 Ubc Zyx'.split()

    return dissociation_stress

def RS025_signature(top=50):
    # Whole fibrotic liver (RS025) markers
    fibr_hsc_markers = ['Gdf10','Lhfp','Hand2','Angptl2','Cxcl12','Dcn','Ntn1','Prelp','Col14a1','Col5a1','Fbln5','Gucy1a3','Fgfr2','Tcf21','Scarf2','Mylk','Pdgfrb','Sept4','Col1a2','Pam','Col3a1','Mustn1','Steap4','Pcolce','Ccdc80','Colec11','Ms4a4d','Reln','Cped1','Fstl1','Des','Sod3','Snhg18','G0s2','Efemp2','Serping1','Rbp1','Lamb1','Lrat','Gdf2','Rbms3','Postn','Frzb','Nkd1','Epha7','Ecm1','Emilin1','Col6a1','Lhx2','Bgn'] 
    fibr_endo_markers = ['Ptprb','Aqp1','Mmrn2','Egfl7','Cd300lg','Tmem2','Rasip1','Esam','Meis2','F8','Sema6a','Fgd5','Cyyr1','Sox18','Erg','Ushbp1','Arhgap31','Tie1','Myct1','Dnase1l3','Plekhg1','Fam167b','Pde2a','4931406P16Rik','Sdpr','Clec4g','Tinagl1','Emcn','Tek','Robo4','Elk3','Cdh5','Adgrl2','Tspan2','Gng11','Snrk','Fam124a','Fabp4','Exoc3l2','Nrp1','Bmp2','Stab1','Mapk12','Adgrg3','Hecw2','Fam171a1','Arhgef15','Fam43a','Rasal2','Tnfaip1']
    fibr_kuppfer_markers = ['Csf1r','C1qc','C1qa','Ccl6','Cybb','C1qb','Mpeg1','Ctss','Lyz2','Unc93b1','Spi1','Adgre1','Cd44','Tbxas1','Sat1','Wfdc17','Sirpa','Clec4a3','Laptm5','Fcgr3','Rab32','Ctsc','Ly86','Clec4a1','Axl','Cmklr1','Lst1','Tyrobp','Cd53','Lilrb4a','Efhd2','Ctsh','Cd300a','Fcgr4','Btk','Tpd52','Ctsa','Nfam1','Pirb','Atp2b1','C3ar1','Lair1','Inpp5d','AF251705','Syk','AI607873','Tnfaip2','Ncf2','Slc15a3','Itgal']
    fibr_chol_markers = ['Spp1','Tmem45a','Clu','Krt18','Sestd1','Tm4sf4','Fxyd3','Ambp','Bicc1','Qpct','Cp','Hnf1b','Krt19','Atp1b1','Cdh1','Serinc2','Sox9','Sorbs2','Tstd1','Fgfr3','Crp','Cystm1','Dsg2','Mal2','Cys1','Alcam','Slc39a4','Snhg18','Cxadr','Sox4','Cgn','Cyr61','Npdc1','Myrf','Plscr1','Shroom3','Lipt2','Spns3','Rnf128','Dsc2','Kcne3','Casp12','Proc','Slc35f5','Bche','Ptprf','Cep19','Anxa4','Rgs5','Tesc']
    fibr_dendrit_markers = ['Fcrla','Rnase6','Siglech','Ccr9','Spib','Bcl11a','Cox6a2','Smim5','Spns3','Sell','Lefty1','St8sia4','Ly6c2','Them6','Ly6d','Arhgap27os2','Lag3','Rpgrip1','Slco4a1','Prkca','Tspan13','Clec10a','Ccdc162','Iglc3','Cd200','Rnf122','Dnajc7','Atp2a3','Atp1b1','Net1','Il21r','Cd7','Cd180','Lrrc16a','Sla2','Scimp','Amica1','Tmem229b','Tbc1d8','Gna15','Blnk','Rgs18','Rabgap1l','Ighm','Nucb2','Ppfia4','Plac8','Irf8','Glcci1','Cxxc5']
    fibr_hep1_markers = ['Chchd10','Gc','Rbp4','Hamp','Hpx','Pcbd1','Cdo1','Apoc4','Aldob','Stard10','Ttr','Gstm1','Gstp1','Ahsg','Krt18','Fgg','Adh1','Ambp','Trf','Ddt','Sord','Apoa2','Cisd1','Serpina1b','Sephs2','Pgrmc1','Hagh','Fgb','Apoc1','Phyh','Apoc3','Fdx1','Dbi','Gsta3','Ttc36','Fam25c','Hmgcs2','Apoa1','Ppa1','Cat','Fabp1','Khk','Hp','Prdx6','Orm1','Serpina1a','Fxyd1','Msrb1','Apom','Hrsp12']
    fibr_hep2_markers = ['Cdo1','Aldob','Hmgcs2','Rbp4','Mat1a','Proc','Sord','Stard10','Ahsg','Hpx','Khk','Fgg','Sephs2','Serpinc1','Ces1c','Cyp2d9','Uox','Gsta3','Tkfc','Pzp','Itih3','Slc27a2','Fgb','Comt','Marc1','Slc38a4','Mug1','Cat','Slc38a3','Insig1','Car3','Apoc4','Ndrg2','Rdh7','C4b','C3','Apob','Cyp2d10','F10','Trf','Tdo2','Agt','Cps1','Ambp','Fasn','Acaa1b','Cp','F2', 'Ttr','Kng1']
    fibr_T_markers = ['Ms4a4b','Ltb','Gimap3','Ptprcap','Skap1','Trbc2','Satb1','Tcf7','Cd2','Txk','Trbc1','Lck','Hcst','Cd3d','Nkg7','Ablim1','Cd3g','Ptpn22','Gimap7','Fam189b','Lat','Thy1','Cd3e','Tbc1d10c','Bcl2','Rac2','Il7r','Coro1a','Il2rb','Dusp2','Sh2d2a','Klk8','Limd2','Cd52','Gimap1','Socs1','Rinl','Jakmip1','Ccl5','Mir142hg','Tnfrsf18','Cd28','Prkcq','Sept6','Arhgap15','Abcb1a','Gm26740','Il18r1', 'B4galnt1','Ikzf3']
    fibr_B1_markers = ['Iglc3','Gpr18','Ly6d','Iglc2','Gm43291','Sell','Mzb1','Vpreb3','Ighm','Cd79b','Ptprcap','Spib','Gm8369','Fam65b','Cytip','Fcrla','H2-Ob','Btla','Gm43603','Igkc','Ltb','Bcl11a','Trim7','Cd2','Iglc1','Cd74','Napsa','Hvcn1','Sbk1','Traf3ip3','Apobec3','Rhoh','H2-Eb1','Blnk','S1pr4','H2-Ab1','Bank1','Ablim1','Coro1a','H2-Aa','Gimap7','Stap1','Gimap1','Tbc1d10c','Jakmip1','Ero1lb','Clec2i','Rabgap1l','Sorl1','Tnfrsf13b']
    fibr_B2_markers = ['Vpreb3','Mzb1','Trp53inp1','Iglc2','Gm43291','Ell2','Iglc3','Txndc11','Gm5547','Igkc','Iglc1','Tram2','St6gal1','Fut8','Prr5','Plpp5','Igha','Fam46c','Chchd10','Ero1lb','Pafah1b3','Tspan13','Atp1b1','Tnfrsf13b','Cd79b','Pqlc3','Pycr1','Blnk','Edem1','2310001H17Rik','Txndc5','Creld2','Dennd5b','Enpp1','Sdc1','Fkbp11','Fam214a','Prdm1','Dnajb11','Sec24a','Edem2','Hook1','Glipr1','Ccpg1','Slamf7','Pim1','Syvn1','Ube2j1','Cytip','Fndc3a']
    fibr_pf_markers = ['Serpinh1','Trio','Gm12840','Ildr2','Flrt2','Dcn','Pbx1','Aebp1','Ier3','Col6a1','Adgrd1','Ptrf','Cav1','Bicc1','Ccdc80','Wt1','Usp53','Selm','Mgp','Tgfb2','Lamb2','Rarres2','Col5a2','Slc16a1','Ptprf','Cryab','Myrf','Cdh3','Tmem98','Eln','Hspb1','Gpm6a','Clu','Lhfp','Saa3','Entpd2','Abcb1b','Timp2','Serping1','Laptm4b','Cd200','Wnt4','Pcolce','Cfb','Stxbp6','Efna5','Gas1','Steap4','Met','Rhoj']
    fibr_pf_spec_markers = ['Cdh3', 'Gm12840', 'Flrt2', 'Adgrd1', 'Mgp', 'Ildr2', 'Wt1', 'Ier3', 'Slc16a1', 'Cav1', 'Trio', 'Ptprf', 'Usp53', 'Myrf']

    loc = locals()
    signature = {
        m: loc[m][:top] for m in loc.keys() if m.endswith('_markers')
    }

    return signature



def han_signature():
    ### Han et al. signatures
    # FC > 1, pct_group > 0.5, pct_rest < 0.5, max 20 genes per cluster
    han_endo = ['Clec4g', 'Plpp3', 'Kdr', 'Nrp1', 'Cyp4b1', 'Socs3', 'Fcgr2b', 'Eng', 'Oit3', 'Aqp1', 'Ehd3', 'Gpihbp1', 'Ptprb', 'Tspan7', 'Sparc', 'Plpp1', 'Egr1', 'Igfbp4', 'Egfl7']
    han_kuppfer = ['Vsig4', 'Fcna', 'Cfp', 'C1qc', 'Ctsc', 'Adgre1', 'Csf1r', 'Lpl', 'Sdc3', 'Folr2', 'Pltp', 'Il18bp', 'Creg1', 'Hpgd']
    han_dendrit = ['Gm2a', 'Ppt1', 'Cytip', 'Plbd1', 'Lsp1', 'Siglech', 'Ccr9', 'Ly6d', 'Klk1', 'Ly6c2', 'Lgals1', 'Irf8', 'Bst2', 'Plac8', 'Xbp1', 'Emp3', 'Ccl4', 'Pld4', 'Mpeg1']
    han_T = ['Trbc2', 'Trac', 'Vps37b', 'Crem', 'Ccl5', 'Gzma', 'Gzmb', 'Xcl1', 'Cd7', 'Il2rb', 'Nkg7', 'Ccl4', 'Ifngr1']
    han_granul = ['S100a9', 'S100a8', 'Clec4d', 'Hdc', 'Slpi', 'Irg1', 'Thbs1', 'Il1rn', 'Il1r2', 'Cd9', 'S100a11', 'S100a6', 'Lilr4b']
    # Clec7a not found?
    han_macro = ['Ccl9', 'Thbs1', 'Lgals3', 'Plac8', 'Clec7a']
    han_hep = ['Fabp1', 'Aldob', 'Wfdc21', 'Car3', 'Chchd10', 'Dbi', 'Bhmt', 'Gstm1', 'Apoa2', 'Apoc3', 'Akr1c6', 'Adh1', 'Serpina1e', 'Ttr', 'Uox', 'Urah', 'Serpina1a', 'Apoc4', 'Gnmt', 'Cdo1', 'Gsta3', 'Hpd', 'Ass1', 'Cyp2e1', 'Mat1a', 'Hmgcs2', 'Serpina1b', 'Rgn', 'Sult2a1', 'Serpina1c', 'Phyh', 'Arg1', 'Aldh1a1', 'Apoa1', 'Fgb', 'Gc', 'Hpx', 'Ahsg', 'Kng1', 'Serpina3k', 'Mup3', 'Rbp4', 'Fgg', 'Apoh', 'Fga', 'Serpinc1', 'Cyp3a11', 'Ces1c', 'H2-Q10', 'Cyp2c29', 'F2',
            'mt-Rnr1', 'mt-Nd5', 'mt-Co1']
    han_epith = ['Spp1', 'Clu', 'Cp', 'Ambp', 'Krt8', 'Tstd1', 'Krt18', 'Tm4sf4', 'Malat1', 'Anxa5', 'Dbi', 'Mgst1', 'Ttr', 'Hes1', 'Ifi27', 'Cd24a', 'Hspb1', 'Gstm1', 'Ndufa6', 'mt-Co1', 'Lcn2', 'Dmbt1', 'Sprr1a', 'Sfn', 'Krt19', 'Wfdc2', 'Ly6d', 'Mmp7', 'Plet1', 'Cldn4', 'Cldn7', 'Tspan8', 'Gsta4', 'S100g', 'Pdzk1ip1', 'Tnfrsf12a']
    han_B = ['Ighm', 'Ccr7', 'Jchain', 'Igha', 'Iglc1', 'Mzb1', 'Iglc2', 'Txndc5', 'Trp53inp1', 'Xbp1', 'Sec11c', 'Serp1', 'Hsp90b1', 'Ssr4', 'H13', 'Manf', 'Herpud1', 'Mtdh', 'Sec61b', 'Krtcap2', 'Spcs2', 'Spcs1']
    han_eryth = ['Hbb-bt', 'Alas2']
    han_stromal = ['Rarres2', 'C3', 'Igfbp6', 'Dcn', 'Serping1', 'Slpi', 'Upk3b', 'Igfbp5', 'Efemp1', 'Crip1', 'Fmod', 'Aebp1', 'Ogn', 'Nkain4', 'Fxyd3', 'Cav1', 'S100a6', 'Fmo2', 'Rnase4', 'Gas6']
    han_neutro = ['Ngp', 'Camp', 'Ltf', 'Retnlg', 'S100a9', 'Lcn2', 'S100a8', 'Chil3', 'Krt83', 'Ifitm6', 'Thbs1', 'Mcemp1', 'Anxa1', 'Dgat1', 'G0s2', 'Hdc', 'Lmnb1', 'Grina', 'Slpi', 'S100a11']

    loc = locals()
    signature = {
        m: loc[m] for m in loc.keys() if m.startswith('han_')
    }

    return signature


def macparland_signature(collapse=True):
    ### MacParland et al. signatures 
    macparland_endo = [['MGP', 'SPARCL1', 'TM4SF1', 'CLEC14A', 'ID1', 'IGFBP7', 'ADIRF', 'CTGF', 'VWF', 'CD9', 'C7', 'SRPX', 'ID3', 'CAV1', 'GNG11', 'AQP1', 'HSPG2', 'EMP1', 'SOX18', 'CLDN5'], ['CCL14', 'CLEC1B', 'FCN2', 'S100A13', 'FCN3', 'CRHBP', 'STAB1', 'GNG11', 'IFI27', 'CLEC4G', 'CLDN5', 'CCL23', 'OIT3', 'RAMP3', 'SGK1', 'DNASE1L3', 'LIFR', 'SPARC', 'ADGRL4', 'EGFL7', 'PCAT19', 'CDKN1C'], ['RAMP3', 'INMT', 'DNASE1L3', 'LIFR', 'PTGDS', 'C7', 'CTGF', 'TIMP3', 'RNASE1', 'ID3', 'ENG', 'MGP', 'PCAT19', 'HSPG2', 'GPM6A', 'PTPRB', 'VWF', 'FAM167B', 'SRPX', 'LTC4S', 'IFI27'],]
    macparland_chol = [['TFF2', 'SCGB3A1', 'FXYD2', 'KRT7', 'DEFB1', 'CD24', 'TFF3', 'LCN2', 'KRT19', 'CXCL1', 'PIGR', 'TFF1', 'CXCL6', 'LGALS2', 'TACSTD2', 'ELF3', 'SPP1', 'MUC5B', 'LGALS4']]
    macparland_hsc = [['ACTA2', 'COL1A1', 'TAGLN', 'COL1A2', 'COL3A1', 'SPARC', 'RBP1', 'DCN', 'MYL9', 'TPM2', 'MEG3', 'BGN', 'IGFBP7', 'IGFBP3', 'CYR61', 'OLFML3', 'IGFBP6', 'CCL2', 'COLEC11', 'CTGF', 'HGF']]
    macparland_macro = [['S100A8', 'LYZ', 'S100A9', 'HLA-DPB1', 'S100A12', 'RP11-1143G9.4', 'EVI2A', 'HLA-DPA1', 'VCAN', 'S100A6', 'CXCL8', 'HLA-DRA', 'MNDA', 'TYR-OBP', 'HLA-DRB1', 'FCN1', 'HLA-DQA1', 'IL18', 'C1QC', 'CD74', 'HLA-DRB5'], ['CD5L', 'MARCO', 'VSIG4', 'CPVL', 'CD163', 'CCDC88A', 'C5AR1', 'LIPA', 'LILRB5', 'MAF', 'CTSB', 'MS4A7', 'VMO1', 'RAB31', 'SLC31A2', 'TTYH3', 'VCAM1', 'KLF4', 'HMOX1', 'AIF1l', 'TMIGD3'], ]
    # Proliferative cluster has been excluded
    macparland_T = [['CD2', 'CD3D', 'TRAC', 'GZMK', 'CCL5', 'CCL4L2', 'PYHIN1', 'TRBC1', 'TRBC2', 'GZMA', 'CD3E', 'JUNB', 'CD69', 'IL7R', 'DUSP2', 'IFNG', 'LTB', 'IL32', 'CD52'], ['GNLY', 'PTGDS', 'GZMB', 'S100B', 'FGFBP2', 'NKG7', 'PRF1', 'KLRF1', 'HOPX', 'CST7', 'KLRD1', 'CTSW', 'SPON2', 'IFITM1', 'GZMA', 'CD247', 'CLIC3', 'CD7', 'ADGRG1', 'CCL5', 'TRDC']]
    macparland_nk = [['CD7', 'CMC1', 'XCL2', 'KLRB1', 'XCL1', 'KLRC1', 'KLRF1', 'IL2RB', 'CD160', 'CCL3', 'KLRD1', 'NKG7', 'TXK', 'ALOX5AP', 'TRDC', 'CD69', 'TMIGD2', 'CLIC3', 'GZMK', 'DUSP2', 'MATK', 'IFITM1', 'CCL4', 'CD247']]
    # CYP2A genes may be to exclude because have different roles in mouse and human (refer to MacParland et al., 2018)
    macparland_hep = [['SCD', 'HMGCS1', 'ACSS2', 'TM7SF2', 'TMEM97', 'CP', 'CRP', 'SLPI', 'C2orf82', 'ACAT2', 'TM4SF5', 'MSMO1', 'LEPR'], ['BCHE', 'G6PC', 'GHR', 'ALDH6A1', 'RCAN1', 'AR', 'RP4-710M16.2', 'LINC00261', 'PLIN1', 'RP11-390F4.3'], ['RPP25L', 'HSD11B1', 'HAMP', 'GHR', 'APOM', 'APOC4-APOC2', 'TKFC', 'G6PC', 'G0S2', 'PON3', 'C1orf53', 'TTC36', 'GOLT1A', 'RCAN1', 'RP4-710M16.2', 'FST', 'MCC', 'AQP9', 'PLIN1'], ['HPR', 'GSTA2', 'AKR1C1', 'MASP2', 'NNT', 'SAA4', 'MRPS18C', 'OCIAD1', 'APOA5', 'TTR'], ['CYP2A7', 'ENTPD5', 'CYP3A7', 'CYP2A6', 'C4B', 'EID2', 'TP53INP2', 'SULT1A1', 'ATIC', 'SERPINH1', 'SAMD5', 'GRB14'], ['SEC16B', 'SLBP', 'RND3', 'ABCD3', 'RHOB', 'EPB41L4B', 'GPAT4', 'SPTBN1', 'SDC2', 'PHLDA1', 'WTAP', 'ACADM'],]
    macparland_B = [['MS4A1', 'LTB', 'CD37', 'CD79B', 'CD52', 'HLA-DQB1', 'TNFRSF13C', 'TCL1A', 'LINC00926', 'STAG3', 'IGHD', 'BANK1', 'IRF8', 'BIRC3', 'P2RX5', 'RP11-693J15.5', 'RP5-887A10.1', 'VPREB3', 'CD22', 'CD74', 'SELL'], ['IGLC2', 'IGHG1', 'IGKC', 'IGHG2', 'IGHG3', 'IGHGP', 'IGLC3', 'JCHAIN', 'IGHA1', 'IGHG4', 'IGHA2', 'IGHM', 'IGLV3-1', 'IGLC7', 'MZB1', 'CD79A', 'SSR4', 'IL16']]

    loc = locals()
    signature = {
        m: loc[m] for m in loc.keys() if m.startswith('macparland_')
    }

    if collapse:
        for k, v in signature.items():
            signature[k] = list(
                set([j for i in signature[k] for j in i])
            )

    return signature


def liver_signatures(target, collapse=True, top=50):
    signatures = {}

    signatures_s = (
        (RS007_signature(top=top), 'mmusculus'),
        (RS025_signature(top=top), 'mmusculus'),
        (misc_mouse_signature(), 'mmusculus'),
        (macparland_signature(collapse=True), 'hsapiens'),
        (han_signature(), 'mmusculus')
    )

    for signature, org in signatures_s:
        if org != target:
            ms = signature.keys()
            signature = {
                m: orthologs(signature[m], org, target) for m in ms
            }
        signatures.update(signature)

    return signatures

def brain_signatures(target, collapse=True):
    signatures = {}

    signatures_s = (
        (abcam_signature(), 'mmusculus'),
        (hippo_signature(), 'mmusculus'),
    )

    for signature, org in signatures_s:
        if org != target:
            ms = signature.keys()
            signature = {
                m: orthologs(signature[m], org, target) for m in ms
            }
        signatures.update(signature)

    return signatures

def prepare_dotplot(signatures):
    var_group_labels = list(signatures.keys())
    i = 0
    var_group_positions = []
    for m in var_group_labels:
        l = len(signatures[m])
        var_group_positions += [(i, i)]
        i += l
    var_names = [m for ml in var_group_labels for m in signatures[ml]]
    return var_names, var_group_labels, var_group_positions

