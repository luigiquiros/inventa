import pandas as pd
import numpy as np
import scipy as sp
from itertools import cycle
from skbio.stats.ordination import pcoa
import cimcb_lite as cb
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

## visualization for the similarity component

def pcoa_2d(filename_col, group_col, title, matrix, data, metric = 'braycurtis'):
    """ Simple 2D PCoA plot of a MEMO matrix using Plotly showing the anomalies

    Args:
        matrix (DataFrame): A DataFrame derived from memo
        data (DataFrame): a dataframe to use as metadata, in this case the SC df
        filename_col (str): Column name in data to match memo_matrix index
        group_col (str): Column name in data to use as groups for plotting
        metric (str, optional): Distance metric to use, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html. Defaults to 'braycurtis'.

    Returns:
        None
    """
    #prepare table and calculate pcoa: 
    matrix= matrix.iloc[: , 1:]
    dist_matrix = sp.spatial.distance.pdist(matrix, metric)

    pcoa_results = pcoa(dist_matrix)
    pc1_perc_var_exp = round(100*pcoa_results.proportion_explained[0], 1)
    pc2_perc_var_exp = round(100*pcoa_results.proportion_explained[1], 1)

    data['anomaly_IF'] = data['anomaly_IF'].astype(str)
    data['anomaly_LOF'] = data['anomaly_LOF'].astype(str)
    data['anomaly_OCSVM'] = data['anomaly_OCSVM'].astype(str)
        
     #figure parameters
    fig = px.scatter(x=pcoa_results.samples['PC1'],
    y=pcoa_results.samples['PC2'],
    color=data[group_col],
    labels={'x': f"PC1 ({pc1_perc_var_exp} %)",
                'y': f"PC2 ({pc2_perc_var_exp} %)",
                'color': group_col
                },
    title=title,
    hover_name=data[filename_col],
    template="simple_white",
            )
    #fig.update_layout({'width':700, 'height':650})
    fig.update_traces(marker=dict(size=10,
                                    line=dict(width=1,
                                    color='DarkSlateGrey'))
        )
    fig.show()

def pcoa_3d(filename_col, group_col, title, matrix, data, metric = 'braycurtis'):
    """ Simple 3D PCoA plot of a MEMO matrix using Plotly

    Args:
        matrix (DataFrame): A DataFrame in the MemoContainer.memo_matrix or MemoContainer.feature_table format
        df_metadata (DataFrame): Metadata of the MEMO matrix samples
        filename_col (str): Column name in df_metadata to match memo_matrix index
        group_col (str): Column name in df_metadata to use as groups for plotting
        metric (str, optional): Distance metric to use, see https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html. Defaults to 'braycurtis'.
        norm (bool, optional): Apply samples normalization. Defaults to False.
        scaling (bool, optional): Apply pareto scaling to MEMO matrix columns. Defaults to False.
        pc_to_plot (list of int, optional): PCs to plot. Defaults to [1,2,3].

    Returns:
        None
    """
    #prepare table and calculate pcoa: 
    matrix= matrix.iloc[:, 1:]
    dist_matrix = sp.spatial.distance.pdist(matrix, metric)
    pcoa_results = pcoa(dist_matrix)

    exp_var_pc1 = round(100*pcoa_results.proportion_explained[0], 1)
    exp_var_pc2 = round(100*pcoa_results.proportion_explained[1], 1)
    exp_var_pc3 = round(100*pcoa_results.proportion_explained[2], 1)

    data['anomaly_IF'] = data['anomaly_IF'].astype(str)
    data['anomaly_LOF'] = data['anomaly_LOF'].astype(str)
    data['anomaly_OCSVM'] = data['anomaly_OCSVM'].astype(str)

    fig = px.scatter_3d(x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'], z=pcoa_results.samples['PC3'],
    color=data[group_col],
    labels={'x': f"PC1 ({exp_var_pc1} %)",
            'y': f"PC2 ({exp_var_pc2} %)",
            'z': f"PC3 ({exp_var_pc3} %)",
            'color': group_col
            },
    title=title,
    hover_name=data[filename_col],
    template="simple_white"
    )
    fig.update_layout(autosize = True, width=800, height=800)
    fig.update_layout(font_family="Times New Roman")
    fig.update_traces(marker=dict(size=8,
                                    line=dict(width=1,
                                    color='DarkSlateGrey')
                                    )
        )
    fig.show()