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
    dist_matrix = sp.spatial.distance.pdist(matrix, 'braycurtis')

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
