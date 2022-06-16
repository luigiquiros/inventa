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
import umap
from plotly.subplots import make_subplots
import ipywidgets as widgets

## visualization for the similarity component

def pcoa_2d(matrix, data, metric, filename_header):
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

    fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('PCoA MEMO aligned (IF)', 'PCoA MEMO aligned (LOF)', 'PCoA MEMO aligned (OCSVM)'))


    #prepare table and calculate pcoa: 
    matrix= matrix.iloc[:, 1:]
    dist_matrix = sp.spatial.distance.pdist(matrix, metric)
    pcoa_results = pcoa(dist_matrix)

    exp_var_pc1 = round(100*pcoa_results.proportion_explained[0], 1)
    exp_var_pc2 = round(100*pcoa_results.proportion_explained[1], 1)
    

    # first plot: Isolation Forest 
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_IF'], 
                    colorscale='YlOrRd',
                    showscale=False,
                ),
            name='[-1]',
            legendgrouptitle_text="Anomaly",
            hovertext=data[filename_header]),
            row=1, col=1)

    # Secon plot: Local factor outlier
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_LOF'], 
                    colorscale='YlOrRd',
                    showscale=False,
                ),
            legendgrouptitle_text="normal",
            showlegend=False,
            hovertext=data[filename_header]),
            row=1, col=2)

    # Third plot: One Class Suport Vector machine
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_OCSVM'], 
                    colorscale='YlOrRd',
                    showscale=False
                ),
            legendgrouptitle_text="Anomaly",
            showlegend=False,
            hovertext=data[filename_header]),
            row=1, col=3)

    fig.update_layout(height=600, width=1500, template = 'simple_white', title_text="Multiple outlier detection results")
    
    # Update xaxis properties
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=1)
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=2)
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=3)

    # Update yaxis properties
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=1)
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=2)
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=3)
    
    fig.update_layout(
        font=dict(
        family="Times New Roman",
        size=16,
        color="black")
    )
    fig.update_traces(marker=dict(size=10)
    )
    fig.update_layout(
    legend=dict(
        orientation="h",
        y= -0.2,
        x = -0.05,
        bordercolor="Black",
        borderwidth=1
        ),
    font=dict(
        family="Times New Roman",
        size=16,
        color="black"
        )    
    )
    fig.update_annotations(font_size=20)
    fig.write_html("../data_out/PCoA_2D.html")    
    fig.show()

def pcoa_3d(matrix, data, metric = 'braycurtis'):#, filename_header):
    
    fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    specs=[[{"type": "scene"}, {'type': "scene"}, {'type': "scene"}]],
                    subplot_titles=('PCoA MEMO aligned (IF)', 'PCoA MEMO aligned (LOF)', 'PCoA MEMO aligned (OCSVM)'))


    #prepare table and calculate pcoa: 
    matrix= matrix.iloc[:, 1:]
    dist_matrix = sp.spatial.distance.pdist(matrix, metric)
    pcoa_results = pcoa(dist_matrix)

    exp_var_pc1 = round(100*pcoa_results.proportion_explained[0], 1)
    exp_var_pc2 = round(100*pcoa_results.proportion_explained[1], 1)
    exp_var_pc3 = round(100*pcoa_results.proportion_explained[2], 1)

    #data['anomaly_IF'] = data['anomaly_IF'].astype(str)
    #data['anomaly_LOF'] = data['anomaly_LOF'].astype(str)
    #data['anomaly_OCSVM'] = data['anomaly_OCSVM'].astype(str)


    # first plot: Isolation Forest 
    fig.add_trace(
        go.Scatter3d(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'], z=pcoa_results.samples['PC3'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_IF'], 
                    colorscale='YlOrRd',
                    showscale=False,
                ),
            name='[-1]',
            legendgrouptitle_text="Anomaly",
            hovertext=data[filename_header]),
            row=1, col=1)

    # Secon plot: Local factor outlier
    fig.add_trace(
        go.Scatter3d(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'], z=pcoa_results.samples['PC3'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_LOF'], 
                    colorscale='YlOrRd',
                    showscale=False,
                ),
            showlegend=False,
            hovertext=data[filename_header]),
            row=1, col=2)
 
    # Third plot: One Class Suport Vector machine
    fig.add_trace(
        go.Scatter3d(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'], z=pcoa_results.samples['PC3'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_OCSVM'], 
                    colorscale='YlOrRd',
                    showscale=False,
                ),
            showlegend=False,
            legendgrouptitle_text="Anomaly",
            hovertext=data[filename_header]),
            row=1, col=3)

    # Update xaxis properties

    fig.update_layout(scene = dict(
                    xaxis_title=f"PC1 ({exp_var_pc1} %)",
                    yaxis_title=f"PC2 ({exp_var_pc2} %)",
                    zaxis_title=f"PC2 ({exp_var_pc2} %)",
                    ), scene2 = dict(
                    xaxis_title=f"PC1 ({exp_var_pc1} %)",
                    yaxis_title=f"PC2 ({exp_var_pc2} %)",
                    zaxis_title=f"PC2 ({exp_var_pc2} %)",
                    ), scene3 = dict(
                    xaxis_title=f"PC1 ({exp_var_pc1} %)",
                    yaxis_title=f"PC2 ({exp_var_pc2} %)",
                    zaxis_title=f"PC2 ({exp_var_pc2} %)",
                    )
    )

    fig.update_layout(height=600, width=1500, template = 'simple_white', title_text="Multiple outlier detection results")
    
    fig.update_layout(
        font=dict(
        family="Times New Roman",
        size=10,
        color="black")
    )

    fig.update_layout(
    legend=dict(
        orientation="h",
        y=-0.15,
        x = 0.0,
        bordercolor="Black",
        borderwidth=1
        ),
    font=dict(
        family="Times New Roman",
        size=16,
        color="black"
        )    
    )
    fig.update_annotations(font_size=20)    
    fig.write_html("../data_out/PCoA_3D.html") 
    fig.show()

def umap_2d(matrix, data, metadata, filename_header):
    
    # 1. Generate df used for plotting 

    memo_aligned_matrix = matrix.iloc[:, 1:]
    memo_aligned_matrix = memo_aligned_matrix.div(memo_aligned_matrix.sum(axis=1), axis=0)


    results_umap = metadata.copy()
    #results_umap = results_umap.dropna()
    metric ='MEMO aligned (Bray-Curtis)'
    matrix = memo_aligned_matrix
    reducer = umap.UMAP(metric='braycurtis')
    embedding = reducer.fit_transform(matrix.values)
    results_umap[metric+'_x'] = embedding[:, 0]
    results_umap[metric+'_y'] = embedding[:, 1]
    colors = 'YlOrRd_r'

    fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('UMAP MEMO aligned: Isolation Forest', 'UMAP MEMO aligned: LOF', 'UMAP MEMO aligned: OCSVM'))

    # first plot: Isolation Forest 
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_IF'], 
                    colorscale=colors,
                    showscale=False,
                ),
            name='Anomaly [-1]',
            #legendgrouptitle_text="Anomaly",
            hovertext=results_umap[filename_header]),
            row=1, col=1)

    # Secon plot: Local factor outlier
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_LOF'], 
                    colorscale=colors,
                    showscale=False,
                ),
            name = 'Normal [1]',
            #legendgrouptitle_text="Normal",
            showlegend=True,
            hovertext=results_umap[filename_header]),
            row=1, col=2)

    # Third plot: One Class Suport Vector machine
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_OCSVM'], 
                    colorscale=colors,
                    showscale=False
                ),
            legendgrouptitle_text="Anomaly",
            showlegend=False,
            hovertext=results_umap[filename_header]),
            row=1, col=3)

    fig.update_layout(height=600, width=1500, template = 'simple_white', title_text="Multiple outlier detection results")
    
    fig.update_layout(
        font=dict(
        family="Times New Roman",
        size=16,
        color="black")
    )
    fig.update_traces(marker=dict(size=10)
    )
    fig.update_layout(
    legend=dict(
        orientation="h",
        y=-0.15,
        x = 0.0,
        bordercolor="Black",
        borderwidth=1
        ),
    font=dict(
        family="Times New Roman",
        size=16,
        color="black"
        )    
    )
    fig.update_annotations(font_size=20)
    fig.write_html("../data_out/UMAP_2D.html")     
    fig.show()

    
def pcoa_umap_2d(matrix, data, metric, filename_header):
    fig = make_subplots(rows=2, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.12,
                    subplot_titles=('PCoA MEMO aligned (IF)', 'PCoA MEMO aligned (LOF)', 'PCoA MEMO aligned (OCSVM)', 'UMAP MEMO aligned (IF)', 'UMAP MEMO aligned (LOF)', 'UMAP MEMO aligned (OCSVM)')
    )

    colors = 'YlOrRd_r'
    #prepare table and calculate pcoa: 
    matrix= matrix.iloc[:, 1:]
    dist_matrix = sp.spatial.distance.pdist(matrix, metric)
    pcoa_results = pcoa(dist_matrix)

    exp_var_pc1 = round(100*pcoa_results.proportion_explained[0], 1)
    exp_var_pc2 = round(100*pcoa_results.proportion_explained[1], 1)
    

    # first plot: Isolation Forest 
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_IF'], 
                    colorscale=colors,
                    showscale=False,
                ),
            showlegend=False,
            #legendgrouptitle_text="Anomaly",
            hovertext=data[filename_header]),
            row=1, col=1)

    # Secon plot: Local factor outlier
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_LOF'], 
                    colorscale=colors,
                    showscale=False,
                ),
            #name='Anomaly [-1]',
            #legendgrouptitle_text="normal",
            showlegend=False,
            hovertext=data[filename_header]),
            row=1, col=2)

    # Third plot: One Class Suport Vector machine
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=pcoa_results.samples['PC1'], y=pcoa_results.samples['PC2'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_OCSVM'], 
                    colorscale=colors,
                    showscale=False
                ),
            #legendgrouptitle_text="Anomaly",
            showlegend=False,
            hovertext=data[filename_header]),
            row=1, col=3)

    fig.update_layout(height=1200, width=1500, template = 'simple_white', title_text="Multiple outlier detection results")
    
    # Update xaxis properties
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=1)
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=2)
    fig.update_xaxes(title_text=f"PC1 ({exp_var_pc1} %)", row=1, col=3)

    # Update yaxis properties
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=1)
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=2)
    fig.update_yaxes(title_text=f"PC2 ({exp_var_pc2} %)", row=1, col=3)
    

    # 1. Generate df used for plotting 

    memo_aligned_matrix = matrix.iloc[:, 1:]
    memo_aligned_matrix = memo_aligned_matrix.div(memo_aligned_matrix.sum(axis=1), axis=0)


    results_umap = data.copy()
    results_umap = results_umap.dropna()
    metric ='MEMO aligned (Bray-Curtis)'
    matrix = memo_aligned_matrix
    reducer = umap.UMAP(metric='braycurtis')
    embedding = reducer.fit_transform(matrix.values)
    results_umap[metric+'_x'] = embedding[:, 0]
    results_umap[metric+'_y'] = embedding[:, 1]


    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_IF'], 
                    colorscale=colors,
                    showscale=False,
                ),
            name='Outlier [-1]',
            #legendgrouptitle_text="Anomaly",
            hovertext=results_umap[filename_header]),
            row=2, col=1)

    # Secon plot: Local factor outlier
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_LOF'], 
                    colorscale=colors,
                    showscale=False,
                ),
                name='Normal [1]',
            #legendgrouptitle_text="normal",
            showlegend=True,
            hovertext=results_umap[filename_header]),
            row=2, col=2)

    # Third plot: One Class Suport Vector machine
    fig.add_trace(
        go.Scatter(
            opacity=0.85,
            mode='markers',
            x=results_umap[metric+'_x'], y=results_umap[metric+'_y'],
            marker=dict(
                    size = 5,
                    line_width=0.5,
                    line_color ='black',
                    color=data['anomaly_OCSVM'], 
                    colorscale=colors,
                    showscale=False
                ),
            legendgrouptitle_text="Anomaly",
            showlegend=False,
            hovertext=results_umap[filename_header]),
            row=2, col=3)

    fig.update_layout(
        font=dict(
        family="Times New Roman",
        size=16,
        color="black")
    )
    fig.update_traces(marker=dict(size=10)
    )
    fig.update_layout(
        legend=dict(
        orientation="h",
        y= -0.05,
        x = -0.05,
        bordercolor="Black",
        borderwidth=1
        ),
    font=dict(
        family="Times New Roman",
        size=16,
        color="black"
        )    
    )
    fig.update_annotations(font_size=20)
    fig.write_html("../data_out/PCoA_UMAP_2D.html")     
    fig.show()

def hist_to_plot(sample, quantitative_data_filename, annotation_df, reduced_df, min_specificity, annotation_preference):

    # 1) recover row ID, m/z, rt nad add annotation status

    df=pd.read_csv(quantitative_data_filename, sep=',', usecols=['row ID','row m/z', 'row retention time'] ,index_col='row ID')
    df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    df.rename(columns = lambda x: x.replace('row retention time', 'retention time (min)'),inplace=True)
    #df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    df.reset_index(inplace=True)
    df= pd.merge(df, annotation_df[['cluster index', 'annotation']], how='left', right_on='cluster index', left_on='row ID').fillna(0)
    df.drop('cluster index', axis=1, inplace=True)

    #2) normalize the filtered table and combine information from specificity and annotation status for each feature
        #normalize row-wise the area of features = relative % of each feature in each sample

    reduced_df_norm = reduced_df.copy()#.transpose()
    reduced_df_norm = reduced_df_norm.div(reduced_df_norm.sum(axis=1), axis=0).fillna(0)
    reduced_df_norm.reset_index(inplace=True)
    reduced_df_norm.rename(columns={'index': 'row ID'}, inplace=True)

    #merge and check status by sample
    df_check = pd.merge(df, reduced_df_norm[['row ID', sample]], how='left', on= 'row ID')
    df_check['status'] = np.where(( (df_check[sample] > min_specificity) & (df_check['annotation'] == annotation_preference)), 'interesting', 'not interesting' )
    df_check['row m/z'].round(decimals = 4)
    #plot
    
    fig = px.histogram(df_check,
                            x='retention time (min)', y=sample,
                            nbins=1000,
                            #histnorm= 'probability density',
                            opacity=0.8, 
                            labels={'retention time (min)':'retention time range (min)'},
                            title='Pseudo chromatogram of sample',
                            template="simple_white",
                            #facet_col="row ID",
                            color ='status',
                            color_discrete_map= {'interesting': '#1E88E5', 'not interesting': '#FFC107'},
                            marginal="rug", # can be 'rug' `box`, `violin`,
                            hover_data=df_check.columns
                            )
                            #, name='row m/z')
    fig.update_layout( # customize font and legend orientation & position
        font_family="Times New Roman",
        legend=dict(
            title=None, orientation="h", y=1, yanchor="bottom", x=0.5, xanchor="center"
        ))
    fig.update_layout(
    autosize=False,
    width=1400,
    height=500)
    
    fig.show()
    
def drop_selection(quant_df):
    
    drop_down = widgets.Dropdown(options=quant_df.columns,
                                    description='Sample to plot',
                                    disabled=False)
    sample = '0'
    def dropdown_handler(change):
        global sample
        print(change.new)
        sample = change.new  # This line isn't working
    drop_down.observe(dropdown_handler, names='value')
    display(drop_down)

def quant_plot(df):
        """ Cleans up the quantitative table to specific format

    Args:
        df = quantitative.csv file, output from MZmine

    Returns:
        None
    """
    df.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    df.rename(columns = lambda x: x.replace('row retention time', 'retention time (min)'),inplace=True)
    df.drop(list(df.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    #df.drop('row m/z', axis=1, inplace=True)
    #df.drop('row retention time', axis=1, inplace=True)
    #df.to_csv('../data_out/quant_df.tsv', sep='\t')
    return df