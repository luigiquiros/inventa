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
            #name='[-1]',
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

def pcoa_3d(matrix, data,filename_header, metric = 'braycurtis'):
    
    fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    specs=[[{"type": "scene"}, {'type': "scene"}, {'type': "scene"}]],
                    subplot_titles=('PCoA MEMO aligned (IF)', 'PCoA MEMO aligned (LOF)', 'PCoA MEMO aligned (OCSVM)'))

    colors = 'YlOrRd_r'
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
                    colorscale=colors,
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
                    colorscale=colors,
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
                    colorscale=colors,
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


    results_umap = data.copy()
    results_umap = results_umap.dropna()
    metric ='MEMO aligned (Bray-Curtis)'
    matrix = memo_aligned_matrix
    reducer = umap.UMAP(metric='braycurtis')
    embedding = reducer.fit_transform(matrix.values)
    results_umap[metric+'_x'] = embedding[:, 0]
    results_umap[metric+'_y'] = embedding[:, 1]

    fig = make_subplots(rows=1, cols=3,
                    shared_xaxes=False,
                    vertical_spacing=0.4,
                    subplot_titles=('UMAP MEMO aligned: Isolation Forest', 'UMAP MEMO aligned: LOF', 'UMAP MEMO aligned: OCSVM'))
    colors = 'YlOrRd_r'
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
  
    df_check.to_csv('../data_out/sample.tsv', sep='\t')
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

def distribution_to_plot(sample, quant_df, reduced_df):
    df1 = quant_df.copy()
    df1.reset_index(inplace=True)
    df1 = df1.replace(0, np.nan)
    df2 = reduced_df.copy()
    df2.reset_index(inplace=True)
    df2 = df2.replace(0, np.nan)

    fig = go.Figure()

    fig.add_trace(go.Violin(y=df1[sample],
                            legendgroup='Non-filtered', scalegroup='Non-filtered', name='Non-filtered',
                            #side='negative',
                            line_color='lightseagreen')
                )
    fig.add_trace(go.Violin(y=df2[sample],
                            legendgroup='Filtered', scalegroup='Filtered', name='Filtered',
                            #side='positive',
                            line_color='mediumpurple')
                )
    # update characteristics shared by all traces
    fig.update_traces(meanline_visible=True,
                    points='all', # show all points
                    jitter=0.05,  # add some jitter on points for better visibility
                    scalemode='count',
                    box_visible=True) #scale violin plot area with total count
    fig.update_layout(
        title_text="Intensity distribution<br><i>before and post filtering",
        violinmode='overlay')
    fig.update_layout(
        autosize=True,
        width=700,
        height=700)
    fig.write_html("../data_out/filtering_plot.html") 
    fig.show()
    
def pseudochromatogram(sample, quantitative_data_filename, annotation_df, metadata_df, reduced_df,  min_specificity, annotation_preference, species_column, organe_column, CC_component, canopus_npc_summary_filename,  min_class_confidence, sirius_annotations, sirius_annotations_filename, min_ConfidenceScore, min_ZodiacScore, use_ion_identity, correlation_groups_df, data_process_origin, filename_header):
    
    if use_ion_identity == True: 
        row_ID_header = 'annotation network number'
    else: 
        row_ID_header = 'row ID'
        
    #get a clean table with the original intensities normalized
    dfq = pd.read_csv(quantitative_data_filename, sep=',')
    dfq.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    dfq.drop(list(dfq.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    dfq.rename(columns ={'row retention time':'retention time (min)'}, inplace=True)
    dfq[sample] = dfq[sample]/dfq[sample].sum()  #normalize to 1 the areas of the particular sample    

    if data_process_origin == 'MZMine3':
        
            if use_ion_identity == True:
        
                dfq = dfq[['row ID', row_ID_header, sample]]
                #complete correlation groups
                dfq[row_ID_header] = dfq[row_ID_header].fillna(dfq['row ID'].apply(str) + 'x')
                dfq.drop('row ID', axis =1, inplace=True)
                dfq = dfq.groupby(row_ID_header, dropna=False).max()

                #recover information from the correlation groups:
                agg_func = {'retention time (min)': 'mean', 'row m/z': 'max',  'adduct (ion identity)': set, 'row ID': set, 'neutral mass (ion identity)': 'max'}
                dfcg = correlation_groups_df.groupby('annotation network number', as_index=False).agg(agg_func)
                dfcg[['adduct (ion identity)', 'row ID']] = dfcg[['adduct (ion identity)', 'row ID']].astype(str)   
                
                #merge with the main data according to sample
                dfq = pd.merge(dfq, dfcg[['annotation network number', 'retention time (min)', 'row m/z', 'row ID', 'adduct (ion identity)','neutral mass (ion identity)']], how ='left', left_on = row_ID_header, right_on=row_ID_header)
                
                #add annotation status 
                df = dfq
                df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)
                df.fillna({'adduct (ion identity)': 'not available', 'neutral mass (ion identity)': 'not available'}, inplace=True)

            else:
                df = dfq[[row_ID_header,'row m/z', 'retention time (min)', sample]]
                df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)
    else:
        df = dfq[[row_ID_header, 'row m/z', 'retention time (min)', sample]]
        df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)

    #2) normalize the filtered table and combine information from specificity and annotation status for each feature
                #normalize row-wise the area of features = relative % of each feature in each sample

    reduced_df_norm = reduced_df.copy()#.transpose()
    reduced_df_norm = reduced_df_norm.div(reduced_df_norm.sum(axis=1), axis=0).fillna(0)
    reduced_df_norm.reset_index(inplace=True)
    reduced_df_norm.rename(columns={'index': row_ID_header}, inplace=True)

    #change the area for the post filter-row wise normalized one
    df_check = df
    df_check[sample] = reduced_df_norm[sample]
    df_check['status'] = np.where(
    ((df_check[sample] > min_specificity) & (df_check['annotation'] == annotation_preference)), 'specific non annotated', 
    (np.where (((df_check[sample] > min_specificity) & (df_check['annotation'] != annotation_preference)),'specific annotated', 'not interesting'))
    )
    df_check.annotation = df_check.annotation.map({1: 'annotated', 0: 'non anotated'})
    #change intensity for the intensity before filtering
    df_check[sample] = dfq[sample]
    #erase all zero values (easy to plot)
    df_check= df_check[df_check[sample] != 0]
    df_check['retention time (min)'] = df_check['retention time (min)'].round(decimals=2)
    df_check = df_check.sort_values(by='retention time (min)', ascending=True)
    #df_check['retention time (min)'] = df_check['retention time (min)'].astype(str)
    df_check['row m/z'] = df_check['row m/z'].round(decimals=4)


    if sirius_annotations == True and use_ion_identity == False:
                
        annot_sirius_df=pd.read_csv(sirius_annotations_filename,
                                        sep='\t', 
                                        usecols =['id','ConfidenceScore','ZodiacScore', 'molecularFormula', 'adduct', 'name'], 
                                        low_memory=False)
        annot_sirius_df['shared name'] = annot_sirius_df['id'].str.split('_').str[-1].astype(int)
        annot_sirius_df.drop('id', axis=1, inplace=True)
        annot_sirius_df.rename(columns={'shared name': 'row ID', 'adduct': 'adduct (sirius)', 'molecularFormula': 'MF (sirius)', 'name': 'Compound name (sirius)'}, inplace=True) 
        annot_sirius_df.drop(annot_sirius_df[(annot_sirius_df.ConfidenceScore >= min_ConfidenceScore) & (annot_sirius_df.ZodiacScore >= min_ZodiacScore)].index, inplace=True) 
        annot_sirius_df.drop('ConfidenceScore', axis=1, inplace=True)         
        annot_sirius_df.drop('ZodiacScore', axis=1, inplace=True) 
        annot_sirius_df['adduct (sirius)'] = annot_sirius_df['adduct (sirius)'].astype(str)

        df_check = pd.merge(df_check, annot_sirius_df, how='left', on= 'row ID')
        df_check.fillna({'MF (sirius)': 'not available', 'adduct (sirius)': 'not available', 'Compound name (sirius)': 'not available'}, inplace=True)
        #df_check.drop('classProbability', axis=1, inplace=True)

    else:
        df_check

    if CC_component == True and use_ion_identity == False:

        canopus_npc_df=pd.read_csv(canopus_npc_summary_filename, sep=',', usecols=['name', 'pathway', 'superclass', 'class', 'classProbability'])
        canopus_npc_df.rename(columns={'name': 'row ID'}, inplace=True)
        #filter by class probability
        canopus_npc_df.drop(canopus_npc_df[canopus_npc_df.classProbability > min_class_confidence].index, inplace=True)
        df_check = pd.merge(df_check, canopus_npc_df, how='left', on= 'row ID')
        df_check.fillna({'pathway': 'not available', 'superclass': 'not available','class': 'not available' }, inplace=True)
        df_check.drop('classProbability', axis=1, inplace=True)
    else:
        df_check    
        df_check.to_csv('../data_out/Interesting_features_for'+sample+'.tsv', sep='\t')
        
    #recover species and organe for the particular sample 
    species=metadata_df.loc[metadata_df[filename_header] == sample, species_column].item()
    organism_part = metadata_df.loc[metadata_df[filename_header] == sample, organe_column].item()

    #plot 
    fig = px.bar(df_check,
                            x='retention time (min)', y=sample,
                            #histnorm= 'probability density',
                            #opacity=0.8, 
                            #labels={'retention time (min)':'retention time range (min)'},
                            title='Pseudo chromatogram for '+sample+'<br>species: <i>'+species+'<i><br>organism part: '+organism_part+'</sup>',

                            template="simple_white",
                            color ='status',
                            color_discrete_map= {'specific non annotated': '#1E88E5', 'not interesting': '#FFC107', 'specific annotated':'#004D40'},
                            hover_data=df_check.columns,
                            facet_col_wrap=2
                            )

    fig.update_layout( # customize font and legend orientation & position
        font_family="Times New Roman",
        legend=dict(
            title=None, orientation="h", y=1, yanchor="bottom", x=0.5, xanchor="center"
        ))
    fig.update_layout(
    autosize=True,
    width=1400,
    height=500,
    barmode ='group',
    bargap=0.15, # gap between bars of adjacent location coordinates.
    bargroupgap=0.1
    )
    fig.update_xaxes(title_text='retention time (min)',showgrid=False, ticks="outside", tickson="boundaries")
    fig.update_yaxes(title_text='relative intensity')
    fig.write_html("../data_out/pseudochromato.html")  
    #fig.write_image("../data_out/pseudochromato.svg")
    fig.show()
    

def chromatogram2D(sample, quantitative_data_filename, annotation_df, metadata_df, reduced_df,  min_specificity, annotation_preference, species_column, organe_column, CC_component, canopus_npc_summary_filename, min_class_confidence, sirius_annotations, sirius_annotations_filename, min_ConfidenceScore, min_ZodiacScore, use_ion_identity, correlation_groups_df, data_process_origin, filename_header):
    
    if use_ion_identity == True: 
        row_ID_header = 'annotation network number'
    else: 
        row_ID_header = 'row ID'
        
    #get a clean table with the original intensities normalized
    dfq = pd.read_csv(quantitative_data_filename, sep=',')
    dfq.rename(columns = lambda x: x.replace(' Peak area', ''),inplace=True)
    dfq.drop(list(dfq.filter(regex = 'Unnamed:')), axis = 1, inplace = True)
    dfq.rename(columns ={'row retention time':'retention time (min)'}, inplace=True)
    dfq[sample] = dfq[sample]/dfq[sample].sum()  #normalize to 1 the areas of the particular sample    

    if data_process_origin == 'MZMine3':
        
            if use_ion_identity == True:
        
                dfq = dfq[['row ID', row_ID_header, sample]]
                #complete correlation groups
                dfq[row_ID_header] = dfq[row_ID_header].fillna(dfq['row ID'].apply(str) + 'x')
                dfq.drop('row ID', axis =1, inplace=True)
                dfq = dfq.groupby(row_ID_header, dropna=False).max()

                #recover information from the correlation groups:
                agg_func = {'retention time (min)': 'mean', 'row m/z': 'max',  'adduct (ion identity)': set, 'row ID': set, 'neutral mass (ion identity)': 'max'}
                dfcg = correlation_groups_df.groupby('annotation network number', as_index=False).agg(agg_func)
                dfcg[['adduct (ion identity)', 'row ID']] = dfcg[['adduct (ion identity)', 'row ID']].astype(str)   
                
                #merge with the main data according to sample
                dfq = pd.merge(dfq, dfcg[['annotation network number', 'retention time (min)', 'row m/z', 'row ID', 'adduct (ion identity)','neutral mass (ion identity)']], how ='left', left_on = row_ID_header, right_on=row_ID_header)
                
                #add annotation status 
                df = dfq
                df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)
                df.fillna({'adduct (ion identity)': 'not available', 'neutral mass (ion identity)': 'not available'}, inplace=True)

            else:
                df = dfq[[row_ID_header,'row m/z', 'retention time (min)', sample]]
                df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)
    else:
        df = dfq[[row_ID_header, 'row m/z', 'retention time (min)', sample]]
        df= pd.merge(df, annotation_df[[row_ID_header, 'annotation']], how='left', on=row_ID_header).fillna(0)


    #2) normalize the filtered table and combine information from specificity and annotation status for each feature
                #normalize row-wise the area of features = relative % of each feature in each sample

    reduced_df_norm = reduced_df.copy()#.transpose()
    reduced_df_norm = reduced_df_norm.div(reduced_df_norm.sum(axis=1), axis=0).fillna(0)
    reduced_df_norm.reset_index(inplace=True)
    reduced_df_norm.rename(columns={'index': row_ID_header}, inplace=True)

    #change the area for the post filter-row wise normalized one
    df_check = df
    df_check[sample] = reduced_df_norm[sample]
    df_check['status'] = np.where(
    ((df_check[sample] > min_specificity) & (df_check['annotation'] == annotation_preference)), 'specific non annotated', 
    (np.where (((df_check[sample] > min_specificity) & (df_check['annotation'] != annotation_preference)),'specific annotated', 'not interesting'))
    )
    df_check.annotation = df_check.annotation.map({1: 'annotated', 0: 'non anotated'})
    #change intensity for the intensity before filtering
    df_check[sample] = dfq[sample]
    #erase all zero values (easy to plot)
    df_check= df_check[df_check[sample] != 0]
    df_check['retention time (min)'] = df_check['retention time (min)'].round(decimals=2)
    df_check = df_check.sort_values(by='retention time (min)', ascending=True)
    #df_check['retention time (min)'] = df_check['retention time (min)'].astype(str)
    #if use_ion_identity == False:
    df_check['row m/z'] = df_check['row m/z'].round(decimals=4)
    #else:
    #    df_check


    if sirius_annotations == True and use_ion_identity == False:
                
        annot_sirius_df=pd.read_csv(sirius_annotations_filename,
                                        sep='\t', 
                                        usecols =['id','ConfidenceScore','ZodiacScore', 'molecularFormula', 'adduct', 'name'], 
                                        low_memory=False)
        annot_sirius_df['shared name'] = annot_sirius_df['id'].str.split('_').str[-1].astype(int)
        annot_sirius_df.drop('id', axis=1, inplace=True)
        annot_sirius_df.rename(columns={'shared name': 'row ID', 'adduct': 'adduct (sirius)', 'molecularFormula': 'MF (sirius)', 'name': 'Compound name (sirius)'}, inplace=True) 
        annot_sirius_df.drop(annot_sirius_df[(annot_sirius_df.ConfidenceScore >= min_ConfidenceScore) & (annot_sirius_df.ZodiacScore >= min_ZodiacScore)].index, inplace=True) 
        annot_sirius_df.drop('ConfidenceScore', axis=1, inplace=True)         
        annot_sirius_df.drop('ZodiacScore', axis=1, inplace=True) 
        annot_sirius_df['adduct (sirius)'] = annot_sirius_df['adduct (sirius)'].astype(str)

        df_check = pd.merge(df_check, annot_sirius_df, how='left', on= 'row ID')
        df_check.fillna({'MF (sirius)': 'not available', 'adduct (sirius)': 'not available', 'Compound name (sirius)': 'not available'}, inplace=True)
        #df_check.drop('classProbability', axis=1, inplace=True)

    else:
        df_check

    if CC_component == True and use_ion_identity == False:

        canopus_npc_df=pd.read_csv(canopus_npc_summary_filename, sep=',', usecols=['name', 'pathway', 'superclass', 'class', 'classProbability'])
        canopus_npc_df.rename(columns={'name': 'row ID'}, inplace=True)
        #filter by class probability
        canopus_npc_df.drop(canopus_npc_df[canopus_npc_df.classProbability > min_class_confidence].index, inplace=True)
        df_check = pd.merge(df_check, canopus_npc_df, how='left', on= 'row ID')
        df_check.fillna({'pathway': 'not available', 'superclass': 'not available','class': 'not available' }, inplace=True)
        df_check.drop('classProbability', axis=1, inplace=True)
    else:
        df_check    
        df_check.to_csv('../data_out/Interesting_features_for'+sample+'.tsv', sep='\t')
        
    #recover species and organe for the particular sample 
    species=metadata_df.loc[metadata_df[filename_header] == sample, species_column].item()
    organism_part = metadata_df.loc[metadata_df[filename_header] == sample, organe_column].item()

    #plot 
    
    fig = px.scatter(df_check,
                            x='retention time (min)', y='row m/z',
                            size= sample,
                            title='2D chromatogram plot for '+sample+'<br>species: <i>'+species+'<i><br>organism part: '+organism_part+'</sup>',
                            template="simple_white",
                            color ='status',
                            color_discrete_map= {'specific non annotated': '#1E88E5', 'not interesting': '#FFC107', 'specific annotated':'#004D40'},
                            hover_data=df_check.columns
                            )

    fig.update_layout( # customize font and legend orientation & position
        font_family="Times New Roman",
        legend=dict(
            title=None, orientation="h", y=1, yanchor="bottom", x=0.5, xanchor="center"
        ))
    fig.update_layout(
    autosize=True,
    width=1400,
    height=500)
    fig.update_xaxes(title_text='retention time (min)',showgrid=False, ticks="outside", tickson="boundaries")
    fig.update_yaxes(title_text='m/z max')
    fig.write_html("../data_out/chromato2D.html") 
    #fig.write_image("../data_out/chromato2D.jpg")
    fig.show()