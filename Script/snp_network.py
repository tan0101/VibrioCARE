# libraries
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys

from matplotlib.lines import Line2D
from sklearn.metrics import pairwise_distances

from matplotlib import cm
import copy
from matplotlib.legend_handler import HandlerBase

class TextHandler(HandlerBase):
    def create_artists(self, legend, orig_handle,xdescent, ydescent,
                        width, height, fontsize,trans):
        h = copy.copy(orig_handle)
        h.set_position((width/2.,height/2.))
        h.set_transform(trans)
        h.set_ha("center");h.set_va("center")
        fp = orig_handle.get_font_properties().copy()
        fp.set_size(fontsize)
        # uncomment the following line, 
        # if legend symbol should have the same size as in the plot
        h.set_font_properties(fp)
        return [h]

def update_progress(progress):
    barLength = 100 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), round(progress*100, 2), status)
    sys.stdout.write(text)
    sys.stdout.flush() 

if __name__ == "__main__":
    mpl.rc('font',family='Arial')
    name_dataset = "Vibrio" 
    folder = "Data"
    results_folder = "Results IDESHI and IEDCR" #"Results IDESHI"
    
    # Load Data Core:
    data_core_df = pd.read_csv(folder+"/"+name_dataset+'_core_genome_snps_data.csv', header = [0], index_col=[0])
    features_name_core = np.array(data_core_df.columns)
    samples_name = np.array(data_core_df.index)
    snp_matrix = np.array(data_core_df)
    snp_matrix[snp_matrix>0] = 1

    # Load Metadata:
    metadata_df = pd.read_csv(folder+"/"+name_dataset+'_metadata.csv', header = [0], index_col=[0])
    samples_metadata = np.array(metadata_df[metadata_df.columns[0]])    

    print(np.array_equal(samples_name,samples_metadata))    
    
    adj_matrix = pairwise_distances(snp_matrix,metric='hamming')
    adj_matrix[np.tril_indices(len(samples_name), 0)] = np.nan
    distribution_adj = adj_matrix.flatten()
    mask = np.where(~np.isnan(distribution_adj))[0]
    distribution_adj = data_core_df.shape[1]*distribution_adj[mask]
    print(np.mean(distribution_adj))
    print(np.amax(distribution_adj))
    print(np.median(distribution_adj))

    # First quartile (Q1) 
    Q1 = np.percentile(distribution_adj, 25, method = 'midpoint') 
    print(Q1)
    
    # Third quartile (Q3) 
    Q3 = np.percentile(distribution_adj, 75, method = 'midpoint') 
    print(Q3)
    
    # Interquaritle range (IQR) 
    IQR = Q3 - Q1 

    adj_matrix = adj_matrix*data_core_df.shape[1]
    adj_matrix[adj_matrix>15] = np.nan

    adj_df = pd.DataFrame(data=adj_matrix,index=samples_name, columns=samples_name)
    lst = adj_df.stack().reset_index()
    lst = lst.rename(columns={lst.columns[0]:"from", lst.columns[1]:"to", lst.columns[2]:"edge_value"})

    T=nx.from_pandas_edgelist(df=lst, source='from', target='to', edge_attr='edge_value', create_using=nx.Graph() )
    
    carac = metadata_df
    carac = carac.set_index("Strain Id")
    sample_carac = carac.index

    carac = carac.reindex(T.nodes())        
    
    # And I need to transform my categorical column in a numerical value: group1->1, group2->2...
    carac['Beast']=pd.Categorical(carac['Beast'])
    my_color_Node_Beast = carac['Beast'].cat.codes

    carac['Year']=pd.Categorical(carac['Year'])
    my_color_node_Date = carac['Year'].cat.codes

    carac['Serology']=pd.Categorical(carac['Serology'])
    my_color_Node_Serology = carac['Serology'].cat.codes

    carac['Location']=pd.Categorical(carac['Location'])
    my_color_Node_Location = carac['Location'].cat.codes

    lst['edge_value'] = pd.Categorical(lst['edge_value'])
    my_color_edge = lst['edge_value'].cat.codes
            
    ColorLegend_Node_Beast = {'BD-1.2': 0,'BD-2': 1}
    ColorLegend_Node_Date = {2015: 0, 2016: 1, 2017: 2, 2018: 3, 2019: 4, 2020: 5, 2021: 6}
    ColorLegend_Node_Serology = {'Ogawa': 0,'Inaba': 1}
    ColorLegend_Node_Location = {'Barisal': 0,'Chittagong': 1,'Dhaka': 2, 'Khulna': 3,
        'Rajshahi': 4, 'Sylhet': 5, "Narayanganj": 6}
    ColorLegend_Edge = {}
    Legend_Edge = {}
    #uni_edge_val = np.unique(lst["edge_value"])

    edge_value_array = []
    for edge in T.edges():
        edge_val = T.get_edge_data(edge[0], edge[1])
        edge_value_array.append(np.ceil(edge_val["edge_value"]))
    uni_edge_val = np.unique(edge_value_array)

    for count, n_gene in enumerate(uni_edge_val):
        ColorLegend_Edge[n_gene] = count
        Legend_Edge[count] = n_gene

    
    values_edge = []
    for edge in T.edges():
        edge_val = T.get_edge_data(edge[0], edge[1])
        values_edge.append(ColorLegend_Edge[np.ceil(edge_val['edge_value'])])

    values_Node_Beast = []
    values_node_Date = []
    values_Node_Serology = []
    values_node_Location = []
    for node  in T.nodes():
        values_Node_Beast.append(ColorLegend_Node_Beast[carac.loc[node,'Beast']])
        values_node_Date.append(ColorLegend_Node_Date[carac.loc[node,'Year']])
        values_Node_Serology.append(ColorLegend_Node_Serology[carac.loc[node,'Serology']])
        values_node_Location.append(ColorLegend_Node_Location[carac.loc[node,'Location']])

    values_Node_Beast = np.array(values_Node_Beast)
    values_node_Date = np.array(values_node_Date)
    values_Node_Serology = np.array(values_Node_Serology)
    values_node_Location = np.array(values_node_Location)
    
    print(len(T.nodes()))
    env_len = len(np.where(values_Node_Beast == 1)[0])
    print("BD-2 = {}".format(env_len))
    chicken_len = len(np.where(values_Node_Beast == 0)[0])
    print("BD-1.2 = {}".format(chicken_len))
   
    # compute maximum value s.t. all colors can be normalised
    maxval_Node_Beast = np.max(values_Node_Beast) 
    maxval_node_Date = np.max(values_node_Date) 
    maxval_Node_Serology = np.max(values_Node_Serology) 
    maxval_node_Location = np.max(values_node_Location) 
    maxval_edge = np.max(values_edge) 
    
    # get colormap
    cmap_Beast=cm.Paired
    cmap_Date=cm.Set1
    cmap_Location=cm.tab10
    cmap_Serology=cm.Paired

    if len(np.unique(values_edge)) < 21:
        cmap_edge=cm.terrain #Accent
    else:
        cmap_edge=cm.gist_rainbow

    pos = nx.nx_agraph.graphviz_layout(T)

    nodes = np.array(T.nodes())
    fig, ax = plt.subplots(nrows = 2, ncols=3, figsize=(16,10))
    ax = ax.ravel()
    
    ## Nodes Source

    nx.draw_networkx_nodes(T, pos, node_size=15, nodelist= nodes, node_color=[cmap_Beast(v/maxval_Node_Beast) for v in values_Node_Beast], 
        node_shape = 'o', ax=ax[0])
    
    nx.draw_networkx_edges(T,pos,alpha = 0.6, edge_color=[cmap_edge(v/maxval_edge) for v in values_edge], edge_cmap=cmap_edge, ax=ax[0], width=0.6)
    legend_elements = []
    for v in set(values_Node_Beast):
        if v == 0:
            label = "BD-1.2"
        elif v == 1:
            label = "BD-2"

        legend_elements.append(Line2D([], [], marker='o', markeredgecolor=cmap_Beast(v/maxval_Node_Beast), label=label,
                        color = 'w', markerfacecolor = cmap_Beast(v/maxval_Node_Beast), markersize=10))
    
    ax[0].legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.03),
        fancybox=True, shadow=True, ncol=5, fontsize=9)
    ax[0].text(0.02, 0.93, "A", transform=ax[0].transAxes, 
                size=15, weight='bold')

    ## Nodes Year

    nx.draw_networkx_nodes(T, pos, node_size=15, nodelist= nodes, node_color=[cmap_Date(v/maxval_node_Date) for v in values_node_Date], 
        node_shape = 'o', ax=ax[1])
    
    #nx.draw_networkx_edges(T,pos,alpha = 0.6, edge_color=[cmap_edge(v/maxval_edge) for v in values_edge], edge_cmap=cmap_edge, ax=ax[1], width=0.6)

    edges = nx.draw_networkx_edges(T,pos,alpha = 0.6, edge_color=[cmap_edge(v/maxval_edge) for v in values_edge], edge_cmap=cmap_edge, ax=ax[1], width=0.6)
    sm = plt.cm.ScalarMappable(cmap=cmap_edge, norm=mpl.colors.Normalize(vmin = 0, vmax = maxval_edge))
    sm._A = []
    plt.colorbar(sm, ax = ax[[2,5]], location='right', shrink=0.6, label="Number of SNPs", alpha = 0.6)
    ax[2].set_visible(False)
    ax[5].set_visible(False)

    legend_elements = []
    for v in set(values_node_Date):
        if v == 0:
            label = "2015"
        elif v == 1:
            label = "2016"
        elif v == 2:
            label = "2017"
        elif v == 3:
            label = "2018"
        elif v == 4:
            label = "2019"
        elif v == 5:
            label = "2020"
        elif v == 6:
            label = "2021"
        legend_elements.append(Line2D([], [], marker='o', markeredgecolor=cmap_Date(v/maxval_node_Date), label=label,
                        color = 'w', markerfacecolor = cmap_Date(v/maxval_node_Date), markersize=10))
    
    ax[1].legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.03),
        fancybox=True, shadow=True, ncol=4, fontsize=9)
    ax[1].text(0.02, 0.93, "B", transform=ax[1].transAxes, 
                size=15, weight='bold')

    ## Nodes Source Type

    nx.draw_networkx_nodes(T, pos, node_size=15, nodelist= nodes, node_color=[cmap_Serology(v/maxval_Node_Serology) for v in values_Node_Serology], 
        node_shape = 'o', ax=ax[3])
    
    nx.draw_networkx_edges(T,pos,alpha = 0.6, edge_color=[cmap_edge(v/maxval_edge) for v in values_edge], edge_cmap=cmap_edge, ax=ax[3], width=0.6)
    legend_elements = []
    for v in set(values_Node_Serology):
        if v == 0:
            label = "Ogawa"
        elif v == 1:
            label = "Inaba"
        
        legend_elements.append(Line2D([], [], marker='o', markeredgecolor=cmap_Serology(v/maxval_Node_Serology), label=label,
                        color = 'w', markerfacecolor = cmap_Serology(v/maxval_Node_Serology), markersize=10))
    
    ax[3].legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.03),
        fancybox=True, shadow=True, ncol=3, fontsize=9)
    ax[3].text(0.02, 0.93, "C", transform=ax[3].transAxes, 
                size=15, weight='bold')

    ## Nodes Location

    nx.draw_networkx_nodes(T, pos, node_size=15, nodelist= nodes, node_color=[cmap_Location(v/maxval_node_Location) for v in values_node_Location], 
        node_shape = 'o', ax=ax[4])
    
    nx.draw_networkx_edges(T,pos,alpha = 0.6, edge_color=[cmap_edge(v/maxval_edge) for v in values_edge], edge_cmap=cmap_edge, ax=ax[4], width=0.6)
    legend_elements = []
    for v in set(values_node_Location):
        if v == 0:
            label = "Barisal"
        elif v == 1:
            label = "Chittagong"
        elif v == 2:
            label = "Dhaka"
        elif v == 3:
            label = "Khulna"
        elif v == 4:
            label = "Rajshahi"
        elif v == 5:
            label = "Sylhet"
        elif v == 6:
            label = "Narayanganj"

        legend_elements.append(Line2D([], [], marker='o', markeredgecolor=cmap_Location(v/maxval_node_Location), label=label,
                        color = 'w', markerfacecolor = cmap_Location(v/maxval_node_Location), markersize=10))
    
    ax[4].legend(handles = legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.03),
        fancybox=True, shadow=True, ncol=3, fontsize=9)
    ax[4].text(0.02, 0.93, "D", transform=ax[4].transAxes, 
                size=15, weight='bold')


    plt.savefig(results_folder+'/SNP_network_'+name_dataset+'.svg', bbox_inches='tight')
        
        
        
    
