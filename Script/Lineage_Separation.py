import numpy as np
import pandas as pd
import sys

from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact

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
    name_dataset = "Vibrio" 
    folder = "Data"
    results_folder = "Results IDESHI and IEDCR"
    
    # Load Data Core:
    data_core_df = pd.read_csv(folder+"/"+name_dataset+'_core_genome_snps_data.csv', header = [0], index_col=[0])
    features_name_core = np.array(data_core_df.columns)
    sample_name_core = np.array(data_core_df.index)
    data_core = np.array(data_core_df)
    data_core_df[data_core_df>0]=1
    
    # Load Data Accessory:
    data_acc_df = pd.read_csv(folder+"/"+name_dataset+'_accessory_genes_data.csv', header = [0], index_col=[0])
    features_name_acc = np.array(data_acc_df.columns)
    sample_name_acc = np.array(data_acc_df.index)
    data_acc = np.array(data_acc_df)
    data_acc_df[data_acc_df>0]=1
    
    # Load Data Intergenic:
    data_inter_df = pd.read_csv(folder+"/"+name_dataset+'_intergenic_snps_data.csv', header = [0], index_col=[0])
    features_name_inter = np.array(data_inter_df.columns)
    sample_name_inter = np.array(data_inter_df.index)
    data_inter = np.array(data_inter_df)
    data_inter_df[data_inter_df>0]=1
    
    print(np.array_equal(sample_name_core, sample_name_acc))
    print(np.array_equal(sample_name_inter, sample_name_acc))
    print(np.array_equal(sample_name_inter, sample_name_core))
    
    # Concatenate Data:
    data_comb_df = pd.concat([data_acc_df, data_core_df, data_inter_df], axis=1)
    features_name_comb = np.array(data_comb_df.columns)
    sample_name_comb = np.array(data_comb_df.index)
    data_comb = np.array(data_comb_df)
    data_comb[data_comb>0]=1

    # Load Antibiotic Data:
    metadata_df = pd.read_csv(folder+"/"+name_dataset+'_metadata.csv', header = [0], index_col=[0])
    samples_metadata = np.array(metadata_df[metadata_df.columns[0]])
    
    print(np.array_equal(sample_name_comb, samples_metadata))
        
    # Change order of samples in the AMR dataframes
    order_id = []
    for count, s_name in enumerate(sample_name_comb):
        idx = np.where(samples_metadata == s_name)[0]
        order_id.append(idx[0])
        
    metadata_df = metadata_df.iloc[order_id, :].reset_index()
    metadata_df.drop(columns="index",axis=1,inplace=True)
    samples_metadata = np.array(metadata_df[metadata_df.columns[0]])
    
    print(np.array_equal(sample_name_comb, samples_metadata))

    data_acc_mean = np.sum(data_acc,axis=1)
    data_core_mean = np.sum(data_core,axis=1)
    data_inter_mean = np.sum(data_inter,axis=1)
    data_df = metadata_df.loc[:,metadata_df.columns[[1,2]]]
    data_df['Year'] = data_df['Year'].astype('int')
    data_df["Accessory Genes"] = data_acc_mean
    data_df["Core Genome SNPs"] = data_core_mean
    data_df["Intergenic Region SNPs"] = data_inter_mean
    data_df["Lineage"] = metadata_df["Beast"]
    data_df["Serotype"] = metadata_df["Serology"]
    
    #Statistics Lineages
    print("Lineage")

    id_bd1 = np.where(data_df["Lineage"] == "BD-1.2")[0]
    id_bd2 = np.where(data_df["Lineage"] == "BD-2")[0]
    
    _, pvalue = mannwhitneyu(data_df.loc[id_bd1,"Accessory Genes"], data_df.loc[id_bd2,"Accessory Genes"])
    print("Pvalue Accessory Genes = {}".format(pvalue))
    print("Mean value Accessory Genes BD-1.2 = {}".format(np.round(np.mean(data_df.loc[id_bd1,"Accessory Genes"]))))
    print("SD value Accessory Genes BD-1.2 = {}".format(np.round(np.std(data_df.loc[id_bd1,"Accessory Genes"]))))
    print("Mean value Accessory Genes BD-2 = {}".format(np.round(np.mean(data_df.loc[id_bd2,"Accessory Genes"]))))
    print("SD value Accessory Genes BD-2 = {}".format(np.round(np.std(data_df.loc[id_bd2,"Accessory Genes"]))))

    _, pvalue = mannwhitneyu(data_df.loc[id_bd1,"Core Genome SNPs"], data_df.loc[id_bd2,"Core Genome SNPs"])
    print("Pvalue Core Genome SNPs = {}".format(pvalue))
    print("Mean value Core Genome SNPs BD-1.2 = {}".format(np.round(np.mean(data_df.loc[id_bd1,"Core Genome SNPs"]))))
    print("SD value Core Genome SNPs BD-1.2 = {}".format(np.round(np.std(data_df.loc[id_bd1,"Core Genome SNPs"]))))
    print("Mean value Core Genome SNPs BD-2 = {}".format(np.round(np.mean(data_df.loc[id_bd2,"Core Genome SNPs"]))))
    print("SD value Core Genome SNPs BD-2 = {}".format(np.round(np.std(data_df.loc[id_bd2,"Core Genome SNPs"]))))

    _, pvalue = mannwhitneyu(data_df.loc[id_bd1,"Intergenic Region SNPs"], data_df.loc[id_bd2,"Intergenic Region SNPs"])
    print("Pvalue Intergenic Region SNPs = {}".format(pvalue))
    print("Mean value Intergenic Region SNPs BD-1.2 = {}".format(np.round(np.mean(data_df.loc[id_bd1,"Intergenic Region SNPs"]))))
    print("SD value Intergenic Region SNPs BD-1.2 = {}".format(np.round(np.std(data_df.loc[id_bd1,"Intergenic Region SNPs"]))))
    print("Mean value Intergenic Region SNPs BD-2 = {}".format(np.round(np.mean(data_df.loc[id_bd2,"Intergenic Region SNPs"]))))
    print("SD value Intergenic Region SNPs BD-2 = {}".format(np.round(np.std(data_df.loc[id_bd2,"Intergenic Region SNPs"]))))

    #Statistics Year

    years_unique = np.unique(data_df['Year'])
    
    df_acc_pvalue = pd.DataFrame()
    for row in years_unique:
        for col in years_unique:
            id_row = np.where(data_df["Year"] == row)[0]
            id_col = np.where(data_df["Year"] == col)[0]

            _, pvalue = mannwhitneyu(data_df.loc[id_row,"Accessory Genes"], data_df.loc[id_col,"Accessory Genes"])
            if pvalue < 0.05:
                df_acc_pvalue.loc[row,col] = True
            else:
                df_acc_pvalue.loc[row,col] = False

    df_core_pvalue = pd.DataFrame()
    for row in years_unique:
        for col in years_unique:
            id_row = np.where(data_df["Year"] == row)[0]
            id_col = np.where(data_df["Year"] == col)[0]

            _, pvalue = mannwhitneyu(data_df.loc[id_row,"Core Genome SNPs"], data_df.loc[id_col,"Core Genome SNPs"])
            if pvalue < 0.05:
                df_core_pvalue.loc[row,col] = True
            else:
                df_core_pvalue.loc[row,col] = False


    df_inter_pvalue = pd.DataFrame()
    for row in years_unique:
        for col in years_unique:
            id_row = np.where(data_df["Year"] == row)[0]
            id_col = np.where(data_df["Year"] == col)[0]

            _, pvalue = mannwhitneyu(data_df.loc[id_row,"Intergenic Region SNPs"], data_df.loc[id_col,"Intergenic Region SNPs"])
            if pvalue < 0.05:
                df_inter_pvalue.loc[row,col] = True
            else:
                df_inter_pvalue.loc[row,col] = False

    print("Accessory")
    print(df_acc_pvalue.head(7))

    print("Core Genome")
    print(df_core_pvalue.head(7))

    print("Intergenic")
    print(df_inter_pvalue.head(7))

    writer = pd.ExcelWriter("Features_lineages.xlsx", engine='xlsxwriter')

    for feat_type in ["Accessory Genes", "Core Genome SNPs", "Intergenic Region SNPs"]:
        if feat_type == "Accessory Genes":
            df = data_acc_df
        elif feat_type == "Core Genome SNPs":
            df = data_core_df
        elif feat_type == "Intergenic Region SNPs":
            df = data_inter_df

        df_res = pd.DataFrame()
        for count, feat in enumerate(df.columns):
            array = np.zeros((2,2))
            array[0,0] = len(np.where(df.iloc[id_bd1,count] == 1)[0])
            array[1,0] = len(np.where(df.iloc[id_bd1,count] == 0)[0])
            array[0,1] = len(np.where(df.iloc[id_bd2,count] == 1)[0])
            array[1,1] = len(np.where(df.iloc[id_bd2,count] == 0)[0])

            df_res.loc[count,"Feature"] = feat
            df_res.loc[count,"Presence BD-1.2"] = 100*array[0,0]/len(id_bd1)
            df_res.loc[count,"Absence BD-1.2"] = 100*array[1,0]/len(id_bd1)
            df_res.loc[count,"Presence BD-2"] = 100*array[0,1]/len(id_bd2)
            df_res.loc[count,"Absence BD-2"] = 100*array[1,1]/len(id_bd2)

            if 100*array[0,0]/len(id_bd1) + 100*array[1,0]/len(id_bd1) < 100 or 100*array[0,1]/len(id_bd2) + 100*array[1,1]/len(id_bd2) < 100:
                print(array)
                print(df_res.loc[count,:])
                print(len(id_bd1))
                print(len(id_bd2))
                print(df.iloc[id_bd1,count])
                print(df.iloc[id_bd2,count])
                input("cont")
            
            res = fisher_exact(array)

            df_res.loc[count,"Fisher Exact"] = res[1]

        df_res.sort_values('Fisher Exact', ignore_index=True, inplace=True)
        df_res.to_excel(writer, sheet_name = feat_type, index = True)

        # Get the dimensions of the dataframe.
        (max_row, max_col) = df_res.shape

        # Get the xlsxwriter workbook and worksheet objects.
        workbook  = writer.book
        worksheet = writer.sheets[feat_type]
        # Apply a conditional format to the required cell range.

        cell_format = workbook.add_format()
        cell_format.set_font_color('red')

        worksheet.conditional_format(1, max_col, max_row, max_col,
                                    {'type':     'cell',
                                        'criteria': 'less than',
                                        'value':     0.05/len(df.columns),
                                        'format':    cell_format})
        
    writer.close()
