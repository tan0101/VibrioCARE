# -*- coding: utf-8 -*-
"""
Created on Fri Sep  24 20:35:00 2021

@author: Alexandre Maciel Guerra
"""

import numpy as np
import pandas as pd
import sys
import os
import pickle

from collections import Counter
from sklearn.model_selection import StratifiedKFold, GridSearchCV,cross_validate
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, ExtraTreesClassifier,GradientBoostingClassifier
from sklearn.metrics import confusion_matrix, accuracy_score, make_scorer, cohen_kappa_score
from imblearn.pipeline import Pipeline
from imblearn.over_sampling import SMOTE
from sklearn.feature_selection import VarianceThreshold
from scipy.stats import fisher_exact
from scipy.stats import chi2_contingency
from scipy.stats.contingency import crosstab

# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category= FutureWarning)
simplefilter(action='ignore', category= UserWarning)
simplefilter(action='ignore', category= DeprecationWarning)

def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]
def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]
def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 1]
scoring = {'tp': make_scorer(tp), 'tn': make_scorer(tn),
           'fp': make_scorer(fp), 'fn': make_scorer(fn),
           'auc': 'roc_auc',
           'acc': make_scorer(accuracy_score),
           'kappa': make_scorer(cohen_kappa_score)}

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
    type_data = "Clinical Classification Paper BD-1.2"
    
    # Load Data Core:
    data_core_df = pd.read_csv(folder+"/"+name_dataset+'_core_genome_snps_data.csv', header = [0], index_col=[0])
    features_name_core = np.array(data_core_df.columns)
    sample_name_core = np.array(data_core_df.index)
    data_core = np.array(data_core_df)
    
    # Load Data Accessory:
    data_acc_df = pd.read_csv(folder+"/"+name_dataset+'_accessory_genes_data.csv', header = [0], index_col=[0])
    features_name_acc = np.array(data_acc_df.columns)
    sample_name_acc = np.array(data_acc_df.index)
    data_acc = np.array(data_acc_df)
    
    # Load Data Intergenic:
    data_inter_df = pd.read_csv(folder+"/"+name_dataset+'_intergenic_snps_data.csv', header = [0], index_col=[0])
    features_name_inter = np.array(data_inter_df.columns)
    sample_name_inter = np.array(data_inter_df.index)
    data_inter = np.array(data_inter_df)
    
    print(np.array_equal(sample_name_core, sample_name_acc))
    print(np.array_equal(sample_name_inter, sample_name_acc))
    print(np.array_equal(sample_name_inter, sample_name_core))
    
    # Concatenate Data:
    data_comb_df = pd.concat([data_acc_df, data_core_df, data_inter_df], axis=1)
    features_name_comb = np.array(data_comb_df.columns)
    sample_name_comb = np.array(data_comb_df.index)

    # Load Clinical Data:
    clinical_data_df = pd.read_csv(folder+"/"+name_dataset+"_clinical_data.csv", header=[0])
    samples_clinical = np.array(clinical_data_df[clinical_data_df.columns[0]])
    
    metadata_df = pd.read_csv(folder+"/"+name_dataset+"_metadata3.csv", header=[0], index_col=[0])
    samples_meta = np.array(metadata_df[metadata_df.columns[0]])

    # Change order of samples in the AMR dataframes
    order_clinical = []
    order_meta = []
    for count, s_name in enumerate(sample_name_comb):
        idx_clinical = np.where(samples_clinical == s_name)[0]
        if len(idx_clinical) == 0:
            continue
        else:
            order_clinical.append(idx_clinical[0])
        
        idx_meta = np.where(samples_meta == s_name)[0]
        order_meta.append(idx_meta[0])
    
    # Clade ids
    metadata_df = metadata_df.iloc[order_meta,:].reset_index()
    metadata_df.drop(columns="index", axis=1, inplace=True)   

    # Clinical Data
    clinical_data_df = clinical_data_df.iloc[order_clinical,:].reset_index()
    clinical_data_df.drop(columns="index", axis=1, inplace=True)
    samples_clinical = np.array(clinical_data_df[clinical_data_df.columns[0]])

    # Check order    
    data_comb_df = data_comb_df.loc[samples_clinical,:]
    
    print(np.array_equal(data_comb_df.index, samples_clinical)) 

    clade_name = "BD-1.2"
    id_clade = np.where(metadata_df["Beast"] == clade_name)[0]

    print(metadata_df.shape)
    print(data_comb_df.shape)
    print(clinical_data_df.shape)
    
    data_comb_df = data_comb_df.loc[samples_clinical[id_clade],:]
    data_comb = np.array(data_comb_df)
    data_comb[data_comb>0]=1

    clinical_data_df = clinical_data_df.iloc[id_clade,:].reset_index()
    clinical_data_df.drop(columns="index", axis=1, inplace=True)

    metadata_df = metadata_df.loc[id_clade,:].reset_index()
    metadata_df.drop(columns="index", axis=1, inplace=True)
    print(np.array_equal(samples_clinical[id_clade], metadata_df[metadata_df.columns[0]]))     
    print(metadata_df.shape)
    print(data_comb_df.shape)
    print(clinical_data_df.shape)
    
    # Nested Cross Validation:
    inner_loop_cv = 3   
    outer_loop_cv = 5
    
    # Number of random trials:
    NUM_TRIALS = 30
    
    # Grid of Parameters:
    C_grid = {"clf__C": [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000]}
    est_grid = {"clf__n_estimators": [2, 4, 8, 16, 32, 64, 128, 256]}
    SVC_grid = {"clf__gamma": [0.0001, 0.0001, 0.001, 0.01, 0.1, 1], "clf__C": [0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000]}
    DT_grid = {"clf__max_depth": [10, 20, 30, 50, 100, 200, 300]}
    XGBoost_grid = {"clf__n_estimators": [2, 4, 8, 16, 32, 64, 128, 256], "clf__learning_rate": [0.0001, 0.001, 0.01, 0.1, 1]}
        
    # Classifiers:
    names = ["Logistic Regression", "Linear SVM", "RBF SVM",
        "Extra Trees", "Random Forest", "AdaBoost", "XGBoost"]

    classifiers = [
        LogisticRegression(),
        LinearSVC(loss='hinge'),
        SVC(),
        ExtraTreesClassifier(),
        RandomForestClassifier(),
        AdaBoostClassifier(),
        GradientBoostingClassifier()
        ]
    
    print(clinical_data_df.columns[5:11])
    for name_clinical in clinical_data_df.columns[6:11]:
        print("Antibiotic: {}".format(name_clinical))
        
        target_orig = np.array(clinical_data_df[name_clinical])
        
        # Check minimum number of samples:
        count_class = Counter(target_orig)
        print(len(target_orig))
        print(count_class)
        print(len(count_class))
        
        if len(count_class) < 2 or len(count_class)>2:
            continue

        count_class = Counter(target_orig)
        most_common = count_class.most_common(1)[0]

        if count_class[0] < 10 or count_class[1] < 10:
            continue 

        data_orig_df = data_comb_df
        data_orig = np.array(data_orig_df)

        features_anti = data_comb_df.columns

        # Remove low variance:
        print("Before removing low variance:{}".format(data_orig.shape))
        selector = VarianceThreshold(threshold=0)
        selector.fit_transform(data_orig)
        cols=selector.get_support(indices=True)
        data_orig = data_orig[:,cols]
        features_anti = features_anti[cols]
        n_features = len(features_anti)
        print("After removing low variance:{}".format(data_orig.shape))

        writer = pd.ExcelWriter(folder+"/"+results_folder+"/"+type_data+"/Confound"+"/Confound_statistics2_"+name_clinical+".xlsx", engine='xlsxwriter')
        df_confound = pd.DataFrame()

        for k, name_confound in enumerate(["Serology", "Year", "Location", "Age", "Sex"]):
            target_confound = np.array(metadata_df[name_confound])

            df_res = pd.DataFrame()
            labels_true = target_confound

            table = crosstab(labels_true, target_orig)
            stat, p, dof, expected = chi2_contingency(table[-1])

            df_confound.loc[k,"Factor"] = name_confound
            df_confound.loc[k,"chi2"] = p

            for count, feat in enumerate(features_anti):
                table = crosstab(labels_true, data_orig[:,count])
                stat, p, dof, expected = chi2_contingency(table[-1])

                df_res.loc[count,"Feature"] = feat
                df_res.loc[count,"chi2"] = p

            df_res.sort_values('chi2', ignore_index=True, inplace=True, ascending=True)
            df_res.to_excel(writer, sheet_name = name_confound, index = True)                        

            # Get the dimensions of the dataframe.
            (max_row, max_col) = df_res.shape

            # Get the xlsxwriter workbook and worksheet objects.
            workbook  = writer.book
            worksheet = writer.sheets[name_confound]
            # Apply a conditional format to the required cell range.

            cell_format = workbook.add_format()
            cell_format.set_font_color('red')

            worksheet.conditional_format(1, max_col, max_row, max_col,
                                        {'type':     'cell',
                                            'criteria': 'less than',
                                            'value':     0.01/df_res.shape[0],
                                            'format':    cell_format})
            

        df_confound.sort_values('chi2', ignore_index=True, inplace=True, ascending=True)
        df_confound.to_excel(writer, sheet_name = "Confound", index = True)

        # Get the xlsxwriter workbook and worksheet objects.
        workbook  = writer.book
        worksheet = writer.sheets["Confound"]
        # Apply a conditional format to the required cell range.

        cell_format = workbook.add_format()
        cell_format.set_font_color('red')

        worksheet.conditional_format(1, max_col, max_row, max_col,
                                    {'type':     'cell',
                                        'criteria': 'less than',
                                        'value':     0.01/df_confound.shape[0],
                                        'format':    cell_format})
            
        writer.close()

        important_factors = np.array(df_confound.loc[np.where(df_confound["chi2"]<0.01/len(df_confound))[0],"Factor"])
        
        features_remove = []
        for factor in important_factors:
            df_factor = pd.read_excel(folder+"/"+results_folder+"/"+type_data+"/Confound"+"/Confound_statistics_"+name_clinical+".xlsx", sheet_name=factor, header=[0], index_col=[0])
            for count, feat in enumerate(df_factor["Feature"]):
                if df_factor.loc[count,"chi2"] < 0.01/df_factor.shape[0]:
                    features_remove.append(feat)

        features_remove = np.unique(features_remove)

        idx_keep = []
        for count, feat in enumerate(features_anti):
            if feat not in features_remove:
                idx_keep.append(count)

        features_anti = features_anti[idx_keep]
        data_orig = data_orig[:,idx_keep]
        n_features = len(features_anti)
        print("After removing confound:{}".format(data_orig.shape))

        features_sel = np.zeros(len(features_anti))
        features_sel_array = np.zeros((len(features_anti),1000))

        update_progress(0)
        for rs in range(0,1000):
            sm = SMOTE(random_state=rs) 
            data, target = sm.fit_resample(data_orig, target_orig)
            data[data>0.5]=1
            data[data<0.5]=0
            
            pvalue_chi2 = np.zeros(n_features,dtype=float)
            for n_feat in range(n_features): 
                id_0 = np.where(target == 0)[0]
                id_1 = np.where(target == 1)[0]
                array = np.zeros((2,2))
                array[0,0] = len(np.where(data[id_0,n_feat] == 1)[0])
                array[1,0] = len(np.where(data[id_0,n_feat] == 0)[0])
                array[0,1] = len(np.where(data[id_1,n_feat] == 1)[0])
                array[1,1] = len(np.where(data[id_1,n_feat] == 0)[0])

                res = fisher_exact(array)
                pvalue_chi2[n_feat] = res[1]

            id_chi2 = np.where(pvalue_chi2 < 0.1)[0]
        
            if len(id_chi2) == 0:
                continue

            features_sel[id_chi2]+=1
            features_sel_array[id_chi2,rs]=1

            update_progress((rs+1)/1000)
    
        id_sel = np.where(features_sel>=750)[0]

        rs_array = []
        update_progress(0)
        for rs in range(0,1000):
            feat_sel = np.where(features_sel_array[:,rs] == 1)[0]
            intersect = np.intersect1d(id_sel,feat_sel)

            if len(intersect) == len(id_sel):
                rs_array.append(rs)

            update_progress((rs+1)/1000)

        print(rs_array)
        print(len(rs_array))

        sm = SMOTE(random_state=rs_array[0])
        data, target = sm.fit_resample(data_orig, target_orig)
        data[data>0.5]=1
        data[data<0.5]=0

        idx = np.where(features_sel_array[:,rs_array[0]] == 1)[0]
        pvalue_sel = np.zeros(len(idx))
        for count_idx, index in enumerate(idx):
            id_0 = np.where(target == 0)[0]
            id_1 = np.where(target == 1)[0]
            array = np.zeros((2,2))
            array[0,0] = len(np.where(data[id_0,index] == 1)[0])
            array[1,0] = len(np.where(data[id_0,index] == 0)[0])
            array[0,1] = len(np.where(data[id_1,index] == 1)[0])
            array[1,1] = len(np.where(data[id_1,index] == 0)[0])

            res = fisher_exact(array)
            pvalue_sel[count_idx] = res[1]

        data = data[:,idx]

        directory = folder+"/"+results_folder+"/"+type_data
        
        if not os.path.exists(directory):
            os.makedirs(directory)

        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(1,1,1)
        ax.hist(features_sel, bins=20, range=(0,1000))
        ax.set_ylabel("Number of Features")
        ax.set_xlabel("Bins")
        ax.set_title(name_clinical)
        fig.tight_layout()
        plt.savefig(directory+"/features_"+name_dataset+"_"+name_clinical+'.svg', dpi=300)
        
        print("Features selected: {}".format(len(idx)))

        #data = data[:,id_sel]
 
        features_anti = features_anti[idx]
            
        with open(directory+"/data_"+name_dataset+"_"+name_clinical+'.pickle', 'wb') as f:
            pickle.dump(data, f)
        
        features_genes_df = pd.DataFrame(features_sel[idx], columns = ["selection times - rs "+str(rs_array[0])])
        features_genes_df["index"] = features_anti
        
        names_f = ["index", "selection times - rs "+str(rs_array[0])]
        features_genes_df = features_genes_df[names_f]
        features_genes_df["pvalue"] = pvalue_sel
        #features_genes_df = pd.DataFrame(features_anti, columns = ["features"])
        features_genes_df.to_csv(directory+"/features_"+name_dataset+"_"+name_clinical+".csv")

        #df_features = pd.DataFrame(concat_array, columns = ["coef"], index=features_anti)
        #df_features.to_csv(directory+"/features_"+name_dataset+"_"+name_clinical+".csv")

        # Initialize Variables:
        scores_auc = np.zeros([NUM_TRIALS,len(classifiers)])
        scores_acc = np.zeros([NUM_TRIALS,len(classifiers)])
        scores_sens = np.zeros([NUM_TRIALS,len(classifiers)])
        scores_spec = np.zeros([NUM_TRIALS,len(classifiers)])
        scores_kappa = np.zeros([NUM_TRIALS,len(classifiers)])
        scores_prec = np.zeros([NUM_TRIALS,len(classifiers)])
        
        # Loop for each trial
        update_progress(0)
        for i in range(NUM_TRIALS):
            #print("Trial = {}".format(i))
        
            inner_cv = StratifiedKFold(n_splits=inner_loop_cv, shuffle=True, random_state=i)
            outer_cv = StratifiedKFold(n_splits=outer_loop_cv, shuffle=True, random_state=i)
        
            k = 0
        
            for name, clf in zip(names, classifiers):
                #print("Classifier = {}".format(name))
                
                #model = Pipeline([('sampling',SMOTEENN(sampling_strategy = "all", random_state=i)),
                #                ('clf', clf)])

                #model = Pipeline([('sampling',SMOTE(random_state=i)),
                #                ('clf', clf)])
                
                model = Pipeline([('clf', clf)])

                if name == "QDA" or name == "LDA" or name == "Naive Bayes":
                    classif = model
                else:
                    if name == "RBF SVM":
                        grid = SVC_grid              
                    elif name == "Random Forest" or name == "AdaBoost" or name == "Extra Trees":
                        grid = est_grid
                    elif name == "Linear SVM":
                        grid = C_grid
                    elif name == "Decision Tree":
                        grid = DT_grid
                    elif name == "XGBoost":
                        grid = XGBoost_grid
                    else:
                        grid = C_grid

                    # Inner Search
                    classif = GridSearchCV(estimator=model, param_grid=grid, cv=inner_cv)
                    classif.fit(data, target)
            
                # Outer Search
                cv_results = cross_validate(classif, data, target, scoring=scoring, cv=outer_cv, return_estimator=True)

                tp = cv_results['test_tp']
                tn = cv_results['test_tn']
                fp = cv_results['test_fp']
                fn = cv_results['test_fn']
                
                sens = np.zeros(outer_loop_cv)
                spec = np.zeros(outer_loop_cv)
                prec = np.zeros(outer_loop_cv)
                
                for j in range(outer_loop_cv):
                    TP = tp[j]
                    TN = tn[j]
                    FP = fp[j]
                    FN = fn[j]
                    
                    # Sensitivity, hit rate, recall, or true positive rate
                    sens[j] = TP/(TP+FN)
                    
                    # Fall out or false positive rate
                    FPR = FP/(FP+TN)
                    spec[j] = 1 - FPR
                    if TP + FP > 0:
                        prec[j] = TP / (TP + FP)
    
                scores_sens[i,k] = sens.mean()
                scores_spec[i,k] = spec.mean()
                scores_prec[i,k] = prec.mean()
                scores_auc[i,k] = cv_results['test_auc'].mean()
                scores_acc[i,k] = cv_results['test_acc'].mean()
                scores_kappa[i,k] = cv_results['test_kappa'].mean()
                
                k = k + 1
                
            update_progress((i+1)/NUM_TRIALS)

        results = np.zeros((12,len(classifiers)))
        scores = [scores_auc, scores_acc, scores_sens, scores_spec, scores_kappa, scores_prec]
        for counter_scr, scr in enumerate(scores):
            results[2*counter_scr,:] = np.mean(scr,axis=0)
            results[2*counter_scr + 1,:] = np.std(scr,axis=0)
            
        names_scr = ["AUC_Mean", "AUC_Std", "Acc_Mean", "Acc_Std", 
            "Sens_Mean", "Sens_Std", "Spec_Mean", "Spec_Std", 
            "Kappa_Mean", "Kappa_Std", "Prec_Mean", "Prec_Std"]

        results_df=pd.DataFrame(results, columns=names, index=names_scr)
        results_df.to_csv(directory+"/SMOTE_results_"+name_dataset+"_"+name_clinical+".csv")

        df_auc = pd.DataFrame(scores_auc, columns=names)
        df_auc.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_auc.csv")
        
        df_acc = pd.DataFrame(scores_acc, columns=names)
        df_acc.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_acc.csv")
        
        df_sens = pd.DataFrame(scores_sens, columns=names)
        df_sens.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_sens.csv")
        
        df_spec = pd.DataFrame(scores_spec, columns=names)
        df_spec.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_spec.csv")
        
        df_kappa = pd.DataFrame(scores_kappa, columns=names)
        df_kappa.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_kappa.csv")
        
        df_prec = pd.DataFrame(scores_prec, columns=names)
        df_prec.to_csv(directory+"/SMOTE_"+name_dataset+"_"+name_clinical+"_prec.csv")
        
        
