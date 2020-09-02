
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from boruta import BorutaPy
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.model_selection import GridSearchCV
import time
import random
import matplotlib.pyplot as plt


# In[2]:


def Boruta_fs(x_train, y_train):
    """Perform feature selection using Boruta
    
    Arguments:
    x_train, y_train
    """
    estimator = RandomForestClassifier(n_jobs=-1, 
                                       random_state=0, 
                                       class_weight='balanced')
    selector = BorutaPy(estimator, 
                        n_estimators='auto', 
                        verbose=2, 
                        random_state=0)  #perc=100, max_iter=100, two_step=True
    selector.fit(x_train.values, y_train.values)
    
    feature_names = x_train.columns.values
    df_rank = pd.DataFrame({'Rank': selector.ranking_, 
                            'Features': feature_names})  #finding ranked list
    confirmed_indices = np.where(selector.ranking_ == 1)  #saving the confirmed features
    confirmed_names = x_train.columns.values[confirmed_indices]
    df_rank_confirmed = pd.DataFrame(confirmed_names)  #print confirmed_names
    df_rank_confirmed.index += 1
    
    return df_rank, df_rank_confirmed


# In[2]:


def RandomForest_fs(x_train, y_train):
    """Perform feature selection using Random Forest
    
    Arguments:
    x_train, y_train
    """
    estimator= RandomForestClassifier(n_jobs=1, 
                                      random_state=0, 
                                      class_weight='balanced') 
    params = {
        'n_estimators': [10, 100, 1000, 10000],  #100 for 2434 genes. 1000 for 1112 genes
        'min_samples_leaf': [1, 2, 3, 4]  #2 for 2434 probes. 2 for 1112 genes                             
    }
    
    CV_rfc = GridSearchCV(estimator, 
                          param_grid=params, 
                          scoring='roc_auc', 
                          cv=5, 
                          verbose=2, 
                          n_jobs=-1)
    CV_rfc.fit(x_train, y_train)
    forest = CV_rfc.best_estimator_
    importances = forest.feature_importances_
    indices = np.argsort(importances)[::-1]
    plt.plot(importances[indices])  #this plot is used to select how many features to select from Random Forest. 200 here. 
    plt.show()
    
    indices_top = indices[0:200] 
    rank_data = x_train.columns.values[indices]
    df_rank = pd.DataFrame(rank_data)
    df_rank.index += 1
    df_rank.columns = ['Features']
    selected_features = pd.DataFrame(x_train.columns.values[indices_top])
    
    return df_rank, selected_features


# In[3]:


def SVC_fs(x_train, y_train):
    """Perform feature selection using Support Vector Machine
    
    Arguments:
    x_train, y_train
    """
    estimator = SVC(kernel = "linear", 
                    probability= True)  #define estimator for feature selection technique
    selector = RFE(estimator)  #define feature selection technique
    params = {
        'step': [1],
        'n_features_to_select': [50, 100, 200, 300, 400, 500],  #200 for 1112 probes
        'estimator__C': [0.01, 0.1, 1.0, 10, 100]  #0.1 for 1112 probes                              
    }
    CV_rfc = GridSearchCV(selector, 
                          param_grid=params, 
                          scoring='roc_auc', 
                          cv=5, 
                          verbose=2, 
                          n_jobs=-1)
    CV_rfc.fit(x_train, y_train)
    
    x_trans_train =  x_train.iloc[:, CV_rfc.best_estimator_.support_]  #slice the data to include only selected features
    final_gene_list = x_trans_train.columns.values
    feature_names = x_train.columns.values 
    df_rank = pd.DataFrame({'Rank': CV_rfc.best_estimator_.ranking_, 
                            'Features': feature_names})
     
    top_indices = np.where(CV_rfc.best_estimator_.ranking_ == 1)  #getting index of all the top features
    top_features = x_train.columns.values[top_indices]  #getting names of the top features in sequence 
    top_ranking_coefs = CV_rfc.best_estimator_.estimator_.coef_  #SVM coefs for the top ranking features
    top_ranking_indices = np.argsort(abs(top_ranking_coefs[0]))[::-1]  #sorted such that the index of highest rank is shown. 
    top_features_ranked = top_features[top_ranking_indices]
    df_rank_top_features = pd.DataFrame(top_features_ranked)
    df_rank_top_features.index += 1
    selected_features = pd.DataFrame(final_gene_list)
    
    return df_rank, df_rank_top_features, selected_features


# In[4]:


def LogisticRegression_Lasso_fs(x_train, y_train):
    """Perform feature selection using Logistic Regression regularized with Lasso
    
    Arguments:
    x_train, y_train
    """
    estimator = LogisticRegression(penalty='l1', 
                                   class_weight='balanced', 
                                   random_state=0)
    params = {
        "C": [0.01, 0.1, 1.0, 10, 100]  #1 for 2,434 probes, 100 for 1112 probes                           
    }
    
    CV_rfc = GridSearchCV(estimator, 
                          param_grid=params, 
                          scoring='roc_auc', 
                          cv=5, 
                          verbose=2, 
                          n_jobs=-1)
    CV_rfc.fit(x_train, y_train)
    
    top_indices = np.nonzero(CV_rfc.best_estimator_.coef_)
    top_indices = np.nonzero(CV_rfc.best_estimator_.coef_)[1]
    top_features = x_train.columns.values[top_indices]  #getting names of the top features in sequence 
    top_ranking_coefs = CV_rfc.best_estimator_.coef_[0][top_indices]  #SVM coefs for the top ranking features
    top_ranking_indices = np.argsort(abs(top_ranking_coefs))[::-1]  #sorted such that the index of highest rank is shown. 
    top_features_ranked = top_features[top_ranking_indices]
    df_rank_top_features = pd.DataFrame(top_features_ranked)
    df_rank_top_features.index += 1
    
    return df_rank_top_features

