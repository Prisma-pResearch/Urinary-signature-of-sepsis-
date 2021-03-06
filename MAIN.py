
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import time
import random
import matplotlib.pyplot as plt
import os

# In[2]:



os.chdir(the directory where your codes are saved)
import preprocessing as pp
import feature_selectors as fs
import voting as vt


# In[3]:


os.chdir(the directory where you want to save your results)


# In[4]:


df =  pd.read_csv("./LimmaTAC_1112probes_177celfiles.txt", sep = "\t", index_col = [0])
df.head()


# In[5]:


limma = pd.read_csv("./Limma_FDR0.01_logFC1.txt", sep = "\t", index_col = [0])
limma.head()


# In[6]:


x_train, y_train = pp.SplitX_y(df)


# In[7]:


ranked_list_boruta, selected_features_boruta = fs.Boruta_fs(x_train, y_train)
#ranked_list_boruta.to_csv("RankData_Boruta_1112probes_default.csv")
#selected_features_boruta.to_csv("SelectedFeatures_Boruta_default_1112probes.csv")


# In[9]:


ranked_list_randomforest, selected_features_randomforest = fs.RandomForest_fs(x_train, y_train)
#ranked_list_randomforest.to_csv("RankData_RandomForest_1112probes.csv")
#selected_features_randomforest.to_csv("SelectedFeatures_RandomForest_1112probes.csv")


# In[10]:


ranked_list_svc, ranked_list_topfeatures_svc, selected_features_svc = fs.SVC_fs(x_train, y_train)
#ranked_list_svc.to_csv("RankData_RFEsvm_1112probes.csv") 
#ranked_list_topfeatures_svc.to_csv("RankData_topfeatures_RFEsvm_1112probes.csv)
#selected_features_svc.to_csv("RankData_RFEsvm_topFeatures_1112probes.csv")


# In[11]:


selected_features_logreg = fs.LogisticRegression_Lasso_fs(x_train, y_train)
#selected_features_logreg.to_csv("RankData_Lasso_1112probes.csv") 


# In[ ]:


voted_probe_df, inter_probe_df, union_probe_df = vt.voting(selected_features_boruta, 
                                                           selected_features_randomforest, 
                                                           selected_features_svc, 
                                                           selected_features_logreg)

