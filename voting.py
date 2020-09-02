
# coding: utf-8

# In[1]:


def voting():
    """Perform majority voting between multiple feature lists"""
    list_boruta = list(boruta.Features)
    list_RFEsvm = list(RFE_svm.Features)
    list_RF = list(RF.Features)
    list_Lasso = list(Lasso_LogReg.Features)
    
    voted_probes = [k for k, v in cnt2.iteritems() if v > 1]
    inter_probes = [k for k, v in cnt2.iteritems() if v == 3]
    union_probes = [k for k, v in cnt2.iteritems() if v > 0]
    
    df_voted_probes = pd.DataFrame(voted_probes)
    df_inter_probes = pd.DataFrame(inter_probes)
    df_union_probes = pd.DataFrame(union_probes)
    
    return df_voted_probes, df_inter_probes, df_union_probes

