
# coding: utf-8

# In[1]:


def SplitX_y(data):
    X = data.drop(['Outcome'], axis=1); y = data['Outcome']
    return X, y


# In[2]:


def Find_GeneDirection(gene_list, limma):
    positive_list, negative_list = [], []
    for count, x in enumerate(gene_list):
        target = positive_list if(np.sign(limma['Fold Change'][x]) == 1) else negative_list
        target.append(x)
    return positive_list, negative_list

