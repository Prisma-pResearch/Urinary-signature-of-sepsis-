
# coding: utf-8

# In[1]:


def SplitX_y(data):
    """Separate X from y"""
    x = data.drop(['Outcome'], axis=1)
    y = data['Outcome']
    
    return x, y


# In[2]:


def Find_GeneDirection(gene_list, limma):
    """Find upregulated and downregulated genes
    
    Arguments:
    gene_list -- Gene list of interest
    limma -- LIMMA analysis on whole genome
    """
    positive_list, negative_list = [], []
    for count, x in enumerate(gene_list):
        target = positive_list if(np.sign(limma['Fold Change'][x]) == 1) else negative_list
        target.append(x)
    
    return positive_list, negative_list

