#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Pandas понадобится нам для чтения денных
import pandas as pd
import numpy as np
import scipy.stats


# In[87]:


from statsmodels.stats.multitest import multipletests


# In[3]:


data_path = "./homework_lecture_5_data.csv"


# In[4]:



# In[5]:




#
# In[16]:


def demonstrate_clt(expressions, sample_size = 100, n_samples = 10000):
    sum_ind_random_sample = []
    for _ in range(n_samples):
        random_sample = np.random.choice(expressions, sample_size)
        mean_of_sample = random_sample.mean()
        sum_ind_random_sample.append(mean_of_sample)
    sum_ind_random_sample = np.array(sum_ind_random_sample)
    return sum_ind_random_sample



def manual_calculate_CI(population, gene, sample, alpha = 0.05): #alpha - допустимый уровень ошибки 1 рода
    tail_alpha = alpha/2
    quantile = scipy.stats.norm.ppf(tail_alpha)
    sigma = population[gene].std()
    sem = sigma / np.sqrt(len(sample))
    left_bound, right_bound = sample.mean() + quantile*sem, sample.mean() - quantile*sem
    return left_bound, right_bound





def check_intervals_intersect(first_ci, second_ci):
    L1, R1 = first_ci
    L2, R2 = second_ci
    if L1 >= R2:
        are_intersect = False
    elif L2 >= R1:
        are_intersect = False
    elif L1 < R2 or L2 > R1:
        are_intersect = True
    return are_intersect # True or False


# In[81]:





# In[92]:


def check_dge_with_ci(one_cell_ds, two_cell_ds, multiple_corrections: int , method = "bonferroni"):
    ci_test_results = []
    alpha =  1 - (1 - 0.05)**len(one_cell_ds.columns[:-1])  # возводим в степень  на кол-во гипотез == генов
    if multiple_corrections:    # если есть multiple_corrections, корректируем alpha для входа в вычисление доверительного интервала
      _, _, _, alpha  = multipletests(pvals = [1]*len(one_cell_ds.columns[:-1]) , alpha=alpha, method=method, is_sorted=False, returnsorted=False)
    for gene in one_cell_ds.columns[:-1]:
        one_cell_sample = np.random.choice(one_cell_ds[gene], 10000)    # сгенерируем случайную выборку для каждого гена
        two_cell_sample = np.random.choice(two_cell_ds[gene], 10000)
        first_ci = manual_calculate_CI(one_cell_ds, gene, one_cell_sample, alpha)
        second_ci = manual_calculate_CI(two_cell_ds, gene, two_cell_sample, alpha)
        res = check_intervals_intersect(first_ci, second_ci)
        ci_test_results.append(res)
    return ci_test_results


# In[93]:





# In[102]:


from statsmodels.stats.weightstats import ztest


# In[103]:


def check_dge_with_ztest(first_table, second_table, multiple_corrections, method = 'bonferroni'):
    z_test_results = []
    p_values = []
    alpha = 1 - (1 - 0.05)**len(first_table.columns[:-1])   #критическое значение ошибки 1 рода для тестирования всех генов
    for gene in first_table.columns[:-1]:
        _, p_value = ztest(first_table[gene], second_table[gene])
        p_values.append(p_value)
    if not multiple_corrections:    # если нет multiple_corrections, сравниваем с не поправленным alpha
        for p_value in p_values:
            if p_value <= alpha:
                z_test_results.append(True)
            else:
                z_test_results.append(False)
    else:
        z_test_results, p_values, _, _ = multipletests(pvals = p_values, alpha=alpha, method=method, is_sorted=False, returnsorted=False)
    return z_test_results, p_values


# In[104]:


def calc_mean_diff(table1, table2):
    genes = table1.columns[:-1]
    mean_diff = []
    for gene in genes:
        mean_b_cell = table1.loc[:, gene].mean()
        mean_nk_cell = table2.loc[:, gene].mean()
        diff = mean_b_cell - mean_nk_cell
        mean_diff.append(diff)
    return mean_diff
    


# In[114]:


def dfg_analyze(data_path, multiple_corrections: int):
    '''

    :param data_path: путь до данных экспрессии двух клеточных типов
    :param multiple_corrections: 1 если есть множественное тестирование, 0 - если нет.
    :return:
    '''
    data = pd.read_csv(data_path, index_col = 0)
    b_cell_data = data.query("Cell_type == 'B_cell'")
    nk_cell_data = data.query("Cell_type == 'NK_cell'")
    ci_test_results = check_dge_with_ci(b_cell_data, nk_cell_data, multiple_corrections=multiple_corrections)
    z_test_results, z_test_p_values = check_dge_with_ztest(b_cell_data, nk_cell_data, multiple_corrections = multiple_corrections)
    mean_diff = calc_mean_diff(b_cell_data, nk_cell_data)
    results = {
    "ci_test_results": ci_test_results,
    "z_test_results": z_test_results,
    "z_test_p_values": z_test_p_values,
    "mean_diff": mean_diff}
    results = pd.DataFrame(results)
    results.to_csv('dfg_analyze.csv', index=False)


# In[115]:




# In[116]:





# In[117]:





# In[63]:





# In[ ]:





#%%

#%%
