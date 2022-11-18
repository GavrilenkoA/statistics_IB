import sys
#%%
from HW3 import dfg_analyze
#%%
import argparse
#%%
#sys.path.append("/home/alexg/IB/stat/multiple_test")

#%%
data_path = input('Enter path to your data: ')
multiple_test = int(input('Enter 1 or 0: '))
#%%
dfg_analyze(data_path, multiple_test)