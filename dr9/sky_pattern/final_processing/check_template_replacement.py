import numpy as np

runs_to_be_replaced = np.array([426, 431, 518, 553, 592, 703, 706, 736])
runs_replaced_into =  np.array([425, 430, 519, 552, 591, 704, 707, 737])

if np.sum(np.abs(runs_to_be_replaced-runs_replaced_into)>1)!=0:
    print("ERROR!!!!")
else:
    print('Fine')
