import sys, os, glob, time, warnings, gc
from multiprocessing import Pool

def f(n):
    print(n)
    for i in range(n):
        for j in range(n):
            1.*1.

if __name__ == '__main__':

    n_processess = 32
    arg = [20000] * n_processess

    print('Start')

    start = time.time()

    with Pool(processes=n_processess) as pool:
        res = pool.map(f, arg)

    end = time.time()
    print('Took {:.1f} seconds'.format(end - start))

    # start = time.time()
    # f(20000)
    # end = time.time()
    # print('Took {:.1f} seconds'.format(end - start))
