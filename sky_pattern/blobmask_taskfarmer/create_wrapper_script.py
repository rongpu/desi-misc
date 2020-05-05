n_task = 15
band = 'g'
for index in range(n_task):
    print('wrapper.sh {} {} {}'.format(band, n_task, index))