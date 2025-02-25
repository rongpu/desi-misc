n_task = 15
band = 'g'
for index in range(n_task):
    print('wrapper.sh {} {} {}'.format(band, n_task, index))


n_task = 10
for index in range(n_task):
    print('wrapper.sh {} {}'.format(n_task, index))


n_task = 50
for index in range(n_task):
    print('wrapper.sh {} {}'.format(n_task, index))
