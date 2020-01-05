#!python2
# >py -2 struct_diff.py -m 33 26 29 27 28 32 31 30 3 4 2 1 5 11 6 7 9 10 8 15 12 13 14 25 23 21 18 17 16 19 22 24 20 -ap --bins 50 --plot_basic --do_dendro --use_ddpv --verbosity 2

import numpy
import pickle
import pylab
from scipy.cluster import hierarchy
import meta

# part1
with open('msd_supp_fig2_1.pickle', 'r') as f:
    [mtr_list, max_mtr, bins, mtr, muscle_dics, n] = pickle.load(f)

pylab.figure()
pylab.hist(mtr_list, bins=numpy.linspace(0, max_mtr, bins), color='k')
pylab.xlim([0, max_mtr])
pylab.xlabel('Difference between muscles, au')
pylab.ylabel('Number of muscle pairs')
pylab.title('Mid len: {}'.format(n))

# heatmap
pylab.figure()
pylab.pcolor(mtr, cmap=pylab.get_cmap('gist_yarg'))
pylab.colorbar()
pylab.xticks([i+0.5 for i in xrange(n)], [i['sMuscle'] for i in muscle_dics], rotation='vertical')
pylab.yticks([i+0.5 for i in xrange(n)], [i['sMuscle'] for i in muscle_dics])
pylab.title('Mid len: {}'.format(n))
pylab.xlim(0, n)
pylab.ylim(0, n)


# part2
with open('msd_supp_fig2_2.pickle', 'r') as f:
    [muscle_dics, n, Z_linkage, fcluster, dendro_link, metamuscle_list] = pickle.load(f)
llf_addlen = 5
def llf(id):
    if id < n:
        return (muscle_dics[id]['sMuscle'] + ' {:d} {:2d}').format(
            muscle_dics[id]['dim'], fcluster[id])
    else:
        return '[%d %1.2f]' % (id, R[n-id,3])
pylab.figure()
R = hierarchy.dendrogram(Z_linkage, orientation='left', leaf_label_func=llf)
pylab.xlabel('{} distance between clusters, au'.format(dendro_link))
pylab.title('Dendrogram')

print('The order of muscles (ids) in the dendrogram:')
print([meta.search_dic(metamuscle_list, sMuscle=i[:-llf_addlen])['idMuscle']
       for i in R['ivl']])
print([i[:-llf_addlen] for i in R['ivl']])
print('And reversed:')
print([meta.search_dic(metamuscle_list, sMuscle=i[:-llf_addlen])['idMuscle']
       for i in reversed(R['ivl'])])
print([i[:-llf_addlen] for i in reversed(R['ivl'])])

pylab.show()
