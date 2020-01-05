#!python2
# >py -2 struct_diff.py -m 1 2 10 23 15 16 20 29 25 24 18 17 13 22 19 12 21 5 4 14 3 11 9 6 7 28 8 27 26 32 33 31 30 -ap --bins 50 --plot_basic --do_dendro --do_pca --verbosity 2
import numpy
import pickle
import pylab
from scipy.cluster import hierarchy
import meta

# part1
with open('msd_fig6_1.pickle', 'r') as f:
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
with open('msd_fig6_2.pickle', 'r') as f:
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


# part 3
with open('msd_fig6_3.pickle', 'r') as f:
    [pcas, spv] = pickle.load(f)
print 'PCA clustering:'
for i in xrange(pcas.n_components_):
    buf = sorted(zip(pcas.components_[i], spv),
                 key=lambda x: -abs(x[0]))
    buf = ' Main comps: ' + ', '.join(
        ['{:10s}: {:+6.4f}'.format(j[1], j[0])
         for j in buf if abs(j[0])>0.2]) + '.'
    print ('\tComponent #{:2d}. Variance explained {:10.8f}.'
           'Total ve {:10.8f}.'
           '{}'
           '\n\t\t[{}].').format(
           i+1, pcas.explained_variance_ratio_[i],
           sum(pcas.explained_variance_ratio_[:i+1]),
           buf,
           ' '.join('{:+6.4f}'.format(i)
                    for i in pcas.components_[i]))

# part 4
with open('msd_fig6_4.pickle', 'r') as f:
    muscle_dics = pickle.load(f)
pylab.figure()
pylab.plot([i['PCA_param_vec'][0] for i in muscle_dics],
           [i['PCA_param_vec'][1] for i in muscle_dics], '.k')
for muscle_dic in muscle_dics:
    pylab.annotate(muscle_dic['sMuscle'],
                   muscle_dic['PCA_param_vec'][:2])
pylab.xlabel('PCA #1')
pylab.ylabel('PCA #2')
pylab.title('PCA scatter of all muscles in 1 and 2')

pylab.show()
