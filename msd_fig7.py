#!python2
# >py -2 struct_diff.py -m all -ap --bins 50 --do_bsbf --do_expbs --verbosity 2
import pickle
import numpy
import pylab

density_weights = lambda x: (numpy.ones_like(x)/float(len(x)))

# part 1
with open('msd_fig7_1.pickle', 'r') as f:
    [bsbf_incluster_dist, bsbf_outcluster_dist, bsbf_in_normt, bsbf_in_normt_p,
     bsbf_ou_normt, bsbf_ou_normt_p, bsbf_u, bsbf_u_nxny, bsbf_u_p,
     bsbf_w, bsbf_w_p, bsbf_bs, bsbf_bs_neg_per, bsbf_bs_pmean, bsbf_bs_t, bsbf_bs_t_p,
     bsbf_bs_z, bsbf_bs_z_p, bsbf_normt, bsbf_normt_p,
     bsbf_bs_w_p, conf, bsbf_bs_w_interv, bsbf_bs_w_norm_interv, bins] = pickle.load(f)
print ('Structural difference between muscles that cross '
       'different dofs')
print ('\tNumber of distances in incluster: {}.').format(
    len(bsbf_incluster_dist))
print ('\tNumber of distances in outcluster: {}.').format(
    len(bsbf_outcluster_dist))
print ('\tIn cluster distance mean {:6.4f} median {:6.4f}'
       ' std {:6.4f}').format(
    numpy.mean(bsbf_incluster_dist),
    numpy.median(bsbf_incluster_dist),
    numpy.std(bsbf_incluster_dist))
print ('\t\tNormality test {}, p-value {}').format(
    bsbf_in_normt, bsbf_in_normt_p)
print ('\tOut cluster distance mean {:6.4f} median {:6.4f}'
       ' std {:6.4f}').format(
    numpy.mean(bsbf_outcluster_dist),
    numpy.median(bsbf_outcluster_dist),
    numpy.std(bsbf_outcluster_dist))
print ('\t\tNormality test {}, p-value {}').format(
    bsbf_ou_normt, bsbf_ou_normt_p)
print ('\tMann-Whitney U {}, nx {}, ny {}, nx*ny {},'
       ' p-value {}.').format(
    bsbf_u, len(bsbf_incluster_dist), len(bsbf_outcluster_dist),
    bsbf_u_nxny, bsbf_u_p)
print ('\tWilcoxon signed-rank test T {}, p-value {}.').format(
    bsbf_w, bsbf_w_p)
print ('\tDifferences mean {} median {} std {}.'
       ' Porion negative {}. Number of elements {}.').format(
    numpy.mean(bsbf_bs), numpy.median(bsbf_bs), numpy.std(bsbf_bs),
    bsbf_bs_neg_per, len(bsbf_bs))
print ('\tT-test for {} mean t {} p-value {}.').format(
    bsbf_bs_pmean, bsbf_bs_t, bsbf_bs_t_p)
print ('\t\tZ-score {}, p {}.').format(
    bsbf_bs_z, bsbf_bs_z_p)
print ('\tDAgostino normality test {} p-value {}'
       ' (low=not normal)').format(
    bsbf_normt, bsbf_normt_p)
print ('\tSign 1-sample test p-value {},'
       ' confidence {} interval {} (assuming normal: {})').format(
    bsbf_bs_w_p, conf, bsbf_bs_w_interv, bsbf_bs_w_norm_interv)

pylab.figure()
bsbf_ax_left = pylab.subplot(2, 2, 1)
pylab.boxplot([bsbf_incluster_dist, bsbf_outcluster_dist],
              vert=False,
              labels=['incluster', 'outcluster'], whis=[0, 100])

bsbf_ax_bottom = pylab.subplot(2, 2, 3, sharex=bsbf_ax_left)
bsbf_bins = numpy.linspace(
    min([min(bsbf_incluster_dist), min(bsbf_outcluster_dist)]),
    max([max(bsbf_incluster_dist), max(bsbf_outcluster_dist)]),
    bins+1)
pylab.hist(
    bsbf_incluster_dist,
    weights=density_weights(bsbf_incluster_dist),
    bins=bsbf_bins, label='incluster', alpha=0.5)
pylab.hist(
    bsbf_outcluster_dist,
    weights=density_weights(bsbf_outcluster_dist),
    bins=bsbf_bins, label='outcluster', alpha=0.5)
pylab.xlim([0, max([max(bsbf_incluster_dist), max(bsbf_outcluster_dist)])])
pylab.xlabel('Distance between vectors')
pylab.ylabel('Probability')
pylab.legend()

bsbf_ax_right = pylab.subplot(2, 2, 2)
pylab.boxplot(bsbf_bs, vert=False, whis=[0, 100])
pylab.axvline(x=0, color='m')
pylab.axvline(x=bsbf_bs_w_interv[1], color='r', linestyle='--')

pylab.subplot(2, 2, 4, sharex=bsbf_ax_right, sharey=bsbf_ax_bottom)
pylab.hist(bsbf_bs, weights=density_weights(bsbf_bs), bins=bins)
pylab.axvline(x=0, color='m')
pylab.axvline(x=bsbf_bs_w_interv[1], color='r', linestyle='--')
pylab.xlabel('Distance between vectors in and out')
pylab.ylabel('Probability')

pylab.suptitle('Crossing same DOF clusters')


# part 2
with open('msd_fig7_2.pickle', 'r') as f:
    [expbs_cl_labels, expbs_numcls, expbs_fbss, expbs_fbss_in,
     expbs_in_tnorm_ps, expbs_fbss_out, expbs_ou_tnorm_ps, expbs_2smpl_zs,
     expbs_2smpl_z_ps, expbs_2smpl_ts, expbs_2smpl_t_ps, expbs_us,
     expbs_u_nxnys, expbs_u_ps, expbs_ws, expbs_w_ps, expbs_fbss,
     expbs_tnorm_ps, expbs_zs, expbs_ts, expbs_t_ps, expbs_z_ps, expbs_1s_w_ps, conf,
     expbs_1s_w_intervs, expbs_1s_w_normintervs, expbs_cl_labelsets, expbs_bs_ps,
     expbs_bs_zs, expbs_bs_z_ps, muscle_dics, expected_clusterss] = pickle.load(f)
print ('Expected functional clusters: (difference between muscles'
       ' that cross the same DOF but belong to different'
       ' clusters)')
for i, expbs_cl_label in enumerate(expbs_cl_labels):
    s = '\t{:2d} clusters. Number of elements: {}'.format(
        expbs_numcls[i], len(expbs_fbss[i]))
    s += ' Labels: {}'.format(expbs_cl_label)

    s += '\n\tTwo sample tests:'
    s += ('\n\t\tIn cluster len {}, mean {}, median {}, std {},'
          ' normality p-value {}').format(
        len(expbs_fbss_in[i]),
        numpy.mean(expbs_fbss_in[i]),
        numpy.median(expbs_fbss_in[i]),
        numpy.std(expbs_fbss_in[i]), expbs_in_tnorm_ps[i])
    s += ('\n\t\tOut cluster len {}, mean {}, median {}, std {},'
          ' normality p-value {}').format(
        len(expbs_fbss_out[i]),
        numpy.mean(expbs_fbss_out[i]),
        numpy.median(expbs_fbss_out[i]),
        numpy.std(expbs_fbss_out[i]), expbs_ou_tnorm_ps[i])
    s += ('\n\t\t2-sample Z-score {}, p-value {}.').format(
        expbs_2smpl_zs[i], expbs_2smpl_z_ps[i])
    s += ('\n\t\t2-sample Welchs T-score {}, p-value {}.').format(
        expbs_2smpl_ts[i], expbs_2smpl_t_ps[i])
    s += ('\n\t\tMann-Whitney U {}, nx*ny {}, p-value {}.').format(
        expbs_us[i], expbs_u_nxnys[i], expbs_u_ps[i])
    s += ('\n\t\tWilcoxon signed-rank test T {},'
          ' p-value {}.').format(
        expbs_ws[i], expbs_w_ps[i])

    s += ('\n\tOne difference sample tests:')
    s += ('\n\t\tDifference len {}, mean {}, median {}, std {},'
          ' normality p-value {}').format(
        len(expbs_fbss[i]),
        numpy.mean(expbs_fbss[i]), numpy.median(expbs_fbss[i]),
        numpy.std(expbs_fbss[i]), expbs_tnorm_ps[i])
    s += ('\n\t\tZ-score from 0: {}, p-value {}').format(
        expbs_zs[i], expbs_z_ps[i])
    s += ('\n\t\tT-score from 0: {}, p-value {}').format(
        expbs_ts[i], expbs_t_ps[i])
    s += ('\n\t\tSign 1-sample test p-value {},'
          ' confidence {} interval {}'
          ' (assuming normal {})').format(
        expbs_1s_w_ps[i], conf, expbs_1s_w_intervs[i],
        expbs_1s_w_normintervs[i])

    print s

    for j in xrange(expbs_numcls[i]):
        if expbs_cl_labelsets[i][j] is None:
            continue
        s = '\t\tCluster {:2d} with {:3d} elements.'.format(
            expbs_cl_labelsets[i][j],
            sum(1 for k in expbs_cl_label
                if k == expbs_cl_labelsets[i][j]))
        s += (' T-test for 0 p-value: {:8e}.').format(
            expbs_bs_ps[i][j])
        s += (' Z-test for 0: {:10f}, p-value: {:8e}.').format(
            expbs_bs_zs[i][j], expbs_bs_z_ps[i][j])
        s += ' Muscles in the cluster: {}'.format(
            ', '.join([muscle_dics[k]['sMuscle']
                       for k, l in enumerate(expbs_cl_label)
                       if l == expbs_cl_labelsets[i][j]]))
        print s

for i in xrange(len(expected_clusterss)):
    pylab.figure()
    try:
        pylab.subplot(2, 2, 1, sharex=bsbf_ax_left)
        expbs_ax_left = bsbf_ax_left
    except NameError:
        expbs_ax_left = pylab.subplot(2, 2, 1)
    pylab.boxplot([expbs_fbss_in[i], expbs_fbss_out[i]],
                  vert=False,
                  labels=['in', 'out'],
                  whis=[0, 100])

    expbs_ax_bottom = pylab.subplot(2, 2, 3, sharex=expbs_ax_left)
    pylab.hist(expbs_fbss_in[i],
               weights=density_weights(expbs_fbss_in[i]),
               bins=bins, label='in', alpha=0.5)
    pylab.hist(expbs_fbss_out[i],
               weights=density_weights(expbs_fbss_out[i]),
               bins=bins, label='out', alpha=0.5)
    expbs_ax_left.set_xlim(xmin=0, xmax=1.5)
    pylab.xlabel('Distance between vectors')
    pylab.ylabel('Probability')
    pylab.legend()

    try:
        pylab.subplot(2, 2, 2, sharex=bsbf_ax_right)
        expbs_ax_right = bsbf_ax_right
    except NameError:
        expbs_ax_right = pylab.subplot(2, 2, 2)
    pylab.boxplot(expbs_fbss[i], vert=False,
                  whis=[0, 100])
    pylab.axvline(x=0, color='m')
    pylab.axvline(x=expbs_1s_w_intervs[i][1],
                  color='r', linestyle='--')

    pylab.subplot(2, 2, 4, sharex=expbs_ax_right,
                  sharey=expbs_ax_bottom)
    pylab.hist(expbs_fbss[i],
               weights=density_weights(expbs_fbss[i]),
               bins=bins)
    pylab.axvline(x=0, color='m')
    pylab.axvline(x=expbs_1s_w_intervs[i][1],
                  color='r', linestyle='--')
    pylab.xlabel('Difference in distance between vectors')
    pylab.ylabel('Probability')

    pylab.suptitle('{} functional clusters on DOF'.format(
        expbs_numcls[i]))

pylab.show()
