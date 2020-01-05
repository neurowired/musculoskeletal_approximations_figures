#!python2
import os
import sys
import csv
import pickle
import time
import argparse
from itertools import izip
import copy
import math
import random

import numpy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm
import scipy.stats
from scipy.cluster import hierarchy
from scipy.spatial import distance
import pylab
import sklearn
import sklearn.cluster
import sklearn.metrics
import sklearn.decomposition

import gepapy

import meta


EXPECTED_CLUSTERS7 = [
    ['BIC_SH', 'BIC_LO', 'SUP'],
    ['PT', 'PQ'],
    ['ECR_LO', 'ECR_BR', 'ECU'],
    ['FCR', 'FCU', 'PL'],
    ['FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5'],
    ['EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERS6 = [
    ['BIC_SH', 'BIC_LO', 'SUP', 'PT', 'PQ'],
    ['ECR_LO', 'ECR_BR', 'ECU'],
    ['FCR', 'FCU', 'PL'],
    ['FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5'],
    ['EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERS5 = [
    ['BIC_SH', 'BIC_LO', 'SUP', 'PT', 'PQ'],
    ['ECR_LO', 'ECR_BR', 'ECU', 'FCR', 'FCU', 'PL'],
    ['FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5'],
    ['EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERS4 = [
    ['BIC_SH', 'BIC_LO', 'SUP', 'PT', 'PQ'],
    ['ECR_LO', 'ECR_BR', 'ECU', 'FCR', 'FCU', 'PL'],
    ['FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5',
     'EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERS3 = [
    ['BIC_SH', 'BIC_LO', 'SUP', 'PT', 'PQ'],
    ['ECR_LO', 'ECR_BR', 'ECU', 'FCR', 'FCU', 'PL',
     'FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5',
     'EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERS2 = [
    ['BIC_SH', 'BIC_LO', 'SUP', 'PT', 'PQ',
     'ECR_LO', 'ECR_BR', 'ECU', 'FCR', 'FCU', 'PL',
     'FDS2', 'FDS3', 'FDS4', 'FDS5', 'FDP2', 'FDP3', 'FDP4', 'FDP5',
     'EDM', 'ED2', 'ED3', 'ED4', 'ED5', 'EIND'],
    ['APL', 'OP', 'APB', 'EPL', 'EPB', 'FPB', 'FPL', 'ADPT'],
]
EXPECTED_CLUSTERSS = [EXPECTED_CLUSTERS2,
                      EXPECTED_CLUSTERS3,
                      EXPECTED_CLUSTERS4,
                      EXPECTED_CLUSTERS5,
                      EXPECTED_CLUSTERS6,
                      EXPECTED_CLUSTERS7
                      ]


def si_muscles(filename1, filename2):
    with open(filename1) as f:
        ma_approx1, len_approx1 = pickle.load(f)
    with open(filename2) as f:
        ma_approx2, len_approx2 = pickle.load(f)

    dim1 = len(ma_approx1)
    dim2 = len(ma_approx2)

    er1 = gepapy.ersatz.Ersatz(dim1, len_approx1['config'],
                               point=len_approx1['point'])
    er2 = gepapy.ersatz.Ersatz(dim2, len_approx2['config'],
                               point=len_approx2['point'])
    return er1.similarity_index(
        er2, agnostic=False, count_reps=False,
        use_point=False,
        id_dofsa=None, id_dofsb=None,
        use_point_difference=False, normalize=True, vsp=False)


def internal_consistency(muscle_ids, plot, bins):
    muscle_ids.sort()
    muscles = []
    with open(meta.METAMUSCLE_FNAME) as f:
        fd = csv.DictReader(f)
        for line in fd:
            if int(line['idMuscle']) in muscle_ids:
                muscles.append(line['sMuscle'])
    print muscles

    sis = []
    for muscle in muscles:
        sis.append(si_muscles(
            meta.apprfname(meta.APPROX_NOT_INTCONS_DIRNAME, muscle),
            meta.apprfname(meta.APPROX_DIRNAME, muscle)))

    sis = [si*100 for si in sis]

    bufnummus = sum(1. for i in sis if i >= 80)

    print 'Differences: {}'.format(
        ', '.join(['{}: {:4.1f}'.format(muscle, s)
                  for muscle, s in izip(muscles, sis)]))
    print 'Mean difference: {}, percent of muscles over 80%: {}, {} total'.format(
        numpy.mean(sis), bufnummus/len(sis)*100., bufnummus)

    if plot:
        pylab.figure()
        bins = numpy.linspace(0, 100, 26)
        pylab.hist(sis, bins=bins, color='k')
        pylab.xlim([0, 100])
        pylab.xlabel('Si between muscles before and after internal '
                     'consistency, %')
        pylab.ylabel('Number of muscles')
        pylab.title('Mids len: {}'.format(len(muscle_ids)))

        with open('msd_fig5_1.pickle', 'w') as f:
            pickle.dump([muscle_ids, sis, bins], f)


def power_support_vector(spv, dim):
    pow_sup_vec = [0]*len(spv)
    for i, e in enumerate(spv):
        if e.count(',') < dim:
            pow_sup_vec[i] = 1
    return pow_sup_vec


def proximity_matrix_from_flat_cluster(n, fcluster):
    X = numpy.zeros((n, n))
    for i in xrange(n):
        for j in xrange(n):
            if fcluster[i] != fcluster[j]:
                X[i][j] = 1.
    return X


def z_score_test(x, val, hypothesis='two-sided'):
    """Hypothesis can be less, greater or two-sided"""
    m = numpy.mean(x)
    std = numpy.std(x)
    if isinstance(val, (int, float)):
        z_score = (val - m) / std
    else:
        m2 = numpy.mean(val)
        std2 = numpy.std(val)
        z_score = (m2 - m) / (std**2/len(x)+std2**2/len(val))**0.5

    if hypothesis == 'two-sided':
        p = min([scipy.stats.norm.cdf(z_score),
                 scipy.stats.norm.cdf(-z_score)])*2
    elif hypothesis == 'greater':
        p = scipy.stats.norm.cdf(-z_score)
    elif hypothesis == 'less':
        p = scipy.stats.norm.cdf(z_score)
    else:
        raise ValueError('Wrong hypothesis specified')
    return (z_score, p)


def sign_test_1sample(x, val, hypothesis='two-sided',
                      zero_method='wilcox', conf=None, conf_norm=False):
    """@http://www.stat.umn.edu/geyer/old03/5102/notes/rank.pdf

    less means that the median is less than value
    greater means that the median is greater than value
    """
    less = sum(1. for i in x if i < val)
    equal = sum(1. for i in x if i == val)
    n = len(x)

    if zero_method == 'wilcox':
        n -= equal
    elif zero_method == 'pratt':
        less += equal
    elif zero_method == 'zsplit':
        less += int(equal/2)
    else:
        raise ValueError('Wrong zero_method specified')

    if hypothesis == 'two-sided':
        p = 2*min([scipy.stats.binom.cdf(less, n, 0.5),
                   scipy.stats.binom.cdf(n-less, n, 0.5)])
    elif hypothesis == 'less':
        p = scipy.stats.binom.cdf(n-less, n, 0.5)
    elif hypothesis == 'greater':
        p = scipy.stats.binom.cdf(less, n, 0.5)
    else:
        raise ValueError('Wrong hypothesis specified')

    if conf is not None:
        sortedx = sorted(x)
        intrv = [-numpy.inf, numpy.inf]
        if hypothesis == 'two-sided':
            for i in sortedx:
                p_ = sign_test_1sample(sortedx, i, hypothesis='less',
                                       zero_method=zero_method)
                if p_ < 1 - conf/2.:
                    intrv[1] = i
                    break
            for i in reversed(sortedx):
                p_ = sign_test_1sample(sortedx, i, hypothesis='greater',
                                       zero_method=zero_method)
                if p_ < 1 - conf/2.:
                    intrv[0] = i
                    break
        elif hypothesis == 'less':
            for i in sortedx[len(x)/2:]:
                p_ = sign_test_1sample(sortedx, i, hypothesis='less',
                                       zero_method=zero_method)
                if p_ < 1 - conf/2.:
                    intrv[1] = i
                    break
        elif hypothesis == 'greater':
            for i in reversed(sortedx)[len(x)/2:]:
                p_ = sign_test_1sample(sortedx, i, hypothesis='greater',
                                       zero_method=zero_method)
                if p_ < 1 - conf/2.:
                    intrv[0] = i
                    break

        if conf_norm:
            # return from normal assumption
            lind = int(n/2.-1.96*float(n)**0.5/2.) - 1
            uind = int(1+n/2.+1.96*float(n)**0.5/2.) - 1
            return (p, intrv, [sortedx[lind], sortedx[uind]])
        return (p, intrv)

    return p


def alltoall(muscle_ids, plot, dendro_link, bins, verbose,
             do_pca, do_dendro,
             do_bsbf, do_expbs, num_expbs, plot_basic, more_plots, use_ddpv):
    approx_cons_dirname = meta.APPROX_DIRNAME
    metamuscle_list = meta.metamuscle_list()
    conf = 0.99

    if verbose > 0:
        print 'Loading data...'
    n = len(muscle_ids)
    muscle_dics = []
    for muscle_id in muscle_ids:
        muscle_dics.append(meta.search_dic(metamuscle_list, idMuscle=muscle_id))
        muscle_dics[-1]['idDOFList'] = [
            meta.DOF_SIMILARITIES[j]
            for j in muscle_dics[-1]['idDOFList']]
        muscle_dics[-1]['dim'] = len(muscle_dics[-1]['idDOFList'])
        approx_fname = meta.apprfname(meta.APPROX_DIRNAME,
                                      muscle_dics[-1]['sMuscle'])
        with open(approx_fname) as f:
            _, len_approx = pickle.load(f)
        muscle_dics[-1]['ersatz'] = gepapy.ersatz.Ersatz(
            muscle_dics[-1]['dim'], len_approx['config'],
            point=len_approx['point'])

    dofs_in_question = set()
    for m in muscle_dics:
        for id in m['idDOFList']:
            dofs_in_question.add(id)
    dofs_in_question = sorted(list(dofs_in_question))

    # expected clustering
    i = [len(i) for i in EXPECTED_CLUSTERSS].index(num_expbs)
    expected_clusterss = [EXPECTED_CLUSTERSS[i]]
    expected_clusters = [next((i for (i, d) in enumerate(expected_clusterss[0])
                          if s['sMuscle'] in d), None)
                         for s in muscle_dics]

    # dimensionality info and report
    avg_dim = numpy.mean([i['dim'] for i in muscle_dics])
    max_dim = max([i['dim'] for i in muscle_dics])
    spv = gepapy.ersatz.parameter_vector_base(max_dim, meta.DEFAULTS['rho'])
    if verbose > 0:
        print '\t', muscle_ids
        print '\t', [i['sMuscle'] for i in muscle_dics]
        print '\tExpected clusters:', expected_clusters
    if verbose > 1:
        print '\tDOFs they cross:', [i['idDOFList'] for i in muscle_dics]
    if verbose > 0:
        print '\tAverage dimensionality of muscles is {}'.format(avg_dim)
        print '\tMaximum dimensionality of muscles is {}'.format(max_dim)
        print '\tDOFs that muscles cross: {}'.format(
            ', '.join(str(i) for i in dofs_in_question))
        print '\tEnclosing supporting vector'
        print '\t', spv

    ######## Generate parameter vectors
    if verbose > 0:
        print 'Generating parameter vectors and proximity matrix...'
    for i in xrange(n):
        muscle_dics[i]['param_vec'] = numpy.array(
            muscle_dics[i]['ersatz'].parameter_vector(spv))
        # check the vector
        if abs(numpy.linalg.norm(muscle_dics[i]['param_vec']) - 1) > 1e-8:
            print 'Muscle {} is not one-vector'.format(
                muscle_dics[i]['sMuscle'])
        if verbose > 1:
            print '\t{:10s}: {}'.format(
                muscle_dics[i]['sMuscle'], ', '.join(
                    ['{:6.4f}'.format(j) for j in muscle_dics[i]['param_vec']]))
    X = [i['param_vec'] for i in muscle_dics]  # for machine learning

    ######## Fill the distance matrix
    mtr = numpy.zeros((n, n))
    for i, imdic in enumerate(muscle_dics):
        for j in xrange(i+1, n):
            jmdic = muscle_dics[j]
            if use_ddpv:
                mtr[i][j] = imdic['ersatz'].structural_difference(
                    jmdic['ersatz'],
                    id_dofsa=imdic['idDOFList'],
                    id_dofsb=jmdic['idDOFList'],
                    normalize=True)*100.
            else:
                mtr[i][j] = numpy.linalg.norm(
                    imdic['param_vec'] - jmdic['param_vec'])

    # Fill other triangle and a list
    mtr_list = [];
    for i in xrange(1, n):
        for j in xrange(i):
            mtr[i][j] = mtr[j][i]
            mtr_list.append(mtr[i][j])

    if verbose > 1:
        print '\tGenerated proximity matrix:'
        print '\n'.join(' '.join(['{:6.4f}'.format(j) for j in i])
                        for i in mtr)

    if use_ddpv:
        max_mtr = 100
    else:
        max_mtr = max(mtr_list)
    pdist_mtr = distance.squareform(mtr)

    ######## Dendrogram clustering
    if do_dendro:
        if verbose > 0:
            print 'Dendrogram generation and clustering...'
        Z_linkage = hierarchy.linkage(
            pdist_mtr, method=dendro_link, optimal_ordering=True)
        cophenetic, _ = scipy.cluster.hierarchy.cophenet(Z_linkage, Y=pdist_mtr)
        R_incons = scipy.cluster.hierarchy.inconsistent(Z_linkage, d=2)
        mean_incons = numpy.mean([i[3] for i in R_incons])
        median_incons = numpy.median([i[3] for i in R_incons])
        std_incons = numpy.std([i[3] for i in R_incons])
        disp_incons = median_incons + std_incons
        using_incons = disp_incons

        fcluster = scipy.cluster.hierarchy.fcluster(
            Z_linkage, using_incons,
            criterion='inconsistent', R=R_incons)
        num_incons_clusters = len(set(fcluster))

    ######## PCA clustering
    if do_pca:
        if verbose > 0:
            print 'Calculating PCA decomposition...'
        pcas = sklearn.decomposition.PCA(svd_solver='full').fit(X)

        # Project points onto cluster centers
        if verbose > 0:
            print 'Projecting points to PCA vectors...'
        for i in xrange(len(muscle_dics)):
            muscle_dics[i]['PCA_param_vec'] = pcas.transform(
                numpy.array(muscle_dics[i]['param_vec']).reshape(1, -1))[0]

    ######## Test on the structure
    if do_bsbf:
        if verbose > 0:
            print 'Analysis of structure...'

        # Compare two distributions
        bsbf_incluster_dist = []
        bsbf_outcluster_dist = []
        for dof in dofs_in_question:
            for m in [_m for _m in muscle_dics if dof in _m['idDOFList']]:
                for mc in muscle_dics:
                    if mc['idMuscle'] == m['idMuscle']:  # skip same muscle
                        continue
                    if dof in mc['idDOFList']:
                        bsbf_incluster_dist.append(
                            numpy.linalg.norm(m['param_vec']-mc['param_vec']))
                    else:
                        bsbf_outcluster_dist.append(
                            numpy.linalg.norm(m['param_vec']-mc['param_vec']))

        bsbf_u, bsbf_u_p = scipy.stats.mannwhitneyu(
            bsbf_incluster_dist, bsbf_outcluster_dist,
            alternative='less')
        bsbf_u_nxny = len(bsbf_incluster_dist)*len(bsbf_outcluster_dist)

        bsbf_w_len = min((len(bsbf_incluster_dist), len(bsbf_outcluster_dist)))
        bsbf_w_in = (bsbf_incluster_dist
                     if len(bsbf_incluster_dist) == bsbf_w_len else
                     random.sample(bsbf_incluster_dist, bsbf_w_len))
        bsbf_w_out = (bsbf_outcluster_dist
                     if len(bsbf_outcluster_dist) == bsbf_w_len else
                     random.sample(bsbf_outcluster_dist, bsbf_w_len))
        bsbf_w, bsbf_w_p = scipy.stats.wilcoxon(
            bsbf_w_in, bsbf_w_out,
            zero_method='wilcox', correction=False)

        bsbf_in_normt, bsbf_in_normt_p = scipy.stats.normaltest(
            bsbf_incluster_dist)
        bsbf_ou_normt, bsbf_ou_normt_p = scipy.stats.normaltest(
            bsbf_outcluster_dist)

        # One distribution of differences
        muscles_in_dof = [[j for j, m in enumerate(muscle_dics)
                           if i in m['idDOFList']]
                          for i in dofs_in_question]
        muscles_out_dof = [[j for j, m in enumerate(muscle_dics)
                            if i not in m['idDOFList']]
                           for i in dofs_in_question]
        bsbf_bs = []
        # all data instead for random sample
        for dofid, dof in enumerate(dofs_in_question):
            for m in muscles_in_dof[dofid]:
                for m_in in muscles_in_dof[dofid]:
                    if m_in == m:  # skip same muscle
                        continue
                    for m_out in muscles_out_dof[dofid]:
                        dist_in = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_in]['param_vec'])
                        dist_out = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_out]['param_vec'])
                        bsbf_bs.append(dist_in - dist_out)

        bsbf_bs_pmean = 0.
        bsbf_bs_t, bsbf_bs_t_p = scipy.stats.ttest_1samp(
            bsbf_bs, bsbf_bs_pmean)
        bsbf_bs_z, bsbf_bs_z_p = z_score_test(bsbf_bs, bsbf_bs_pmean,
                                              hypothesis='greater')
        bsbf_bs_neg_per = float(sum(1 for i in bsbf_bs if i < 0))/len(bsbf_bs)

        my_ttest = lambda a, m: (numpy.mean(a)-m)/(
            numpy.std(a)/numpy.sqrt(len(a)))
        bsbf_bs_my_t = my_ttest(bsbf_bs, bsbf_bs_pmean)
        bsbf_bs_my_t_p = scipy.stats.t.cdf(bsbf_bs_my_t, df=len(bsbf_bs)-1)

        bsbf_normt, bsbf_normt_p = scipy.stats.normaltest(bsbf_bs)
        bsbf_bs_w_p, bsbf_bs_w_interv, bsbf_bs_w_norm_interv = sign_test_1sample(
            bsbf_bs, bsbf_bs_pmean,
            hypothesis='less', zero_method='wilcox', conf=conf, conf_norm=True)

    ######## Test on the function
    if do_expbs:
        if verbose > 0:
            print 'Analysis of expected functional clusters...'
        expbs_cl_labels = [[next((i for (i, d) in enumerate(ec)
                                  if s['sMuscle'] in d), None)
                          for s in muscle_dics]
                         for ec in expected_clusterss]
        expbs_cl_labelsets = [sorted(list(set(i))) for i in expbs_cl_labels]
        expbs_numcls = [len(i) for i in expbs_cl_labelsets]
        expbs_p_matrs = [proximity_matrix_from_flat_cluster(n, i)
                         for i in expbs_cl_labels]

        if verbose > 0:
            print '\tInout cluster distances'
        expbs_inouts = []
        # [i][j][k]: i - expected clustering, j - cluster, k - measurement
        expbs_bss = []
        expbs_bs_ps = []
        expbs_bs_zs = []
        expbs_bs_z_ps = []
        for expbs_numcl in expbs_numcls:
            expbs_inouts.append([])
            expbs_bss.append([])
            expbs_bs_ps.append([])
            expbs_bs_zs.append([])
            expbs_bs_z_ps.append([])
            for _ in xrange(expbs_numcl):
                expbs_inouts[-1].append({'in': [], 'out': []})
                expbs_bss[-1].append([])
                expbs_bs_ps[-1].append([])
                expbs_bs_zs[-1].append([])
                expbs_bs_z_ps[-1].append([])

        for i, expbs_cl_label in enumerate(expbs_cl_labels):
            clusters = expbs_cl_labelsets[i]
            muscles_in_cluster = [[j for j in xrange(len(muscle_dics))
                                   if expbs_cl_label[j]==k]
                                  for k in clusters]
            muscles_out_cluster = [[j for j in xrange(len(muscle_dics))
                                    if expbs_cl_label[j]!=k]
                                   for k in clusters]
            for cluster_id, cluster in enumerate(clusters):
                for m in muscles_in_cluster[cluster_id]:
                    for m_in in muscles_in_cluster[cluster_id]:
                        if m_in == m:  # skip same muscle
                            continue
                        dist_in = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_in]['param_vec'])
                        expbs_inouts[i][cluster_id]['in'].append(dist_in)

                        for m_out in muscles_out_cluster[cluster_id]:
                            dist_out = numpy.linalg.norm(
                                muscle_dics[m]['param_vec'] -
                                muscle_dics[m_out]['param_vec'])
                            expbs_bss[i][cluster_id].append(dist_in-dist_out)

                    for m_out in muscles_out_cluster[cluster_id]:
                        dist_out = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_out]['param_vec'])
                        expbs_inouts[i][cluster_id]['out'].append(dist_out)

            # test the above for being zero
            for cluster_id in xrange(expbs_numcls[i]):
                _, expbs_bs_ps[i][cluster_id] = scipy.stats.ttest_1samp(
                    expbs_bss[i][cluster_id], 0.)
                (expbs_bs_zs[i][cluster_id],
                 expbs_bs_z_ps[i][cluster_id]) = z_score_test(
                    expbs_bss[i][cluster_id], 0, hypothesis='greater')


        # Functional testing histogram
        if verbose > 0:
            print '\tGenerating functional bootstrapping histogram'
        expbs_diq = dofs_in_question
        if 10 in expbs_diq:
            expbs_diq.pop(expbs_diq.index(10))
            if verbose > 0:
                print '\t\tRemoving DIP from dofs in question'
        expbs_fbss = []
        expbs_fbss_in = []
        expbs_fbss_out = []
        expbs_us = [0]*len(expected_clusterss)
        expbs_u_ps = [0]*len(expected_clusterss)
        expbs_u_nxnys = [0]*len(expected_clusterss)
        expbs_ws = [0]*len(expected_clusterss)
        expbs_w_ps = [0]*len(expected_clusterss)
        expbs_zs = [0]*len(expected_clusterss)
        expbs_z_ps = [0]*len(expected_clusterss)
        expbs_ts = [0]*len(expected_clusterss)
        expbs_t_ps = [0]*len(expected_clusterss)
        expbs_1s_w_ps = [0]*len(expected_clusterss)
        expbs_1s_w_intervs = [0]*len(expected_clusterss)
        expbs_1s_w_normintervs = [0]*len(expected_clusterss)
        expbs_2smpl_ts = [0]*len(expected_clusterss)
        expbs_2smpl_t_ps = [0]*len(expected_clusterss)
        expbs_2smpl_zs = [0]*len(expected_clusterss)
        expbs_2smpl_z_ps = [0]*len(expected_clusterss)
        expbs_tnorms = [0]*len(expected_clusterss)
        expbs_tnorm_ps = [0]*len(expected_clusterss)
        expbs_in_tnorms = [0]*len(expected_clusterss)
        expbs_in_tnorm_ps = [0]*len(expected_clusterss)
        expbs_ou_tnorms = [0]*len(expected_clusterss)
        expbs_ou_tnorm_ps = [0]*len(expected_clusterss)
        for c in xrange(len(expected_clusterss)):
            expbs_fbss.append([])
            expbs_fbss_in.append([])
            expbs_fbss_out.append([])
            # generate subcluster groups
            muscles_in_dof = [[j for j, m in enumerate(muscle_dics)
                               if i in m['idDOFList']]
                              for i in expbs_diq]

            for dofind, dof in enumerate(expbs_diq):
                for m in muscles_in_dof[dofind]:
                    # muscles from same cluster crossing the DOF
                    mfsc = [i for i in muscles_in_dof[dofind]
                            if expbs_cl_labels[c][i] == expbs_cl_labels[c][m]]
                    # muscles from other clusters crossing the DOF
                    mfoc = [i for i in muscles_in_dof[dofind]
                            if expbs_cl_labels[c][i] != expbs_cl_labels[c][m]]
                    for m_in in mfsc:
                        if m == m_in:
                            continue
                        dist_in = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_in]['param_vec'])
                        expbs_fbss_in[-1].append(dist_in)
                        for m_out in mfoc:
                            dist_out = numpy.linalg.norm(
                                muscle_dics[m]['param_vec'] -
                                muscle_dics[m_out]['param_vec'])

                            expbs_fbss[-1].append(dist_in - dist_out)

                    for m_out in mfoc:
                        dist_out = numpy.linalg.norm(
                            muscle_dics[m]['param_vec'] -
                            muscle_dics[m_out]['param_vec'])

                        expbs_fbss_out[-1].append(dist_out)

                if verbose > 2:
                    print ('DOF #{:2d}, muscle {:>8s}, cluster {},'
                           ' muscle in cluster {:>8s}, out cluster {:>8s}'
                           '(cluster {}).').format(
                        dof, muscle_dics[m]['sMuscle'], expbs_cl_labels[c][m],
                        muscle_dics[m_in]['sMuscle'],
                        muscle_dics[m_out]['sMuscle'],
                        expbs_cl_labels[c][m_out])

            # probability testing time
            expbs_in_tnorms[c], expbs_in_tnorm_ps[c] = scipy.stats.normaltest(
                expbs_fbss[-1])
            expbs_ou_tnorms[c], expbs_ou_tnorm_ps[c] = scipy.stats.normaltest(
                expbs_fbss[-1])

            expbs_us[c], expbs_u_ps[c] = scipy.stats.mannwhitneyu(
                expbs_fbss_in[-1], expbs_fbss_out[-1],
                alternative='less')
            expbs_u_nxnys[c] = len(expbs_fbss_in[-1])*len(expbs_fbss_out[-1])

            expbs_w_len = min((len(expbs_fbss_in[-1]), len(expbs_fbss_out[-1])))
            expbs_w_in = (expbs_fbss_in[-1]
                         if len(expbs_fbss_in[-1]) == expbs_w_len else
                         random.sample(expbs_fbss_in[-1], expbs_w_len))
            expbs_w_out = (expbs_fbss_out[-1]
                         if len(expbs_fbss_out[-1]) == expbs_w_len else
                         random.sample(expbs_fbss_out[-1], expbs_w_len))
            expbs_ws[c], expbs_w_ps[c] = scipy.stats.wilcoxon(
                expbs_w_in, expbs_w_out,
                zero_method='wilcox', correction=False)

            expbs_2smpl_zs[c], expbs_2smpl_z_ps[c] = z_score_test(
                expbs_fbss_in[-1], expbs_fbss_out[-1],
                hypothesis='greater')

            expbs_2smpl_ts[c], expbs_2smpl_t_ps[c] = scipy.stats.ttest_ind(
                expbs_fbss_in[-1], expbs_fbss_out[-1],
                equal_var=False)

            expbs_zs[c], expbs_z_ps[c] = z_score_test(
                expbs_fbss[-1], 0, hypothesis='greater')

            expbs_ts[c], expbs_t_ps[c] = scipy.stats.ttest_1samp(
                expbs_fbss[-1], 0.)

            expbs_tnorms[c], expbs_tnorm_ps[c] = scipy.stats.normaltest(
                expbs_fbss[-1])
            expbs_1s_w_ps[c], expbs_1s_w_intervs[c], expbs_1s_w_normintervs[c]  = sign_test_1sample(
                expbs_fbss[-1], 0, hypothesis='less', zero_method='wilcox',
                conf=conf, conf_norm=True)

    ######## Print
    if verbose > 0:
        print '\t~~~~REPORT~~~~'

    # Report on inconsistency dendro clustering
    if do_dendro:
        if verbose > 0:
            print 'Linkage using {}:'.format(dendro_link)
        if verbose > 1:
            print '\n'.join([
                'Cluster {:2d} contains {:2d} and {:2d} clusters. '
                'Height of link: {:10.6f}, std {:10.6f}, contains {:2d} elements. '
                'Inconsistency numlinks {:2d}, coefficient {:10.6f}'.format(
                    i+n, int(e[0]), int(e[1]), e[2], c[1], int(e[3]),
                    int(c[2]), c[3])
                for i, (e, c) in enumerate(izip(Z_linkage, R_incons))])
        if verbose > 0:
            print ('Mean inconsistency coefficient {}, median {}, std {},'
                   ' disp {}').format(
                mean_incons, median_incons, std_incons, disp_incons)
            print 'Reduced tree has {} clusters, varlen: {}:'.format(
                len(set(fcluster)), len(fcluster))
            print fcluster
            print 'Grouped hierarchical incons clusters:'
            for i in xrange(max(fcluster)):
                print '\tCluster', i+1, [muscle_dics[j]['sMuscle']
                                         for j, fc in enumerate(fcluster)
                                         if fc==i+1]
            print 'Cophenetic correlation distance:', cophenetic

    # Report on PCA
    if do_pca:
        if verbose > 0:
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
            with open('msd_fig6_3.pickle', 'w') as f:
                pickle.dump([pcas, spv], f)

    # Report on structural difference
    if do_bsbf:
        if verbose > 0:
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
            with open('msd_fig7_1.pickle', 'w') as f:
                pickle.dump([
                    bsbf_incluster_dist, bsbf_outcluster_dist, bsbf_in_normt, bsbf_in_normt_p,
                    bsbf_ou_normt, bsbf_ou_normt_p, bsbf_u, bsbf_u_nxny, bsbf_u_p,
                    bsbf_w, bsbf_w_p, bsbf_bs, bsbf_bs_neg_per, bsbf_bs_pmean, bsbf_bs_t, bsbf_bs_t_p,
                    bsbf_bs_z, bsbf_bs_z_p, bsbf_normt, bsbf_normt_p,
                    bsbf_bs_w_p, conf, bsbf_bs_w_interv, bsbf_bs_w_norm_interv, bins], f)

    # Report on analysis of expected functional clusters
    if do_expbs:
        if verbose > 0:
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

            with open('msd_fig7_2.pickle', 'w') as f:
                pickle.dump([expbs_cl_labels, expbs_numcls, expbs_fbss, expbs_fbss_in,
                             expbs_in_tnorm_ps, expbs_fbss_out, expbs_ou_tnorm_ps, expbs_2smpl_zs,
                             expbs_2smpl_z_ps, expbs_2smpl_ts, expbs_2smpl_t_ps, expbs_us,
                             expbs_u_nxnys, expbs_u_ps, expbs_ws, expbs_w_ps, expbs_fbss,
                             expbs_tnorm_ps, expbs_zs, expbs_z_ps, expbs_ts, expbs_t_ps,
                             expbs_1s_w_ps, conf,
                             expbs_1s_w_intervs, expbs_1s_w_normintervs, expbs_cl_labelsets,
                             expbs_bs_ps,
                             expbs_bs_zs, expbs_bs_z_ps, muscle_dics, expected_clusterss], f)


    ######## Plot
    if plot:
        density_weights = lambda x: (numpy.ones_like(x)/float(len(x)))

        if plot_basic:
            # histogram
            pylab.figure()
            pylab.hist(mtr_list, bins=numpy.linspace(0, max_mtr, bins), color='k')
            pylab.xlim([0, max_mtr])
            pylab.xlabel('Difference between muscles, au')
            pylab.ylabel('Number of muscle pairs')
            pylab.title('Mid len: {}'.format(n))

            # heatmap
            pylab.figure()
            if use_ddpv:
                pylab.pcolor(mtr, cmap=pylab.get_cmap('gist_yarg'),
                             vmin=0, vmax=100)
            else:
                pylab.pcolor(mtr, cmap=pylab.get_cmap('gist_yarg'))
            pylab.colorbar()
            pylab.xticks([i+0.5 for i in xrange(n)],
                         [i['sMuscle'] for i in muscle_dics],
                         rotation='vertical')
            pylab.yticks([i+0.5 for i in xrange(n)],
                         [i['sMuscle'] for i in muscle_dics])
            pylab.title('Mid len: {}'.format(n))
            pylab.xlim(0, n)
            pylab.ylim(0, n)

            if use_ddpv:
                with open('msd_supp_fig2_1.pickle', 'w') as f:
                    pickle.dump([mtr_list, max_mtr, bins, mtr, muscle_dics, n], f)
            else:
                with open('msd_fig6_1.pickle', 'w') as f:
                    pickle.dump([mtr_list, max_mtr, bins, mtr, muscle_dics, n], f)

        # dendrogram
        if do_dendro:
            llf_addlen = 5
            def llf(id):
                if id < n:
                    return (muscle_dics[id]['sMuscle'] +
                            ' {:d} {:2d}').format(
                        muscle_dics[id]['dim'], fcluster[id])
                else:
                    return '[%d %d %1.2f]' % (id, count, R[n-id,3])
            pylab.figure()
            R = hierarchy.dendrogram(
                Z_linkage, orientation='left',
                leaf_label_func=llf,
                )
            pylab.xlabel('{} distance between clusters, au'.format(dendro_link))
            pylab.title('Dendrogram')
            if use_ddpv:
                pylab.xlim([100, 0])

            if verbose > 1:
                print('The order of muscles (ids) in the dendrogram:')
                print([meta.search_dic(metamuscle_list,
                                       sMuscle=i[:-llf_addlen])['idMuscle']
                       for i in R['ivl']])
                print([i[:-llf_addlen] for i in R['ivl']])
                print('And reversed:')
                print([meta.search_dic(metamuscle_list,
                                       sMuscle=i[:-llf_addlen])['idMuscle']
                       for i in reversed(R['ivl'])])
                print([i[:-llf_addlen] for i in reversed(R['ivl'])])

            if use_ddpv:
                with open('msd_supp_fig2_2.pickle', 'w') as f:
                    pickle.dump([muscle_dics, n, Z_linkage, fcluster, dendro_link, metamuscle_list], f)
            else:
                with open('msd_fig6_2.pickle', 'w') as f:
                    pickle.dump([muscle_dics, n, Z_linkage, fcluster, dendro_link, metamuscle_list], f)

        # scatter plot in two main PCs
        if do_pca:
            pylab.figure()
            pylab.plot([i['PCA_param_vec'][0] for i in muscle_dics],
                       [i['PCA_param_vec'][1] for i in muscle_dics], '.k')
            for muscle_dic in muscle_dics:
                pylab.annotate(muscle_dic['sMuscle'],
                               muscle_dic['PCA_param_vec'][:2])
            pylab.xlabel('PCA #1')
            pylab.ylabel('PCA #2')
            pylab.title('PCA scatter of all muscles in 1 and 2')

            if more_plots:
                pylab.figure()
                pylab.plot([i['PCA_param_vec'][1] for i in muscle_dics],
                           [i['PCA_param_vec'][2] for i in muscle_dics], '.k')
                for muscle_dic in muscle_dics:
                    pylab.annotate(muscle_dic['sMuscle'],
                                   [muscle_dic['PCA_param_vec'][1],
                                    muscle_dic['PCA_param_vec'][2]])
                pylab.xlabel('PCA #2')
                pylab.ylabel('PCA #3')
                pylab.title('PCA scatter of all muscles in 2 and 3')

                pylab.figure()
                pylab.plot([i['PCA_param_vec'][0] for i in muscle_dics],
                           [i['PCA_param_vec'][2] for i in muscle_dics], '.k')
                for muscle_dic in muscle_dics:
                    pylab.annotate(muscle_dic['sMuscle'],
                                   [muscle_dic['PCA_param_vec'][0],
                                    muscle_dic['PCA_param_vec'][2]])
                pylab.xlabel('PCA #1')
                pylab.ylabel('PCA #3')
                pylab.title('PCA scatter of all muscles in 1 and 3')
            with open('msd_fig6_4.pickle', 'w') as f:
                pickle.dump(muscle_dics, f)

        # Structural difference
        if do_bsbf:
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
            pylab.xlim([0, max([max(bsbf_incluster_dist),
                                max(bsbf_outcluster_dist)])])
            pylab.xlabel('Distance between vectors')
            pylab.ylabel('Probability')
            pylab.legend()

            bsbf_ax_right = pylab.subplot(2, 2, 2)
            pylab.boxplot(bsbf_bs, vert=False, whis=[0, 100])
            pylab.axvline(x=0, color='m')
            pylab.axvline(x=bsbf_bs_w_interv[1], color='r', linestyle='--')

            pylab.subplot(2, 2, 4, sharex=bsbf_ax_right, sharey=bsbf_ax_bottom)
            pylab.hist(bsbf_bs, weights=density_weights(bsbf_bs),
                       bins=bins)
            pylab.axvline(x=0, color='m')
            pylab.axvline(x=bsbf_bs_w_interv[1], color='r', linestyle='--')
            pylab.xlabel('Distance between vectors in and out')
            pylab.ylabel('Probability')

            pylab.suptitle('Crossing same DOF clusters')

        # Functional difference
        if do_expbs:
            if more_plots:
                for i, expbs_cl_label in enumerate(expbs_cl_labels):
                    pylab.figure()
                    ax1 = pylab.subplot(2, 2, 1)
                    boxplot = []
                    boxplot_labels = []
                    for cl, lbl in izip(expbs_inouts[i], expbs_cl_labelsets[i]):
                        boxplot.append(cl['in'])
                        boxplot.append(cl['out'])
                        boxplot_labels += [str(lbl) + ' in', str(lbl) + ' out']

                    pylab.boxplot(boxplot, vert=False,
                                  labels=boxplot_labels,
                                  whis=[0, 100],
                                  )

                    pylab.subplot(2, 2, 3, sharex=ax1)
                    for h, l in izip(boxplot, boxplot_labels):
                        pylab.hist(h, weights=density_weights(h),
                                   bins=bins, label=l, alpha=0.5)
                    pylab.xlabel('Distance between vectors')
                    pylab.ylabel('Probability')
                    pylab.legend()

                    exbps_ax2 = pylab.subplot(2, 2, 2)
                    pylab.boxplot(expbs_bss[i], vert=False,
                                  labels=[str(j) for j in expbs_cl_labelsets[i]],
                                  whis=[0, 100])
                    pylab.axvline(x=0, color='m')

                    pylab.subplot(2, 2, 4, sharex=exbps_ax2)
                    for h, l in izip(expbs_bss[i], expbs_cl_labelsets[i]):
                        pylab.hist(h, weights=density_weights(h),
                                   bins=bins, label=str(l), alpha=0.5)
                    pylab.axvline(x=0, color='m')
                    pylab.xlabel('Difference in distance between vectors')
                    pylab.ylabel('Probability')
                    pylab.legend()

                    pylab.suptitle('Expected {} clusters'.format(
                        len(set(expbs_cl_label))))

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


def complexity(muscle_ids, plot, bins):
    metamuscle_list = meta.metamuscle_list()
    approx_cons_dirname = meta.APPROX_DIRNAME

    n = len(muscle_ids)
    muscle_dics = []
    for muscle_id in muscle_ids:
        muscle_dics.append(meta.search_dic(metamuscle_list, idMuscle=muscle_id))

    for i in range(len(muscle_dics)):
        with open(meta.apprfname(meta.APPROX_DIRNAME, muscle_dics[i]['sMuscle'])) as f:
            mom_approx, len_approx = pickle.load(f)
        dim = len(muscle_dics[i]['idDOFList'])
        er = gepapy.ersatz.Ersatz(dim, len_approx['config'])
        muscle_dics[i]['complexity'] = er.complexity(meta.DEFAULTS['rho'])*100
        muscle_dics[i]['abs_complexity'] = len(er)

    muscle_dics_resorted = sorted(muscle_dics, key=lambda j: j['complexity'])
    print 'List of complexities:'
    print '\n'.join(['{}\t dim {}: {}'.format(
        i['sMuscle'], len(i['idDOFList']), i['complexity'])
        for i in muscle_dics_resorted])

    if plot:
        y_subplots = 2
        # Histogram
        pylab.figure()
        pylab.subplot(y_subplots, 1, 1)
        l_bins = numpy.linspace(0, 100, bins/2)
        pylab.hist([i['complexity'] for i in muscle_dics], l_bins, color='k')
        pylab.xlim(0, 100)
        pylab.xlabel('Relative complexity of polynomials')
        pylab.ylabel('Number of muscles')

        pylab.subplot(y_subplots, 1, 2)
        lp_bins = numpy.linspace(0, 60, bins/2)
        pylab.hist([i['abs_complexity'] for i in muscle_dics], lp_bins, color='k')
        pylab.xlabel('Length of polynomial')
        pylab.ylabel('Number of muscles')
        pylab.xlim([0, 60])

        # Scatter plot
        pylab.figure()
        pylab.subplot(y_subplots, 1, 1)
        x = [len(i['idDOFList']) for i in muscle_dics]
        y = [i['complexity'] for i in muscle_dics]
        pylab.plot(x, y, '.k')
        grps = {}
        for md in muscle_dics:
            crd = (len(md['idDOFList']), md['complexity'])
            if crd not in grps.keys():
                grps[crd] = []
            grps[crd].append(md['sMuscle'])
        for key, item in grps.iteritems():
            pylab.annotate(','.join(item), key)

        comp_linreg = scipy.stats.linregress(x, y)
        comp_poly1d = lambda x: comp_linreg[0]*x+comp_linreg[1]
        print ('Complexity linreg: slope {}, intercept {}, r_value {},'
               ' r_sq_value {}, p_value {}, std_err {}').format(
            comp_linreg[0], comp_linreg[1], comp_linreg[2], comp_linreg[2]**2,
            comp_linreg[3], comp_linreg[4])
        pylab.plot([min(x), max(x)],
                   [comp_poly1d(min(x)), comp_poly1d(max(x))], 'k')
        pylab.ylim([0, 100])
        pylab.xlabel('Muscle dimensionality')
        pylab.ylabel('Normalized polynomial complexity')

        pylab.subplot(y_subplots, 1, 2)
        y = [i['abs_complexity'] for i in muscle_dics]
        pylab.plot(x, y, '.k')
        grps = {}
        for md in muscle_dics:
            crd = (len(md['idDOFList']), md['abs_complexity'])
            if crd not in grps.keys():
                grps[crd] = []
            grps[crd].append(md['sMuscle'])
        for key, item in grps.iteritems():
            pylab.annotate(','.join(item), key)

        acomp_linreg = scipy.stats.linregress(x, y)
        acomp_poly1d = lambda x: acomp_linreg[0]*x+acomp_linreg[1]
        print ('Absolute omplexity linreg: slope {}, intercept {}, r_value {},'
               ' r_sq_value {}, p_value {}, std_err {}').format(
            acomp_linreg[0], acomp_linreg[1], acomp_linreg[2], acomp_linreg[2]**2,
            acomp_linreg[3], acomp_linreg[4])
        pylab.plot([min(x), max(x)],
                   [acomp_poly1d(min(x)), acomp_poly1d(max(x))], 'k')
        pylab.ylim([0, 60])
        pylab.xlabel('Muscle dimensionality')
        pylab.ylabel('Polynomial length')

        f_out = {
            'bins': bins,
            'muscle_dics': muscle_dics,
            'comp_linreg': comp_linreg,
            'acomp_linreg': acomp_linreg,
        }
        with open('msd_fig5_2.pickle', 'w') as f:
            pickle.dump(f_out, f)


def main(muscle_ids,
         sef=False, plot=False, al=False,
         show=False, bins=10,
         dendro_link='average',
         do_complexity=False,
         verbosity=1,
         do_pca=False, do_dendro=False,
         do_bsbf=False, do_expbs=False, num_expbs=7, plot_basic=False,
         more_plots=False, use_ddpv=False):
    muscle_ids = meta.get_mids_from_args(muscle_ids)

    if sef:
        internal_consistency(muscle_ids, plot, bins)

    if al:
        alltoall(muscle_ids, plot, dendro_link, bins, verbosity,
                 do_pca, do_dendro,
                 do_bsbf, do_expbs, num_expbs,  plot_basic,
                 more_plots, use_ddpv)

    if do_complexity:
        complexity(muscle_ids, plot, bins)

    if show:
        pylab.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Compare structures of muscles.')

    parser.add_argument('-s', '--self',
                        dest='s',
                        action='store_true',
                        help='Compare structure before and after internal '
                        'consistency enforce.')
    parser.add_argument('--bins',
                        dest='bins', type=int,
                        action='store', default=10,
                        help='Number of bins for histograms')
    parser.add_argument('-a', '--all',
                        dest='a',
                        action='store_true',
                        help='Compare structure of all muscles.')
    parser.add_argument('--dendro_link',
                        dest='dendro_link', type=str,
                        action='store', default='average',
                        help='Modifier for -a. Linkage method for dendrogram')
    parser.add_argument('-p', '--plot',
                        dest='p',
                        action='store_true',
                        help='Make plots.')
    parser.add_argument('-m', '--muscle-ids',
                        dest='m', default=[],
                        action='store', nargs='+',
                        help='muscle ids for analysis. You can either specify'
                        ' an ID, or a preset code word: `all\' or `hand\' '
                        'or `special\'.')
    parser.add_argument('-c', '--complexity',
                        dest='c',
                        action='store_true',
                        help='Make a histopgram of polynomial complexity.')
    parser.add_argument('--verbosity',
                        dest='verbosity',
                        type=int, default=1,
                        help='Level of verbosity.')
    parser.add_argument('--do_pca',
                        dest='do_pca',
                        action='store_true',
                        help='Do a PCA clustering analysis.')
    parser.add_argument('--do_dendro',
                        dest='do_dendro',
                        action='store_true',
                        help='Do a dendrogram clustering analysis.')
    parser.add_argument('--do_bsbf',
                        dest='do_bsbf',
                        action='store_true',
                        help='Bruteforce bootstrapping to test if function'
                        ' is within the structure.')
    parser.add_argument('--do_expbs',
                        dest='do_expbs',
                        action='store_true',
                        help='Do bootstrapping analysis for expected '
                        'functional clusters.')
    parser.add_argument('--num_expbs',
                        dest='num_expbs', type=int,
                        action='store', default=7,
                        help='Set number of size of expected cluster. If 1'
                        ' will go through all of them.')
    parser.add_argument('--plot_basic',
                        dest='plot_basic',
                        action='store_true',
                        help='Plots dendrogram and heatmap.')
    parser.add_argument('--more_plots',
                        dest='more_plots',
                        action='store_true',
                        help='More plots.')
    parser.add_argument('--use_ddpv',
                        dest='use_ddpv',
                        action='store_true',
                        help='Use DOF-DEPENDENT polynomial vectors.')

    args = parser.parse_args()

    main(args.m, sef=args.s, plot=args.p, al=args.a,
         show=args.p, bins=args.bins,
         dendro_link=args.dendro_link,
         do_complexity=args.c,
         verbosity=args.verbosity,
         do_pca=args.do_pca,
         do_dendro=args.do_dendro,
         do_bsbf=args.do_bsbf, do_expbs=args.do_expbs, num_expbs=args.num_expbs,
         plot_basic=args.plot_basic, more_plots=args.more_plots,
         use_ddpv=args.use_ddpv)
