#!python2
import os
import sys
import pickle
import argparse

import pylab

import gepapy
import meta


def get_labels(dim, dump):
    labels = []
    prev_ers = gepapy.ersatz.Ersatz(dim, {})
    for d in dump:
        best_ers = gepapy.ersatz.Ersatz(dim, d['config'])
        diff_ers = best_ers - prev_ers
        lblparts = []
        for key in diff_ers.sorted_keys():
            if key == 'a':
                lblparts.append('a')
            elif key == 'M':
                for el in diff_ers.config_dic[key]:
                    lblparts.append('K' + ','.join(str(j+1) for j in el))
        labels.append(';'.join(lblparts))
        prev_ers = best_ers

    return labels


def analyze_one(muscle_dic, filename, dof_labels):
    with open(filename) as f:
        mom_arm_approxs, length_approx, dump, dump_m = pickle.load(f)
    dim = len(mom_arm_approxs)

    labels = get_labels(dim, dump)
    print labels

    fig = pylab.figure()
    ax = pylab.gca()
    pylab.semilogy(range(len(labels)-1), [d['metric'] for d in dump[:-1]], 'k')

    for j, d in enumerate(dump_m):
        x = range(len(labels))
        y = [di['metric'] for di in d]
        i = len(labels)-1
        while y[i] == y[i-1]:
            del y[i]
            del x[i]
            i -= 1
        print 'Moment arm around dof {} has final precision of {}'.format(
            dof_labels[j], y[-1])
        pylab.semilogy(x, y)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels)
    pylab.xlim([-0.5, len(labels)-1.5])
    pylab.ylim([3e-4, 10])
    pylab.xlabel('Iteration, Added terms')
    pylab.ylabel('Precision of fit, au')
    pylab.legend(['l'] + dof_labels)
    pylab.title(muscle_dic['sMuscle'])


def main(muscle_ids, show):
    muscle_ids = meta.get_mids_from_args(muscle_ids)
    metamuscle = meta.metamuscle_list()
    metadof = meta.metadof_list()

    for mid in muscle_ids:
        muscle_dic = meta.search_dic(metamuscle, idMuscle=mid)
        filename = meta.apprfname(meta.APPROX_PROGRESS_DIRNAME,
                                  muscle_dic['sMuscle'])
        dof_labels = [meta.search_dic(metadof, idDOF=i)['sDOF']
                      for i in muscle_dic['idDOFList']]
        analyze_one(muscle_dic, filename, dof_labels)

    if show:
        pylab.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Displays progress in muscle optimization.')

    parser.add_argument('-s', '--show',
                        dest='show',
                        action='store_false',
                        help='Run pylab show.')

    parser.add_argument('-m', '--muscle_ids',
                        dest='m', default=[],
                        action='store', nargs='+',
                        help='muscle ids for analysis. You can either specify'
                        ' an ID, or a preset code word: `all\' or `hand\' '
                        'or `special\'.')

    args = parser.parse_args()

    main(args.m, args.show)
