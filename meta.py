#!/usr/bin/env python2
import os
import csv
import time
from scipy.io import loadmat


METADOF_FNAME = 'arm_metaDOF.csv'
METAMUSCLE_FNAME = 'arm_metaMuscle.csv'

MUSCLE_IDS = {
    'hand': range(1, 34),
    'fullhand': range(1, 34),
    'all': range(1, 34),
    'hand_nothumb': range(1, 26),
    'hand_noepl': range(1, 26) + range(27, 34)
}

def get_mids_from_args(arg):
    if len(arg) == 0:
        return []
    elif arg[0] in MUSCLE_IDS.keys():
        return MUSCLE_IDS[arg[0]]
    else:
        return [int(i) for i in arg]

DEFAULTS = {
    'epow': 'mixed',
    'lytype': 'l3',
    'mytype': 'm1',
    'rho': 5,
    'pps': 6,
    'mid': None,
    'dimid': None,
    'overwrite': False,
    'dump': 0,
}

METRIC_TYPE_DIC = {
    'm1': 'rms_n2am',
    'm2': 'rms',
    'l1': 'nrms',
    'l2': 'rms_n2mv',
    'l3': 'rms_n2rom'
}

DOF_SIMILARITIES = range(19)
# 3rd finger         = 4rd finger           = 5th finger           = 2nd finger
DOF_SIMILARITIES[10] = DOF_SIMILARITIES[13] = DOF_SIMILARITIES[16] = 7
DOF_SIMILARITIES[11] = DOF_SIMILARITIES[14] = DOF_SIMILARITIES[17] = 8
DOF_SIMILARITIES[12] = DOF_SIMILARITIES[15] = DOF_SIMILARITIES[18] = 9

DATA_DIRNAME = 'DATA 9 Point'
APPROX_DIRNAME = 'approximations'
APPROX_NOT_INTCONS_DIRNAME = 'approximations_not_intcons'
APPROX_PROGRESS_DIRNAME = 'approximations_progress'
FUNCTION_DIRNAME = 'out'

datafname = lambda data_dirname, muscle: os.path.join(
    data_dirname, muscle + 'MomentArm.mat')
apprfname = lambda approx_dirname, muscle: os.path.join(
    approx_dirname, muscle + 'MomentArm_approxs.pcl')

MOMENTARM_OUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'momarmfunc_c.c')
MUSCLELENGTH_OUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'lenfunc_c.c')
MOMENTARM_NOT_INTCONS_OUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'momarmfunc_not_intcons_c.c')
MUSCLELENGTH_NOT_INTCONS_OUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'lenfunc_not_intcons_c.c')

MOMENTARM_MOUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'momarmfunc.m')
MUSCLELENGTH_MOUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'lenfunc.m')
MOMENTARM_NOT_INTCONS_MOUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'momarmfunc_not_intcons.m')
MUSCLELENGTH_NOT_INTCONS_MOUT_FUNCTION_NAME = os.path.join(FUNCTION_DIRNAME, 'lenfunc_not_intcons.m')

QUALITY_THRESHOLD = None
QUALITY_SDEF_THRESHOLD = 1e-4
QUALITY_CHECK_THRESHOLD = 0.05

def time_stamp(gmtime=None, sep=':'):
    """Formats time string."""
    if gmtime is None:
        gmtime = time.gmtime()
    timestring = '%Y-%m-%d_%H'+sep+'%M'+sep+'%S'
    return time.strftime(timestring, gmtime)


#------------------------------------------------------------------------------
#---------------------- Meta
#------------------------------------------------------------------------------
def metadof_list(metadof_fname=None):
    if metadof_fname is None:
        metadof_fname = METADOF_FNAME
    dof_list = []
    with open(metadof_fname) as f:
        fd = csv.DictReader(f)
        for line in fd:
            dof_list.append({
                'idDOF': int(line['idDOF']),
                'sDOF': line['sDOF'],
                'nRangeMin': float(line['nRangeMin']),
                'nRangeMax': float(line['nRangeMax']),
                'nIC': float(line['nIC']),
                })

    return dof_list


def metamuscle_list(metamuscle_fname=None):
    if metamuscle_fname is None:
        metamuscle_fname = METAMUSCLE_FNAME
    muscle_list = []
    with open(metamuscle_fname) as f:
        fd = csv.DictReader(f)
        for line in fd:
            muscle_list.append({
                'idMuscle': int(line['idMuscle']),
                'sMuscle': line['sMuscle'],
                'idDOFList': [int(i) for i in line['idDOFList'].split()],
                })
    return muscle_list


def search_dic(dic, **argv):
    """Searches a list of dics for what would satisfy argv

    Example:
        search_dic(muscle_list, idMuscle=34)

    Returns a full dictionary from the list
    """
    if len(argv) == 0:
        raise ValueError('Need to provide args what to search for.')
    what = argv.items()[0]
    for i in xrange(len(dic)):
        if dic[i][what[0]] == what[1]:
            return dic[i]
    raise ValueError('Element not found in list.')
