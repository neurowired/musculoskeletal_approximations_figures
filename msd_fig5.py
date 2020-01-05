#!python2
# >py -2 struct_diff.py -scp -m all --bins 40
import pickle
import numpy
import pylab


with open('msd_fig5_1.pickle', 'r') as f:
    [muscle_ids, sis, bins] = pickle.load(f)

pylab.figure()
bins = numpy.linspace(0, 100, 26)
pylab.hist(sis, bins=bins, color='k')
pylab.xlim([0, 100])
pylab.xlabel('Si between muscles before and after internal '
             'consistency, %')
pylab.ylabel('Number of muscles')
pylab.title('Mids len: {}'.format(len(muscle_ids)))

with open('msd_fig5_2.pickle', 'r') as f:
    f_in = pickle.load(f)

bins = f_in['bins']
muscle_dics = f_in['muscle_dics']
comp_linreg = f_in['comp_linreg']
acomp_linreg = f_in['acomp_linreg']

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

pylab.show()
