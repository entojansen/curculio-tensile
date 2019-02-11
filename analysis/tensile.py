# -*- coding: utf-8 -*-
'''
Created on Tue Aug  1 13:59:50 2017

@author: Flashcock
'''

import os
from scipy.stats import linregress as linreg
from itertools import takewhile
from numpy import array, poly1d, sqrt, absolute, arcsin, pi
from scipy.optimize import curve_fit, minimize_scalar
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from numpy import log as ln

class specimen:
    '''
    Args should be formatted as follows:
        asuid (str or int), species (str), length (str), exo (str), endo (str),
        total (str), N (list [float]), mm (list [float]), Sec (list [float]).
    These are converted to their respective types below and used to calculate
        true and engineering stress and strain, young's modulus, and UTS.
    '''
    def __init__(self, asuid, species, length, exo, endo, total, N, mm, Sec):
        self.__asuid__ = asuid
        self.__species__ = species.lower()

        self.__length__ = float(length)
        self.__exo__ = float(exo)
        self.__endo__ = float(endo)
        self.__total__ = float(total)

        self.__N__ = array(N)
        self.__mm__ = array(mm) - min(mm)
        self.__Sec__ = array(Sec) - Sec[0]

        self.__sigmaN__ = self.__N__ / self.__total__ #P/A0
        self.__epsilonN__ = (self.__mm__ / self.__length__) #dL/L0
        self.__sigmaT__ = self.__sigmaN__ * (1 + self.__epsilonN__) #s(1+e)
        self.__epsilonT__ = ln(1 + self.__epsilonN__) #ln(1+e)

        self.__fmax__ = max(N)
        self.__uts__ = max(self.__sigmaN__) #evaluated from engineering stress

        self.__dmax__ = max(self.__mm__)
        self.__ufs__ = max(self.__epsilonN__) #evaluated from engineering strain


        self.__Eindex__ = []
        lower = self.__ufs__ * 0
        upper = self.__ufs__ * 0.33
        self.__Eindex__.append(len([x for x in takewhile(lambda x: x[1] <= lower, enumerate(self.__epsilonN__))]))
        self.__Eindex__.append(len([y for y in takewhile(lambda y: y[1] <= upper, enumerate(self.__epsilonN__))]))
        xvals = self.__epsilonN__[self.__Eindex__[0]:self.__Eindex__[1]]
        yvals = self.__sigmaN__[self.__Eindex__[0]:self.__Eindex__[1]]
        self.__E__ = linreg(xvals, yvals) #engineering stress and strain via Hooke's Law

        self.__Etanl__ = self.__E__[0]
        self.__Etanh__ = 0

        t_upper = 0
        t_lower = 0

        for i in range(len(self.__epsilonN__)-1):
            ub = self.__epsilonN__[i+1]
            lb = self.__epsilonN__[i]
            uh = self.__sigmaN__[i+1]
            lh = self.__sigmaN__[i]
            t_upper += (ub - lb) * uh
            t_lower += (ub - lb) * lh

        self.__U__ = (t_upper + t_lower) / 2

    def asuid(self):
        return self.__asuid__

    def species(self):
        return self.__species__

    def length(self):
        return self.__length__

    def exo(self):
        return self.__exo__

    def endo(self):
        return self.__endo__

    def total(self):
        return self.__total__

    def N(self):
        return self.__N__

    def mm(self):
        return self.__mm__

    def Sec(self):
        return self.__Sec__

    def sigmaN(self):
        return self.__sigmaN__

    def epsilonN(self):
        return self.__epsilonN__

    def sigmaT(self):
        return self.__sigmaT__

    def epsilonT(self):
        return self.__epsilonT__

    def Eindex(self, index):
        return self.__Eindex__[index]

    def recalcE(self, lower=0.00, upper=0.015, graphical=False, normal=False, limits='E_limits.txt'):
        self.__Eindex__ = []
        if normal:
            eps = self.__epsilonT__
            sig = self.__sigmaT__
        else:
            eps = self.__epsilonN__
            sig = self.__sigmaN__
        if graphical:
            file = open(limits, 'r')
            for line in file:
                asuid, species, lims = line[:-1].split('|', 2)
                if asuid == self.__asuid__:
                    x, y = lims.split(',')
                    lowerg = float(x[1:])
                    upperg = float(y[:-1])
                else:
                    pass
            self.__Eindex__.append(len([x for x in takewhile(lambda x: x[1] <= lowerg, enumerate(self.__Sec__))]))
            self.__Eindex__.append(len([y for y in takewhile(lambda y: y[1] <= upperg, enumerate(self.__Sec__))]))

        else:
            self.__Eindex__=[0, 1]
            low = False
            high = False
            while (not low) and (not high):
                for i, c in enumerate(self.__epsilonN__):
                    if (c >=lower) and (not low):
                        self.__Eindex__[0] = i
                        low = True
                    elif (c >= upper) and (not high):
                        self.__Eindex__[1] = i
                        high = True
                    else:
                        pass
        xvals = eps[self.__Eindex__[0] : self.__Eindex__[1]]
        yvals = sig[self.__Eindex__[0] : self.__Eindex__[1]]

        self.__E__ = linreg(xvals, yvals) #engineering stress and strain via Hooke's Law

    def E(self, index = 0):
        return self.__E__[index]

    def Etanh(self):
        return self.__Etanh__

    def Etanl(self):
        return self.__Etanl__

    def fmax(self):
        return self.__fmax__

    def uts(self):
        return self.__uts__

    def dmax(self):
        return self.__dmax__

    def ufs(self):
        return self.__ufs__

    def U(self):
        return self.__U__

def as_si(x, ndp):
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))

def get_files(path='.\data'):
    '''
    Retrieves filenames in directory.
    '''
    filelist = os.listdir(path)
    return filelist

def spec_data(files, dims='imagedata.csv'):
    '''
    Returns a set of specimens.
    '''
    spec_set = set()

    file1 = open(dims, 'r')
    for line in file1:
        if line[0] == 'A':
            pass
        else:
            asuid, species, length, exo, endo, total = line[:-1].split(',', 5)
            N = []
            mm = []
            Sec = []
            normal_order = True
            for name in files:
                if name == '{}-{}-edit.dat'.format(asuid, species):
                    file2=open('.\data\\' + name,'r')
                    for line in file2:
                        if line[0] == 'S':
                            normal_order = False
                        elif line[0] in '1234567890':
                            x, y, z = line[:-1].split('\t', 2)
                            if normal_order:
                                N.append(float(x))
                                mm.append(float(y))
                                Sec.append(float(z))
                            else:
                                N.append(float(z))
                                mm.append(float(y))
                                Sec.append(float(x))
                        else:
                            pass
                    file2.close()
                    identifier, etc= name[:-9].split('-')
                    spec = specimen(asuid, species, length, exo, endo, total, N, mm, Sec)
                    spec_set.add(spec)
                else:
                    pass
    file1.close()
    return spec_set

def stress_strain(spec_set, save=True, i=0, j=-1):
    for spec in spec_set:
        ssplot, (ax1, ax3) = plt.subplots(2, sharex=False)

        ax1.plot(spec.Sec()[i:j], spec.N()[i:j], color='b')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Axial Force (N)', color='b')
        ax1.tick_params('y', colors='b')
        ax1.set_title('Tensile Test Raw Data - ASUHIC_{} - Curculio {}'.format(spec.asuid(), spec.species()))

        ax2 = ax1.twinx()
        ax2.plot(spec.Sec()[i:j], spec.mm()[i:j], color='r')
        ax2.set_ylabel('Axial Displacement (mm)', color='r')
        ax2.tick_params('y', colors='r')

        ax3.plot(spec.epsilonT()[i:j], spec.sigmaT()[i:j], color = 'gray')
        ax3.plot(spec.epsilonN()[i:j], spec.sigmaN()[i:j], color = 'black')
        ax3.set_ylabel('Stress (MPa)')
        ax3.set_xlabel('Strain (ln(L/L0) or dL/L0)')

        print(spec.asuid()+' '+spec.species())
        new_l = spec.ufs() * 0.00
        new_u = spec.ufs() * 0.33
        spec.recalcE(lower = new_l, upper = new_u, graphical=False, normal=False)
        ax3.plot(spec.epsilonN()[spec.Eindex(0):spec.Eindex(1)],
                 spec.E()*spec.epsilonN()[spec.Eindex(0):spec.Eindex(1)] + spec.E(1),
                 color='r')

        new_l = spec.ufs() * 0.67
        new_u = spec.ufs() * 1.00
        spec.recalcE(lower = new_l, upper = new_u, graphical=False, normal=False)
        ax3.plot(spec.epsilonN()[spec.Eindex(0):spec.Eindex(1)],
                 spec.E()*spec.epsilonN()[spec.Eindex(0):spec.Eindex(1)] + spec.E(1),
                 color='r')

        true = mpatches.Patch(color='gray', label='True')
        eng = mpatches.Patch(color='black', label='Nominal')
        E1 = mpatches.Patch(color='r', label='E1 = {}MPa'.format(int(round(spec.Etanl()))))
        E2 = mpatches.Patch(color='r', label='E2 = {}MPa'.format(int(round(spec.E()))))
        ax3.legend(handles=[true, eng, E1, E2], loc=2, fontsize = 6)
        #ax3.set_title('Stress vs. Nominal Strain')
        #ax3.set_xlim(0, 1)
        #ax3.set_ylim(0, 300)
        if save:
            ssplot.savefig('plots/rawdata/{}-{}.pdf'.format(spec.asuid(), spec.species()))
            #plt.close()
        else:
            plt.show()
            #plt.close()
            print('Plot not saved!!')

def xyplot(spec_set, save = True):
    spl = []
    regx = []
    regy = []

    fig = plt.figure()
    ax = fig.add_subplot(111)

    #plt.ylim(0,0.04)
    #plt.xlim(0,8)
    xlab = 'test (units)'
    ylab = 'test (units)'
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title('{}_vs_{}.pdf'.format(xlab, ylab))

    car = mpatches.Patch(color='r', label='caryae')
    pro = mpatches.Patch(color='lightcoral', label='proboscideus')
    uni = mpatches.Patch(color='b', label='uniformis')
    hum = mpatches.Patch(color='skyblue', label='humeralis')
    sul = mpatches.Patch(color='dodgerblue', label='sulcatulus')
    vic = mpatches.Patch([], [], color='navy', label='victoriensis')
    ax.legend(handles=[car, pro, hum, sul, uni, vic], loc=4, fontsize = 7)

    for spec in spec_set:
        new_l = spec.ufs() * 0.67
        new_u = spec.ufs() * 1.00
        spec.recalcE(lower = new_l, upper = new_u, graphical=False, normal=False)
        x = ln(spec.endo() / spec.exo())
        y = ln(spec.ufs())
        #spec.uts() / spec.ufs()

        if spec.species() == 'caryae':
            ax.plot(x, y, color='r', marker ='o')
        elif spec.species() == 'proboscideus':
            ax.plot(x, y, color='lightcoral', marker ='o')
        elif spec.species() == 'uniformis':
            ax.plot(x, y, color='b', marker ='o')
        elif spec.species() == 'humeralis':
            ax.plot(x, y, color='skyblue', marker ='o')
        elif spec.species() == 'sulcatulus':
            ax.plot(x, y, color='dodgerblue', marker ='o')
        elif spec.species() == 'victoriensis':
            ax.plot(x, y, color='navy', marker ='o')
        else:
            pass

        regx.append(x)
        regy.append(y)
        spl.append(str(x)+' '+str(y)+' '+spec.asuid()+' '+spec.species())

    x = array(regx)
    y = array(regy)
    regline = linreg(x, y)
    ax.plot(x, regline[0]*x+regline[1], color='black')
    ax.text(0.03,
            3,
            r'$y={:.3f}x+{:.3f}$'.format(regline[0], regline[1])
            +'\n'
            +r'$R={:.3f}$'.format(regline[2])
            +'\n'
            #+r'$p={:.3f}$'.format(regline[3]),
            +r'$p={0:s}$'.format(as_si(regline[3], 2)),
            fontsize=10)

    for  i in sorted(spl):
        print(i)

    plt.show()

    if save:
        fig.savefig('plots/xyplots/{}_vs_{}.pdf'.format(xlab, ylab))
        #plt.close()
    else:
        #plt.close()
        print('Plot not saved!!')

    print('m = {}'.format(regline[0]),
          'b = {}'.format(regline[1]),
          'r = {}'.format(regline[2]),
          'p = {}'.format(regline[3]),
          'std_err = {}'.format(regline[4]),
          sep='\n')

def xyzplot(spec_set):
    spl = []
    regx = []
    regy = []
    regz = []
    xyzfig = plt.figure()
    ax = Axes3D(xyzfig)
    for spec in spec_set:
        spec.recalcE(graphical=False, normal=False)
        x = (spec.length())
        y = (spec.U())
        z = (spec.exo() / spec.endo())
        if spec.species() == 'caryae':
            ax.scatter(x, y, z, c='r')
        elif spec.species() == 'proboscideus':
            ax.scatter(x, y, z, c='lightcoral')
        elif spec.species() == 'uniformis':
            ax.scatter(x, y, z, c='b')
        elif spec.species() == 'humeralis':
            ax.scatter(x, y, z, c='navy')
        elif spec.species() == 'sulcatulus':
            ax.scatter(x, y, z, c='dodgerblue')
        elif spec.species() == 'victoriensis':
            ax.scatter(x, y, z, c='skyblue')
        else:
            pass
        regx.append(x)
        regy.append(y)
        regz.append(z)
        spl.append(str(x)+' '+str(y)+' '+str(z)+' '+spec.asuid()+' '+spec.species())
    x = array(regx)
    y = array(regy)
    z = array(regz)
    ax.set_xlabel('Length')
    ax.set_ylabel('Strain')
    ax.set_zlabel('Stress')
    #regline = linreg(x, y, z)
    #plt.plot(x, regline[0]*x+regline[1], color='black')
    #plt.ylim(0,0.04)
    #plt.xlim(0,8)
    for  i in sorted(spl):
        print(i)
    plt.show()
    #print('m = {}'.format(regline[0]), 'b = {}'.format(regline[1]), 'r = {}'.format(regline[2]), 'p = {}'.format(regline[3]), 'std_err = {}'.format(regline[4]), sep='\n')

def write_data(specimens, filename = 'tensile.csv'):
    outfile=open(filename,'w')
    print('asuid,species,specimen,length,exo,endo,total,fmax,uts,dmax,ufs,Etanl,Etanh,Esec,Ut', file=outfile)
    humeralis = 0
    caryae = 0
    proboscideus = 0
    sulcatulus = 0
    uniformis = 0
    victoriensis = 0
    current_sp = 0
    for spec in specimens:
        if spec.species() == 'humeralis':
            humeralis += 1
            current_sp = humeralis
        elif spec.species() == 'caryae':
            caryae += 1
            current_sp = caryae
        elif spec.species() == 'proboscideus':
            proboscideus += 1
            current_sp = proboscideus
        elif spec.species() == 'sulcatulus':
            sulcatulus += 1
            current_sp = sulcatulus
        elif spec.species() == 'uniformis':
            uniformis += 1
            current_sp = uniformis
        elif spec.species() == 'victoriensis':
            victoriensis += 1
            current_sp = victoriensis
        else:
            print('unrecognized species')
        print(spec.asuid(),
              spec.species(),
              spec.species()+str(current_sp),
              spec.length(),
              spec.exo(),
              spec.endo(),
              spec.total(),
              spec.fmax(),
              spec.uts(),
              spec.dmax(),
              spec.ufs(),
              spec.Etanl(),
              spec.E(),
              spec.uts() / spec.ufs(),
              spec.U(),
              sep=',',
              file=outfile)
    outfile.close()

def main():
    specimens = spec_data(get_files())
    #stress_strain(specimens, save=True)
    xyplot(specimens, save=False)
    #xyzplot(specimens)
    write_data(specimens)
    return specimens

testlib = main()