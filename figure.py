'''
Cindy Olvera Perez 
May 16, 2020
figure.py
plotting functions
'''
from astropy.table import Column
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def label(n):
    ''' returns a list of column names for velocities'''
    col = [] # column name list
    for i in range(n):
        iname = str("initial velocity " + str(i+1)) # initial velocity column name
        fname =  str("final velocity " + str(i+1))  # and final velocity column name
        col.append(iname)
        col.append(fname)
    return col

def sLabels(n):
    ''' returns a list of column names for solutions'''
    col = [] # column name list
    for i in range(n):
        name = str("solution " + str(i+1)) # initial velocity column name
        col.append(name)
    return col

def Vcolumn(tableName, idata, fdata, columnNamei, columnNamef):
    '''
    adds velocity columns to astropy table
    '''
    tableName.add_column(Column(data=idata, name=columnNamei))
    tableName.add_column(Column(data=fdata, name=columnNamef))

def Scolumn(tableName, data, columnName):
    '''
    adds solution columns to astropy table
    '''
    tableName.add_column(Column(data=data, name=columnName))
    
def histogram(v_initial, v_final,
                    title1, xlable, name,
                    binw=0.2, binStart=5.4, binStop=15.7,
                    label_vi = "initial velocity", label_vf = "final velocity", 
                    ledgend_loc='upper right'):
    '''
    histogram code
    '''
    # setting plotting options
    grid_style =   {     'alpha' : '0.75',
                     'linestyle' : ':' }
    legend_style = {  'fontsize' : '10' }
    font_syle =    {      'size' : '14' }
    mpl.rc(  'font', **font_syle)
    mpl.rc(  'grid', **grid_style)
    mpl.rc('legend', **legend_style)
    
    bins = np.arange(binStart, binStop, binw)
    fig, ax=plt.subplots(1,1, figsize=(8,5))

    ax.hist(v_initial, bins, label= label_vi, color="red", density=1.0, linewidth=1.3, alpha=0.5)
    ax.hist(v_final, bins, label= label_vf, color="black", density=1.0, linewidth=1.3, alpha=0.5)
    
    ax.grid()
    ax.legend(bbox_to_anchor=(1.05, 1), loc=ledgend_loc, shadow=True)
    ax.set_facecolor('whitesmoke')
    ax.patch.set_alpha(0.1)
    ax.set_title(title1)
    ax.set_xlabel(xlable)
    fig.tight_layout()
    fig.savefig(name)

def histograms(textfile,
                    binw=0.2, binStart=5.4, binStop=15.7,
                    label_vi = "initial velocity", label_vf = "final velocity", velHist='', vel='',
                    ledgend_loc='upper right'):
    '''
    reads a data textfile and creates histograms
    uses function galaxyHistogram from h_plotting.py
    '''
    plt.close('all')
    galaxyData = ascii.read(textfile, format='commented_header')
    
#    for i in range(len(galaxyData.columns)/2):
#        name = str("figure " + str(i))
    i=1
    for colName in galaxyData.columns:
        name = str("figure " + str(i))
        if i%2!=0: #checking for odd number columns
            ivel=galaxyData[colName]
        else:      #checking for even number columns
            fvel=galaxyData[colName]
            histogram(ivel, fvel, velHist, vel, name, binw, binStart, binStop,
                                 label_vi, label_vf,
                                 ledgend_loc)
        i+=1

