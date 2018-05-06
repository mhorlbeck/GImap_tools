import numpy as np
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
from scipy import stats

####	Analysis pipeline	####

#version that only filters based on cycledCol as many sgRNAs have median rep of 0 at the end
def calcLog2e_cycledonly(t0Col, cycledCol, doublesTable, filterThreshold = 1.0, pseudocount = 1.0, doublingDifferences = 1.0):
    meanCounts = pd.concat((cycledCol.groupby(doublesTable['name_a']).agg(np.median),cycledCol.groupby(doublesTable['name_b']).agg(np.median)),axis=1, keys=['a','b'])
    
    sgsToFilter = set(meanCounts.loc[meanCounts.loc[:,'b'] < filterThreshold].index).union(set(meanCounts.loc[meanCounts.loc[:,'a'] < filterThreshold].index))
    doublesTable_filt = doublesTable.loc[doublesTable.apply(lambda row: row['name_a'] not in sgsToFilter and row['name_b'] not in sgsToFilter, axis=1)]
    print str(len(doublesTable_filt)) + ' pairs of ' + str(len(doublesTable)) + ' passing filter'
    
    countsRatio = (t0Col.loc[doublesTable_filt.index] + pseudocount).sum()*1.0/(cycledCol.loc[doublesTable_filt.index] + pseudocount).sum()
    log2es = np.log2((cycledCol.loc[doublesTable_filt.index] + pseudocount)/(t0Col.loc[doublesTable_filt.index] + pseudocount)/countsRatio)

    doubleNegatives = doublesTable.apply(lambda row: row['gene_a'] == 'negative' and row['gene_b'] == 'negative', axis=1)

    log2es -= log2es.loc[doubleNegatives].median()

    log2es /= doublingDifferences
    
    return log2es

#for a specified variable position and sgRNA, get single phenotypes, double phenotypes, and optionally single phenotype std dev.
def getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition, returnXerr=False):
    if not returnXerr:
        return singlePhenotypes[variablePosition+'.mean'], \
            phenotypeMatrix.loc[sgRNA,:] if fixedPosition == 'a' else phenotypeMatrix.loc[:,sgRNA], \
            singlePhenotypes.loc[sgRNA, fixedPosition +'.mean']
    else:
        return singlePhenotypes[variablePosition+'.mean'], \
            phenotypeMatrix.loc[sgRNA,:] if fixedPosition == 'a' else phenotypeMatrix.loc[:,sgRNA], \
            singlePhenotypes.loc[sgRNA, fixedPosition +'.mean'], singlePhenotypes[variablePosition+'.std']
        
#convert phenotypes into square matrix
def generatePhenotypeMatrix(phenotypes):
    numSingles = int(np.sqrt(len(phenotypes)))
    phenotypeMatrix = np.zeros((numSingles,numSingles))
    singlesTable = []
    for i, (sgPair, counts) in enumerate(phenotypes.sort_index().iteritems()):
        phenotypeMatrix[i/numSingles, i%numSingles] = counts
        if i%numSingles == 0:
            singlesTable.append(sgPair.split('++')[0])

    phenotypeMatrix = pd.DataFrame(phenotypeMatrix, index=singlesTable, columns=singlesTable)
    singlesTable = pd.DataFrame([s.split('_')[0] for s in singlesTable], index=singlesTable, columns=['gene'])
    
    singlePhenotypes = pd.concat((phenotypeMatrix.loc[singlesTable['gene'] == 'negative',:].apply(np.nanmean, axis=0), 
                                  phenotypeMatrix.loc[singlesTable['gene'] == 'negative',:].apply(np.nanstd, axis=0), 
                                  phenotypeMatrix.loc[:, singlesTable['gene'] == 'negative'].apply(np.nanmean, axis=1),
                                 phenotypeMatrix.loc[:, singlesTable['gene'] == 'negative'].apply(np.nanstd, axis=1)), 
                                 axis=1, keys=['b.mean','b.std','a.mean','a.std'])
    
    return phenotypeMatrix, singlesTable, singlePhenotypes

def abbaAveragePhenotypes(phenotypeMatrix, singlesTable):
	phenotypeMatrix_abba = (phenotypeMatrix + phenotypeMatrix.T) / 2

	singlePhenotypes_abba = pd.concat((phenotypeMatrix_abba.loc[singlesTable['gene'] == 'negative',:].apply(np.nanmean, axis=0), 
                                  phenotypeMatrix_abba.loc[singlesTable['gene'] == 'negative',:].apply(np.nanstd, axis=0), 
                                  phenotypeMatrix_abba.loc[:, singlesTable['gene'] == 'negative'].apply(np.nanmean, axis=1),
                                  phenotypeMatrix_abba.loc[:, singlesTable['gene'] == 'negative'].apply(np.nanstd, axis=1)), 
                                 axis=1, keys=['b.mean','b.std','a.mean','a.std'])

	return phenotypeMatrix_abba, singlePhenotypes_abba


#calculate epistasis interactions, optionally z-standardizing based on negative controls
def calculateInteractions(phenotypeMatrix, singlePhenotypes, singlesTable, fitFunction, zstandardize=True):
    emap1 = pd.DataFrame(np.zeros(phenotypeMatrix.shape), index=phenotypeMatrix.index, columns=phenotypeMatrix.columns)
    variablePosition, fixedPosition = 'a','b'
    for i, sgRNA in enumerate(phenotypeMatrix.index):
        xdata, ydata, bdata = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition)
        
        fit = fitFunction(xdata, ydata, bdata)
        epistasis = ydata - fit(xdata)

        if zstandardize:
	        emap1.loc[sgRNA,:] = epistasis / epistasis.loc[singlesTable['gene'] == 'negative'].std()
       	else:
	        emap1.loc[sgRNA,:] = epistasis 

    emap2 = pd.DataFrame(np.zeros(phenotypeMatrix.shape), index=phenotypeMatrix.index, columns=phenotypeMatrix.columns)
    variablePosition, fixedPosition = 'b','a'
    for i, sgRNA in enumerate(phenotypeMatrix.index):
        xdata, ydata, bdata = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition)
        
        fit = fitFunction(xdata, ydata, bdata)
        epistasis = ydata - fit(xdata)

        if zstandardize:
	        emap2.loc[sgRNA,:] = epistasis / epistasis.loc[singlesTable['gene'] == 'negative'].std()
       	else:
	        emap2.loc[sgRNA,:] = epistasis 

    emap12 = (emap1+emap2)/2
    
    emap_ave = (emap12 + emap12.T) / 2
    
    return emap1, emap2, emap_ave

#calculate all pairwise intra-sgRNA or intra-gene correlations
def calculateCorrelationMatrix(matrix, diagNull=True):
    corrMatrix = pd.DataFrame(np.corrcoef(matrix), index=matrix.index, columns=matrix.columns)
    
    if diagNull:
        for i in range(len(corrMatrix)):
            corrMatrix.iloc[i,i] = np.nan
            
    return corrMatrix

#find correlations between sgRNAs targeting the same gene and negative controls
def calculateIntrageneCorrelation(sgCorrMatrix, singlePhenotypes, singlesTable):
    sameGeneCorrTups = []
    negCorrTups = []
    for gene, sgs in singlesTable.groupby('gene'):
        for i, (sg1, row) in enumerate(sgCorrMatrix.loc[sgs.index, sgs.index].iterrows()):
            for j, (sg2, val) in enumerate(row.iteritems()):
                if i>j:
                    if gene != 'negative':
                        sameGeneCorrTups.append((sg1, sg2, 
                                                 singlePhenotypes.loc[sg1,['a.mean','b.mean']].mean(), 
                                                 singlePhenotypes.loc[sg2,['a.mean','b.mean']].mean(),
                                                val))
                    else:
                        negCorrTups.append((sg1, sg2, 
                                                 singlePhenotypes.loc[sg1,['a.mean','b.mean']].mean(), 
                                                 singlePhenotypes.loc[sg2,['a.mean','b.mean']].mean(),
                                                val))
                        
    return sameGeneCorrTups, negCorrTups


#generate a gene map by averaging sgRNA epistasis
def generateGeneMap(emap_sgRNA, singlesTable):
    emap_gene = pd.DataFrame(np.zeros((len(set(singlesTable['gene'])),len(set(singlesTable['gene'])))), index = sorted(set(singlesTable['gene'])), columns = sorted(set(singlesTable['gene'])))
    for gene_a, rowgroup in emap_sgRNA.groupby(singlesTable['gene']):
        for gene_b, colgroup in rowgroup.groupby(singlesTable['gene'], axis=1):
            emap_gene.loc[gene_a, gene_b] = colgroup.sum().sum() / (colgroup.shape[0] * colgroup.shape[1])
            
    return emap_gene

### fit functions for calculating interactions and plotting
def linearFitForceIntercept(xdata, ydata, bdata):
    m1 = optimize.fmin(lambda m, x, y: ((m*x + bdata - y)**2).sum(), x0=0.1, args=(xdata, ydata), disp=0)[0]
    
    return lambda x1: m1*np.array(x1) + bdata

def quadFitForceIntercept(xdata, ydata, bdata):
    m1 = optimize.fmin(lambda m, x, y: ((m[0]*(x**2) + m[1]*x + bdata - y)**2).sum(), x0=[0.1,0.1], args=(xdata, ydata), disp=0)
    
    return lambda x1: m1[0]*(np.array(x1)**2) + m1[1]*np.array(x1) + bdata

####	Plotting functions	####

cdict = {'red':((0.0,0.0,0.0),
                (0.5,0.0,0.0),
                (1.0,1.0,1.0)),
        'green':((0.0,0.0,0.0),
                (0.5,0.0,0.0),
                (1.0,1.0,1.0)),
        'blue': ((0.0,1.0,1.0),
                (0.5,0.0,0.0),
                (1.0,0.0,0.0))}
blue_yellow = matplotlib.colors.LinearSegmentedColormap('BlueYellow',cdict)
plt.register_cmap(cmap=blue_yellow)

almost_black = '#262626'

#plot single vs double sgRNA phenotypes showing both a and b positions
def plotSingleVsDouble(sgRNA, phenotypeMatrix, singlePhenotypes, emap, fitFunction=None, showXerr=True):
    fig, axes = plt.subplots(1,3,figsize=(8,4), gridspec_kw={'width_ratios':[6,6,1]})
    
    variablePosition, fixedPosition = 'a','b'
    xdata1, ydata1, bdata1, xerr1 = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition, True)
    variablePosition, fixedPosition = 'b','a'
    xdata2, ydata2, bdata2, xerr2 = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition, True)
    
    minVal = np.nanmin([np.nanmin(xdata1), np.nanmin(ydata1), np.nanmin(xdata2), np.nanmin(ydata2)])
    maxVal = np.nanmax([np.nanmax(xdata1), np.nanmax(ydata1), np.nanmax(xdata2), np.nanmax(ydata2)])
    
    minVal *= 1.2
    maxVal *= 1.2
    
    axis = axes[0]
    
    if showXerr:
        axis.errorbar(xdata1, ydata1, xerr=xerr1, fmt='none', ecolor=almost_black, alpha=.1, zorder=1)
        
    axis.scatter(xdata1, ydata1, c = emap.loc[:,sgRNA], edgecolor = almost_black, cmap='BlueYellow')
    
    
    if fitFunction:
        fit = fitFunction(xdata1, ydata1, bdata1)
        fitrange = np.linspace(minVal,maxVal, 100)
        axis.plot(fitrange, fit(fitrange), color=almost_black)

    axis.plot((minVal,maxVal), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (0,0), color='#BFBFBF',lw=.5)
    axis.plot((0,0), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (bdata1,bdata1), color=almost_black, lw=.5)
    axis.set_xlim((minVal, maxVal))
    axis.set_ylim((minVal, maxVal))
    
    axis.set_xlabel('sgRNA x negative control', fontsize=10)
    axis.set_ylabel('sgRNA x ' + sgRNA, fontsize=10)

    axis.xaxis.set_tick_params(labelsize=10)
    axis.yaxis.set_tick_params(labelsize=10)
    axis.yaxis.tick_left()
    axis.xaxis.tick_bottom()

    axis = axes[1]
    
    if showXerr:
        axis.errorbar(xdata2, ydata2, xerr=xerr2, fmt='none', ecolor=almost_black, alpha=.1, zorder=1)
    
    result = axis.scatter(xdata2, ydata2, c = emap.loc[:,sgRNA], edgecolor = almost_black, cmap='BlueYellow')
    
    if fitFunction:
        fit = fitFunction(xdata2, ydata2, bdata2)
        fitrange = np.linspace(minVal,maxVal, 100)
        axis.plot(fitrange, fit(fitrange), color=almost_black)

    axis.plot((minVal,maxVal), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (0,0), color='#BFBFBF',lw=.5)
    axis.plot((0,0), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (bdata2,bdata2), color=almost_black, lw=.5)
    axis.set_xlim((minVal, maxVal))
    axis.set_ylim((minVal, maxVal))
    
    axis.set_xlabel('negative control x sgRNA', fontsize=10)
    axis.set_ylabel(sgRNA + ' x sgRNA', fontsize=10)
    
    axis.xaxis.set_tick_params(labelsize=10)
    axis.yaxis.set_tick_params(labelsize=10)
    axis.yaxis.tick_left()
    axis.xaxis.tick_bottom()
    
    axis = axes[2]
    axis.spines['top'].set_visible(False)
    axis.spines['bottom'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_xticks([])
    axis.set_yticks([])
    
    cbar = fig.colorbar(result, ax=axis)
    cbar.ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    
    return fig


#plot single vs double sgRNA phenotypes showing only a positions (assuming phenotype matrix has been abba averaged)
def plotSingleVsDouble_abba(sgRNA, phenotypeMatrix, singlePhenotypes, emap, fitFunction=None, showXerr=True):
    fig, axes = plt.subplots(1,2,figsize=(5,4), gridspec_kw={'width_ratios':[6,1]})
    
    variablePosition, fixedPosition = 'a','b'
    xdata1, ydata1, bdata1, xerr1 = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition, True)
    variablePosition, fixedPosition = 'b','a'
    xdata2, ydata2, bdata2, xerr2 = getXYB(sgRNA, singlePhenotypes, phenotypeMatrix, variablePosition, fixedPosition, True)
    
    minVal = np.nanmin([np.nanmin(xdata1), np.nanmin(ydata1), np.nanmin(xdata2), np.nanmin(ydata2)])
    maxVal = np.nanmax([np.nanmax(xdata1), np.nanmax(ydata1), np.nanmax(xdata2), np.nanmax(ydata2)])
    
    minVal *= 1.1
    maxVal *= 1.5
    
    axis = axes[0]
    
    if showXerr:
        axis.errorbar(xdata1, ydata1, xerr=xerr1, fmt='none', ecolor=almost_black, alpha=.1, lw=.5, capsize=2, zorder=1)
        
    colorRange = np.percentile(np.abs(np.reshape(emap.values, (len(emap)**2,))), 90)
    result = axis.scatter(xdata1, ydata1, c = emap.loc[:,sgRNA], s=8, alpha=.75, edgecolor = almost_black, cmap='BlueYellow', vmin=-1*colorRange, vmax=colorRange)
    
    
    if fitFunction:
        fit = fitFunction(xdata1, ydata1, bdata1)
        fitrange = np.linspace(minVal,maxVal, 100)
        axis.plot(fitrange, fit(fitrange), color='#BFBFBF')

    axis.plot((minVal,maxVal), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (0,0), color='#BFBFBF',lw=.5)
    axis.plot((0,0), (minVal,maxVal), color='#BFBFBF',lw=.5)
    axis.plot((minVal,maxVal), (bdata1,bdata1), color=almost_black, lw=.5)
    axis.set_xlim((minVal, maxVal))
    axis.set_ylim((minVal, maxVal))
    
    axis.set_xlabel('sgRNA single phenotype', fontsize=8)
    axis.set_ylabel('sgRNA double phenotype with ' + sgRNA, fontsize=8)

    axis.xaxis.set_tick_params(labelsize=8)
    axis.yaxis.set_tick_params(labelsize=8)
    axis.yaxis.tick_left()
    axis.xaxis.tick_bottom()
    
    axis = axes[1]
    axis.spines['top'].set_visible(False)
    axis.spines['bottom'].set_visible(False)
    axis.spines['left'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.set_xticks([])
    axis.set_yticks([])
    
    cbar = fig.colorbar(result, ax=axis)
    cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    
    return fig

def plotContourFromScatter(x, y, densityLevels = [1,5,10,25,50,75,100], xlabel = '', ylabel = ''):
	phenGrid, phenExtents, phenDensities = fast_kde(x, y, sample=True)

	fig, axis = plt.subplots(figsize=(3,3))

	axis.set_aspect('equal')

	axis.xaxis.tick_bottom()
	axis.yaxis.tick_left()
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)

	axis.scatter(x, y, s=1, c=almost_black, alpha=.5, rasterized=True)
	cf = axis.contour(phenGrid, extent=phenExtents, origin='lower', levels = np.percentile(phenDensities, densityLevels), colors=almost_black, linewidths=.5)
	axis.contourf(phenGrid, extent=phenExtents, origin='lower', levels = np.percentile(phenDensities, densityLevels), colors='w')

	axis.plot((0,0), (-25,10), color='#BFBFBF', lw=.5)
	axis.plot((-25,10), (0,0), color='#BFBFBF', lw=.5)

	axis.set_xlim((np.floor(min(phenExtents)),np.ceil(max(phenExtents))))
	axis.set_ylim((np.floor(min(phenExtents)),np.ceil(max(phenExtents))))
	axis.xaxis.set_tick_params(labelsize=8)
	axis.yaxis.set_tick_params(labelsize=8)

	axis.set_xlabel(xlabel, fontsize=8)
	axis.set_ylabel(ylabel, fontsize=8)

	return fig, axis

def plotSingleVsSelf(phenotypeMatrix, singlePhenotypes):
    fig, axes = plt.subplots(1,2, figsize=(7,3.5))

    diagonal = [phenotypeMatrix.iloc[i,i] for i in range(len(phenotypeMatrix))]
    
    
    minVal = np.nanmin((np.nanmin(singlePhenotypes['a.mean'] - singlePhenotypes['a.std']), np.nanmin(singlePhenotypes['b.mean'] - singlePhenotypes['b.std']), diagonal))
    maxVal = np.nanmax((np.nanmax(singlePhenotypes['a.mean'] + singlePhenotypes['a.std']), np.nanmax(singlePhenotypes['b.mean'] + singlePhenotypes['b.std']), diagonal))

    minVal *= 1.1
    maxVal *= 1.1

    axis = axes[0]
    axis.yaxis.tick_left()
    axis.xaxis.tick_bottom()
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.errorbar(singlePhenotypes['a.mean'], diagonal, xerr=singlePhenotypes['a.std'], fmt='none', ecolor=almost_black, alpha=.1, lw=.5, capsize=2, zorder=1)
    axis.scatter(singlePhenotypes['a.mean'], diagonal, c=almost_black, s=8)

    axis.plot((0,0),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(0,0),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)

    axis.set_xlim((minVal,maxVal))
    axis.set_ylim((minVal,maxVal))

    axis.xaxis.set_tick_params(labelsize=8)
    axis.yaxis.set_tick_params(labelsize=8)

    axis.set_xlabel('sgRNA x negative control gamma', fontsize=8)
    axis.set_ylabel('sgRNA x sgRNA gamma', fontsize=8)

    axis = axes[1]
    axis.yaxis.tick_left()
    axis.xaxis.tick_bottom()
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    axis.errorbar(singlePhenotypes['b.mean'], diagonal, xerr=singlePhenotypes['b.std'], fmt='none', ecolor=almost_black, alpha=.1, lw=.5, capsize=2, zorder=1)
    axis.scatter(singlePhenotypes['b.mean'], diagonal, c=almost_black, s=8)


    axis.plot((0,0),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(0,0),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)

    axis.set_xlim((minVal,maxVal))
    axis.set_ylim((minVal,maxVal))
    
    axis.set_yticklabels([])

    axis.xaxis.set_tick_params(labelsize=8)
    axis.yaxis.set_tick_params(labelsize=8)

    axis.set_xlabel('negative control x sgRNA gamma', fontsize=8)
#     axis.set_ylabel('sgRNA x sgRNA')

    plt.tight_layout()
    
    return fig

def plotSingleAvsB(singlePhenotypes, singlesTable):
    fig, axis = plt.subplots(figsize=(3.5,3.5))

    axis.set_aspect('equal')

    axis.xaxis.tick_bottom()
    axis.yaxis.tick_left()
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)

    minVal = np.nanmin((np.nanmin(singlePhenotypes['a.mean'] - singlePhenotypes['a.std']), np.nanmin(singlePhenotypes['b.mean'] - singlePhenotypes['b.std'])))
    maxVal = np.nanmax((np.nanmax(singlePhenotypes['a.mean'] + singlePhenotypes['a.std']), np.nanmax(singlePhenotypes['b.mean'] + singlePhenotypes['b.std'])))

    minVal *= 1.1
    maxVal *= 1.1

    axis.errorbar(singlePhenotypes['a.mean'], 
                  singlePhenotypes['b.mean'],
                  xerr = singlePhenotypes['a.std'], 
                  yerr = singlePhenotypes['b.std'],
                fmt='none', ecolor=almost_black, alpha=.1, lw=.5, capsize=2, zorder=1)

    axis.scatter(singlePhenotypes['a.mean'], 
                  singlePhenotypes['b.mean'], c=almost_black, s=6, alpha=.9)

    axis.scatter(singlePhenotypes.loc[singlesTable['gene'] == 'negative', 'a.mean'], 
                  singlePhenotypes.loc[singlesTable['gene'] == 'negative', 'b.mean'], c='r', s=6, alpha=.9)

    linfit = np.poly1d(np.polyfit(singlePhenotypes['a.mean'], singlePhenotypes['b.mean'], 1))

    axis.plot((-10,3), linfit((-10,3)), c='#BFBFBF', lw=1)

    axis.plot((0,0),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(0,0),'#BFBFBF',alpha=.5, lw=.5)
    axis.plot((minVal,maxVal),(minVal,maxVal),'#BFBFBF',alpha=.5, lw=.5)

    axis.set_xlim((minVal,maxVal))
    axis.set_ylim((minVal,maxVal))

    axis.text(minVal,.1,'y=%0.2fx+%0.2f' % tuple(np.polyfit(singlePhenotypes['a.mean'], singlePhenotypes['b.mean'], 1)), fontsize=6)
    axis.text(minVal,0,'R=%0.2f;P=%0.2f' % stats.pearsonr(singlePhenotypes['a.mean'], singlePhenotypes['b.mean']), fontsize=6)

    # axis.set_title('mean of all single phenotypes', fon)
    axis.set_xlabel('targeting sgRNA in A position', fontsize=8)
    axis.set_ylabel('targeting sgRNA in B position', fontsize=8)

    axis.xaxis.set_tick_params(labelsize=8)
    axis.yaxis.set_tick_params(labelsize=8)

    plt.tight_layout()

    return fig

#### Contour utility functions ####
#### from https://gist.github.com/joferkington/d95101a61a02e0ba63e5

import scipy.sparse
import scipy.ndimage
import scipy.stats
import scipy.signal
def fast_kde(x, y, gridsize=(400, 400), extents=None, weights=None,
             sample=False):
    """
    Performs a gaussian kernel density estimate over a regular grid using a
    convolution of the gaussian kernel with a 2D histogram of the data.
    This function is typically several orders of magnitude faster than
    scipy.stats.kde.gaussian_kde for large (>1e7) numbers of points and
    produces an essentially identical result.
    Input:
        x: array-like
            The x-coords of the input data points
        y: array-like
            The y-coords of the input data points
        gridsize: tuple, optional
            An (nx,ny) tuple of the size of the output
            grid. Defaults to (400, 400).
        extents: tuple, optional
            A (xmin, xmax, ymin, ymax) tuple of the extents of output grid.
            Defaults to min/max of x & y input.
        weights: array-like or None, optional
            An array of the same shape as x & y that weighs each sample (x_i,
            y_i) by each value in weights (w_i).  Defaults to an array of ones
            the same size as x & y.
        sample: boolean
            Whether or not to return the estimated density at each location.
            Defaults to False
    Output:
        density : 2D array of shape *gridsize*
            The estimated probability distribution function on a regular grid
        extents : tuple
            xmin, xmax, ymin, ymax
        sampled_density : 1D array of len(*x*)
            Only returned if *sample* is True.  The estimated density at each
            point.
    """
    #---- Setup --------------------------------------------------------------
    x, y = np.atleast_1d([x, y])
    x, y = x.reshape(-1), y.reshape(-1)

    if x.size != y.size:
        raise ValueError('Input x & y arrays must be the same size!')

    nx, ny = gridsize
    n = x.size

    if weights is None:
        # Default: Weight all points equally
        weights = np.ones(n)
    else:
        weights = np.squeeze(np.asarray(weights))
        if weights.size != x.size:
            raise ValueError('Input weights must be an array of the same size'
                    ' as input x & y arrays!')

    # Default extents are the extent of the data
    if extents is None:
        xmin, xmax = x.min(), x.max()
        ymin, ymax = y.min(), y.max()
    else:
        xmin, xmax, ymin, ymax = map(float, extents)
    extents = xmin, xmax, ymin, ymax
    dx = (xmax - xmin) / (nx - 1)
    dy = (ymax - ymin) / (ny - 1)

    #---- Preliminary Calculations -------------------------------------------

    # Most of this is a hack to re-implment np.histogram2d using `coo_matrix`
    # for better memory/speed performance with huge numbers of points.

    # First convert x & y over to pixel coordinates
    # (Avoiding np.digitize due to excessive memory usage!)
    ij = np.column_stack((y, x))
    ij -= [ymin, xmin]
    ij /= [dy, dx]
    ij = np.floor(ij, ij).T

    # Next, make a 2D histogram of x & y
    # Avoiding np.histogram2d due to excessive memory usage with many points
    grid = scipy.sparse.coo_matrix((weights, ij), shape=(ny, nx)).toarray()

    # Calculate the covariance matrix (in pixel coords)
    cov = image_cov(grid)

    # Scaling factor for bandwidth
    scotts_factor = np.power(n, -1.0 / 6) # For 2D

    #---- Make the gaussian kernel -------------------------------------------

    # First, determine how big the kernel needs to be
    std_devs = np.diag(np.sqrt(cov))
    kern_nx, kern_ny = np.round(scotts_factor * 2 * np.pi * std_devs)

    # Determine the bandwidth to use for the gaussian kernel
    inv_cov = np.linalg.inv(cov * scotts_factor**2)

    # x & y (pixel) coords of the kernel grid, with <x,y> = <0,0> in center
    xx = np.arange(kern_nx, dtype=np.float) - kern_nx / 2.0
    yy = np.arange(kern_ny, dtype=np.float) - kern_ny / 2.0
    xx, yy = np.meshgrid(xx, yy)

    # Then evaluate the gaussian function on the kernel grid
    kernel = np.vstack((xx.flatten(), yy.flatten()))
    kernel = np.dot(inv_cov, kernel) * kernel
    kernel = np.sum(kernel, axis=0) / 2.0
    kernel = np.exp(-kernel)
    kernel = kernel.reshape((int(kern_ny), int(kern_nx)))

    #---- Produce the kernel density estimate --------------------------------

    # Convolve the gaussian kernel with the 2D histogram, producing a gaussian
    # kernel density estimate on a regular grid

    # Big kernel, use fft...
    if kern_nx * kern_ny > np.product(gridsize) / 4.0:
        grid = scipy.signal.fftconvolve(grid, kernel, mode='same')
    # Small kernel, use ndimage
    else:
        grid = scipy.ndimage.convolve(grid, kernel, mode='constant', cval=0)

    # Normalization factor to divide result by so that units are in the same
    # units as scipy.stats.kde.gaussian_kde's output.
    norm_factor = 2 * np.pi * cov * scotts_factor**2
    norm_factor = np.linalg.det(norm_factor)
    norm_factor = n * dx * dy * np.sqrt(norm_factor)

    # Normalize the result
    grid /= norm_factor

    if sample:
        i, j = ij.astype(int)
        return grid, extents, grid[i, j]
    else:
        return grid, extents
    
def image_cov(data):
    """Efficiently calculate the cov matrix of an image."""
    def raw_moment(data, ix, iy, iord, jord):
        data = data * ix**iord * iy**jord
        return data.sum()

    ni, nj = data.shape
    iy, ix = np.mgrid[:ni, :nj]
    data_sum = data.sum()

    m10 = raw_moment(data, ix, iy, 1, 0)
    m01 = raw_moment(data, ix, iy, 0, 1)
    x_bar = m10 / data_sum
    y_bar = m01 / data_sum

    u11 = (raw_moment(data, ix, iy, 1, 1) - x_bar * m01) / data_sum
    u20 = (raw_moment(data, ix, iy, 2, 0) - x_bar * m10) / data_sum
    u02 = (raw_moment(data, ix, iy, 0, 2) - y_bar * m01) / data_sum

    cov = np.array([[u20, u11], [u11, u02]])
    return cov
