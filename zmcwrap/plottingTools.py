"""
Functions related to plotting the results from an MCMC run
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sf

def chainPlot(parsByWalker, labels):
    """
    Make a plot of the values of each walker, for each parameter,
    as a function of step in the chain.  Helpful for diagnosing
    the convergence of the MCMC chain.
    
    Takes:
    parsByWalker - samples from the MCMC chain.
             An array with dimensions (nParams, nSteps, nWalkers)
    labels - names to be used for the parameters, as axis labels.
             A list of strings with length nParams

    Returns:
    fig - a matplotlib figure containing the plot of the chain
          Use plt.show() to display the figure, or use the save method
          of the figure to save it: fig.savefig('corner.pdf')
    """
    fmt_str = '{x:.6g}'
    npars = len(labels)
    fig, axs = plt.subplots(npars, 1, sharex=True,
        gridspec_kw={'hspace':0, 'wspace':0, 'top':0.97, 'bottom':0.08})
    for i in range(npars):
        axs[i].set_ylabel('{:}'.format(labels[i]))
        axs[i].plot(parsByWalker[i,:,:], 'k', alpha=0.1)
        axs[i].tick_params(axis='both', direction='in', top=True, right=True)
        axs[i].yaxis.set_major_formatter(fmt_str)
    axs[-1].set_xlabel('Step number')
    return fig


def cornerPlotLight(samples, labels, bins=20, contour_sigmas=[1.0, 2.0],
                    hist_sigmas=None, trim_range=True, figsize=(8, 8)):
    """
    Make a corner plot.
    Mostly useful for plotting the results of an MCMC chain.
    Similar to, and inspired by, corner.py from Foreman-Mackey (2016).

    Takes:
    samples - samples from the MCMC chain.
              An array with dimensions (nSampes, nParams)
    labels - names to be used for the parameters, as axis labels.
             A list of strings with length nParams
    bins - the number of bins to be used in the 2D and 1D histograms (integer)
    contour_sigmas - A list of sigma levels for drawing contours on top of the
                     2D histograms. A list of floats.  Or, if this is set to
                     None, no contours are drawn.
    hist_sigmas - A list of sigma levels for drawing lines on top of the 
                  1D histograms.  A list of floats.  Or, if this is set to None,
                  no lines are drawn.
    trim_range - a boolean flag for whether to a slightly restricted range (True)
                 or the full range (False) of the values in samples.
    figsize - the size of the matplotlib figure to generate.

    Returns:
    fig - a matplotlib figure containing the corner plot
          Use plt.show() to display the figure, or use the save method
          of the figure to save it: fig.savefig('corner.pdf')
    """

    nsamples = samples.shape[0]
    nparam = samples.shape[-1]
    minSamp = np.amin(samples, axis=0)
    maxSamp = np.amax(samples, axis=0)
    if trim_range:
        minPer = np.percentile(samples, 1.0, axis=0)
        maxPer = np.percentile(samples, 99.0, axis=0)
        minPlot = minPer - (maxPer-minPer)*0.1
        maxPlot = maxPer + (maxPer-minPer)*0.1
    else:
        minPlot = minSamp
        maxPlot = maxSamp

    # This follows the contour 'sigma' levels from corner.py
    # (Foreman-Mackey, 2016, JOSS, 1(2), 24
    # In general, a two dimensional Gaussian would be
    # f(x,y) = A*exp(-(((x - x0)^2)/(2*sigma_x^2) + ((y - y0)^2)/(2*sigma_y^2)))
    # for that to integrate to unity, A = 1/(2*pi*sigma_x*sigma_y).
    #
    # In the symmetric case where sigma = sigma_x = sigma_y and r^2 = (x - x0)^2 + (y - y0)^2
    # [Or in an elliptical case where sigma = sigma_x/a = sigma_y/b and r^2 = ((x - x0)/a)^2 + ((y - y0)/b)^2.]
    # f(x,y) = 1/(2*pi*sigma^2) * exp(-(r^2)/(2*sigma^2))
    #
    # then the cumulative probability distribution function,
    # using polar coordinates to make life easy in this symmetric approximation, is
    # cdf = integral_{from 0 to 2pi} integral_{from 0 to x} 1/(2*pi*sigma^2) * exp(-0.5*(r/sigma)^2) r dr dtheta
    #     = integral_{from 0 to x} r/sigma^2 * exp(-0.5*(r/sigma)^2) dr
    #     = 1 - exp(-0.5*(r/sigma)^2)     [at r = 0 cdf = 0]
    #
    # This leads to different fractions for 1, 2, and 3 sigma values from
    # what one would expect in the 1D Gaussian case.
    # 1 sigma is ~0.39347 of the distribution, not ~0.68269,
    # and 2 sigma is ~0.86466 of the distribution, not 0.95450.
    
    #Get contour levels
    if contour_sigmas is not None:
        cdf_lvls = 1.0 - np.exp(-0.5 * np.array(contour_sigmas)**2)

    fig = plt.figure(figsize=figsize)

    fmt_str = '{x:.5g}'

    nframe = 0
    for i in range(nparam): #loop over y positions in grid of plots
        j=-1
        for j in range(i): #loop over x positions in grid of plots
            nframe = i*(nparam)+j+1
            if i == nparam-1 and j==0:
                ax = plt.subplot(nparam, nparam, nframe, 
                                 xlabel=labels[j], ylabel=labels[i])
                ax.xaxis.set_major_formatter(fmt_str)
                ax.yaxis.set_major_formatter(fmt_str)
            elif i == nparam-1:
                ax = plt.subplot(nparam, nparam, nframe, 
                                 xlabel=labels[j], yticklabels=[])
                ax.xaxis.set_major_formatter(fmt_str)
            elif j == 0:
                ax = plt.subplot(nparam, nparam, nframe, 
                                 ylabel=labels[i], xticklabels=[])
                ax.yaxis.set_major_formatter(fmt_str)
            else:
                ax = plt.subplot(nparam, nparam, nframe, 
                                 xticklabels=[], yticklabels=[])
             
            #H, X, Y = np.histogram2d(samples[:,j], samples[:,i], bins=bins)
            H, X, Y, img = ax.hist2d(samples[:,j], samples[:,i],
                                     bins=bins, cmap='binary') 
            Xc = (X[1:] + X[:-1])*0.5
            Yc = (Y[1:] + Y[:-1])*0.5
            
            # generate a cumulative histogram for the pixels/histogram in H
            if contour_sigmas is not None:
                Hsort = np.sort(H, axis=None, kind='stable')
                Hcumdist = np.cumsum(Hsort)/float(nsamples)
                # and set the contour levels based on the cumulative histogram
                contours = np.zeros_like(cdf_lvls)
                for ilvl, lvl in enumerate(cdf_lvls):
                    contours[ilvl] = Hsort[Hcumdist > (1.0 - lvl)][0]
                contours.sort()                
                ax.contour(Xc, Yc, H.T, levels=contours, colors='black')

                # One could smooth the 2D histogram to get smother contours
                #import scipy.ndimage as simg
                #Hsmooth = simg.gaussian_filter(H, sigma=1)
                #ax.contour(Xc, Yc, Hsmooth.T, levels=contours, colors='green')
            
            # Set limits and tick marks
            ax.set_xlim((minPlot[j],maxPlot[j]))
            ax.set_ylim((minPlot[i],maxPlot[i]))
            ax.tick_params(axis='both', direction='in', top=True, right=True,
                           labelrotation=45)
            ax.locator_params(axis='both', nbins=4, min_n_ticks=2, prune='lower')
    
        #Add histograms at the top of the columns of plots
        if (nframe+1 == (nparam)**2):
            ax = plt.subplot(nparam, nparam, nframe+1,
                             yticklabels=[], xlabel=labels[j+1])
            ax.xaxis.set_major_formatter(fmt_str)
        else:
            ax = plt.subplot(nparam, nparam, nframe+1,
                             yticklabels=[], xticklabels=[])

        hval, bbin, img = ax.hist(samples[:,j+1], 
                                  histtype='step', color='k', bins=bins)

        # Optionally, plot vertical lines on the histograms for specific
        # confidence levels
        if hist_sigmas is not None:
            # Gaussian cumulative probability distribution,
            # evaluated for for + and - sigma levels
            sigmas = np.array(hist_sigmas)
            sigmas = np.concatenate((-sigmas, [0.0], sigmas))
            cdf_lvls_h = 0.5*(1.0 + sf.erf(sigmas/np.sqrt(2.0)))
            
            percentiles = np.percentile(samples[:,j+1], 100.*cdf_lvls_h)
            for per in percentiles:
                if per >= bbin[-1]:
                    ymax = hval[-1]
                else:
                    ymax = hval[np.nonzero(bbin <= per)[0][-1]]
                ax.plot([per, per], [0., ymax], ':k', lw=1.0)

        # set up ticks and limits
        ax.set_xlim((minPlot[j+1],maxPlot[j+1]))
        ax.tick_params(axis='both', direction='in', left=False, labelrotation=45)
        plt.locator_params(axis='both', nbins=4, min_n_ticks=2)
    
    plt.tight_layout()
    fig.subplots_adjust(wspace=0.0, hspace=0.0)
    return fig
