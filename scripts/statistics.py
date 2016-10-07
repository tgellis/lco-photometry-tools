# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:11:59 2016

@author: rstreet
"""

import numpy as np

def calc_weighted_mean(data, sigma):
    """Function to calculate the mean of a data vector, weighted by the
    given uncertainties"""

    idx = np.isfinite(data)
    if True in idx:
        sigma_sq = 1.0/(sigma**2)
        wmean = sum(data[idx]*sigma_sq[idx])/sum(sigma_sq[idx])
        sigma_wmean = 1.0/sum(sigma_sq[idx])

    else:
        wmean = np.NaN
        sigma_wmean = np.NaN

    return wmean, sigma_wmean

def calc_weighted_sigma(data, sigma, wmean):
    """Function to calculate the standard deviation of a data vector, weighted by the
    uncertainties"""

    idx = np.isfinite(data)
    if True in idx:
        weights = 1.0/(sigma[idx]**2)
        ddata = data[idx] - wmean
        wsigma = np.sqrt(sum(ddata**2*weights)/sum(weights))

    else:
        wsigma = np.NaN

    return wsigma

def calc_mad(data):
    """Function to calculate the MAD of a data vector and uncertainties"""

    idx = np.isfinite(data)
    if True in idx:
        return np.median(abs(data[idx]-np.mean(data[idx])))
    else:
        return np.NaN

if __name__ == '__main__':
    _data = np.random.normal(12.0, 0.2, 100)
    _sigma = np.random.normal(0.2, 0.05, 100)
    MAD = calc_mad(_data)
    print "MAD: %0.3f"%MAD
    