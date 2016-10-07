# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 16:42:38 2016

@author: rstreet
"""
from sys import argv
from datetime import datetime
from astropy.coordinates import SkyCoord
from astropy import units as u
from pyslalib import slalib

def calc_hjd(date_obs, ra, dec, debug=False):
    """Function to convert a timestamp in %Y-%m-%dT%H:%M:%S UTC format
    to Heliocentric Julian Date for a given RA, Dec of target"""

    # Convert RA, Dec to radians:
    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Convert the timestamp into a DateTime object:
    if 'T' in date_obs:
        try:
            dt = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S")
        except ValueError:
            dt = datetime.strptime(date_obs, "%Y-%m-%dT%H:%M:%S.%f")
    else:
        try:
            dt = datetime.strptime(date_obs, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            dt = datetime.strptime(date_obs, "%Y-%m-%d %H:%M:%S.%f")

    # Calculate the MJD (UTC) timestamp:
    mjd_utc = datetime_to_mjd_utc(dt)

    # Correct the MJD to TT:
    mjd_tt = mjd_utc_to_mjd_tt(mjd_utc)

    # Calculating MJD of 1st January that year:
    mjd_jan1 = slalib.sla_cldj(dt.year, 1, 1)[0]

    # Calculating the MJD difference between the DateObs and Jan 1 of the same year:
    tdiff = mjd_tt - mjd_jan1

    # Calculating the RV and time corrections to the Sun:
    tcorr = slalib.sla_ecor(coord.ra.radian, coord.dec.radian,\
                                dt.year, int(tdiff), (tdiff-int(tdiff)))[1]
    # Calculating the HJD:
    hjd = mjd_tt+(tcorr/86400.0)+2400000.5

    if debug:
        print "Dec: %s, Dec [rads]: %s"%(dec, coord.dec.radian)
        print "RA: %s, Dec. RA: %s"%(ra, coord.ra.radian)
        print "MCD_UTC: %s"%mjd_utc
        print "MCD_TT: %s"%mjd_tt
        print "MCD of 1 Jan %s: %s"%(dt.year, mjd_jan1)
        print "Time difference from 1 Jan to Obs. Date, %s: %s%(dt.year, tdiff"
        print "Time correctiopn to the Sun: %s"%tcorr
        print "HJD: %s"%hjd

    return hjd

def datetime_to_mjd_utc(utc_date):
    """Function to calculate MJD for a given UTC"""

    mjd, status = slalib.sla_cldj(utc_date.year, utc_date.month, utc_date.day)

    if status != 0:
        return None
    fday, status = slalib.sla_dtf2d(utc_date.hour, utc_date.minute,
                                    utc_date.second+(utc_date.microsecond/1e6))
    if status != 0:
        return None

    return mjd + fday

def mjd_utc_to_mjd_tt(mjd_utc, debug=False):
    """Converts a MJD in UTC (MJD_UTC) to a MJD in TT (Terrestial Time) which is
    needed for any position/ephemeris-based calculations.UTC
    UTC->TT consists of: UTC->TAI = 10s offset + 24 leapseconds (last one 2009 Jan 1.)
    TAI->TT  = 32.184s fixed offset"""

    tt_utc = slalib.sla_dtt(mjd_utc)
    mjd_tt = mjd_utc + (tt_utc/86400.0)

    if debug:
        print "TT-UTC(s): %s"%tt_utc
        print "MJD(TT): %s"%mjd_tt

    return mjd_tt


if __name__ == '__main__':

    if len(argv) == 1:
        print 'Call sequence: python calctime.py DateObsStr RA Dec'
        exit(0)
    else:
        print calc_hjd(*argv[1:], debug=True)
