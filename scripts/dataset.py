# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 10:22:41 2016

@author: rstreet
"""
from os import path
import glob
import numpy as np

from matplotlib import use as useBackend
from matplotlib import pyplot

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
import time_functions
import statistics
useBackend('Agg')

class Star(object):
    """Class simply keeps the individual stars in the analysis separate and
    easy to access for coordinates and IDs"""
    def __init__(self, idstr, rastr, decstr):
        self.id = idstr

        if ':' in str(rastr) and ':' in str(decstr):
            self.coord = SkyCoord(rastr, decstr, unit=(u.hourangle, u.deg))

        else:
            self.coord = SkyCoord(ra=rastr*u.deg, dec=decstr*u.deg)

class DataSet(object):
    """DataSet contains all of the directories, images, and stars used in the
    differential photometry. This class also carries out the analysis."""

    def __init__(self, params=None):
        self.input_dir = None
        self.output_dir = None
        self.input_frames = None
        self.valid_frames = []
        self.data_cube = np.zeros(1)
        self.residuals_cube = np.zeros(1)
        self.ensemble_flux = np.zeros(1)
        self.ensemble_fluxerr = np.zeros(1)
        self.qc = []
        self.nstars = 0
        self.nframes = 0
        self.star_list = []
        self.star_ids = []

        if params:
            for key, value in params.items():
                setattr(self, key, value)

    def load_star_list(self, params):
        """Reads in the star list text file and creates Stars for each item"""
        self.star_list = []
        self.star_ids = []

        entries = open(params['star_file'], 'r').readlines()

        for line in entries:
            line.replace('\n', '')
            if line and (not line.startswith("#")):
                #idstr, rastr, decstr = line.split()
                star = Star(*line.split())
                self.star_list.append(star)
                self.star_ids.append(star.id)

    def get_star_idx(self, star):
        """Gets the id of a star"""
        #Why is this here? It doesn't seem to really save any typing or add clarity
        return self.star_ids.index(star.id)

    def print_star_list(self):
        """Prints the list of stars and coordinates to the terminal"""
        for star in self.star_list:
            print "%s %s"%(str(star.id), star.coord.to_string('hmsdms'))

    def load_data(self, params):
        """Reads in the data images and creates the appropriate table"""
        self.nstars = len(self.star_list)

        if self.input_dir is None:
            self.input_dir = params['input_dir']

        self.input_frames = glob.glob(path.join(self.input_dir, '*.fits'))
        self.input_frames.sort()
        self.nframes = len(self.input_frames)

        self.data_cube = np.zeros([self.nstars, self.nframes, 7])
        self.residuals_cube = np.zeros([self.nstars, self.nframes, 2])
        self.data_cube.fill(np.NaN)
        self.residuals_cube.fill(np.NaN)
        self.ensemble_flux = np.zeros([self.nframes])
        self.ensemble_flux.fill(np.NaN)
        self.ensemble_fluxerr = np.zeros([self.nframes])
        self.ensemble_fluxerr.fill(np.NaN)

        for i in xrange(0, len(self.input_frames)):
            print "Reading data for %s"%path.basename(self.input_frames[i])
            table = BanzaiTable(self.input_frames[i])

            if table.got_data:
                self.valid_frames.append(table.data_file)
                stars_data = extract_star_data(self.star_list, table)
                ts = time_functions.calc_hjd(table.date_obs,
                                             stars_data[:, 2], stars_data[:, 3])
                self.data_cube[:, i, 0] = ts - 2450000.0
                for k in xrange(0, stars_data.shape[1]):
                    self.data_cube[:, i, k+1] = stars_data[:, k]

    def diff_photometry(self):
        """Method to compute differential photometry for the selected stars"""

        self.calc_residuals()
        print "Calculated residuals"
        self.plot_lightcurves(self.star_list, suffix='norm')
        print "Plotted normalized lightcurves"
        self.calc_ensemble_lightcurve()
        print "Calculated the ensemble lightcurve"
        self.calc_differential_lightcurve()
        print "Completed differential photometry"

    def calc_residuals(self, debug=False):
        """Method to calculate the residual lightcurve for all stars in the
        data cube by dividing by the Median Absolute Deviation from each
        lightcurve, weighted by the photometric errors"""

        for j in range(0, self.nstars):
            flux = self.data_cube[j, :, 5]
            fluxerr = self.data_cube[j, :, 6]
            wmean = statistics.calc_weighted_mean(flux, fluxerr)[0]
            idx = np.isfinite(flux)
            self.residuals_cube[j, idx, 0] = flux[idx]/wmean
            self.residuals_cube[j, idx, 1] = fluxerr[idx]/wmean

            if debug:
                print "Mean flux of star %s: %0.3f"%(self.star_ids[j], wmean)

    def calc_ensemble_lightcurve(self):
        """Method to calculate the ensemble lightcurve for differential
        photometry"""

        residuals = self.residuals_cube[1:, :, 0]
        sigmas = self.residuals_cube[1:, :, 1]
        for i in xrange(0, self.nframes):
            self.ensemble_flux[i] = statistics.calc_weighted_mean(residuals[:, i],
                                                                  sigmas[:, i])[0]
            self.ensemble_fluxerr[i] = statistics.calc_weighted_sigma(residuals[:, i],
                                                                      sigmas[:, i],
                                                                      self.ensemble_flux[i])


    def calc_differential_lightcurve(self):
        """Method to calculate the differential lightcurves of all stars
        in the data cube"""

        ensemble_sigma_sq = 1.0/self.ensemble_fluxerr**2
        for j in xrange(0, self.nstars):
            self.residuals_cube[j, :, 0] = self.residuals_cube[j, :, 0]/self.ensemble_flux
            sigmas = self.residuals_cube[j, :, 1]
            weights = (1.0/sigmas**2)+ensemble_sigma_sq
            self.residuals_cube[j, :, 1] = np.sqrt(1.0/weights)

    def quality_control(self):
        """Method to identify suspect frames in the dataset from the scatter
        in the photometry of the comparison stars"""

        idx = np.arange(0, self.nframes)
        flags = []
        for star in self.star_list[1:]:
            j = self.get_star_idx(star)
            flux = self.residuals_cube[j, :, 0]
            kdx = np.isfinite(flux)

            med = np.median(flux[kdx])
            mad = statistics.calc_mad(flux[kdx])

            jdx = np.where(abs(flux[kdx]-med) > mad)[0].tolist()
            flags += idx[kdx][jdx].tolist()

        for i in xrange(0, self.nframes):
            if flags.count(i) == self.nstars-1:
                self.qc.append(i)

    def plot_lightcurves(self, selected_stars, suffix=None):
        """Method to plot lightcurve files of a selected range of stars
        from the star list.  Expects a list of Star objects"""

        for star in selected_stars:
            j = self.get_star_idx(star)
            ts = self.data_cube[j, :, 0]
            flux = self.residuals_cube[j, :, 0]
            flux_err = self.residuals_cube[j, :, 1]
            idx = np.isfinite(flux)
            if len(idx) > 0:
                pyplot.errorbar(ts[idx], flux[idx], yerr=flux_err[idx], fmt='k.',
                                mfc='k', mec='k', ms=2, capsize=1)
                pyplot.plot(ts[self.qc], flux[self.qc], 'rx', ms=4)
                med = np.median(flux[idx])
                mad = statistics.calc_mad(flux[idx])
#                wmean = statistics.calc_weighted_mean(flux[idx], flux_err[idx])[0]
#                wsig = statistics.calc_weighted_sigma(flux[idx], flux_err[idx], wmean)

                pyplot.plot([ts[0], ts[-1]], [med, med], 'r-')
                pyplot.plot([ts[0], ts[-1]], [med+mad, med+mad], 'r-.')
                pyplot.plot([ts[0], ts[-1]], [med-mad, med-mad], 'r-.')
                pyplot.xlabel('HJD-2450000.0')
                pyplot.ylabel('Residual flux')
                pyplot.title("Lightcurve of %s"%star.id)

                if suffix:
                    plt_file = path.join(self.output_dir,
                                         "lightcurve_%s_%s.png"%(suffix, j))

                else:
                    plt_file = path.join(self.output_dir, "lightcurve_%s.png"%j)

                pyplot.savefig(plt_file)
                pyplot.close(1)
            else:
                print "Warning: Insufficient data for lightcurve: %s"%star.id

    def plot_ensemble_lightcurve(self):
        """Method to plot the ensemble lightcurve"""

        pyplot.figure(1)
        ts = self.data_cube[0, :, 0]
        flux = self.ensemble_flux
        flux_err = self.ensemble_fluxerr
        idx = np.isfinite(flux)
        pyplot.errorbar(ts[idx], flux[idx], yerr=flux_err[idx], fmt='k.', mfc='k',
                        mec='k', ms=2, capsize=1)
        pyplot.xlabel('HJD-2450000.0')
        pyplot.ylabel('Residual flux')
        pyplot.title('Ensemble lightcurve')
        plt_file = path.join(self.output_dir, 'ensemble_lightcurve.png')
        pyplot.savefig(plt_file)
        pyplot.close(1)

    def output_lightcurves(self, selected_stars):

        def flux_to_mag(flux):
            mag = 2.5*np.log(flux)
            return mag

        def calc_mag(flux, flux_err):
            mag = flux_to_mag(flux)
            df1 = flux - flux_err
            m1 = flux_to_mag(df1)
            df2 = flux + flux_err
            m2 = flux_to_mag(df2)
            mag_err = (m2-m1)/2.0
            return mag, mag_err

        for star in selected_stars:
            j = self.get_star_idx(star)
            out_file = path.join(self.output_dir, "lightcurve_%s.txt"%j)
            head_str = "HJD-2450000.0    X[pix]    Y[pix]    Flux    Flux_err    \
                        Res_flux    Res_flux_err    Res_mag Res_mag_err QC_flags"
            frames = self.input_frames
            ts = np.array(self.data_cube[j, :, 0])
            x = np.array(self.data_cube[j, :, 1])
            y = np.array(self.data_cube[j, :, 2])
            flux = np.array(self.data_cube[j, :, 5])
            flux_err = np.array(self.data_cube[j, :, 6])
            rflux = np.array(self.residuals_cube[j, :, 0])
            rflux_err = np.array(self.residuals_cube[j, :, 1])
            mag, mag_err = calc_mag(rflux, rflux_err)

            idx = np.where(flux > -99.0)
            q_flags = np.zeros(len(frames))
            q_flags[self.qc] == 1
            out_array = np.array([ts[idx], x[idx], y[idx], flux[idx], flux_err[idx],
                                  rflux[idx], rflux_err[idx], mag[idx], mag_err[idx],
                                  q_flags[idx]])

            np.savetxt(out_file, np.transpose(out_array), fmt="%s", header=head_str)


class BanzaiTable(object):
    """I think this is a big data table thing that reads in FITS images?"""
    def __init__(self, file_path=None):
        self.data_file = None
        self.date_obs = None
        self.exptime = None
        self.obs_lat = None
        self.obs_long = None
        self.obs_height = None
        self.x = None
        self.y = None
        self.coord = None
        self.flux = None
        self.flux_err = None
        self.got_data = False

        if file_path:
            self.read_banzai_table(file_path)

    def check_header(self, header):
        """Method to verify that a Banzai table has the required header
        parameters"""

        status = True
        required_params = ['DATE-OBS', 'EXPTIME', 'LATITUDE', 'LONGITUD', 'HEIGHT']
        for par in required_params:
            if par not in header.keys():
                status = False
        return status

    def check_table(self, table):
        """Method to verify that a Banzai photometry table contains the
        required columns"""

        status = True
        required_columns = ['X', 'Y', 'RA', 'DEC', 'FLUX', 'FLUXERR']
        for col in required_columns:
            if col not in table.dtype.fields:
                status = False
        return status

    def read_banzai_table(self, file_path):
        """Method to read the photometry table from a single reduced BANZAI
        FITS data product"""

        hdu_list = fits.open(file_path)
        if len(hdu_list) == 3:
            table = hdu_list[1].data
            header_ok = self.check_header(hdu_list[0].header)
            table_ok = self.check_table(table)

            if header_ok and table_ok:
                self.data_file = path.basename(file_path)
                self.date_obs = hdu_list[0].header['DATE-OBS']
                self.exptime = float(hdu_list[0].header['EXPTIME'])
                self.obs_lat = hdu_list[0].header['LATITUDE']
                self.obs_long = hdu_list[0].header['LONGITUD']
                self.obs_height = hdu_list[0].header['HEIGHT']

                self.x = table.field('X')
                self.y = table.field('Y')
                self.coord = SkyCoord(ra=(table.field('RA')*u.degree),
                                      dec=(table.field('DEC')*u.degree))
                self.flux = table.field('FLUX')
                self.flux_err = table.field('FLUXERR')
                self.got_data = True

            else:
                print "ERROR: Data table for %s incomplete"%path.basename(file_path)
                self.got_data = False

        else:
            print "ERROR: %s has too few FITS extensions"%path.basename(file_path)
            self.got_data = False

def extract_star_data(star_list, table):
    """Function to extract the data for a given star from a table array
    star should be a SkyCoord object
    """
    data = np.zeros([len(star_list), 6])
    for j in range(0, len(star_list)):
        sep = table.coord.separation(star_list[j].coord)
        idx = np.where(sep == sep.min())[0]
        data[j] = table.x[idx], table.y[idx], table.coord[idx].ra.value,\
            table.coord[idx].dec.value, table.flux[idx], table.flux_err[idx]
#        data[j, 0] = table.x[idx]
#        data[j, 1] = table.y[idx]
#        data[j, 2] = table.coord[idx].ra.value
#        data[j, 3] = table.coord[idx].dec.value
#        data[j, 4] = table.flux[idx]
#        data[j, 5] = table.flux_err[idx]

    return data
