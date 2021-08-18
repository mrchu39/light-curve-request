import base64
import matplotlib.pyplot as plt
import numpy as np
import os
import shlex
import shutil
import sncosmo
import subprocess
import sys
import time

from astropy.io import ascii
from astropy.table import QTable
from astropy.cosmology import FlatLambdaCDM
from func import *

# These numbers come from running model fits on ~500 Type Ia supernovae
z = 0.060574239946858476
z_std = 0.023994157056121096
x1 = -0.14238796934437334
x1_std = 1.4557579021314682
c = 0.08928354223298558
c_std = 0.15670291093588692
x0 = 0.0007648532623426458
x0_std = 0.0004363803462578883

with open('info.info', 'r') as f:
    SNID_loc = f.read().split('\n')[0].split(':')[1].strip()

def get_photometry(ztfname, format='flux'):

    ''' Info : Retrieves photometry data for a source from Fritz and filters out Nonetype points
        Input : Source name and brightness format ("flux" or "mag")
        Returns : Astropy QTable with data that feeds into sncosmo.fit_lc
    '''

    url = BASEURL+'api/sources/'+ztfname+'/photometry' # Access photometry

    if format == 'flux':
        data = {"format": "flux"}
    elif format == 'mag':
        data = {'format': 'mag'}

    response = api('GET', url, params=data, timeout=10)

    if format == 'flux':

        flux = []
        fluxerr = []
        band = []
        mjd = []
        zpsys = []
        zp = []

        for d in response['data']:
            if d['flux'] != None and (d['filter'] == 'ztfg' or d['filter'] == 'ztfr' or d['filter'] == 'ztfi'):
                flux.append(d['flux'])
                fluxerr.append(d['fluxerr'])
                band.append(d['filter'])
                mjd.append(d['mjd'])
                zpsys.append(d['magsys'])
                zp.append(d['zp'])

        return QTable([mjd, band, flux, fluxerr, zp, zpsys], names=('mjd', 'filter', 'flux','fluxerr', 'zp', 'zpsys'))

    elif format == 'mag':

        mag = []
        magerr = []
        band = []
        mjd = []
        zpsys = []

        for d in response['data']:
            if d['mag'] != None and (d['filter'] == 'ztfg' or d['filter'] == 'ztfr' or d['filter'] == 'ztfi'):
                mag.append(d['mag'])
                magerr.append(d['magerr'])
                band.append(d['filter'])
                mjd.append(d['mjd'])
                zpsys.append(d['magsys'])

        return QTable([mjd, band, mag, magerr, zpsys], names=('mjd', 'filter', 'mag','magerr', 'zpsys'))

def model_lc(source):

    ''' Info : Fits photometry data to light curve using sncosmo.
        Input : source
        Returns : photometry data, fitted parameters, plottable model
    '''

    data = get_photometry(source)

    with open('test_table.txt', 'w') as f:
        f.write(str(data))

    red, red_err = get_redshift(source, True)
    model = sncosmo.Model(source='salt2')

    if red != 'No redshift found':

        if red_err != 'No redshift error found': # If both redshift and error present, use error as bounds to fit redshift
            result, fitted_model = sncosmo.fit_lc(
                data, model,
                ['z', 't0', 'x0', 'x1', 'c'],  # parameters of model to vary
                bounds={'z':(red-red_err, red+red_err)}, minsnr=5)  # bounds on parameters (if any)
        else: # If no redshift error, don't determine redshift
            model.set(z=red)

            result, fitted_model = sncosmo.fit_lc(
                data, model,
                ['t0', 'x0', 'x1', 'c'],
                guess_z=False, minsnr=5)

    else:
        result, fitted_model = sncosmo.fit_lc(
            data, model,
            ['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z':(0,0.3)}, minsnr=5)

    return data, result, fitted_model

def get_sigmas(result):
    x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
    c_nstds = np.round(np.abs((result.parameters[4]-c))/c_std, 1)

    return x1_nstds, c_nstds

def get_peak_absmag(z, x0):
    peak_mag = -2.5*np.log10(x0) + 10.635
    cosmo = FlatLambdaCDM(H0=70,Om0=0.3)
    mu = cosmo.distmod(z).value
    absmag = peak_mag -mu

    return absmag

def post_lc(source):
    data = get_photometry(source)
    comment_infos = get_source_api(source)['comments']

    for i in range (len(get_source_api(source)['comments'])):

        comment_info = comment_infos[i]
        comment = comment_info['text']

        if 'sncosmo light curve fit' in comment:
            if int(comment[int(comment.index('n='))+2:].split(',')[0]) != len(data) or 'x1_nstds' not in comment or 'c_nstds' not in comment or 'M_peak' not in comment:

                try:
                    dfit, result, fitted_model = model_lc(source)
                except RuntimeError:
                    return

                x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
                c_nstds = np.round(np.abs((result.parameters[4]-c))/c_std, 1)

                sncosmo.plot_lc(dfit, model=fitted_model)
                plt.savefig('temp.png')

                resp = edit_comment(source, comment_info['id'], comment_info['author_id'], 'sncosmo light curve fit n='+str(len(data))+', M_peak = '+str(np.round(get_peak_absmag(result.parameters[0], result.parameters[2]),1))+
                    ', x1_nstds = '+str(x1_nstds)+', c_nstds = '+str(c_nstds), 'temp.png', source+'_sncosmo_lc.png')

                plt.close('all')

                return resp
            else:
                return False

    try:
        dfit, result, fitted_model = model_lc(source)
    except RuntimeError:
        print(bcolors.FAIL + 'sncosmo encountered runtime error. Skipping...' + bcolors.ENDC)
        return

    x1_nstds = np.round(np.abs((result.parameters[3]-x1)/x1_std), 1)
    c_nstds = np.round(np.abs((result.parameters[4]-c)/c_std), 1)

    sncosmo.plot_lc(dfit, model=fitted_model)
    plt.savefig('temp.png')

    resp = resp = post_comment(source, 'sncosmo light curve fit n='+str(len(data))+', M_peak = '+str(np.round(get_peak_absmag(result.parameters[0], result.parameters[2]),1))+
        ', x1_nstds = '+str(x1_nstds)+', c_nstds = '+str(c_nstds), 'temp.png', source+'_sncosmo_lc.png')

    plt.close('all')

    return resp
