import datetime
from astropy.time import Time as astropyTime
import warnings
import numpy as np

import libpySat as pySat

def UTC2MJD(epoch: datetime.datetime):
    """
    Converts UTC epoch to MJD
    :return: mjd
    """
    return astropyTime(epoch, format='datetime', scale='utc').mjd

def UT12MJD(epoch: datetime.datetime):
    """
    Converts UTC epoch to MJD
    :return: mjd
    """
    return astropyTime(epoch, format='datetime', scale='ut1').mjd

def getdTJ2000(epoch:datetime.datetime,scale: str):
    """

    :param epoch - time in UTC
    :param scale - string corresponding to astopy's format ('tt', 'utc' ...)
    :return: difference w.r.t to epoch J2000 in julian days in scale
    """
    # 1) Define the J2000 epoch and given epoch in UTC,TT,UT1
    j2000 = astropyTime('2000-01-01T12:00:00', format='isot', scale=scale)
    date = astropyTime(epoch, format='datetime', scale='utc')

    return date.jd1-j2000.jd1+date.jd1-j2000.jd2

def Sod(epoch: datetime.datetime):
    """
    Calculates seconds from midnight
    :return: seconds of the day
    """
    return epoch.hour*3600+epoch.minute*60+epoch.second+epoch.microsecond*1e-6

def UTC2JD1JD2(epoch: datetime.datetime):
    """

    :param epoch:
    :return:
    """
    return astropyTime(epoch, format='datetime', scale='ut1').jd1, astropyTime(epoch, format='datetime', scale='ut1').jd2

def getUT1(epoch: datetime.datetime):
    """
    Converts gets UT1 at UTC epoch
    :return: ut1
    """
    iers_tab = pySat.GetIersTable()
    ut1=epoch
    mjd = UTC2MJD(epoch)
    try:
        ut1utc = np.interp(mjd, iers_tab['MJD'], iers_tab['UT1_UTC'])
        ut1=epoch+datetime.timedelta(seconds=ut1utc)
    except:
        warnings.warn('Cannot calculate UT1: using UTC...')

    return ut1

def getUT1UTC(epoch: datetime.datetime):
    """
    Converts gets UT1 at UTC epoch
    :return: ut1
    """
    iers_tab = pySat.GetIersTable()
    ut1utc = epoch
    mjd = UTC2MJD(epoch)
    try:
        ut1utc = np.interp(mjd, iers_tab['MJD'], iers_tab['UT1_UTC'])
        ut1utc = epoch + datetime.timedelta(seconds=ut1utc)
    except:
        warnings.warn('Cannot calculate UT1UTC: using UTC...')

    return ut1utc

