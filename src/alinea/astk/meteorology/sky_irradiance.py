# -*- python -*-
#
#       Copyright 2016 INRIA - CIRAD - INRA
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       WebSite : https://github.com/openalea-incubator/astk
#
#       File author(s): Christian Fournier <Christian.Fournier@supagro.inra.fr>
#
# ==============================================================================

""" Equation for determining global horizontal irradiance (GHI),
direct normal irradiance (DNI) and diffuse horizontal irradiance under clearsky
condition or estimate them from meteorological data

This module is mainly a collection of syntactic sugar to pvlib clearsky and
irradiances packages.
"""
import numpy
import pandas
from alinea.astk.meteorology.sun_position import sun_position, \
    sun_extraradiation

try:
    import pvlib
except ImportError as e:
    raise ImportError(
        '{0}\npvlib not found on your system, you may use sun_position_astk '
        'instead OR install ephem and use sun_position_ephem OR install pvlib '
        '(recommended)'.format(e))

# default location and dates
_daydate = '2000-06-21'
_timezone = 'Europe/Paris'
_longitude = 3.52
_latitude = 43.36
_altitude = 56


def horizontal_irradiance(normal_irradiance, elevation):
    """ irradiance measured on an horizontal surface from a source
    with known elevation (degrees) and known normal irradiance
    """
    return normal_irradiance * numpy.sin(numpy.radians(elevation))


def normal_irradiance(horizontal_irradiance, elevation):
    """ irradiance measured on an surface perpendicular
    to a source with known elevation (degrees) and horizontal irradiance
    """
    return horizontal_irradiance / numpy.sin(numpy.radians(elevation))


def air_mass(zenith, altitude=0):
    """Estimate the pressure-corrected air mass
    (optical path length relative to zenital path at a location)

    Args:
        zenith : an array-like object of zenital directions (degrees)
        altitude : (float)
    """
    airmass = pvlib.atmosphere.get_relative_airmass(zenith)
    pressure = pvlib.atmosphere.alt2pres(altitude)
    am = pvlib.atmosphere.get_absolute_airmass(airmass, pressure)
    return am


def clearness(dni, dhi, sun_zenith):
    """Perez formula for clearness index

    Args:
        dni:
        dhi:
        sun_zenith:

    Returns:

    """
    z = numpy.radians(sun_zenith)
    return ((dhi + dni) / dhi + 1.041 * z**3) / (1 + 1.041 * z**3)


def brightness(air_mass, dhi, dni_extra):
    """perez formula for brightness index

    Args:
        air_mass:
        dhi:
        dni_extra:

    Returns:

    """
    return air_mass * dhi / dni_extra


def clear_sky_irradiances(dates=None, daydate=_daydate, longitude=_longitude,
                          latitude=_latitude, altitude=_altitude,
                          timezone=_timezone):
    """ Estimate component of sky irradiance for clear sky conditions

    Args:
        dates: A pandas datetime index (as generated by pandas.date_range). If
            None, daydate is used.
        daydate: (str) yyyy-mm-dd (not used if dates is not None).
        longitude: (float) in degrees
        latitude: (float) in degrees
        altitude: (float) in meter
        timezone:(str) the time zone (not used if dates are already localised)

    Returns:
        a pandas dataframe with global horizontal irradiance, direct normal
        irradiance and diffuse horizontal irradiance.

    Details:
        P. Ineichen and R. Perez, "A New airmass independent formulation for
        the Linke turbidity coefficient", Solar Energy, vol 73, pp. 151-157,
        2002
    """

    df = sun_position(dates=dates, daydate=daydate, latitude=latitude,
                      longitude=longitude, altitude=altitude,
                      timezone=timezone)

    tl = pvlib.clearsky.lookup_linke_turbidity(df.index, latitude,
                                               longitude)
    am = air_mass(df['zenith'], altitude)
    dni_extra = sun_extraradiation(df.index)
    clearsky = pvlib.clearsky.ineichen(df['zenith'], am, tl,
                                       dni_extra=dni_extra,
                                       altitude=altitude)
    clearsky = pandas.concat([df, clearsky], axis=1)

    return clearsky.loc[:, ['ghi', 'dni', 'dhi']]


def actual_sky_irradiances(dates=None, daydate=_daydate, ghi=None,
                           attenuation=None,
                           pressure=101325, temp_dew=None, longitude=_longitude,
                           latitude=_latitude, altitude=_altitude,
                           timezone=_timezone):
    """ Estimate component of sky irradiances from measured actual global
    horizontal irradiance or attenuated clearsky conditions.

    Args:
        dates: A pandas datetime index (as generated by pandas.date_range). If
            None, daydate is used.
        daydate: (str) yyyy-mm-dd (not used if dates is not None).
        ghi: (array_like) : global horizontal irradiance (W. m-2).If None
         (default) clear_sky irradiance are used
        attenuation: (float) a attenuation factor for ghi (actual_ghi =
         attenuation * ghi). If None (default), no attenuation is applied
        pressure: the site pressure (Pa) (for dirint model)
        temp_dew: the dew point temperature (dirint model)
        longitude: (float) in degrees
        latitude: (float) in degrees
        altitude: (float) in meter
        timezone:(str) the time zone (not used if dates are already localised)

    Returns:
        a pandas dataframe with global horizontal irradiance, direct normal
        irradiance and diffuse horizontal irradiance.

    Details:
        Perez, R., P. Ineichen, E. Maxwell, R. Seals and A. Zelenka, (1992).
        Dynamic Global-to-Direct Irradiance Conversion Models.
        ASHRAE Transactions-Research Series, pp. 354-369
    """

    df = sun_position(dates=dates, daydate=daydate, latitude=latitude,
                      longitude=longitude, altitude=altitude,
                      timezone=timezone)

    if ghi is None:
        cs = clear_sky_irradiances(dates=df.index, latitude=latitude,
                                   longitude=longitude, altitude=altitude,
                                   timezone=timezone)
        ghi = cs['ghi']

    df['ghi'] = ghi

    if attenuation is not None:
        df.ghi *= attenuation

    df['dni'] = pvlib.irradiance.dirint(df.ghi, 90 - df.elevation, df.index,
                                        pressure=pressure, temp_dew=temp_dew)
    df['dhi'] = df.ghi - horizontal_irradiance(df.dni, df.elevation)

    return df.loc[:, ('ghi', 'dhi', 'dni')]


def sky_irradiances(dates=None, daydate=_daydate, ghi=None, dhi=None,
                           attenuation=None,
                           pressure=101325, temp_dew=None, longitude=_longitude,
                           latitude=_latitude, altitude=_altitude,
                           timezone=_timezone):
    """ Estimate variables related to sky irradiance.

    Args:
        dates: A pandas datetime index (as generated by pandas.date_range). If
            None, daydate is used.
        daydate: (str) yyyy-mm-dd (not used if dates is not None).
        ghi: (array_like) : global horizontal irradiance (W. m-2).If None
         (default) clear_sky irradiance are used
        dhi: (array-like): diffuse horizontal irradiance
        attenuation: (float) a attenuation factor for ghi (actual_ghi =
         attenuation * ghi). If None (default), no attenuation is applied. If
         dhi is not None, this parameter is not taken into account.
        pressure: the site pressure (Pa) (for dirint model)
        temp_dew: the dew point temperature (dirint model)
        longitude: (float) in degrees
        latitude: (float) in degrees
        altitude: (float) in meter
        timezone:(str) the time zone (not used if dates are already localised)

    Returns:
        a pandas dataframe with azimuth, zenital and elevation angle of the sun,
        clearness and brightness indices, global horizontal irradiance, direct
        normal irradiance and diffuse horizontal irradiance of the sky.
    """

    df = sun_position(dates=dates, daydate=daydate, latitude=latitude,
                      longitude=longitude, altitude=altitude,
                      timezone=timezone)

    if ghi is None or dhi is None:
        irr = actual_sky_irradiances(dates=df.index, ghi=ghi,
                                     attenuation=attenuation, pressure=pressure,
                                     temp_dew=temp_dew, latitude=latitude,
                                     longitude=longitude, altitude=altitude,
                                     timezone=timezone)
        df = pandas.concat([df, irr], axis=1)
    else:
        df['ghi'] = ghi
        df['dhi'] = dhi
        df['dni'] = normal_irradiance(numpy.array(ghi) - numpy.array(dhi),
                                      df.elevation)

    am = air_mass(df['zenith'], altitude)
    dni_extra = sun_extraradiation(df.index)
    df['brightness'] = brightness(am, df['dhi'], dni_extra)
    df['clearness'] = clearness(df['dni'], df['dhi'], df['zenith'])

    # twilight conditions (sun_el < 0, ghi > 0)
    if len(df) < 1 and ghi is not None:
        df = sun_position(dates=dates, daydate=daydate, latitude=latitude,
                          longitude=longitude, altitude=altitude,
                          timezone=timezone, filter_night=False)
        df['ghi'] = ghi
        df['dhi'] = ghi
        df['dni'] = 0
        df['clearness'] = None
        df['brightness'] = None
        df = df.loc[df.ghi > 0, :]

    return df.loc[:,
           ['azimuth', 'zenith', 'elevation', 'clearness', 'brightness', 'ghi',
            'dni', 'dhi']]
