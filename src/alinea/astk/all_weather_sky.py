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

""" Implementation of the all weather directional sky luminance model.

Perez R, Seals R, Michalsky J (1993) All-weather model for sky luminance
distribution - preliminary configuration and validation Solar energy 50: 235-245.
"""
import numpy
import pandas
import pvlib

# some defaults
_day = '2000-06-21'
_dates = pandas.date_range(_day, periods=24, freq='H')
_timezone = 'Europe/Paris'
_longitude = 3.87
_latitude = 43.61
_altitude = 45


def sun_position(dates=_dates, latitude=_latitude, longitude=_longitude, altitude=_altitude, timezone=_timezone):
    """ Sun position

    Args:
        dates: a pandas.DatetimeIndex specifying the dates at which sun position
        is required
        latitude: float
        longitude: float
        altitude: (float) altitude in m
        timezone: a pytz identifier of the timezone associated to dates
        (default UTC)

    Returns:
        a pandas dataframe with sun position at requested dates indexed by
        localised dates
    """
    times = dates.tz_localize(timezone)
    df = pvlib.solarposition.get_solarposition(times, latitude, longitude,
                                                 altitude)
    df['UTC'] = times.tz_convert('UTC')

    return df.loc[df['apparent_elevation'] > 0,:]


def air_mass(sunpos, altitude=_altitude):
    """Estimate the pressure-corrected air mass
    (optical path length relative to zenital path path at the location)

    Args:
        sunpos : a pandas dataframe as returned by sun_position
        altitude : (float)
    """

    airmass = pvlib.atmosphere.relativeairmass(sunpos['apparent_zenith'])
    pressure = pvlib.atmosphere.alt2pres(altitude)
    return pvlib.atmosphere.absoluteairmass(airmass, pressure)


def clear_sky_irradiances(dates=_dates, longitude=_longitude, latitude=_latitude, altitude=_altitude, timezone=_timezone):
    df = sun_position(dates, longitude, latitude, altitude, timezone)
    df['am'] = air_mass(df, altitude)
    tl = pvlib.clearsky.lookup_linke_turbidity(df.index, latitude, longitude)
    df['dni_extra'] = pvlib.irradiance.extraradiation(df.index)
    df = pandas.concat([df, pvlib.clearsky.ineichen(df['apparent_zenith'], df['am'], tl, dni_extra = df['dni_extra'], altitude = altitude)], axis=1)
    return df.loc[:,['UTC', 'azimuth', 'apparent_zenith', 'apparent_elevation', 'am', 'dni_extra', 'ghi', 'dni', 'dhi' ]]


def actual_sky_irradiances(weather, dates):
    """ derive a sky irradiance dataframe from actual weather data"""
    loc = weather.localisation
    data = weather.data.loc[dates,:]
    df = sun_position(dates, loc['longitude'], loc['latitude'], loc['altitude'], loc['timezone'])
    df['am'] = air_mass(df, loc['altitude'])
    df['dni_extra'] = pvlib.irradiance.extraradiation(df.index)
    df['ghi'] = data['global']
    if 'dhi' not in data.columns:
        pressure = pvlib.atmosphere.alt2pres(altitude)
        df['dni'] = pvlib.irradiance.dirint(df['ghi'], df['apparent_zenith'], df.index, pressure=pressure, temp_dew=None)
    else:
        df['dhi'] = data['dhi']
        df['dni'] = (df['ghi'] - df['dhi']) / numpy.sin(df['apparent_elevation'])
    return df.loc[:,['UTC', 'azimuth', 'apparent_zenith', 'apparent_elevation', 'am', 'dni_extra', 'ghi', 'dni', 'dhi' ]]


def all_weather_sky(irradiances):
    """ Return parameters of the all_weather_sky_model

        Args:
            irradiances: a pandas dataframe similar to clear_sky_irradiance output

    """

    df = irradiances.loc[:,['UTC', 'azimuth', 'apparent_zenith', 'apparent_elevation']]
    dni = irradiances['dni']
    dhi = irradiances['dhi']
    z = numpy.radians(irradiances['apparent_zenith'])
    df['clearness'] = ((dhi + dni) / dhi + 1.041 * z**3) / (1 + 1.041 * z**3)
    df['brightness'] = irradiances['am'] * dhi / irradiances['dni_extra']
    elevation = numpy.radians(df['apparent_elevation'])
    df['parameters'] = map(lambda x: aw_parameters(*x), zip(elevation, df['clearness'],df['brightness']))

    return df


def angle(a, b):
    """ angle (rad) between vector a and b"""
    return numpy.arctan2(numpy.linalg.norm(numpy.cross(a, b)), numpy.dot(a, b))



def aw_parametric_fit(p1, p2, p3, p4, z, brightness):
    """ Estimate the parameter fitting function at (z, brightness)

    Args:
        p1 to p4: parameter 1 to 4 of the function
        z: solar zenith angle
        brightness: sky brightness

    Returns:
        the value of the parameter fitting function

    """
    p1 = numpy.array(p1)
    p2 = numpy.array(p2)
    p3 = numpy.array(p3)
    p4 = numpy.array(p4)

    return p1 + p2 * z + brightness * (p3 + p4 * z)


def aw_parameters(sun_elevation, clearness, brightness):
    """ parameters of the relative luminance equation in the z, cleaness,
    brihghtness space

    Args:

    Returns: the five coefficient of the relative luminance equation (a, b, c, d
    and e)

    """

    bins = [1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2]
    a1 = (1.3525, -1.2219, -1.1000, -0.5484, -0.6000, -1.0156, -1.0000, -1.0500)
    a2 = (-0.2576, -0.7730, -0.2215, -0.6654, -0.3566, -0.3670, 0.0211, 0.0289)
    a3 = (-0.2690, 1.4148, 0.8952, -0.2672, -2.5000, 1.0078, 0.5025, 0.4260)
    a4 = (-1.4366, 1.1016, 0.0156, 0.7117, 2.3250, 1.4051, -0.5119, 0.3590)
    b1 = (-0.7670, -0.2054, 0.2782, 0.7234, 0.2937, 0.2875, -0.3000, -0.3250)
    b2 = (0.0007, 0.0367, -0.1812, -0.6219, 0.0496, -0.5328, 0.1922, 0.1156)
    b3 = (1.2734, -3.9128, -4.5000, -5.6812, -5.6812, -3.8500, 0.7023, 0.7781)
    b4 = (-0.1233, 0.9156, 1.1766, 2.6297, 1.8415, 3.3750, -1.6317, 0.0025)
    c1 = (2.8000, 6.9750, 24.7219, 33.3389, 21.0000, 14.0000, 19.0000, 31.0625)
    c2 = (
    0.6004, 0.1774, -13.0812, -18.3000, -4.7656, -0.9999, -5.0000, -14.5000)
    c3 = (
    1.2375, 6.4477, -37.7000, -62.2500, -21.5906, -7.1406, 1.2438, -46.1148)
    c4 = (1.0000, -0.1239, 34.8438, 52.0781, 7.2492, 7.5469, -1.9094, 55.3750)
    d1 = (1.8734, -1.5798, -5.0000, -3.5000, -3.5000, -3.4000, -4.0000, -7.2312)
    d2 = (0.6297, -0.5081, 1.5218, 0.0016, -0.1554, -0.1078, 0.0250, 0.4050)
    d3 = (0.9738, -1.7812, 3.9229, 1.1477, 1.4062, -1.0750, 0.3844, 13.3500)
    d4 = (0.2809, 0.1080, -2.6204, 0.1062, 0.3988, 1.5702, 0.2656, 0.6234)
    e1 = (0.0356, 0.2624, -0.0156, 0.4659, 0.0032, -0.0672, 1.0468, 1.5000)
    e2 = (-0.1246, 0.0672, 0.1597, -0.3296, 0.0766, 0.4016, -0.3788, -0.6426)
    e3 = (-0.5718, -0.2190, 0.4199, -0.0876, -0.0656, 0.3017, -2.4517, 1.8564)
    e4 = (0.9938, -0.4285, -0.5562, -0.0329, -0.1294, -0.4844, 1.4656, 0.5636)

    z = numpy.pi / 2. - sun_elevation
    index = max(0, numpy.searchsorted(bins, clearness) - 1)
    a = aw_parametric_fit(a1, a2, a3, a4, z, brightness)[index]
    b = aw_parametric_fit(b1, b2, b3, b4, z, brightness)[index]
    c = aw_parametric_fit(c1, c2, c3, c4, z, brightness)[index]
    d = aw_parametric_fit(d1, d2, d3, d4, z, brightness)[index]
    e = aw_parametric_fit(e1, e2, e3, e4, z, brightness)[index]

    if clearness <= 1.065:
        c = numpy.exp(numpy.power(brightness * (c1[0] + c2[0] * z), c3[0])) - \
            c4[0]
        d = -numpy.exp(brightness * (d1[0] + d2[0] * z)) + \
            d3[0] + d4[0] * brightness

    return a, b, c, d, e


def relative_luminance(theta=0, phi=0, aw_sky):
    """ Relative luminance of a sky element

    Args:
        theta: zenith angle (rad) of the sky element
        gamma: the angle (rad) between the sky element and the sun
        aw_sky : a pandas Dataframe with all weather sky model parameters and sun position

    Returns:
        the relative luminance

    """
    a, b, c, d, e = map(numpy.array, zip(*aw_sky['parameters']))
    az = numpy.radians(aw_sky['azimuth'])
    el = numpy.radians(aw_sky['apparent_elevation'])
    sun = zip(numpy.cos(az), numpy.sin(az), numpy.sin(el))
    element =(numpy.cos(phi), numpy.sin(phi), numpy.cos(theta))
    gamma = map(lambda x: abs(angle(element, x)), sun)

    return (1 + a * numpy.exp(b / numpy.cos(theta))) * (
    1 + c * numpy.exp(d * gamma) + e * numpy.cos(gamma) ** 2)


