# -*- python -*-
#
#       Copyright 2016 INRIA - CIRAD - INRA
#
#       Distributed under the Cecill-C License.
#       See accompanying file LICENSE.txt or copy at
#           http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
#
#       WebSite : https://github.com/openalea-incubator/astk/meteorology
#
#       File author(s): Christian Fournier <Christian.Fournier@supagro.inra.fr>
#
# ==============================================================================

""" Directional sky luminance models from CIE standard and Perez all weather model

CIE, 2002, Spatial distribution of daylight CIE standard general sky,
    CIE standard, CIE Central Bureau, Vienna
Perez R, Seals R, Michalsky J (1993) All-weather model for sky luminance
distribution - preliminary configuration and validation Solar energy 50: 235-245.
"""

import numpy

# CIE relative luminance distribution

def cie_luminance_gradation(theta, a, b):
    """ function giving the dependance of the relative luminance of a sky element
    to its elevation angle

    CIE, 2002, Spatial distribution of daylight CIE standard general sky,
    CIE standard, CIE Central Bureau, Vienna

    theta: zenital angle of the sky element
    a, b : coefficient for the type of sky
    """
    z = numpy.array(theta)
    phi_0 = 1 + a * numpy.exp(b)
    phi_z = numpy.where(numpy.cos(z) == 0, 1,
                        1 + a * numpy.exp(b / numpy.cos(z)))
    return phi_z / phi_0


def cie_scattering_indicatrix(theta, phi, sun_zenith, sun_azimuth,
                                c, d, e):
    """ function giving the dependance of the relative luminance of a sky element
    to its azimuth distance to the sun

    CIE, 2002, Spatial distribution of daylight CIE standard general sky,
    CIE standard, CIE Central Bureau, Vienna

    elevation : elevation angle of the sky element (rad)
    d, e : coefficient for the type of sky
    """
    z = numpy.array(theta)
    zs = numpy.radians(numpy.array(sun_zenith))
    alpha = numpy.array(phi)
    alpha_s = numpy.radians(numpy.array(sun_azimuth))
    ksi = numpy.arccos(
        numpy.cos(zs) * numpy.cos(z) + numpy.sin(zs) * numpy.sin(z) * numpy.cos(
            numpy.abs(alpha - alpha_s)))

    f_ksi = 1 + c * (
        numpy.exp(d * ksi) - numpy.exp(d * numpy.pi / 2)) + e * numpy.power(
        numpy.cos(ksi), 2)
    f_zs = 1 + c * (
        numpy.exp(d * zs) - numpy.exp(d * numpy.pi / 2)) + e * numpy.power(
        numpy.cos(zs), 2)

    return f_ksi / f_zs


def cie_relative_luminance(theta, phi=None, sun_zenith=None,
                           sun_azimuth=None, type='soc'):
    """ cie relative luminance of a sky element relative to the luminance
    at zenith

    angle in radians
    type is one of 'soc' (standard overcast sky), 'uoc' (uniform luminance)
    , 'clear_sky' (standard clear sky low turbidity) or 'blended' (mixture of overcast and clearsky)
    for blended CIE skies, mixing function is after Mardaljevic
    """

    if type == 'clear_sky' and (
                        sun_zenith is None or sun_azimuth is None or phi is None):
        raise ValueError, 'Clear sky requires sun position'

    if type == 'soc':
        return cie_luminance_gradation(theta, 4, -0.7)
    elif type == 'uoc':
        return cie_luminance_gradation(theta, 0, -1)
    elif type == 'clear_sky':
        return cie_luminance_gradation(theta, -1,
                                       -0.32) * cie_scattering_indicatrix(
            theta, phi, sun_zenith, sun_azimuth, 10, -3,
            0.45)
    elif type == 'blended':
        # f_clear = min(1, (clearness - 1) / (1.41 - 1))
        raise ValueError('not yet implemented')
    else:
        raise ValueError, 'Unknown sky type'


# all weather sky model


def aw_parameters(z, clearness, brightness):
    """ parameters of the all weather sky relative luminance equation in the z, cleaness,
    brihghtness space

    Args:
        z:
        clearness:
        brightness:

    Returns:

    """

    def _awfit(p1, p2, p3, p4, zen, br):
        p1 = numpy.array(p1)
        p2 = numpy.array(p2)
        p3 = numpy.array(p3)
        p4 = numpy.array(p4)
        return p1 + p2 * zen + br * (p3 + p4 * zen)

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

    index = max(0, numpy.searchsorted(bins, clearness) - 1)
    a = _awfit(a1, a2, a3, a4, z, brightness)[index]
    b = _awfit(b1, b2, b3, b4, z, brightness)[index]
    c = _awfit(c1, c2, c3, c4, z, brightness)[index]
    d = _awfit(d1, d2, d3, d4, z, brightness)[index]
    e = _awfit(e1, e2, e3, e4, z, brightness)[index]

    if clearness <= 1.065:
        c = numpy.exp(numpy.power(brightness * (c1[0] + c2[0] * z), c3[0])) - \
            c4[0]
        d = -numpy.exp(brightness * (d1[0] + d2[0] * z)) + d3[0] + d4[
                                                                       0] * brightness

    return a, b, c, d, e


def angle(a, b):
    """ angle (rad) between vector a and b"""
    return numpy.arctan2(numpy.linalg.norm(numpy.cross(a, b)), numpy.linalg.norm(numpy.dot(a, b)))


def cartesian(theta, phi):
    return numpy.sin(theta) * numpy.cos(phi), numpy.sin(theta) * numpy.sin(phi), numpy.cos(theta)

def aw_relative_luminance(theta, phi, sun_zenith, sun_azimuth, clearness, brightness):
    """

    Args:
        theta: (float)
        phi: (float)
        sun_elevation: (array)
        sun_azimuth: (array)
        clearness: (array)
        brightness: (array)

    Returns:

    """
    # to do normalise by zenith for consistency with cie
    z = numpy.radians(sun_zenith)
    az = numpy.radians(sun_azimuth)
    pars = numpy.frompyfunc(aw_parameters, 3, 5)
    a, b, c, d, e = map(lambda x: numpy.array(map(float, x)), pars(z, clearness, brightness))
    sun = zip(*cartesian(z, az))
    element = cartesian(theta,phi)
    gamma = map(lambda x: abs(angle(element, x)), sun)

    def _lum(t, g):
        return (1 + a * numpy.exp(b / numpy.cos(t))) * (
            1 + c * numpy.exp(d * g) + e * numpy.cos(g) ** 2)

    return _lum(theta, gamma) / _lum(0, numpy.abs(z))


def directional_irradiances(theta, phi, ghi=1, model='cie_soc', sun_zenith=None, sun_azimuth=None, clearness=None, brightness=None):
    """

    Args:
        theta: (array)
        phi: (array)
        ghi:
        model:
        sun_zenith:
        sun_azimuth:
        clearness:
        brightness:

    Returns:

    """
    theta = numpy.array(theta)
    phi = numpy.array(phi)

    if model == 'all_weather':
        if any(map(lambda x: x is None, [ghi, sun_zenith, sun_azimuth, clearness, brightness])):
            raise ValueError('All weather model needs all inputs to be set')
        # relatives_irradiances(t) for all directions
        r_irrt = [aw_relative_luminance(t,p,sun_zenith, sun_azimuth,clearness,brightness) * numpy.cos(t) for t, p in zip(theta,phi)]
    elif model == 'cie_soc':
        r_irrt = [[cie_relative_luminance(t) * numpy.cos(t)] * len(ghi) for t in theta]
    else:
        raise ValueError('unknown model')
    # sky integrals_of relative_irradiances(t)
    sky_irrt = map(sum, zip(*r_irrt))
    # irradiances
    irradiances = [sum(r / sky_irrt * ghi) for r in r_irrt]
    return irradiances


