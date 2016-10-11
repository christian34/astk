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

""" Name, description, units,synonyms and simple converters between  meterological variables
"""

import numpy


# converters

# PAR / global
# 1 WattsPAR.m-2 = 4.6 micromol.m-2.sec-1, 1 Wglobal = 0.48 WattsPAR
# 1 W = 1 J.s-1


def PAR_to_global(PAR):
    return PAR / 0.48


def PPFD_to_global(PPFD):
    return PPFD / 4.6 / 0.48


def diffuse_PAR_to_diffuse(PAR, diffuse_PAR):
    beam_fraction = 1 - diffuse_PAR / PAR
    return diffuse_PAR / (0.48 - 0.48 * beam_fraction ** 4)


def global_to_PAR(global_radiation):
    return global_radiation * 0.48


def global_to_PPFD(global_radiation):
    return global_radiation * 0.48 * 4.6


def PPFD_to_PAR(PPFD):
    return PPFD / 4.6


def PAR_to_PPFD(PAR):
    return PAR * 4.6


def diffuse_to_diffuse_PAR(global_radiation, diffuse_radiation):
    beam_fraction = 1 - diffuse_radiation / global_radiation
    return diffuse_radiation * (0.48 - 0.48 * beam_fraction ** 4)


def Psat(T):
    """ Saturating water vapor pressure (kPa) at temperature T (Celcius) with Tetens formula
    """
    return 0.6112 * numpy.exp(17.62 * T / (243.12 + T))


def humidity_to_vapor_pressure(humidity, Tair):
    """ Convert the relative humidity (%) in water vapor pressure (kPa)
    """
    return humidity / 100. * Psat(Tair)


def humidity_to_Tdew(humidity, Tair):
    """ dew point temperature"""
    e = humidity / 100. * Psat(Tair)
    ep = e / 0.6112
    return 243.12 * numpy.log(ep) / (17.62 - numpy.log(ep))


variables = {
    # temperature
    'air_temperature': {'unit': 'degree Celcius',
                        'desc': 'Air temperature at screen height',
                        'synonym': ['Tair']},
    # radiation
    'global_radiation': {'unit': 'W.m-2',
                         'desc': 'hemispherical irradiance on horizontal surface',
                         'synonym': ['global_horizontal_irradiance', 'GHI'],
                         'convert': {'PAR': PAR_to_global,
                                     'PPFD': PPFD_to_global}},
    'diffuse_radiation': {'unit': 'W.m-2',
                          'desc': 'diffuse horizontal irradiance from the sky',
                          'synonym': ['diffuse_horizontal_irradiance', 'DHI'],
                          'convert': {
                              ('PAR', 'diffuse_PAR'): diffuse_PAR_to_diffuse},
                          'PAR': {'unit': 'W.m-2',
                                  'desc': 'Photosynthetically active radiation (direct and diffuse'},
                          'convert': {'global_radiation': global_to_PAR,
                                      'PPFD': PPFD_to_PAR}},
    'PPFD': {'unit': 'micromol.m-2.s-1',
             'desc': 'Photosynthetic Photon Flux Density (direct and diffuse)',
             'convert': {'global_radiation': global_to_PPFD,
                         'PAR': PAR_to_PPFD}},
    'diffuse_PAR': {'unit': 'W.m-2',
                    'desc': 'Photosynthetically active radiation of the sky',
                    'convert': {('global_radiation', 'diffuse_radiation'): diffuse_to_diffuse_PAR}},
    # humidity
    'relative_humidity': {'unit': 'percent',
                          'desc': 'relative humidity of the air at screen height',
                          'synonym': ['HR']},
    'vapor_pressure': {'unit': 'kPa', 'desc': 'vapor pressure of the air',
                       'synonym': [],
                       'convert': {('HR', 'Tair'): humidity_to_vapor_pressure}},
    'dew_point_temperature': {'unit': 'degrees Celcius',
                              'desc': 'dew point temperature of the air',
                              'synonym': ['Tdew'],
                              'convert': {('HR', 'Tair'): humidity_to_Tdew}},
    # wind
    'wind_speed': {'unit': 'm.s-1', 'desc': 'velocity of wind'},
    # rain
    'rain': {'unit': 'mm', 'desc': 'precipitations'}, }


def meteorological_variables():
    """ return an expanded (one entry for all synonims) dict of variables"""

    expanded = variables.copy()
    for k in variables:
        if 'synonym' in variables[k]:
            for s in variables[k]['synonym']:
                d = {kk: vv for kk, vv in variables[k].iteritems()}
                d['synonym'] = [k] + [ss for ss in variables[k]['synonym'] if
                                      ss != s]
                expanded[s] = d

    return expanded
