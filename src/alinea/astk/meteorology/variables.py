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

""" Name, description, units and synonyms of meterological variables
"""
variables = {
    # radiation
    'global_radiation': {'unit': 'W.m-2',
        'desc': 'hemispherical irradiance on horizontal surface',
        'synonym': ['global_horizontal_irradiance', 'GHI']},
    'direct_solar_radiation': {'unit': 'W.m-2',
        'desc': 'direct normal irradiance of the sun',
        'synonym': ['direct_normal_irradiance', 'DNI']},
    'diffuse_solar_radiation': {'unit': 'W.m-2',
        'desc': 'diffuse horizontal irradiance from the sky',
        'synonym': ['diffuse_horizontal_irradiance', 'DHI']},
    'PAR': {'unit':'W.m-2',
        'desc': 'Photosynthetically active radiation'},
    'PPFD': {'unit':'micromol.m-2.s-1',
        'desc': 'Photosynthetic Photon Flux Density'},

}
