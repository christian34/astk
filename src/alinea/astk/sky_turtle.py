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
""" A class interface to 'turtle' sky representation
"""

from alinea.astk.icosphere import turtle_dome, sample_faces, display
from alinea.astk.all_weather_sky import clear_sky_irradiances, all_weather_sky, relative_irradiance
from alinea.astk.colormap import jet_colors

class SkyTurtle(object):
    """ A class interface to 'turtle' sky representation
    """

    def __init__(self, discretisation_level=3):
        self.discretisation_level = discretisation_level
        self.vertices, self.faces = turtle_dome(discretisation_level)

    def face_samples(self, iter=None):
        """ return a list of theta, phi direction tuples sampling the turtle  faces"""
        return sample_faces(self.vertices, self.faces, iter=iter, spheric=True)

    def sample_irradiances(self, iter=None, sky_irradiances=None):
        if sky_irradiances is None:
            sky_irradiances = clear_sky_irradiances().iloc[[7],:]
        sky = all_weather_sky(sky_irradiances)
        directions, tags = self.face_samples(iter)
        irr_t = [relative_irradiance(t, p, sky) for t, p in directions]
        total_irr_t = map(sum, zip(*irr_t))
        irradiances = [sum(it / total_irr_t * sky['ghi']) for it in irr_t]
        return directions, irradiances, tags

    def display(self, property=None):
        colors = None
        if property is not None:
            colors = jet_colors(property)
        display(self.vertices, self.faces, colors)
