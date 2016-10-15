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
from alinea.astk.meteorology.sky_irradiance import clear_sky_irradiances
from alinea.astk.meteorology.sky_luminance import directional_irradiances
from alinea.astk.colormap import jet_colors
import pandas
import numpy


class SkyTurtle(object):
    """ A class interface to 'turtle' sky representation
    """

    def __init__(self, discretisation_level=3):
        self.discretisation_level = discretisation_level
        self.vertices, self.faces = turtle_dome(discretisation_level)

    def face_samples(self, iter=None):
        """ return a list of theta, phi direction tuples sampling the turtle faces"""
        return sample_faces(self.vertices, self.faces, iter=iter, spheric=True)

    def directional_irradiances(self, model='all_weather', iter=None, add_sun=False,
                                sky_irradiances=None, ghi=(1,), zenith=(0,),
                                azimuth=(0,), clearness=(1,), brightness=(0.2,)):
        """

        Args:
            iter: number of face spliting iteration to perform to subsample faces
            add_sun : shoul sun directions (zenith/azimuth) be added to the turtle face sampling set ?
            sky_irradiances: a pandas dataframe with columns or a dict of list
            model:

        Returns:

        """
        if sky_irradiances is None:
            sky = {'ghi':ghi, 'azimuth': azimuth, 'zenith':zenith, 'clearness':clearness, 'brightness':brightness}
        else:
            sky = sky_irradiances
        directions, tags = self.face_samples(iter)
        irradiances = directional_irradiances(sky['ghi'], directions,
                                              model=model,
                                              sun_zenith=sky['zenith'],
                                              sun_azimuth=sky['azimuth'],
                                              clearness=sky['clearness'],
                                              brightness=sky['brightness'])
        df = pandas.DataFrame({'tag':tags, 'irradiance':irradiances})
        irr_faces = df.groupby('tag').agg(numpy.sum)
        # to do : choose argmax dir
        order = numpy.argsort(tags)
        return directions, irradiances,irr_faces, tags, order

    def sample_luminances(self, iter=None, sky_irradiances=None):
        directions, irradiances, tags, order = self.sample_irradiances(iter, sky_irradiances)
        theta, phi = zip(*directions)
        cost = numpy.where(numpy.cos(theta) <= 0, 0, numpy.cos(theta))
        return directions, irradiances * cost, tags, order


    def display(self, property=None):
        colors = None
        if property is not None:
            colors = jet_colors(property)
        display(self.vertices, self.faces, colors)
