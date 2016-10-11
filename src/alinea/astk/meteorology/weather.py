# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 14:29:15 2013

@author: lepse
"""

import pandas
import numpy
import pytz
from datetime import datetime, timedelta
from math import exp

from alinea.astk.TimeControl import *
import alinea.astk.sun_and_sky as sunsky

from alinea.astk.meteorology.localisations import localisations
from openalea.deploy.shared_data import shared_data
import alinea.septo3d
sample_meteo = shared_data(alinea.septo3d, 'meteo00-01.txt')
sample_location = localisations['Grignon']

def septo3d_reader(data_file):
    """ reader for septo3D meteo files """

    data = pandas.read_csv(data_file,sep='\t')
    year = data['An']
    dayofyear = data['Jour']
    hour = data['hhmm'] / 100
    # 24h -> 0h
    last_day = max(dayofyear)
    dayofyear = dayofyear + numpy.where(hour == 24, 1, 0)
    dayofyear = numpy.where(dayofyear > last_day, last_day, dayofyear)
    hour = hour - numpy.where(hour == 24, 24, 0)

    def _date(y, d, h):
        return datetime.strptime('%d %d %d' % (y, d, h), '%Y %j %H')
    date = map(lambda x: _date(*x), zip(year, dayofyear, hour))
    data = data.rename(columns={'PAR':'PPFD', 'Vent': 'wind_speed',
                                'Pluie': 'rain'})
    data.index = date
    data.index = data.index.tz_localize('UTC')
    return data.drop(['An', 'Jour', 'hhmm'],1)



def linear_degree_days(data, start_date=None, base_temp=0., max_temp=35.):
    df = data['temperature_air'].copy()
    if start_date is None:
        start_date = data.index[0]
    df[df < base_temp] = 0.
    df[df > max_temp] = 0.
    dd = numpy.zeros(len(df))
    if isinstance(start_date, str):
        start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
    ind_start = len(df.ix[:start_date + timedelta(0, 60)])
    seq = pandas.date_range(start=df.index[0],
                            end=start_date - timedelta(0, 60), freq='H')
    seq = seq.order(ascending=False)
    dd[ind_start:] = numpy.cumsum(
        (df.ix[start_date + timedelta(0, 60):].values - base_temp) / 24.)
    dd[:len(seq)] = -numpy.cumsum((df.ix[seq].values[::-1] - base_temp) / 24.)[
                     ::-1]
    return dd


def diffuse_fraction(data, localisation):
    """ Estimate the diffuse to global fraction 
    """
    heureTU = data.index.hour + data.index.minute / 60.
    DOY = data.index.dayofyear
    Rg = data['global_radiation']
    rdrs = sunsky.diffuse_fraction(Rg, heureTU, DOY,
                                   longitude=localisation['longitude'],
                                   latitude=localisation['latitude'])
    return rdrs


class Weather(object):
    """ Class compliying echap local_microclimate model protocol (meteo_reader).
        expected variables of the data_file are:
            - 'An'
            - 'Jour'
            - 'hhmm' : hour and minutes (universal time, UTC)
            - 'PAR' : Quantum PAR (ppfd) in micromol.m-2.sec-1
            - 'Pluie' : Precipitation (mm)
            - 'Tair' : Temperature of air (Celcius)
            - 'HR': Humidity of air (%)
            - 'Vent' : Wind speed (m.s-1)
        - localisation is a {'name':city, 'lontitude':lont, 'latitude':lat} dict
        - timezone indicates the standard timezone name (see pytz infos) to be used for interpreting the date (default 'UTC')
    """

    def __init__(self, data_file=sample_meteo, reader=septo3d_reader, wind_screen=2,
                 temperature_screen=2,
                 localisation=sample_location):
        self.data_path = data_file
        self.localisation = localisation
        self.latitude = localisation['latitude']
        self.longitude = localisation['longitude']
        self.timezone = localisation['timezone']

        if data_file is '':
            self.data = None
        else:
            self.data = reader(data_file)
            self.data.index = self.data.index.tz_convert(self.timezone)
        self.wind_screen = wind_screen
        self.temperature_screen = temperature_screen


    def date_range_index(self, start, end=None, by=24):
        """ return a (list of) time sequence that allow indexing one or several time intervals between start and end every 'by' hours
        if end is None, only one time interval of 'by' hours is returned
        
        start and end are expected in local time
        """
        if end is None:
            seq = pandas.date_range(start=start, periods=by, freq='H',
                                    tz=self.timezone.zone)
            return seq.tz_convert('UTC')
        else:
            seq = pandas.date_range(start=start, end=end, freq='H',
                                    tz=self.timezone.zone)
            seq = seq.tz_convert('UTC')
            bins = pandas.date_range(start=start, end=end, freq=str(by) + 'H',
                                     tz=self.timezone.zone)
            bins = bins.tz_convert('UTC')
            return [seq[(seq >= bins[i]) & (seq < bins[i + 1])] for i in
                    range(len(bins) - 1)]

    def get_weather(self, time_sequence):
        """ Return weather data for a given time sequence
        """
        return self.data.truncate(before=time_sequence[0],
                                  after=time_sequence[-1])

    def get_weather_start(self, time_sequence):
        """ Return weather data at start of timesequence
        """
        return self.data.truncate(before=time_sequence[0],
                                  after=time_sequence[0])

    def get_variable(self, what, time_sequence):
        """
        return values of what at date specified in time sequence
        """
        return self.data[what][time_sequence]

    def check(self, varnames=[], models={}, args={}):
        """ Check if varnames are in data and try to create them if absent using defaults models or models provided in arg.
        Return a bool list with True if the variable is present or has been succesfully created, False otherwise.
        
        Parameters: 
        
        - varnames : a list of name of variable to check
        - models a dict (name: model) of models to use to generate the data. models receive data as argument
        """

        models.update(self.models)

        check = []

        for v in varnames:
            if v in self.data.columns:
                check.append(True)
            else:
                if v in models.keys():
                    values = models[v](self.data, **args.get(v, {}))
                    self.data[v] = values
                    check.append(True)
                else:
                    check.append(False)
        return check

    def split_weather(self, time_step, t_deb, n_steps):

        """ return a list of sub-part of the meteo data, each corresponding to one time-step"""
        tdeb = pandas.date_range(t_deb, periods=1, freq='H')[0]
        tstep = [tdeb + i * timedelta(hours=time_step) for i in range(n_steps)]
        return [self.data.truncate(before=t,
                                   after=t + timedelta(hours=time_step - 1)) for
                t in tstep]

    def sun_path(self, seq, azimuth_origin='North'):
        """ Return position of the sun corresponing to a sequence of date
        """
        data = self.get_weather(seq)
        hUTC = data.index.hour + data.index.minute / 60.
        dayofyear = data.index.dayofyear
        latitude = self.localisation['latitude']
        longitude = self.localisation['longitude']
        data['sun_elevation'] = sunsky.sun_elevation(hUTC, dayofyear, longitude,
                                                     latitude)
        data['sun_azimuth'] = sunsky.sun_azimuth(hUTC, dayofyear, longitude,
                                                 latitude,
                                                 origin=azimuth_origin)
        data['sun_irradiance'] = sunsky.sun_irradiance(hUTC, dayofyear,
                                                       longitude, latitude)
        data['sun_ground_irradiance'] = sun_clear_sky_direct_normal_irradiance(
            hUTC, dayofyear, longitude, latitude, 'horizontal')
        return data

    def light_sources(self, seq, what='global_radiation', sky='turtle46',
                      azimuth_origin='North', irradiance='horizontal',
                      scale=1e-6):
        """ return direct and diffuse ligh sources representing the sky and the sun
         for a given time period indicated by seq
         Irradiance are accumulated over the whole time period and multiplied by the duration of the period (second) and by scale
        """
        self.check([what, 'diffuse_fraction'], args={
            'diffuse_fraction': {'localisation': self.localisation}})
        data = self.get_weather(seq)
        latitude = self.localisation['latitude']
        longitude = self.localisation['longitude']
        # expect the first interval to be like the second
        dt = numpy.diff((data.index - data.index[0]).seconds).astype(
            int).tolist()
        dt = numpy.array([dt[0]] + dt)

        hUTC = data.index.hour + data.index.minute / 60.
        dayofyear = data.index.dayofyear
        sun_elevation = sunsky.sun_elevation(hUTC, dayofyear, longitude,
                                             latitude)
        sun_azimuth = sunsky.sun_azimuth(hUTC, dayofyear, longitude, latitude,
                                         origin=azimuth_origin)
        sun_irradiance = data[what] * (
        1 - data['diffuse_fraction']) * dt * scale

        if irradiance == 'normal':
            sun_irradiance = sunsky.normal_irradiance(sun_irradiance,
                                                      sun_elevation)

        sky_elevation, sky_azimuth, sky_fraction = sunsky.sky_discretisation(
            type=sky)
        # to do : compute clear sky / diffuse sky depending on Rd/Rs
        # irradiance is supposed to be horizontal in meteo file
        sky_irradiance = sunsky.diffuse_light_irradiance(sky_elevation,
                                                         sky_azimuth,
                                                         sky_fraction,
                                                         sky_type='soc',
                                                         irradiance=irradiance)
        sky_irradiance *= (data[what] * data[
            'diffuse_fraction']).sum() * dt.sum() * scale

        sunny = numpy.where((sun_elevation > 0) & (sun_irradiance > 0))
        sun = {'elevation': sun_elevation[sunny], 'azimuth': sun_azimuth[sunny],
               'irradiance': numpy.array(sun_irradiance)[sunny]}
        sky = {'elevation': sky_elevation, 'azimuth': sky_azimuth,
               'irradiance': sky_irradiance}

        return sun, sky

    def daylength(self, seq):
        """
        """
        return sunsky.day_length(self.localisation['latitude'], seq.dayofyear)


def weather_node(weather_path):
    return Weather(weather_path)


def weather_check_node(weather, vars, models):
    ok = weather.check(vars, models)
    if not numpy.all(ok):
        print "weather_check: warning, missing  variables!!!"
    return weather


def weather_data_node(weather):
    return weather.data


def weather_start_node(timesequence, weather):
    return weather.get_weather_start(timesequence),


def date_range_node(start, end, periods, freq, tz, normalize,
                    name):  # nodemodule = pandas in wralea result in import errors
    return pandas.date_range(start, end, periods, freq, tz, normalize, name)


def sample_weather(periods=24):
    """ provides a sample weather instance for testing other modules
    """
    from openalea.deploy.shared_data import shared_data
    import alinea.septo3d

    meteo_path = shared_data(alinea.septo3d, 'meteo00-01.txt')
    t_deb = "2000-10-01 01:00:00"
    seq = pandas.date_range(start="2000-10-02", periods=periods, freq='H')
    weather = Weather(data_file=meteo_path)
    weather.check(
        ['temperature_air', 'PPFD', 'relative_humidity', 'wind_speed', 'rain',
         'global_radiation', 'vapor_pressure'])
    return seq, weather


def sample_weather_with_rain():
    seq, weather = sample_weather()
    every_rain = rain_filter(seq, weather)
    rain_timing = IterWithDelays(*time_control(seq, every_rain, weather.data))
    return rain_timing.next().value


def climate_todict(x):
    if isinstance(x, pandas.DataFrame):
        return x.to_dict('list')
    elif isinstance(x, pandas.Series):
        return x.to_dict()
    else:
        return x



        # def add_global_radiation(self):
        # """ Add the column 'global_radiation' to the data frame.
        # """
        # data = self.data
        # global_radiation = self.PPFD_to_global(data['PPFD'])
        # data = data.join(global_radiation)

        # def add_vapor_pressure(self, globalclimate):
        # """ Add the column 'global_radiation' to the data frame.
        # """
        # vapor_pressure = self.humidity_to_vapor_pressure(globalclimate['relative_humidity'], globalclimate['temperature_air'])
        # globalclimate = globalclimate.join(vapor_pressure)
        # mean_vapor_pressure = globalclimate['vapor_pressure'].mean()
        # return mean_vapor_pressure, globalclimate

        # def fill_data_frame(self):
        # """ Add all possible variables.

        # For instance, call the method 'add_global_radiation'.
        # """
        # self.add_global_radiation()

        # def next_date(self, timestep, t_deb):
        # """ Return the new t_deb after the timestep 
        # """
        # return t_deb + timedelta(hours=timestep)

#
# To do /add (pour ratp): 
# file meteo exemples
# add RdRs (ratio diffus /global)
# add NIR = RG - PAR
# add Ratmos = epsilon sigma Tair^4, epsilon = 0.7 clear sky, eps = 1 overcast sky
# add CO2
#
# peut etre aussi conversion hUTC -> time zone 'euroopean' 

##
# sinon faire des generateur pour tous les fichiers ratp
#
