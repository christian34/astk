""" Example script for reading / manipulating meteo data of mangosim
"""

import pandas
import openalea.plantgl.all as pgl
from alinea.astk.sun_and_sky import sun_sky_sources
from alinea.caribu.CaribuScene import CaribuScene
from alinea.caribu.light import light_sources


from vplants.fractalysis.light.directLight import scene_irradiance

def reader(data_file='rayostpierre2002.csv'):
    """ reader for mango meteo files """
    


    data = pandas.read_csv(data_file, parse_dates=['Date'],
                               delimiter = ';',
                               usecols=['Date','Rayonnement','Temperature_Air','HR'], dayfirst=True)
    data = data.rename(columns={'Date':'date',
                                 'Rayonnement':'global_radiation',
                                 'Temperature_Air':'temperature_air',
                                 'HR':'relative_humidity'})
    # convert J.cm2.h-1 to W.m-2
    data['global_radiation'] *= (10000. / 3600)
    index = pandas.DatetimeIndex(data['date']).tz_localize('Indian/Reunion')
    data = data.set_index(index)
    return data
  
# a strange mango tree
mango = pgl.Scene('fruitstructure.bgeom')
pgl.Viewer.display(mango)
localisation={'latitude':-21.32, 'longitude':55.5, 'timezone': 'Indian/Reunion'}

meteo = reader()

dates = pandas.date_range(start='2002-09-02', end = '2002-09-15', freq='H',tz='Indian/Reunion')
ghi = meteo.loc[dates,'global_radiation']

sun, sky = sun_sky_sources(ghi=ghi, dates=dates, **localisation)

# Caribu
import time
t = time.time()
light = light_sources(*sun) + light_sources(*sky)
cs = CaribuScene(mango, light=light, scene_unit='cm')
raw, agg = cs.run(direct=True, simplify=True)
res = pandas.DataFrame(agg)
print 'made in', time.time() - t

# #muslim
t = time.time()
# permute az, el, irr to el, az, irr
sun_m = sun[1], sun[0], sun[2]
sky_m = sky[1], sky[0], sky[2]
directions = zip(*sun_m) + zip(*sky_m)
resm =scene_irradiance(mango, directions, horizontal=True, scene_unit='cm')
print 'made in', time.time() - t


#ratp


# compare
res = res.rename(columns={'area': 'area_c'})
df = res.join(resm)
df.plot('irradiance', 'Ei',xlim=(0,35000), ylim=(0,35000), style='p')
df.plot('area', 'area_c',xlim=(0,0.2), ylim=(0,0.2), style='p')






