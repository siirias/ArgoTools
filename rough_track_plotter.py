import pandas as pd
import folium
from datetime import datetime
from geopy.distance import geodesic
import requests
import geopandas as gpd
from shapely.geometry import Point
import os
import zipfile

# === Raw Data (Same as before) ===
raw_data = """
f8906.000.20250609T071932.science_log.csv:GPS,20250609T071758,58.995140,20.547571,4
f8906.000.20250609T084554.science_log.csv:GPS,20250609T081941,58.994282,20.544962,4
f8906.000.20250609T092956.science_log.csv:GPS,20250609T092834,58.993828,20.543348,7
f8906.000.20250609T094036.science_log.csv:GPS,20250609T093904,58.993713,20.543554,6
f8906.000.20250609T095118.science_log.csv:GPS,20250609T094943,58.993698,20.543512,6
f8906.000.20250609T100432.science_log.csv:GPS,20250609T100028,58.993607,20.543394,6
f8906.000.20250609T102314.science_log.csv:GPS,20250609T101545,58.993515,20.543579,6
f8906.000.20250609T103710.science_log.csv:GPS,20250609T103428,58.993336,20.543974,5
f8906.000.20250609T104740.science_log.csv:GPS,20250609T104601,58.993343,20.544489,4
f8906.000.20250609T110536.science_log.csv:GPS,20250609T110031,58.993244,20.545095,9
f8906.000.20250609T111628.science_log.csv:GPS,20250609T111519,58.993191,20.546087,9
f8906.000.20250609T114758.science_log.csv:GPS,20250609T114418,58.992973,20.548231,9
f8906.000.20250609T120418.science_log.csv:GPS,20250609T120224,58.992886,20.550064,8
f8906.000.20250609T141020.science_log.csv:GPS,20250609T140801,58.991405,20.565422,8
f8906.000.20250609T161518.science_log.csv:GPS,20250609T161306,58.991814,20.581919,5
f8906.000.20250609T184942.science_log.csv:GPS,20250609T181901,58.994259,20.598503,11
f8906.000.20250609T205714.science_log.csv:GPS,20250609T205525,58.998875,20.612928,4
f8906.000.20250609T230320.science_log.csv:GPS,20250609T230131,59.004341,20.624578,11
f8906.000.20250610T010824.science_log.csv:GPS,20250610T010611,59.009457,20.637844,6
f8906.000.20250610T031452.science_log.csv:GPS,20250610T031204,59.013027,20.652676,9
f8906.000.20250610T052100.science_log.csv:GPS,20250610T051852,59.014885,20.669643,5
f8906.000.20250610T072730.science_log.csv:GPS,20250610T072428,59.015831,20.683870,7
f8906.000.20250610T093528.science_log.csv:GPS,20250610T093322,59.019611,20.699446,4
f8906.000.20250610T114350.science_log.csv:GPS,20250610T114025,59.025356,20.717291,9
f8906.000.20250610T134844.science_log.csv:GPS,20250610T134653,59.032372,20.734741,4
f8906.000.20250610T155450.science_log.csv:GPS,20250610T155255,59.041504,20.749443,10
f8906.000.20250610T180518.science_log.csv:GPS,20250610T180321,59.053123,20.763247,5
f8906.000.20250610T201440.science_log.csv:GPS,20250610T201251,59.064575,20.777031,10
f8906.000.20250610T222022.science_log.csv:GPS,20250610T221756,59.075489,20.789154,5
f8906.000.20250611T002702.science_log.csv:GPS,20250611T002501,59.087444,20.800179,7
f8906.000.20250611T025426.science_log.csv:GPS,20250611T023232,59.098511,20.811762,10
"""  # replace this with full text or keep as before

# === Parse data ===
lines = [line for line in raw_data.strip().split('\n') if line]
data = []
for line in lines:
    parts = line.split(',')
    timestamp = datetime.strptime(parts[1], '%Y%m%dT%H%M%S')
    lat, lon, hdop = float(parts[2]), float(parts[3]), int(parts[4])
    data.append((timestamp, lat, lon, hdop))

df = pd.DataFrame(data, columns=['time', 'lat', 'lon', 'hdop'])
df = df.sort_values('time').reset_index(drop=True)

# === Speed calculations ===
distances_km = []
time_deltas_h = []
for i in range(1, len(df)):
    p1 = (df.lat[i - 1], df.lon[i - 1])
    p2 = (df.lat[i], df.lon[i])
    distance_km = geodesic(p1, p2).km
    delta_h = (df.time[i] - df.time[i - 1]).total_seconds() / 3600.0
    distances_km.append(distance_km)
    time_deltas_h.append(delta_h)

# Handle any divide-by-zero
speeds_knots = [d / t * 0.539957 if t > 0 else 0 for d, t in zip(distances_km, time_deltas_h)]

average_speed = sum(speeds_knots) / len(speeds_knots)
latest_speed = speeds_knots[-1] if speeds_knots else 0

print(f"Average speed: {average_speed:.2f} knots")
print(f"Latest segment speed: {latest_speed:.2f} knots")

# === Create map ===
m = folium.Map(location=[df.lat.iloc[0], df.lon.iloc[0]], zoom_start=6, tiles="OpenStreetMap")

# Add polyline
folium.PolyLine(df[['lat', 'lon']].values, color='blue', weight=3).add_to(m)

# Add markers
for i, row in df.iterrows():
    folium.CircleMarker(
        location=(row['lat'], row['lon']),
        radius=1,
        color='red',
        fill=True,
        fill_opacity=0.7,
        popup=row['time'].strftime("%Y-%m-%d %H:%M")
    ).add_to(m)


folium.raster_layers.WmsTileLayer(
    url='https://geo.vliz.be/geoserver/MarineRegions/wms',
    layers='MarineRegions:eez',
    name='EEZ (Marine Regions)',
    fmt='image/png',
    transparent=True,
    overlay=True,
    control=True
).add_to(m)

folium.LayerControl().add_to(m)

# Save map
m.save("float_track_with_eez.html")
print("Map saved to float_track_with_eez.html")
