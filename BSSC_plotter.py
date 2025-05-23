import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
from cartopy.geodesic import Geodesic
from datetime import datetime
from geopy.distance import geodesic
from cmocean import cm

CLOSE_FIGURES = False  # Set to False to keep figures open after saving, True to close.
PLOT_MAPS = True
PLOT_PROFILES = True
def load_profiles_within_radius(directory, target_lat, target_lon, radius_nm, variable):
    """
    Load NetCDF profiles from files in `directory` within `radius_nm` nautical miles from
    (target_lat, target_lon), and extract data for the specified variable.

    Returns a DataFrame with columns: time, depth, value, profile_id
    """
    radius_km = radius_nm * 1.852  # Convert nautical miles to kilometers
    files = glob.glob(os.path.join(directory, "*.nc"))
    records = []

    for file in files:
        try:
            ds = xr.open_dataset(file)
            lat = ds['LATITUDE'].values
            lon = ds['LONGITUDE'].values
            time = ds['TIME'].values
            pres = ds['PRES'].values
            var_data = ds[variable].values

            for i in range(len(time)):
                point = (lat[i], lon[i])
                dist = geodesic((target_lat, target_lon), point).km
                if dist <= radius_km:
                    for j in range(pres.shape[1]):
                        pval = pres[i, j]
                        vval = var_data[i, j]
                        if np.isfinite(pval) and np.isfinite(vval):
                            records.append({
                                'time': pd.to_datetime(str(time[i])),
                                'depth': pval,
                                'value': vval,
                                'profile_id': f"{os.path.basename(file)}_{i}"
                            })
            ds.close()
        except Exception as e:
            print(f"Skipping {file}: {e}")

    df = pd.DataFrame(records)
    if df.empty or "time" not in df.columns:
        return None    
    return df.sort_values("time")

def pick_colormap(variable):
        # Choose a colormap based on the variable name
    if "TEMP" in variable.upper():
        cmap = cmocean.cm.thermal
    elif "PSAL" in variable.upper():
        cmap = cmocean.cm.haline
    elif "DOX2" in variable.upper():
        cmap = cmocean.cm.oxy
    else:
        cmap = 'viridis'  # fallback
    return cmap

def plot_profiles(df, variable, save_path=None, depth_bins = 400):
    """
    Plot profiles from dataframe (time, depth, value) using pcolormesh-style plotting.
    Masks values beyond the deepest actual measurement in each profile.
    Optionally saves the figure if `save_path` is given.
    """
    depths = np.linspace(df['depth'].min(), df['depth'].max(), depth_bins)
    times = pd.to_datetime(sorted(df['time'].unique()))

    time_labels = []
    value_grid = []

    for time in times:
        profile = df[df['time'] == time].sort_values("depth")
        if len(profile) < 5:
            continue

        profile_depths = profile['depth'].values
        profile_values = profile['value'].values

        # Interpolate only within actual profile range
        interp_values = np.full_like(depths, np.nan, dtype=np.float32)
        valid = (depths >= profile_depths.min()) & (depths <= profile_depths.max())
        interp_values[valid] = np.interp(depths[valid], profile_depths, profile_values)

        value_grid.append(interp_values)
        time_labels.append(time)

    value_grid = np.array(value_grid).T  # shape: depth x time

    plt.figure(figsize=(12, 6))
    cmap = pick_colormap(variable)
    im = plt.pcolormesh(time_labels, depths, value_grid, shading='auto', cmap=cmap)
    plt.colorbar(im, label=variable)
    plt.gca().invert_yaxis()
    plt.xlabel("Time")
    plt.ylabel("Depth [dbar]")
    plt.title(f"{variable} profiles over time")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Profile plot saved to: {save_path}")
    else:
        plt.show()

    if CLOSE_FIGURES:
        plt.close()


def plot_profile_locations_map(directory, target_lat, target_lon, radius_nm, save_path=None):
    """
    Plots a map of all profile locations. Highlights those within radius.
    Optionally saves the figure if `save_path` is given.
    """
    radius_km = radius_nm * 1.852
    files = glob.glob(os.path.join(directory, "*.nc"))
    locations = []
    accepted = []

    for file in files:
        try:
            ds = xr.open_dataset(file)
            lat = ds['LATITUDE'].values
            lon = ds['LONGITUDE'].values
            for i in range(len(lat)):
                pt = (lat[i], lon[i])
                dist = geodesic((target_lat, target_lon), pt).km
                locations.append(pt)
                accepted.append(dist <= radius_km)
            ds.close()
        except Exception as e:
            print(f"Skipping {file}: {e}")

    if not locations:
        print("No locations found.")
        return

    lats, lons = zip(*locations)
    proj = ccrs.AzimuthalEquidistant(central_longitude=target_lon, central_latitude=target_lat)

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=proj)

    margin = 0.5
    extent = [min(lons)-margin, max(lons)+margin, min(lats)-margin, max(lats)+margin]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.gridlines(draw_labels=True)

    for (lat, lon), ok in zip(locations, accepted):
        color = 'lightblue' if ok else 'darkred'
        ax.plot(lon, lat, marker='o', color=color, markersize=2, transform=ccrs.PlateCarree())

    geod = Geodesic()
    circle = geod.circle(lon=target_lon, lat=target_lat, radius=radius_km * 1000, n_samples=180)
    circle_lons, circle_lats = zip(*circle)
    ax.plot(circle_lons, circle_lats, color='black', linewidth=1, transform=ccrs.PlateCarree(), label=f"{radius_nm} NM radius")
    ax.plot(target_lon, target_lat, marker='x', color='black', markersize=8, transform=ccrs.PlateCarree())

    ax.set_title(f"Profile locations near ({target_lat:.2f}°N, {target_lon:.2f}°E)")
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Map plot saved to: {save_path}")
    else:
        plt.show()
    if CLOSE_FIGURES:
        plt.close()

class PlotType:
    def __init__(self, directory, lat, lon, radius, variable, fnamebase):
        self.directory = directory
        self.lat = lat
        self.lon = lon
        self.radius = radius
        self.variable = variable
        self.filename_base = f"{fnamebase}_{variable}_r{radius}"

def main():
    save_directory = r"C:\Data\ArgoData\Figures\BSSC2025\\"
    radiis = [5.0, 10.0, 15.0, 20.0]  #nautical miles
    for the_radius in radiis:
        plot_sets = []
        parameter_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOX2_ADJUSTED']
        for param in parameter_list:
            plot_sets.append(PlotType(r"C:\Data\ArgoData\BSSC2025\GotlandDeep\\",
                                 57.3, 20.0,
                                 the_radius,
                                 param,
                                 'GotlandDeep'))
        
            plot_sets.append(PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianBay\\",
                                 64.9, 23.3,
                                 the_radius,
                                 param,
                                 'BothnianBay'))
    
            plot_sets.append(PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianSea\\",
                                 61.3, 20.05,
                                 the_radius,
                                 param,
                                 'BothnianBay1'))
            plot_sets.append(PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianSea\\",
                                 62.15, 19.9,
                                 the_radius,
                                 param,
                                 'BothnianBay2'))
            plot_sets.append(PlotType(r"C:\Data\ArgoData\BSSC2025\NBP\\",
                                     58.9, 20.3,
                                     the_radius,
                                     param,
                                     'NBP'))
        
        os.makedirs(save_directory, exist_ok=True)
        for the_plot in plot_sets:
            df = load_profiles_within_radius(the_plot.directory, 
                                             the_plot.lat, 
                                             the_plot.lon, 
                                             the_plot.radius,
                                             the_plot.variable)
            if df is not None and not df.empty:
                profile_savefile = os.path.join(save_directory, f"{the_plot.filename_base}_profile.png")
                map_savefile = os.path.join(save_directory, f"{the_plot.filename_base}_map.png")
                if PLOT_PROFILES:
                    plot_profiles(df, the_plot.variable, profile_savefile)
                if PLOT_MAPS:
                    plot_profile_locations_map(the_plot.directory,
                                           the_plot.lat, 
                                           the_plot.lon, 
                                           the_plot.radius,
                                           map_savefile)
            else:
                print("No profiles found within radius.")

if __name__ == "__main__":
    main()
