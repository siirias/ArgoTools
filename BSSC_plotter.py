import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
from cartopy.geodesic import Geodesic
from datetime import datetime
from geopy.distance import geodesic
from cartopy.feature import NaturalEarthFeature
from matplotlib.colors import ListedColormap, BoundaryNorm, LinearSegmentedColormap
from datetime import timedelta
import time
from colorama import Fore, Style, init
init(autoreset=True)  # ensures colors don’t bleed into next prints


DATE_RANGE_MIN = datetime(2000, 1, 1)
DATE_RANGE_MAX = datetime(2025, 12, 31)
CLOSE_FIGURES = True
PLOT_MAPS = True
PLOT_PROFILES = True

def load_profiles_within_radius(directory, target_lat, target_lon, radius_nm, variable):
    radius_km = radius_nm * 1.852
    files = glob.glob(os.path.join(directory, "*.nc"))
    records = []
    included_files = []
    skipped_files = []
    for file in files:
        try:
            ds = xr.open_dataset(file)
            lat = ds['LATITUDE'].values
            lon = ds['LONGITUDE'].values
            time = ds['TIME'].values
            pres = ds['PRES'].values
            var_data = ds[variable].values
            unit = ds[variable].attrs.get('units', '')

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
            included_files.append(os.path.basename(file))
        except Exception as e:
            skipped_files.append(os.path.basename(file))
    print(f"From {directory}")
    print(f"{Fore.GREEN}included:{', '.join(included_files)}\n\n")
    print(f"{Fore.RED}excluded:{', '.join(skipped_files)}\n\n")
    df = pd.DataFrame(records)
    if df.empty or "time" not in df.columns:
        return None, None    
    return df.sort_values("time"), unit

def interpolate_profiles_to_grid(df, depth_grid):
    times = sorted(pd.to_datetime(df['time'].unique()))
    interp_profiles = []
    valid_times = []

    for t in times:
        prof = df[df['time'] == t].sort_values("depth")
        if len(prof) < 5:
            continue
        d = prof['depth'].values
        v = prof['value'].values
        try:
            interp = np.interp(depth_grid, d, v, left=np.nan, right=np.nan)
            interp_profiles.append(interp)
            valid_times.append(t)
        except Exception as e:
            print(f"Interpolation failed at {t}: {e}")

    return np.array(interp_profiles), valid_times

def pick_colormap(variable):
    if "TEMP" in variable.upper():
        return cmocean.cm.thermal
    elif "PSAL" in variable.upper():
        return cmocean.cm.haline
    elif "DOX2" in variable.upper():
        return cmocean.cm.matter
    return 'viridis'

def plot_profiles(interp_profiles, interp_times, depth_grid, variable, 
                  save_path=None, unit='', max_gap_days = None):
    # Convert interp_times to list of pandas Timestamps
    interp_times = pd.to_datetime(interp_times)
    value_grid = np.array(interp_profiles)
    if max_gap_days is not None:
        new_profiles = []
        new_times = []

        for i in range(len(interp_times)):
            new_profiles.append(value_grid[i])
            new_times.append(interp_times[i])

            if i < len(interp_times) - 1:
                time_diff = (interp_times[i + 1] - interp_times[i]).days
                if time_diff > max_gap_days:
                    new_profiles.append([np.nan] * value_grid.shape[1])
                    new_times.append(interp_times[i] + timedelta(days=np.min([time_diff/2.0,max_gap_days/2.1])))
                    new_profiles.append([np.nan] * value_grid.shape[1])
                    new_times.append(interp_times[i+1] - timedelta(days=np.min([time_diff/2.0,max_gap_days/2.1])))

        interp_profiles = np.array(new_profiles)
        interp_times = np.array(new_times)


    value_grid = np.array(interp_profiles).T
    plt.figure(figsize=(12, 6))
    cmap = pick_colormap(variable)
    im = plt.pcolormesh(interp_times, depth_grid, value_grid, shading='nearest', cmap=cmap)
    plt.colorbar(im, label=f"{variable} ({unit})")
    plt.gca().invert_yaxis()
    plt.xlabel("Time")
    plt.ylabel("Depth [dbar]")
    plt.title(f"{variable} profiles over time")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"saved: {save_path}")
    if CLOSE_FIGURES:
        plt.close()


def plot_profile_cloud(df, variable, save_path=None, unit='', colortype="time_range",
                        interp_profiles=None, interp_times=None, depth_grid=None):
    if interp_profiles is None or interp_times is None or depth_grid is None:
        print("Missing interpolated data.")
        return

    fig, ax = plt.subplots(figsize=(8, 10))
    interp_times = np.array(pd.to_datetime(interp_times))

    if colortype == "months":
        seasonal_colors = LinearSegmentedColormap.from_list(
            "seasonal_months",
            ['#0055FF', '#00D4FF', '#00FF88', '#88FF00', '#FFFF00', '#FFC000',
             '#FF8000', '#FF4000', '#FF0000', '#D00070', '#8800FF', '#4400FF'], N=12)
        month_colors = seasonal_colors(np.linspace(0, 1, 12))
        months = np.array([pd.Timestamp(t).month for t in interp_times])
        a_value = 0.1 
        if(len(interp_profiles)>200): #too many profiles, the cloud gets too dense.
            a_value = 0.03
        for i, profile in enumerate(interp_profiles):
            ax.plot(profile, depth_grid, color=month_colors[months[i]-1], alpha=a_value, linewidth=1)

        for month in range(1, 13):
            mask = months == month
            if not np.any(mask):
                continue
            mean_profile = np.nanmean(interp_profiles[mask], axis=0)
            ax.plot(mean_profile, depth_grid, color=month_colors[month - 1], linewidth=2.0)

        cmap = ListedColormap(month_colors)
        bounds = np.arange(13) - 0.5
        norm = BoundaryNorm(bounds, cmap.N)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, ticks=np.arange(12))
        cbar.set_ticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                             'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
        cbar.set_label("Month")

    else:
        time_nums = mdates.date2num(interp_times)
        norm = plt.Normalize(time_nums.min(), time_nums.max())
        cmap = plt.colormaps["plasma"]

        for i, profile in enumerate(interp_profiles):
            ax.plot(profile, depth_grid, color=cmap(norm(time_nums[i])), alpha=0.1, linewidth=1)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax)
        ticks = np.linspace(time_nums.min(), time_nums.max(), num=6)
        tick_labels = [mdates.num2date(tick).strftime('%Y-%m') for tick in ticks]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(tick_labels)
        cbar.set_label('Time (old → new)')

    ax.invert_yaxis()
    ax.set_xlabel(f"{variable} ({unit})")
    ax.set_ylabel("Depth [dbar]")
    ax.set_title(f"{variable} profile cloud")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"saved: {save_path}")
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

#    ax.coastlines(resolution='10m')
    #ax.add_feature(cfeature.LAND, facecolor='lightgray')
    land = NaturalEarthFeature('physical', 'land', scale='50m',\
                          edgecolor='face', facecolor='lightgray')    
    ax.add_feature(land)
    
    ax.gridlines(draw_labels=True)

    for (lat, lon), ok in zip(locations, accepted):
        color = 'blue' if ok else 'darkred'
        ax.plot(lon, lat, marker='o', color=color, markersize=1, transform=ccrs.PlateCarree())

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
        self.set_name = fnamebase

def main():
    start_time = time.time()
    save_directory = r"C:\Data\ArgoData\Figures\BSSC2025\\"
    radiis = [5, 10, 15, 20, 2000]  # nautical miles
    parameter_list = ['TEMP_ADJUSTED', 'PSAL_ADJUSTED', 'DOX2_ADJUSTED']

    for the_radius in radiis:
        plot_sets = \
            [PlotType(r"C:\Data\ArgoData\BSSC2025\GotlandDeep\\", 57.3, 20.0, the_radius, param, 'GotlandDeep')
             for param in parameter_list]\
        + [PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianBay\\", 64.9, 23.3, the_radius, param, 'BothnianBay')
             for param in parameter_list]\
        + [PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianSea\\", 61.3, 20.05, the_radius, param, 'BothnianSea1')
            for param in parameter_list]\
        + [PlotType(r"C:\Data\ArgoData\BSSC2025\BothnianSea\\", 62.15, 19.9, the_radius, param, 'BothnianSea2')
             for param in parameter_list]\
        + [PlotType(r"C:\Data\ArgoData\BSSC2025\NBP\\", 58.9, 20.3, the_radius, param, 'NBP')
             for param in parameter_list]

        for the_plot in plot_sets:
            the_save_directory = os.path.join(save_directory, the_plot.set_name, the_plot.variable)
            os.makedirs(the_save_directory, exist_ok=True)

            df, unit = load_profiles_within_radius(the_plot.directory, the_plot.lat, the_plot.lon,
                                                   the_plot.radius, the_plot.variable)
            if df is not None and not df.empty:
                depth_grid = np.linspace(df['depth'].min(), df['depth'].max(), 400)
                interp_profiles, interp_times = interpolate_profiles_to_grid(df, depth_grid)

                map_savefile = os.path.join(the_save_directory, f"{the_plot.filename_base}_map.png")
                profile_savefile = os.path.join(the_save_directory, f"{the_plot.filename_base}_profile.png")
                profile_gap_savefile = os.path.join(the_save_directory, f"{the_plot.filename_base}_profile_gap.png")
                cloud_savefile_tr = os.path.join(the_save_directory, f"{the_plot.filename_base}_cloud_tr.png")
                cloud_savefile_m = os.path.join(the_save_directory, f"{the_plot.filename_base}_cloud_m.png")

                if PLOT_PROFILES:
                    plot_profiles(interp_profiles, interp_times, depth_grid, 
                                  the_plot.variable, profile_savefile, unit)
                    plot_profiles(interp_profiles, interp_times, depth_grid, 
                                  the_plot.variable, profile_gap_savefile, unit, 20)
                    plot_profile_cloud(df, the_plot.variable, cloud_savefile_tr, unit, "time_range",
                                        interp_profiles, interp_times, depth_grid)
                    plot_profile_cloud(df, the_plot.variable, cloud_savefile_m, unit, "months",
                                        interp_profiles, interp_times, depth_grid)
                if PLOT_MAPS:
                    plot_profile_locations_map(the_plot.directory, the_plot.lat, the_plot.lon, the_plot.radius, map_savefile)
    end_time = time.time()
    elapsed = end_time - start_time
    print(f"Total plotting time: {(elapsed/60.0):.2f} minutes")
    
if __name__ == "__main__":
    main()

