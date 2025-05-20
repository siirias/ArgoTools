import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from datetime import datetime
from geopy.distance import geodesic

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
    return df.sort_values("time")

def plot_profiles(df, variable):
    """
    Plot profiles from dataframe (time, depth, value) using pcolormesh-style plotting.
    """
    # Pivot to wide format with time-depth grid
    time_vals = pd.to_datetime(sorted(df['time'].unique()))
    depths = np.linspace(df['depth'].min(), df['depth'].max(), 200)
    times = pd.to_datetime(df['time'].unique())
    times = sorted(times)

    # Interpolate values to a common depth grid for each time
    time_labels = []
    value_grid = []
    for time in times:
        profile = df[df['time'] == time].sort_values("depth")
        if len(profile) < 5:
            continue
        interp = np.interp(depths, profile['depth'], profile['value'])
        value_grid.append(interp)
        time_labels.append(time)

    value_grid = np.array(value_grid).T  # shape: depth x time

    plt.figure(figsize=(12, 6))
    im = plt.pcolormesh(time_labels, depths, value_grid, shading='auto')
    plt.colorbar(im, label=variable)
    plt.gca().invert_yaxis()
    plt.xlabel("Time")
    plt.ylabel("Depth [dbar]")
    plt.title(f"{variable} profiles over time")
    plt.tight_layout()
    plt.show()

def main():
#    directory = r"C:\Data\ArgoData\BSSC2025\NBP"  # Adjust this to your actual path
    directory = r"C:\Data\ArgoData\BSSC2025\GotlandDeep"  # Adjust this to your actual path
    lat = 57.3
    lon = 20.0
    radius_nm = 5  # nautical miles
    variable = "TEMP"

    df = load_profiles_within_radius(directory, lat, lon, radius_nm, variable)
    if not df.empty:
        plot_profiles(df, variable)
    else:
        print("No profiles found within radius.")

if __name__ == "__main__":
    main()
