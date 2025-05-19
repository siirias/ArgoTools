import re
from datetime import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

data = input("Enter the data (press Enter when done):\n")
entries = data.strip().split('\n\n')

latitudes = []
longitudes = []
timestamps = []

for entry in entries:
    lines = entry.strip().split('\n')
    coordinates = []

    for line in lines:
        match = re.search(r'Fix:\s+([\d.]+)\s+([\d.]+)\s+(\d{2}/\d{2}/\d{4})\s+(\d{2})(\d{2})(\d{2})', line)
        if match:
            longitude = float(match.group(1))
            latitude = float(match.group(2))
            longitudes.append(longitude)
            latitudes.append(latitude)
            date_str = match.group(3)
            hour = int(match.group(4))
            minute = int(match.group(5))
            second = int(match.group(6))

            timestamp = datetime.strptime(date_str, "%m/%d/%Y").replace(hour=hour, minute=minute, second=second)
            timestamps.append(timestamp)

            coordinates.append((longitude, latitude))

    if coordinates:
        # Create figure and axes with Cartopy projection
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())

        # Add map features
        ax.add_feature(cfeature.LAND, color='white')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

        # Plot track
        lon, lat = zip(*coordinates)
        alpha_values = np.linspace(0.1, 1.0, len(coordinates))
        label_skip = 2  # Set the frequency of plotted labels (e.g., every 2nd point)

        for i in range(len(coordinates) - 1):
            ax.plot([lon[i], lon[i+1]], [lat[i], lat[i+1]], linewidth=2, color='red', alpha=alpha_values[i],
                    transform=ccrs.PlateCarree())
            if i % label_skip == 0:
                ax.text(coordinates[i][0], coordinates[i][1], timestamps[i].strftime("%H:%M:%S UTC"),
                        fontsize=8, color='black', alpha=0.7, transform=ccrs.PlateCarree())

        # Always plot the latest point with label
        ax.plot(lon[-1], lat[-1], 'ro', markersize=6, transform=ccrs.PlateCarree())
        ax.text(coordinates[-1][0], coordinates[-1][1], timestamps[-1].strftime("%H:%M:%S UTC"),
                fontsize=8, color='black', alpha=0.7, transform=ccrs.PlateCarree())

        # Set plot extent
        lon_min, lon_max = min(lon), max(lon)
        lat_min, lat_max = min(lat), max(lat)
        ax.set_extent([lon_min-1, lon_max+1, lat_min-1, lat_max+1], crs=ccrs.PlateCarree())

        plt.title("GPS Track")
        plt.show()

# Print list of coordinates with timestamps
print("Coordinates:")
for i in range(0, len(coordinates), label_skip):
    print("Latitude: {:.3f}, Longitude: {:.3f}, Time (UTC): {}".format(
        latitudes[i], longitudes[i], timestamps[i].strftime("%H:%M:%S")
    ))
