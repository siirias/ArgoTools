# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 10:39:27 2024

@author: siirias
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the provided file
file_name = 'Descent profile CTD Aanderaa Message.csv'
file_path = 'C:/Users/siirias/Downloads/FloatData/'
encoding = 'ISO-8859-1'
temperature = 'CTD - Temperature (°C)'
salinity = 'CTD - Salinity (PSU)'
pressure = 'CTD - Pressure (dbar)'

# Clean the data by removing rows where pressure is exactly 1000.0
def clean_data(data):
    filtered_data = data[data[pressure] != 1000.0]
    return filtered_data


# Plotting the data with pressure as an inverted y-axis
def plot_data(data):
    # Assuming 'Pressure' and 'Temperature' are column names; adjust as needed
    plt.figure(figsize=(10, 6))
    
    # Plot Temperature vs Pressure
    plt.subplot(1, 2, 1)
    # + np.arange(len(data)) * 0.1
    plt.plot(data[temperature], data[pressure], marker='o', linestyle='-')
    plt.gca().invert_yaxis()
    plt.xlabel('Temperature (°C)')
    plt.ylabel('Pressure (dbar)')
    plt.title('Temperature vs Pressure')

    # Plot Salinity vs Pressure 
    plt.subplot(1, 2, 2)
    plt.plot(data[salinity], data[pressure], marker='o', linestyle='-')
    plt.gca().invert_yaxis()
    plt.xlabel('Salinity (PSU)')
    plt.ylabel('Pressure (dbar)')
    plt.title('Salinity vs Pressure')
    
    plt.tight_layout()
    plt.show()

def plot_float_pressure(file_path):
    """
    Reads hydraulic message file and plots the pressure over time.
    
    Parameters:
    file_path (str): Path to the CSV file.
    """
    # Read the CSV file
    data = pd.read_csv(file_path, delimiter=';', encoding=encoding)

    # Display the first few rows of the data to understand its structure
    print(data.head())

    # Clean column names and fill forward the missing values
    data.columns = data.columns.str.strip().str.replace('"', '')
    data = data.ffill()

    # Create a datetime-like column for plotting
    data['Time'] = data['Cycle day (day)'] * 24 * 60 + data['Cycle hour (minute)'] + data['Hour (minute)']

    # Convert time to hours
    data['Time'] = data['Time'] / 60.0

    # Filter out irrelevant rows (if needed)
    data = data[data['Pressure (dbar)'] > 0]

    # Plotting the data
    plt.figure(figsize=(12, 6))
    plt.plot(data['Time'], data['Pressure (dbar)'], marker='o', linestyle='-')
    plt.gca().invert_yaxis()
    plt.xlabel('Time (hours)')
    plt.ylabel('Pressure (dbar)')
    plt.title('Pressure vs Time')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

data = pd.read_csv(file_path + file_name, delimiter = ';', encoding=encoding)

data = clean_data(data)

# Call the plotting function
plot_data(data)




# Example usage
hydraulic_file = file_path + 'Hydraulic Message.csv'
plot_float_pressure(hydraulic_file)

