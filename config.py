# config.py
"""
Configuration file for satellite constellation simulation.
Edit this file to change satellites, stations, map, and simulation settings.
"""

from datetime import datetime, timedelta
from skyfield.api import utc, load
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- Simulation mode ---
# Choose 'MODIS', 'WorldView-3', or 'Both'
simulation_mode = "GeoEye-1"

# --- Earth parameters ---
R_earth_km = 6371
mu_km3_s2 = 398600.4418

# --- Satellite definitions (add more as needed) ---
satellites = {
    "MODIS_Terra": {
        "altitude_km": 681,
        "inclination_deg": 98,
        "fov_deg": 110,
        "color": "blue",
        "tle_line1": "1 25994U 99068A   25141.50000000  .00000714  00000-0  39998-4 0  9991",
        "tle_line2": "2 25994  98.2041 100.6060 0001200  80.5935 279.5400 14.57107008270001"
    },
    "GeoEye-1": {
        "altitude_km": 681,
        "inclination_deg": 98,
        "fov_deg": 1.28, 
        "color": "red",
        "tle_line1": "1 40115U 14045A   25141.50000000  .00000714  00000-0  39998-4 0  9992",
        "tle_line2": "2 40115  97.8743  47.2925 0001529  71.0323  289.1058 14.83390144560002",
        "resolutions_m_px": { # New: Define land and ocean resolutions
            "land": 1.0,
            "ocean": 10.0
        }
    }
}

# --- Time configuration ---
start_datetime = datetime(2025, 5, 21, 0, 0, 0, tzinfo=utc)
simulation_duration_hours = 6  # 24h analysis
time_step_seconds = 60


# --- Reception Stations ---
reception_stations = {
    'Utqiagvik_USA': {'lat': 71.27, 'lon': -156.81, 'marker': 'o', 'color': 'red'},
    'Dundee_Scotland': {'lat': 56.40, 'lon': -3.18, 'marker': 's', 'color': 'blue'},
    'Chitose_Japan': {'lat': 42.77, 'lon': 141.62, 'marker': '^', 'color': 'green'},
    'Mojave_USA': {'lat': 35.05, 'lon': -118.15, 'marker': 'D', 'color': 'purple'},
    'Dubai_UAE': {'lat': 24.94, 'lon': 55.35, 'marker': '*', 'color': 'orange'},
    'Paumalu_USA': {'lat': 21.67, 'lon': -158.03, 'marker': 'v', 'color': 'cyan'},
    'Harmon_Guam': {'lat': 13.51, 'lon': 144.82, 'marker': '<', 'color': 'magenta'},
    'Mwulire_Rwanda': {'lat': -1.96, 'lon': 30.39, 'marker': '>', 'color': 'teal'},
    'Tahiti_French_Polynesia': {'lat': -17.63, 'lon': -149.60, 'marker': 'h', 'color': 'brown'},
    'Awarua_New_Zealand': {'lat': -46.52, 'lon': 168.38, 'marker': 'p', 'color': 'pink'},
    'Ojebyn_Sweden': {'lat': 65.33, 'lon': 21.42, 'marker': 'X', 'color': 'gray'},
    'North_Pole_USA': {'lat': 64.79, 'lon': -147.53, 'marker': '8', 'color': 'olive'},
    'Guildford_Pole_UK': {'lat': 51.24, 'lon': -0.61, 'marker': 'P', 'color': 'navy'},
    'Obihiro_Japan': {'lat': 42.59, 'lon': 143.45, 'marker': 'H', 'color': 'maroon'},
    'Pendergrass_USA': {'lat': 34.17, 'lon': -83.67, 'marker': 'd', 'color': 'lime'},
    'Accra_Ghana': {'lat': 5.74, 'lon': -0.30, 'marker': '1', 'color': 'crimson'},
    'Pretoria_South_Africa': {'lat': -25.88, 'lon': 27.70, 'marker': '2', 'color': 'gold'},
    'Cordoba_Argentina': {'lat': -31.52, 'lon': -64.46, 'marker': '3', 'color': 'slateblue'},
    'Ushuaia_Argentina': {'lat': -54.51, 'lon': -67.11, 'marker': '4', 'color': 'coral'},
    'Sodankyla_Finland': {'lat': 67.36, 'lon': 26.63, 'marker': '+', 'color': 'black'},
    'Alice_Springs_Australia': {'lat': -23.75, 'lon': 133.88, 'marker': 'x', 'color': 'chocolate'},
    'Mingenew_Australia': {'lat': -29.01, 'lon': 115.34, 'marker': '.', 'color': 'orchid'},
    'Stockholm_Sweden': {'lat': 59.33, 'lon': 18.06, 'marker': ',', 'color': 'skyblue'},
    'Dublin_Ireland': {'lat': 53.41, 'lon': 8.24, 'marker': '|', 'color': 'darkgreen'},
    'Portland_USA': {'lat': 45.52, 'lon': -122.67, 'marker': '_', 'color': 'khaki'},
    'Columbus_USA': {'lat': 39.96, 'lon': -83.0, 'marker': '1', 'color': 'plum'},
    'Seoul_South_Korea': {'lat': 37.55, 'lon': 126.99, 'marker': '2', 'color': 'turquoise'},
    'Deadhorse_USA': {'lat': 70.20, 'lon': -148.45, 'marker': '3', 'color': 'salmon'},
    'Zallaq_Bahrain': {'lat': 26.04, 'lon': 50.48, 'marker': '4', 'color': 'indigo'},
    'Kapolei_USA': {'lat': 21.33, 'lon': -158.08, 'marker': '+', 'color': 'steelblue'},
    'Singapore': {'lat': 1.35, 'lon': 103.81, 'marker': 'x', 'color': 'green'},
    'Cape_Town_South_Africa': {'lat': -33.92, 'lon': 18.42, 'marker': 'o', 'color': 'navy'},
    'Dubbo_Australia': {'lat': -32.24, 'lon': 148.61, 'marker': 's', 'color': 'orange'},
    'Punta_Arenas_Chile': {'lat': -53.16, 'lon': -70.90, 'marker': '^', 'color': 'purple'},
}


# --- Map configuration ---
use_custom_map = True  # Set to True to use your custom map
custom_map_path = 'custom_map.png'  # Place your PNG/GeoTIFF in the same folder
custom_map_extent = [-180, 180, -90, 90]  # [lon_min, lon_max, lat_min, lat_max]

# --- Animation parameters ---
animation_fps = 10
animation_bitrate = 1800
fov_history_length = 5  # Number of previous FOVs to show with transparency
output_video_filename = "dual_satellite_animation.mp4"

# --- Coverage analysis ---
latitude_limit = 66.73  # Exclude coverage above/below +/- this latitude (e.g., 80 deg)

# Tryb zliczania zdjęć i obliczania objętości danych
PHOTO_COUNTING_MODE = "DAY_ONLY" # "DAY_ONLY" or "DAY_AND_NIGHT"

SATELLITE_PHOTO_SIZES_MB = {
    "geoeye_land": 302.96,  # GeoEye 1m/px resolution (over land)
    "geoeye_ocean": 3.03,   # GeoEye 10m/px resolution (over ocean)    
    "modis": 86.60     # Modis 250m/px resolution
}
