# orbit.py
"""
Orbit propagation and FOV calculation functions for satellite simulation.
"""
import numpy as np
from shapely.geometry import box
from shapely.affinity import rotate
from skyfield.api import load, EarthSatellite
from config import R_earth_km, start_datetime, time_step_seconds, simulation_duration_hours, satellites as SATELLITE_CONFIG # Import SATELLITE_CONFIG
from datetime import timedelta

def propagate_orbits(satellites_from_main): # Renamed input to avoid conflict
    ts = load.timescale()
    total_seconds = int(simulation_duration_hours * 3600)
    num_steps = total_seconds // time_step_seconds
    sky_datetime_objects = [start_datetime + timedelta(seconds=i * time_step_seconds) for i in range(num_steps)]
    sky_times = ts.utc(sky_datetime_objects)
    
    for name, params in satellites_from_main.items(): # Use renamed input
        satellite_obj = EarthSatellite(params["tle_line1"], params["tle_line2"], name, ts)
        geocentric = satellite_obj.at(sky_times)
        params['latitudes'] = geocentric.subpoint().latitude.degrees
        params['longitudes'] = geocentric.subpoint().longitude.degrees
        params['altitudes_actual_km'] = geocentric.subpoint().elevation.km
        lats = params['latitudes']
        lons = params['longitudes']
        alts_actual = params['altitudes_actual_km']
        fov_polygons_shapely = []
        h_actual_avg = np.mean(alts_actual)

        # Get the base FOV from the config for this satellite
        base_fov_deg = SATELLITE_CONFIG[name].get('fov_deg', 1.0) # Default if not found

        # --- DYNAMIC FOV CALCULATION BASED ON RESOLUTION (EXAMPLE FOR GEOEYE) --
        # This section is a placeholder for how you might adjust FOV if resolution changes affect it.
        # For now, we assume the fov_deg in config is for the 1m/px resolution.
        # If 10m/px implies a different angular FOV for the sensor, that logic would go here.
        # For this task, we are primarily changing the *data size* based on resolution,
        # not necessarily the ground footprint size if the sensor's angular FOV remains constant.
        # However, if the "resolution change" means switching to a different sensor mode with a different angular FOV,
        # then params['fov_deg'] would need to be updated here based on terrain.
        # For simplicity, we'll use the configured fov_deg for geometric calculations
        # and handle data volume changes in the visualization/counting part.

        current_fov_deg = base_fov_deg
        # --- END DYNAMIC FOV ---

        alpha_sensor_rad = np.deg2rad(current_fov_deg / 2.0)

        val_for_arcsin = ((R_earth_km + h_actual_avg) / R_earth_km) * np.sin(alpha_sensor_rad)

        if val_for_arcsin >= 1.0:
            earth_central_half_angle_rad = np.arccos(R_earth_km / (R_earth_km + h_actual_avg))
        elif val_for_arcsin <= -1.0: 
            earth_central_half_angle_rad = -np.arccos(R_earth_km / (R_earth_km + h_actual_avg))
        else:
            angle_beta_rad = np.arcsin(val_for_arcsin) 
            earth_central_half_angle_rad = angle_beta_rad - alpha_sensor_rad
        
        earth_central_half_angle_rad = max(0, earth_central_half_angle_rad)
        delta_rad_fov_width = 2 * earth_central_half_angle_rad 
        delta_deg_fov_width = np.rad2deg(delta_rad_fov_width)
        
        # Store the calculated ground FOV width for potential use in data volume calculation
        params['ground_fov_width_deg'] = delta_deg_fov_width
        params['current_fov_deg'] = current_fov_deg # Store the FOV used for this step's geometry

        for i in range(len(lats)):
            cx, cy = lons[i], lats[i]
            angle_deg = 0
            if 0 < i < len(lats) - 1:
                dx = lons[i+1] - lons[i-1]; dy = lats[i+1] - lats[i-1]
                if abs(dx) > 180: dx = np.sign(dx) * (360 - abs(dx)) if dx !=0 else 0
                angle_deg = np.rad2deg(np.arctan2(dy, dx))
            elif i > 0 :
                dx = lons[i] - lons[i-1]; dy = lats[i] - lats[i-1]
                if abs(dx) > 180: dx = np.sign(dx) * (360 - abs(dx)) if dx !=0 else 0
                angle_deg = np.rad2deg(np.arctan2(dy, dx))
            fov_side_deg = delta_deg_fov_width
            current_half_lat_span = fov_side_deg / 2.0
            current_half_lon_span = (fov_side_deg / 2.0) / np.cos(np.deg2rad(cy)) if np.cos(np.deg2rad(cy)) > 0.05 else (fov_side_deg / 2.0)
            current_half_lon_span = min(current_half_lon_span, 60)
            try:
                rect_shapely = box(cx - current_half_lon_span, cy - current_half_lat_span,
                                   cx + current_half_lon_span, cy + current_half_lat_span)
                rotated_fov = rotate(rect_shapely, angle_deg, origin=(cx, cy), use_radians=False)
                fov_polygons_shapely.append(rotated_fov)
            except Exception:
                fov_polygons_shapely.append(None)
        params['fov_polygons_shapely'] = fov_polygons_shapely
        d_lon_diff = np.diff(params['longitudes'])
        params['jump_mask'] = np.concatenate(([False], np.abs(d_lon_diff) > 300))
    return satellites_from_main, sky_datetime_objects # Return renamed input

def filter_satellites(satellites, mode):
    if mode.lower() == 'modis':
        return {k: v for k, v in satellites.items() if 'modis' in k.lower()}
    elif mode.lower() == 'geoeye-1': # Jeśli simulation_mode to "GeoEye-1"
        return {k: v for k, v in satellites.items() if 'geoeye-1' in k.lower()} # Poprawiony błąd: filtruje po 'geoeye-1'
    else:
        return satellites