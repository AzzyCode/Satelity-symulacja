# visualization.py
"""
Visualization functions for satellite simulation: static plot and animation.
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature.nightshade import Nightshade
from shapely.geometry import box, MultiPolygon
from shapely.ops import unary_union
from matplotlib.patheffects import withStroke
from skyfield.api import load, wgs84
from config import (
    reception_stations,
    use_custom_map,
    custom_map_path,
    custom_map_extent,
    animation_fps,
    animation_bitrate,
    fov_history_length,
    output_video_filename,
    latitude_limit,
    PHOTO_COUNTING_MODE,     
    SATELLITE_PHOTO_SIZES_MB,
    satellites as SATELLITE_CONFIG  # Import SATELLITE_CONFIG
)

def plot_static(satellites, sky_datetime_objects):
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    if use_custom_map:
        try:
            custom_map_image = plt.imread(custom_map_path)
            ax.imshow(
                custom_map_image,
                origin='upper',
                extent=custom_map_extent,
                transform=ccrs.PlateCarree(),
                zorder=0,
            )
        except FileNotFoundError:
            print(f"Custom map '{custom_map_path}' not found. Using stock image.")
            ax.stock_img(zorder=0)
    else:
        ax.stock_img(zorder=0)

    ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='gray', zorder=1)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=2)
    map_boundary_box = box(-180, -90, 180, 90)

    for name, params in satellites.items():
        lons = params['longitudes']
        lats = params['latitudes']
        jump_mask = params['jump_mask']

        # plot ground tracks
        for i in range(len(lons) - 1):
            if not jump_mask[i+1]:
                ax.plot(
                    [lons[i], lons[i+1]],
                    [lats[i], lats[i+1]],
                    color=params['color'],
                    linewidth=1.5,
                    transform=ccrs.Geodetic(),
                    zorder=3,
                )

        # sample FOVs
        N_fov_samples = min(20, len(params['fov_polygons_shapely']))
        indices = list(range(0, len(params['fov_polygons_shapely']), max(1, len(params['fov_polygons_shapely']) // N_fov_samples)))
        for idx in indices:
            poly = params['fov_polygons_shapely'][idx]
            if poly and poly.is_valid and not poly.is_empty:
                vis = poly.intersection(map_boundary_box)
                if vis.is_empty:
                    continue
                if isinstance(vis, MultiPolygon):
                    for part in vis.geoms:
                        if part.is_valid and not part.is_empty:
                            x, y = part.exterior.xy
                            ax.plot(
                                x,
                                y,
                                color=params['color'],
                                alpha=0.3,
                                linewidth=1.0,
                                transform=ccrs.Geodetic(),
                                zorder=2,
                            )
                else:
                    x, y = vis.exterior.xy
                    ax.plot(
                        x,
                        y,
                        color=params['color'],
                        alpha=0.3,
                        linewidth=1.0,
                        transform=ccrs.Geodetic(),
                        zorder=2,
                    )

    # plot reception stations
    for st_name, st_info in reception_stations.items():
        ax.plot(
            st_info['lon'],
            st_info['lat'],
            marker=st_info['marker'],
            color=st_info['color'],
            markersize=10,
            transform=ccrs.Geodetic(),
            linestyle='',
            zorder=4,
        )
        ax.text(
            st_info['lon'] + 0.5,
            st_info['lat'] + 0.5,
            st_name,
            transform=ccrs.Geodetic(),
            color=st_info['color'],
            fontsize=9,
            zorder=5,
            path_effects=[withStroke(linewidth=2, foreground='white')],
        )

    plt.title("Satellite Ground Tracks & Sample FOVs")
    plt.tight_layout()
    plt.show()


def plot_coverage_overlay(satellites):
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_global()

    if use_custom_map:
        try:
            img = plt.imread(custom_map_path)
            ax.imshow(
                img,
                origin='upper',
                extent=custom_map_extent,
                transform=ccrs.PlateCarree(),
                zorder=0,
            )
        except FileNotFoundError:
            ax.stock_img(zorder=0)
    else:
        ax.stock_img(zorder=0)

    ax.add_feature(cfeature.COASTLINE, linewidth=0.7, edgecolor='gray', zorder=1)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False, zorder=2)

    # union of all FOVs clipped
    clip_box = box(-180, -latitude_limit, 180, latitude_limit)
    all_polys = []
    for params in satellites.values():
        for poly in params.get('fov_polygons_shapely', []):
            if poly and poly.is_valid and not poly.is_empty:
                clipped = poly.intersection(clip_box)
                if clipped.is_valid and not clipped.is_empty:
                    all_polys.append(clipped)

    if all_polys:
        coverage_union = unary_union(all_polys)
        if isinstance(coverage_union, MultiPolygon):
            geoms = coverage_union.geoms
        else:
            geoms = [coverage_union]

        for p in geoms:
            if p.is_valid and not p.is_empty:
                x, y = p.exterior.xy
                ax.fill(x, y, color='yellow', alpha=0.3, zorder=3, label='Covered Area')

    for st_name, st_info in reception_stations.items():
        ax.plot(
            st_info['lon'], st_info['lat'],
            marker=st_info['marker'], color=st_info['color'],
            markersize=10, transform=ccrs.Geodetic(), linestyle='', zorder=4
        )
        ax.text(
            st_info['lon'] + 0.5, st_info['lat'] + 0.5,
            st_name,
            transform=ccrs.Geodetic(), color=st_info['color'], fontsize=9, zorder=5,
            path_effects=[withStroke(linewidth=2, foreground='white')]
        )

    plt.title(f"Total Coverage Overlay ({latitude_limit}° to -{latitude_limit}°)")
    plt.tight_layout()
    plt.show()


def animate_coverage(satellites, sky_datetime_objects):
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    map_box = box(-180, -90, 180, 90)
    clip_box = box(-180, -latitude_limit, 180, latitude_limit)
    num_steps = len(sky_datetime_objects)

    # Counters: day land/ocean only
    day_land = 0
    day_ocean = 0

    cumulative_photo_counts = {name: 0 for name in satellites.keys()}
    cumulative_data_volumes_mb = {name: 0.0 for name in satellites.keys()}

    # classify overlap with land
    # Pre-fetch land geometries once to optimize
    try:
        LAND_GEOMETRIES = list(cfeature.LAND.geometries())
    except Exception as e:
        print(f"Error fetching land geometries: {e}. Assuming all ocean for classification.")
        LAND_GEOMETRIES = []

    def classify_surface(poly):
        if not LAND_GEOMETRIES: # If land geometries failed to load
            return 'ocean' # Default to ocean
        # Check for intersection with any land geometry
        # For performance, consider simplifying poly or land_geom if they are very complex
        # or using a spatial index if there are many land_geom features.
        for land_geom in LAND_GEOMETRIES:
            if poly.intersects(land_geom):
                return 'land'
        return 'ocean'

    ts = load.timescale()
    eph = load('de421.bsp')  # Load the planetary ephemeris once

    def animate(frame):
        nonlocal day_land, day_ocean
        ax.clear()
        ax.set_global()

        current_time_utc = sky_datetime_objects[frame]
        t_sky = ts.utc(current_time_utc) # Use the passed ts object

        # build nightshade feature - keep for visual representation
        ns = Nightshade(current_time_utc, alpha=0) # Use current_time_utc
        night_polys_visual = list(ns.geometries()) # For visual display only
        ax.add_feature(Nightshade(current_time_utc, alpha=0.3), zorder=1) # Use current_time_utc

        # background
        if use_custom_map:
            try:
                img = plt.imread(custom_map_path)
                ax.imshow(img, origin='upper', extent=custom_map_extent,
                          transform=ccrs.PlateCarree(), zorder=0)
            except FileNotFoundError:
                ax.stock_img(zorder=0)
        else:
            ax.stock_img(zorder=0)

        ax.add_feature(cfeature.COASTLINE, linewidth=0.5, zorder=2)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, zorder=2)
        ax.gridlines(draw_labels=False, color='lightgray', alpha=0.5, zorder=2)

        # cumulative coverage (optional)
        polys = []
        for params in satellites.values():
            for poly in params['fov_polygons_shapely'][:frame+1]:
                if poly and poly.is_valid and not poly.is_empty:
                    cp = poly.intersection(clip_box)
                    if cp.is_valid and not cp.is_empty:
                        polys.append(cp)
        if polys:
            cov = unary_union(polys)
            geoms = cov.geoms if isinstance(cov, MultiPolygon) else [cov]
            for g in geoms:
                x, y = g.exterior.xy
                ax.fill(x, y, color='yellow', alpha=0.2, zorder=3)

        # plot tracks & FOV count day photos
        for name, params in satellites.items():
            lons = params['longitudes']
            lats = params['latitudes']
            mask = params['jump_mask']
            # ground track
            for i in range(frame):
                if not mask[i+1]:
                    ax.plot([lons[i], lons[i+1]], [lats[i], lats[i+1]],
                            color=params['color'], transform=ccrs.Geodetic(), zorder=4)
            # current FOV
            curr = params['fov_polygons_shapely'][frame]
            if curr and curr.is_valid and not curr.is_empty:
                vis = curr.intersection(map_box)
                parts = vis.geoms if isinstance(vis, MultiPolygon) else [vis]
                for part in parts:
                    x, y = part.exterior.xy
                    ax.plot(x, y, color=params['color'], alpha=0.7,
                            linewidth=1.5, transform=ccrs.Geodetic(), zorder=5)

                # NOWA LOGIKA: Zliczanie zdjęć i objętości danych
                sat_type_key = None
                # Determine the base satellite type (e.g., 'geoeye', 'modis')
                # This assumes SATELLITE_PHOTO_SIZES_MB keys are like 'geoeye_land', 'modis'
                base_sat_type = ""
                if "geoeye" in name.lower():
                    base_sat_type = "geoeye"
                elif "modis" in name.lower():
                    base_sat_type = "modis"
                # Add other satellite types if necessary

                photo_size_mb = 0.0
                current_resolution_m_px = None # For display

                if base_sat_type == "geoeye" and 'resolutions_m_px' in SATELLITE_CONFIG[name]:
                    surface_type = classify_surface(curr) # curr is the current FOV polygon
                    if surface_type == 'land':
                        sat_type_key = "geoeye_land"
                        current_resolution_m_px = SATELLITE_CONFIG[name]['resolutions_m_px']['land']
                    else: # ocean
                        sat_type_key = "geoeye_ocean"
                        current_resolution_m_px = SATELLITE_CONFIG[name]['resolutions_m_px']['ocean']
                    photo_size_mb = SATELLITE_PHOTO_SIZES_MB.get(sat_type_key, 0.0)
                elif base_sat_type == "modis": # MODIS or other satellites without dynamic resolution
                    sat_type_key = base_sat_type # Assumes key like "modis" in SATELLITE_PHOTO_SIZES_MB
                    photo_size_mb = SATELLITE_PHOTO_SIZES_MB.get(sat_type_key, 0.0)
                else: # Fallback for other satellites or if config is missing
                    for key_pattern in SATELLITE_PHOTO_SIZES_MB.keys():
                        if key_pattern in name.lower():
                            sat_type_key = key_pattern
                            photo_size_mb = SATELLITE_PHOTO_SIZES_MB.get(sat_type_key, 0.0)
                            break
                
                should_count_photo_this_frame = False
                if PHOTO_COUNTING_MODE == 'DAY_AND_NIGHT':
                    should_count_photo_this_frame = True
                elif PHOTO_COUNTING_MODE == 'DAY_ONLY':
                    # --- NOWA LOGIKA OKREŚLANIA DNIA/NOCY ---
                    # Pobierz współrzędne środka FOV lub bezpośrednio subpointu satelity
                    # Dla uproszczenia użyjemy subpoint satelity, który już masz
                    sat_lat = lats[frame]
                    sat_lon = lons[frame]
                    
                    # Utwórz obiekt GeographicPosition dla subpointu satelity
                    # wgs84 jest już importowane na górze pliku
                    subpoint = wgs84.latlon(sat_lat, sat_lon)
                    
                    # Oblicz pozycję Słońca względem Ziemi
                    earth_observer = eph['earth'] + subpoint 
                    
                    # Obserwuj Słońce z subpointu na Ziemi
                    # eph jest przekazywane do animate_coverage
                    sun_astrometric = earth_observer.at(t_sky).observe(eph['sun'])
                    alt, az, distance = sun_astrometric.apparent().altaz() # Domyślne parametry dla refrakcji są OK

                    is_day_at_fov = alt.degrees > 0 # Słońce jest nad horyzontem
                    # --- KONIEC NOWEJ LOGIKI OKREŚLANIA DNIA/NOCY ---
                    
                    if is_day_at_fov:
                        should_count_photo_this_frame = True
                
                if should_count_photo_this_frame:
                    cumulative_photo_counts[name] += 1
                    cumulative_data_volumes_mb[name] += photo_size_mb

                # Oryginalna logika zliczania zdjęć dziennych nad lądem/oceanem
                # Możemy ją również zaktualizować, aby korzystała z nowej metody, jeśli chcesz
                # Na razie zostawmy ją, aby porównać wyniki, lub zaktualizujmy:
                if is_day_at_fov: # Użyj nowego warunku
                    surf = classify_surface(curr)
                    if surf == 'land': 
                        day_land += 1
                    else:           
                        day_ocean += 1
                # Jeśli nie chcesz już starych liczników day_land/day_ocean opartych na Nightshade,
                # możesz usunąć poniższy blok:
                # cen = curr.centroid 
                # in_night = any(n.contains(cen) for n in night_polys_visual) # Użyj night_polys_visual
                # if not in_night:
                #     surf = classify_surface(curr)
                #     if surf == 'land': 
                #         day_land += 1
                #     else:           
                #         day_ocean += 1

        # reception stations
        for st in reception_stations.values():
            ax.plot(st['lon'], st['lat'], marker=st['marker'], color=st['color'],
                    markersize=8, transform=ccrs.Geodetic(), zorder=6)
        for st in reception_stations.values():
            ax.text(st['lon']+0.5, st['lat']+0.5, st['marker'],
                    transform=ccrs.Geodetic(), color=st['color'], fontsize=8,
                    zorder=7, path_effects=[withStroke(linewidth=1.5, foreground='white')])

        # NOWE WYŚWIETLANIE LICZNIKÓW
        display_texts = []
        total_sim_photos = sum(cumulative_photo_counts.values())
        # Konwersja sumarycznej objętości na GB
        total_sim_volume_gb = sum(cumulative_data_volumes_mb.values()) / 1024.0

        display_texts.append(f"Photo Counting Mode: {PHOTO_COUNTING_MODE}")
        # Wyświetlanie w GB
        display_texts.append(f"Total: {total_sim_photos} photos, {total_sim_volume_gb:.2f} GB")
        for sat_name_display in satellites.keys():
            data_volume_gb = cumulative_data_volumes_mb[sat_name_display] / 1024.0
            res_text = ""
            if "geoeye" in sat_name_display.lower() and SATELLITE_CONFIG[sat_name_display].get('resolutions_m_px'):
                # Attempt to show current resolution; this is tricky as resolution can change per frame
                # For simplicity, we might just indicate it has variable resolution
                # or try to get the *last* determined resolution for this satellite in this frame (if available)
                # This part needs more robust state passing if we want to display live resolution accurately.
                # For now, let's just add a marker if it's GeoEye
                res_text = " (varies)"
            
            display_texts.append(
                f"{sat_name_display}{res_text}: {cumulative_photo_counts[sat_name_display]} photos, {data_volume_gb:.2f} GB"
            )
        # Dodanie istniejących liczników zdjęć dziennych
        display_texts.append(f"Day Photos (Land): {day_land}, Day Photos (Ocean): {day_ocean}")

        text_y_start = -0.02  # Startowa pozycja Y dla tekstu (blisko dolnej krawędzi)
        text_x_pos = 0.01    # Startowa pozycja X dla tekstu (blisko lewej krawędzi)
        line_height = 0.03   # Odstęp między liniami
        background_alpha = 0.5 # Przezroczystość tła dla tekstu

        for i, line in enumerate(display_texts):
            ax.text(
                text_x_pos, text_y_start - (i * line_height),
                line,
                transform=ax.transAxes, 
                ha='left', 
                va='top', 
                fontsize=8, # Zmniejszona czcionka dla większej ilości tekstu
                color='black', # Kolor tekstu
                path_effects=[withStroke(linewidth=1.5, foreground='white')], # Obrys dla czytelności
                bbox=dict(boxstyle='round,pad=0.2', fc='lightgray', alpha=background_alpha, ec='none') # Tło dla tekstu
            )
        
        ax.set_title(f"Satellite Simulation: {current_time_utc.strftime('%Y-%m-%d %H:%M:%S')} UTC", fontsize=10)
        return []

    anim = FuncAnimation(fig, animate, frames=num_steps, interval=1000/animation_fps, blit=False)
    writer = FFMpegWriter(fps=animation_fps, metadata=dict(artist='Satellite Simulator'), bitrate=animation_bitrate)
    anim.save(output_video_filename, writer=writer)
    plt.close(fig)


def animate_modis_only(satellites, sky_datetime_objects, output_filename="modis_animation.mp4"):
    """
    Animate only the MODIS satellite: blue ground track, blue FOV, red square for satellite position.
    """
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    map_boundary_box = box(-180, -90, 180, 90)
    num_steps = len(sky_datetime_objects)
    # Find MODIS satellite key
    modis_key = None
    for k in satellites:
        if "modis" in k.lower():
            modis_key = k
            break
    if not modis_key:
        raise ValueError("MODIS satellite not found in satellites dictionary.")
    params = satellites[modis_key]
    def animate(frame_num):
        ax.clear()
        ax.set_global()
        if use_custom_map:
            try:
                custom_map_image = plt.imread(custom_map_path)
                ax.imshow(custom_map_image, origin='upper', extent=custom_map_extent, transform=ccrs.PlateCarree(), zorder=0)
            except FileNotFoundError:
                ax.stock_img(zorder=0)
        else:
            ax.stock_img(zorder=0)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5, edgecolor='gray', zorder=1)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='gray', zorder=1)
        ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False, zorder=2, color='lightgray', alpha=0.5)
        current_sim_time = sky_datetime_objects[frame_num]
        # Plot ground track up to current frame
        lons_hist = params['longitudes'][:frame_num+1]
        lats_hist = params['latitudes'][:frame_num+1]
        jump_mask_hist = params['jump_mask'][:frame_num+1]
        for i in range(len(lons_hist) - 1):
            if not jump_mask_hist[i+1]:
                ax.plot([lons_hist[i], lons_hist[i+1]], [lats_hist[i], lats_hist[i+1]], color='blue', linewidth=2.0, transform=ccrs.Geodetic(), zorder=3)
        # Plot FOV (current frame)
        current_fov = params['fov_polygons_shapely'][frame_num]
        if current_fov and current_fov.is_valid and not current_fov.is_empty:
            visible_poly = current_fov.intersection(map_boundary_box)
            if not visible_poly.is_empty:
                if isinstance(visible_poly, MultiPolygon):
                    for p_part in visible_poly.geoms:
                        if p_part.is_valid and not p_part.is_empty:
                            x, y = p_part.exterior.xy
                            ax.fill(x, y, color='blue', alpha=0.3, zorder=4)
                else:
                    x, y = visible_poly.exterior.xy
                    ax.fill(x, y, color='blue', alpha=0.3, zorder=4)
        # Plot current satellite position as a red square
        sat_lon = params['longitudes'][frame_num]
        sat_lat = params['latitudes'][frame_num]
        ax.plot(sat_lon, sat_lat, marker='s', color='red', markersize=14, transform=ccrs.Geodetic(), linestyle='', zorder=5)
        # Plot reception stations (static)
        for st_name, st_info in reception_stations.items():
            ax.plot(st_info["lon"], st_info["lat"], marker=st_info["marker"], color=st_info["color"], markersize=8, transform=ccrs.Geodetic(), linestyle='', zorder=6)
            ax.text(st_info["lon"] + 0.5, st_info["lat"] + 0.5, st_name, transform=ccrs.Geodetic(), color=st_info["color"], fontsize=8, zorder=7, path_effects=[withStroke(linewidth=1.5, foreground='white')])
        # Update title with current time
        ax.set_title(f"MODIS Satellite Simulation: {current_sim_time.strftime('%Y-%m-%d %H:%M:%S')} UTC", fontsize=12)
        return []
    anim = FuncAnimation(fig, animate, frames=num_steps, interval=1000/animation_fps, blit=False)
    writer = FFMpegWriter(fps=animation_fps, metadata=dict(artist='Satellite Simulator'), bitrate=animation_bitrate)
    try:
        anim.save(output_filename, writer=writer)
        print(f"\nAnimation saved as {output_filename}")
    except Exception as e:
        print(f"\nError saving animation: {e}")
        print("Please ensure FFMpeg is installed and in your system's PATH.")
    plt.close(fig)