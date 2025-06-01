# main.py
"""
Main script for satellite constellation simulation and visualization.
Run this script to generate both a static plot and an animation video.
"""
from orbit import propagate_orbits
from visualization import animate_coverage, plot_coverage_overlay
from config import simulation_mode, satellites

def filter_satellites(satellites_dict, mode_string):
    """Filtruje słownik satelitów na podstawie podanego trybu."""
    mode_lower = mode_string.lower()

    # Obsługa trybów specjalnych
    if mode_lower == 'both' or mode_lower == 'all':
        return satellites_dict

    # Sprawdzenie, czy mode_string jest bezpośrednim kluczem w satellites_dict
    if mode_string in satellites_dict:
        return {mode_string: satellites_dict[mode_string]}

    # Filtrowanie na podstawie części nazwy (np. "geoeye" w "GeoEye-1")
    # To pozwala na elastyczność, np. simulation_mode = "GeoEye" powinno znaleźć "GeoEye-1"
    # simulation_mode = "MODIS" powinno znaleźć "MODIS_Terra"
    filtered_sats = {}
    # Sprawdzamy, czy mode_lower (np. "geoeye") jest częścią klucza satelity (np. "geoeye-1")
    # lub czy klucz satelity (np. "geoeye-1") zawiera mode_lower.
    # Dla pewności sprawdzamy obie strony, jeśli np. mode_string to "GeoEye" a klucz to "GeoEye-1"
    for sat_key, sat_params in satellites_dict.items():
        if mode_lower in sat_key.lower():
            filtered_sats[sat_key] = sat_params
            
    if filtered_sats:
        return filtered_sats
    else:
        # Jeśli nic nie pasuje, można zwrócić wszystkie z ostrzeżeniem lub pusty słownik
        print(f"Ostrzeżenie: Tryb symulacji '{mode_string}' nie pasuje do żadnego zdefiniowanego satelity ani słowa kluczowego ('both', 'all', typ satelity). Symuluję wszystkie dostępne satelity.")
        return satellites_dict

if __name__ == "__main__":
    # Filter satellites based on simulation_mode
    selected_sats = filter_satellites(satellites, simulation_mode)
    # Propagate orbits and calculate FOVs
    satellites_with_tracks, sky_datetime_objects = propagate_orbits(selected_sats)
    # Generate and save animation video (with coverage overlay)
    animate_coverage(satellites_with_tracks, sky_datetime_objects)
    # Plot static coverage overlay at the end
    plot_coverage_overlay(satellites_with_tracks)
    print("\nSimulation and visualization complete. Check the output video file and coverage plot in your folder.")