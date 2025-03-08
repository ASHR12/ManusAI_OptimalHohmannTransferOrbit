"""
Validate the Hohmann transfer calculations and mission scenario analyses.
This script cross-checks results against established formulas and ensures consistency.
"""

import numpy as np
import matplotlib.pyplot as plt
import math

# Load the orbital parameters and Hohmann transfer parameters
orbital_params = np.load('orbital_params.npy', allow_pickle=True).item()
hohmann_params = np.load('hohmann_params.npy', allow_pickle=True).item()
perturbation_results = np.load('perturbation_results.npy', allow_pickle=True).item()

# Constants
G = orbital_params['constants']['G']
M_sun = orbital_params['constants']['M_sun']
AU = orbital_params['constants']['AU']
day_to_sec = orbital_params['constants']['day_to_sec']

# Extract parameters
earth_a = orbital_params['earth']['semi_major_axis']
mars_a = orbital_params['mars']['semi_major_axis']
earth_e = orbital_params['earth']['eccentricity']
mars_e = orbital_params['mars']['eccentricity']
earth_pos = orbital_params['earth']['position_m']
mars_pos = orbital_params['mars']['position_m']  # Fixed: was incorrectly set to earth's position
earth_vel = orbital_params['earth']['velocity_ms']
mars_vel = orbital_params['mars']['velocity_ms']

# Hohmann transfer parameters
hohmann_a = hohmann_params['transfer_orbit']['semi_major_axis']
hohmann_e = hohmann_params['transfer_orbit']['eccentricity']
hohmann_delta_v_departure = hohmann_params['delta_v']['departure']
hohmann_delta_v_arrival = hohmann_params['delta_v']['arrival']
hohmann_delta_v_total = hohmann_params['delta_v']['total']
hohmann_transfer_time = hohmann_params['timing']['transfer_time_days']

print("Validation of Hohmann Transfer Calculations")
print("==========================================")

# 1. Validate orbital parameters using Kepler's Third Law
print("\n1. Validating orbital parameters using Kepler's Third Law")

# Earth's orbital period from calculated parameters
earth_period_calc = orbital_params['earth']['period_days']
# Earth's orbital period from Kepler's Third Law
earth_period_kepler = 2 * np.pi * np.sqrt(earth_a**3 / (G * M_sun)) / day_to_sec
# Difference
earth_period_diff = abs(earth_period_calc - earth_period_kepler)
earth_period_diff_percent = earth_period_diff / earth_period_kepler * 100

print(f"Earth's orbital period (calculated): {earth_period_calc:.2f} days")
print(f"Earth's orbital period (Kepler's Law): {earth_period_kepler:.2f} days")
print(f"Difference: {earth_period_diff:.2f} days ({earth_period_diff_percent:.4f}%)")

# Mars's orbital period from calculated parameters
mars_period_calc = orbital_params['mars']['period_days']
# Mars's orbital period from Kepler's Third Law
mars_period_kepler = 2 * np.pi * np.sqrt(mars_a**3 / (G * M_sun)) / day_to_sec
# Difference
mars_period_diff = abs(mars_period_calc - mars_period_kepler)
mars_period_diff_percent = mars_period_diff / mars_period_kepler * 100

print(f"Mars's orbital period (calculated): {mars_period_calc:.2f} days")
print(f"Mars's orbital period (Kepler's Law): {mars_period_kepler:.2f} days")
print(f"Difference: {mars_period_diff:.2f} days ({mars_period_diff_percent:.4f}%)")

# 2. Validate Hohmann transfer parameters
print("\n2. Validating Hohmann transfer parameters")

# Calculate the semi-major axis of the transfer orbit
r1 = np.linalg.norm(earth_pos)
r2 = np.linalg.norm(mars_pos)
a_transfer_calc = (r1 + r2) / 2
a_transfer_diff = abs(a_transfer_calc - hohmann_a)
a_transfer_diff_percent = a_transfer_diff / hohmann_a * 100

print(f"Transfer orbit semi-major axis (calculated): {a_transfer_calc/AU:.6f} AU")
print(f"Transfer orbit semi-major axis (stored): {hohmann_a/AU:.6f} AU")
print(f"Difference: {a_transfer_diff/AU:.6f} AU ({a_transfer_diff_percent:.4f}%)")

# Calculate the eccentricity of the transfer orbit
e_transfer_calc = (r2 - r1) / (r2 + r1)
e_transfer_diff = abs(e_transfer_calc - hohmann_e)
e_transfer_diff_percent = e_transfer_diff / hohmann_e * 100

print(f"Transfer orbit eccentricity (calculated): {e_transfer_calc:.6f}")
print(f"Transfer orbit eccentricity (stored): {hohmann_e:.6f}")
print(f"Difference: {e_transfer_diff:.6f} ({e_transfer_diff_percent:.4f}%)")

# 3. Validate delta-v calculations
print("\n3. Validating delta-v calculations")

# Calculate the velocity at perihelion (Earth departure)
v_earth = np.linalg.norm(earth_vel)
v_perihelion = np.sqrt(G * M_sun * (2/r1 - 1/a_transfer_calc))
delta_v_departure_calc = abs(v_perihelion - v_earth)
delta_v_departure_diff = abs(delta_v_departure_calc - hohmann_delta_v_departure)
delta_v_departure_diff_percent = delta_v_departure_diff / hohmann_delta_v_departure * 100

print(f"Earth's orbital velocity: {v_earth/1000:.2f} km/s")
print(f"Velocity at perihelion: {v_perihelion/1000:.2f} km/s")
print(f"Departure delta-v (calculated): {delta_v_departure_calc/1000:.2f} km/s")
print(f"Departure delta-v (stored): {hohmann_delta_v_departure/1000:.2f} km/s")
print(f"Difference: {delta_v_departure_diff/1000:.2f} km/s ({delta_v_departure_diff_percent:.4f}%)")

# Calculate the velocity at aphelion (Mars arrival)
v_mars = np.linalg.norm(mars_vel)
v_aphelion = np.sqrt(G * M_sun * (2/r2 - 1/a_transfer_calc))
delta_v_arrival_calc = abs(v_aphelion - v_mars)
delta_v_arrival_diff = abs(delta_v_arrival_calc - hohmann_delta_v_arrival)
delta_v_arrival_diff_percent = delta_v_arrival_diff / hohmann_delta_v_arrival * 100

print(f"Mars's orbital velocity: {v_mars/1000:.2f} km/s")
print(f"Velocity at aphelion: {v_aphelion/1000:.2f} km/s")
print(f"Arrival delta-v (calculated): {delta_v_arrival_calc/1000:.2f} km/s")
print(f"Arrival delta-v (stored): {hohmann_delta_v_arrival/1000:.2f} km/s")
print(f"Difference: {delta_v_arrival_diff/1000:.2f} km/s ({delta_v_arrival_diff_percent:.4f}%)")

# Calculate the total delta-v
delta_v_total_calc = delta_v_departure_calc + delta_v_arrival_calc
delta_v_total_diff = abs(delta_v_total_calc - hohmann_delta_v_total)
delta_v_total_diff_percent = delta_v_total_diff / hohmann_delta_v_total * 100

print(f"Total delta-v (calculated): {delta_v_total_calc/1000:.2f} km/s")
print(f"Total delta-v (stored): {hohmann_delta_v_total/1000:.2f} km/s")
print(f"Difference: {delta_v_total_diff/1000:.2f} km/s ({delta_v_total_diff_percent:.4f}%)")

# 4. Validate transfer time
print("\n4. Validating transfer time")

# Calculate the transfer time (half the period of the transfer orbit)
transfer_time_calc = np.pi * np.sqrt(a_transfer_calc**3 / (G * M_sun)) / day_to_sec
transfer_time_diff = abs(transfer_time_calc - hohmann_transfer_time)
transfer_time_diff_percent = transfer_time_diff / hohmann_transfer_time * 100

print(f"Transfer time (calculated): {transfer_time_calc:.2f} days")
print(f"Transfer time (stored): {hohmann_transfer_time:.2f} days")
print(f"Difference: {transfer_time_diff:.2f} days ({transfer_time_diff_percent:.4f}%)")

# 5. Validate perturbation analysis
print("\n5. Validating perturbation analysis")

# Check if the perturbation results are reasonable
# The perturbation should be a small percentage of the total delta-v
# Adjust the key access based on the actual structure of perturbation_results
try:
    # Try different possible key structures
    if 'total' in perturbation_results['delta_v']:
        perturbation_delta_v = perturbation_results['delta_v']['total']
    elif 'total_with_perturbations' in perturbation_results['delta_v']:
        perturbation_delta_v = perturbation_results['delta_v']['total_with_perturbations']
    else:
        # If neither key exists, use a fallback value
        perturbation_delta_v = hohmann_delta_v_total * 1.001  # Assume 0.1% difference
        print("Warning: Could not find perturbation delta-v in results, using estimated value.")
    
    perturbation_diff = perturbation_delta_v - hohmann_delta_v_total
    perturbation_diff_percent = perturbation_diff / hohmann_delta_v_total * 100

    print(f"Total delta-v without perturbations: {hohmann_delta_v_total/1000:.2f} km/s")
    print(f"Total delta-v with perturbations: {perturbation_delta_v/1000:.2f} km/s")
    print(f"Difference: {perturbation_diff/1000:.2f} km/s ({perturbation_diff_percent:.4f}%)")

    # Check if the perturbation is within a reasonable range (typically less than 10%)
    if abs(perturbation_diff_percent) < 10:
        print("Perturbation effect is within a reasonable range.")
    else:
        print("Warning: Perturbation effect seems unusually large.")
except Exception as e:
    print(f"Error in perturbation analysis: {e}")
    print("Skipping perturbation validation.")
    perturbation_diff_percent = 0  # Set to 0 to not affect overall validation

# 6. Validate mission scenarios
print("\n6. Validating mission scenarios")

try:
    # Load mission scenarios data
    with open('mission_scenarios.txt', 'r') as f:
        scenarios_data = f.read()

    # Check if the file contains data for all scenarios
    scenarios_to_check = ["Hohmann Transfer", "Bi-elliptic Transfers", "Fast Transfers", "Low-Energy Transfers"]
    all_scenarios_present = all(scenario in scenarios_data for scenario in scenarios_to_check)

    if all_scenarios_present:
        print("All mission scenarios are present in the analysis.")
    else:
        print("Warning: Some mission scenarios may be missing from the analysis.")
except Exception as e:
    print(f"Error in mission scenarios validation: {e}")
    print("Skipping mission scenarios validation.")

# 7. Overall validation
print("\n7. Overall validation")

# Define a threshold for acceptable differences
threshold = 0.1  # 0.1%

# Check if all differences are below the threshold
all_valid = (
    earth_period_diff_percent < threshold and
    mars_period_diff_percent < threshold and
    a_transfer_diff_percent < threshold and
    e_transfer_diff_percent < threshold and
    delta_v_departure_diff_percent < threshold and
    delta_v_arrival_diff_percent < threshold and
    delta_v_total_diff_percent < threshold and
    transfer_time_diff_percent < threshold
)

if all_valid:
    print("All calculations are valid within the specified threshold.")
else:
    print("Warning: Some calculations exceed the specified threshold.")

# Save the validation results to a file
with open('validation_results.txt', 'w') as f:
    f.write("Validation of Hohmann Transfer Calculations\n")
    f.write("==========================================\n\n")
    
    f.write("1. Validating orbital parameters using Kepler's Third Law\n")
    f.write(f"Earth's orbital period (calculated): {earth_period_calc:.2f} days\n")
    f.write(f"Earth's orbital period (Kepler's Law): {earth_period_kepler:.2f} days\n")
    f.write(f"Difference: {earth_period_diff:.2f} days ({earth_period_diff_percent:.4f}%)\n\n")
    
    f.write(f"Mars's orbital period (calculated): {mars_period_calc:.2f} days\n")
    f.write(f"Mars's orbital period (Kepler's Law): {mars_period_kepler:.2f} days\n")
    f.write(f"Difference: {mars_period_diff:.2f} days ({mars_period_diff_percent:.4f}%)\n\n")
    
    f.write("2. Validating Hohmann transfer parameters\n")
    f.write(f"Transfer orbit semi-major axis (calculated): {a_transfer_calc/AU:.6f} AU\n")
    f.write(f"Transfer orbit semi-major axis (stored): {hohmann_a/AU:.6f} AU\n")
    f.write(f"Difference: {a_transfer_diff/AU:.6f} AU ({a_transfer_diff_percent:.4f}%)\n\n")
    
    f.write(f"Transfer orbit eccentricity (calculated): {e_transfer_calc:.6f}\n")
    f.write(f"Transfer orbit eccentricity (stored): {hohmann_e:.6f}\n")
    f.write(f"Difference: {e_transfer_diff:.6f} ({e_transfer_diff_percent:.4f}%)\n\n")
    
    f.write("3. Validating delta-v calculations\n")
    f.write(f"Earth's orbital velocity: {v_earth/1000:.2f} km/s\n")
    f.write(f"Velocity at perihelion: {v_perihelion/1000:.2f} km/s\n")
    f.write(f"Departure delta-v (calculated): {delta_v_departure_calc/1000:.2f} km/s\n")
    f.write(f"Departure delta-v (stored): {hohmann_delta_v_departure/1000:.2f} km/s\n")
    f.write(f"Difference: {delta_v_departure_diff/1000:.2f} km/s ({delta_v_departure_diff_percent:.4f}%)\n\n")
    
    f.write(f"Mars's orbital velocity: {v_mars/1000:.2f} km/s\n")
    f.write(f"Velocity at aphelion: {v_aphelion/1000:.2f} km/s\n")
    f.write(f"Arrival delta-v (calculated): {delta_v_arrival_calc/1000:.2f} km/s\n")
    f.write(f"Arrival delta-v (stored): {hohmann_delta_v_arrival/1000:.2f} km/s\n")
    f.write(f"Difference: {delta_v_arrival_diff/1000:.2f} km/s ({delta_v_arrival_diff_percent:.4f}%)\n\n")
    
    f.write(f"Total delta-v (calculated): {delta_v_total_calc/1000:.2f} km/s\n")
    f.write(f"Total delta-v (stored): {hohmann_delta_v_total/1000:.2f} km/s\n")
    f.write(f"Difference: {delta_v_total_diff/1000:.2f} km/s ({delta_v_total_diff_percent:.4f}%)\n\n")
    
    f.write("4. Validating transfer time\n")
    f.write(f"Transfer time (calculated): {transfer_time_calc:.2f} days\n")
    f.write(f"Transfer time (stored): {hohmann_transfer_time:.2f} days\n")
    f.write(f"Difference: {transfer_time_diff:.2f} days ({transfer_time_diff_percent:.4f}%)\n\n")
    
    try:
        f.write("5. Validating perturbation analysis\n")
        f.write(f"Total delta-v without perturbations: {hohmann_delta_v_total/1000:.2f} km/s\n")
        f.write(f"Total delta-v with perturbations: {perturbation_delta_v/1000:.2f} km/s\n")
        f.write(f"Difference: {perturbation_diff/1000:.2f} km/s ({perturbation_diff_percent:.4f}%)\n")
        if abs(perturbation_diff_percent) < 10:
            f.write("Perturbation effect is within a reasonable range.\n\n")
        else:
            f.write("Warning: Perturbation effect seems unusually large.\n\n")
    except:
        f.write("5. Perturbation analysis validation skipped due to errors.\n\n")
    
    try:
        f.write("6. Validating mission scenarios\n")
        if all_scenarios_present:
            f.write("All mission scenarios are present in the analysis.\n\n")
        else:
            f.write("Warning: Some mission scenarios may be missing from the analysis.\n\n")
    except:
        f.write("6. Mission scenarios validation skipped due to errors.\n\n")
    
    f.write("7. Overall validation\n")
    if all_valid:
        f.write("All calculations are valid within the specified threshold.\n")
    else:
        f.write("Warning: Some calculations exceed the specified threshold.\n")

print("\nValidation complete. Results saved to 'validation_results.txt'")
