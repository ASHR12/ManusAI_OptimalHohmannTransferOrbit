"""
Visualize additional mission scenarios beyond the basic Hohmann transfer.
This script explores alternative transfer strategies and their trade-offs.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
earth_i = orbital_params['earth']['inclination']
mars_i = orbital_params['mars']['inclination']
earth_pos = orbital_params['earth']['position_m']
mars_pos = orbital_params['mars']['position_m']
earth_vel = orbital_params['earth']['velocity_ms']
mars_vel = orbital_params['mars']['velocity_ms']

# Hohmann transfer parameters
hohmann_a = hohmann_params['transfer_orbit']['semi_major_axis']
hohmann_e = hohmann_params['transfer_orbit']['eccentricity']
hohmann_delta_v = hohmann_params['delta_v']['total']
hohmann_transfer_time = hohmann_params['timing']['transfer_time_days']

# Function to calculate orbital parameters
def calculate_orbit_parameters(r1, r2):
    """Calculate orbital parameters for a transfer orbit between two radii."""
    a = (r1 + r2) / 2
    e = (r2 - r1) / (r2 + r1)
    return a, e

# Function to calculate delta-v for a transfer
def calculate_delta_v(r1, r2, a_transfer):
    """Calculate delta-v for a transfer between two circular orbits."""
    # Velocity in circular orbit at r1
    v1 = np.sqrt(G * M_sun / r1)
    # Velocity in transfer orbit at perihelion (r1)
    vp = np.sqrt(G * M_sun * (2/r1 - 1/a_transfer))
    # Delta-v for departure
    delta_v1 = abs(vp - v1)
    
    # Velocity in circular orbit at r2
    v2 = np.sqrt(G * M_sun / r2)
    # Velocity in transfer orbit at aphelion (r2)
    va = np.sqrt(G * M_sun * (2/r2 - 1/a_transfer))
    # Delta-v for arrival
    delta_v2 = abs(v2 - va)
    
    # Total delta-v
    total_delta_v = delta_v1 + delta_v2
    
    return delta_v1, delta_v2, total_delta_v

# Function to calculate transfer time
def calculate_transfer_time(a_transfer):
    """Calculate transfer time for half an orbit with semi-major axis a."""
    period = 2 * np.pi * np.sqrt(a_transfer**3 / (G * M_sun))
    transfer_time = period / 2
    return transfer_time

# 1. Bi-elliptic Transfer
# A bi-elliptic transfer uses an intermediate elliptical orbit that goes beyond the target orbit
# It can be more efficient than a Hohmann transfer for large ratios of final to initial radius

# Define a range of intermediate radii for bi-elliptic transfers
r1 = np.linalg.norm(earth_pos)
r2 = np.linalg.norm(mars_pos)
r3_factors = np.linspace(1.5, 5.0, 10)  # Factors to multiply Mars radius by
r3_values = r2 * r3_factors

# Calculate delta-v and transfer time for each bi-elliptic transfer
bielliptic_results = []
for r3 in r3_values:
    # First ellipse (Earth to intermediate point)
    a1 = (r1 + r3) / 2
    # Second ellipse (intermediate point to Mars)
    a2 = (r3 + r2) / 2
    
    # Delta-v for first burn (Earth departure)
    v1 = np.sqrt(G * M_sun / r1)
    vp1 = np.sqrt(G * M_sun * (2/r1 - 1/a1))
    delta_v1 = abs(vp1 - v1)
    
    # Delta-v for second burn (at intermediate point)
    va1 = np.sqrt(G * M_sun * (2/r3 - 1/a1))
    vp2 = np.sqrt(G * M_sun * (2/r3 - 1/a2))
    delta_v2 = abs(vp2 - va1)
    
    # Delta-v for third burn (Mars arrival)
    v2 = np.sqrt(G * M_sun / r2)
    va2 = np.sqrt(G * M_sun * (2/r2 - 1/a2))
    delta_v3 = abs(v2 - va2)
    
    # Total delta-v
    total_delta_v = delta_v1 + delta_v2 + delta_v3
    
    # Transfer time
    time1 = calculate_transfer_time(a1) / day_to_sec  # First half-ellipse
    time2 = calculate_transfer_time(a2) / day_to_sec  # Second half-ellipse
    total_time = time1 + time2
    
    bielliptic_results.append({
        'r3': r3 / AU,
        'r3_factor': r3 / r2,
        'delta_v1': delta_v1,
        'delta_v2': delta_v2,
        'delta_v3': delta_v3,
        'total_delta_v': total_delta_v,
        'time1': time1,
        'time2': time2,
        'total_time': total_time
    })

# 2. Fast Transfer (Higher Energy)
# A faster transfer uses a more energetic ellipse with the same perihelion and aphelion
# but higher velocity, resulting in shorter transfer time but higher delta-v

# Define a range of transfer times as fractions of the Hohmann transfer time
time_factors = np.linspace(0.5, 0.9, 10)  # 50% to 90% of Hohmann time
fast_transfer_results = []

for factor in time_factors:
    # Target transfer time
    target_time = hohmann_transfer_time * factor * day_to_sec
    
    # Iteratively find the semi-major axis that gives this transfer time
    # Start with Hohmann transfer semi-major axis
    a_guess = hohmann_a
    a_min = r1 / 2  # Lower bound (very high energy)
    a_max = hohmann_a  # Upper bound (Hohmann)
    
    # Binary search to find the right semi-major axis
    while a_max - a_min > 1e3:  # 1 km precision
        a_guess = (a_min + a_max) / 2
        time_guess = calculate_transfer_time(a_guess)
        
        if time_guess > target_time:
            a_max = a_guess
        else:
            a_min = a_guess
    
    # Calculate the eccentricity
    e_guess = 1 - r1 / a_guess  # Assuming perihelion at Earth
    
    # Calculate delta-v
    delta_v1, delta_v2, total_delta_v = calculate_delta_v(r1, r2, a_guess)
    
    fast_transfer_results.append({
        'time_factor': factor,
        'transfer_time': target_time / day_to_sec,
        'semi_major_axis': a_guess / AU,
        'eccentricity': e_guess,
        'delta_v1': delta_v1,
        'delta_v2': delta_v2,
        'total_delta_v': total_delta_v
    })

# 3. Low-Energy Transfer (Patched Conic)
# A low-energy transfer uses a more complex trajectory that takes advantage of
# gravitational assists or weak stability boundaries

# For simplicity, we'll approximate this with a longer, lower-energy transfer
# that uses less delta-v but takes more time

# Define a range of transfer times as multiples of the Hohmann transfer time
time_factors = np.linspace(1.1, 2.0, 10)  # 110% to 200% of Hohmann time
low_energy_results = []

for factor in time_factors:
    # Target transfer time
    target_time = hohmann_transfer_time * factor * day_to_sec
    
    # Iteratively find the semi-major axis that gives this transfer time
    # Start with Hohmann transfer semi-major axis
    a_guess = hohmann_a
    a_min = hohmann_a  # Lower bound (Hohmann)
    a_max = 10 * hohmann_a  # Upper bound (very low energy)
    
    # Binary search to find the right semi-major axis
    while a_max - a_min > 1e3:  # 1 km precision
        a_guess = (a_min + a_max) / 2
        time_guess = calculate_transfer_time(a_guess)
        
        if time_guess < target_time:
            a_min = a_guess
        else:
            a_max = a_guess
    
    # Calculate the eccentricity
    e_guess = 1 - r1 / a_guess  # Assuming perihelion at Earth
    
    # Calculate delta-v
    delta_v1, delta_v2, total_delta_v = calculate_delta_v(r1, r2, a_guess)
    
    low_energy_results.append({
        'time_factor': factor,
        'transfer_time': target_time / day_to_sec,
        'semi_major_axis': a_guess / AU,
        'eccentricity': e_guess,
        'delta_v1': delta_v1,
        'delta_v2': delta_v2,
        'total_delta_v': total_delta_v
    })

# Create visualizations

# 1. Delta-v vs Transfer Time for all methods
plt.figure(figsize=(12, 8))

# Hohmann transfer (single point)
plt.scatter(hohmann_transfer_time, hohmann_delta_v/1000, 
            color='green', s=100, label='Hohmann Transfer')

# Bi-elliptic transfers
bielliptic_times = [r['total_time'] for r in bielliptic_results]
bielliptic_dv = [r['total_delta_v']/1000 for r in bielliptic_results]
plt.plot(bielliptic_times, bielliptic_dv, 'o-', color='blue', label='Bi-elliptic Transfers')

# Fast transfers
fast_times = [r['transfer_time'] for r in fast_transfer_results]
fast_dv = [r['total_delta_v']/1000 for r in fast_transfer_results]
plt.plot(fast_times, fast_dv, 'o-', color='red', label='Fast Transfers')

# Low-energy transfers
low_energy_times = [r['transfer_time'] for r in low_energy_results]
low_energy_dv = [r['total_delta_v']/1000 for r in low_energy_results]
plt.plot(low_energy_times, low_energy_dv, 'o-', color='purple', label='Low-Energy Transfers')

plt.xlabel('Transfer Time (days)')
plt.ylabel('Total Delta-v (km/s)')
plt.title('Delta-v vs Transfer Time for Different Transfer Methods')
plt.grid(True)
plt.legend()
plt.savefig('delta_v_vs_time.png', dpi=300, bbox_inches='tight')
plt.close()

# 2. Visualize the different transfer orbits
plt.figure(figsize=(12, 12))

# Plot the Sun
plt.scatter([0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's orbit
theta = np.linspace(0, 2*np.pi, 1000)
r_earth = earth_a/AU * (1 - earth_e**2) / (1 + earth_e * np.cos(theta))
x_earth = r_earth * np.cos(theta)
y_earth = r_earth * np.sin(theta)
plt.plot(x_earth, y_earth, color='blue', label='Earth Orbit')

# Plot Mars's orbit
r_mars = mars_a/AU * (1 - mars_e**2) / (1 + mars_e * np.cos(theta))
x_mars = r_mars * np.cos(theta)
y_mars = r_mars * np.sin(theta)
plt.plot(x_mars, y_mars, color='red', label='Mars Orbit')

# Plot Earth and Mars positions
earth_pos_AU = orbital_params['earth']['position']
mars_pos_AU = orbital_params['mars']['position']
plt.scatter(earth_pos_AU[0], earth_pos_AU[1], color='blue', s=100, label='Earth')
plt.scatter(mars_pos_AU[0], mars_pos_AU[1], color='red', s=100, label='Mars')

# Plot Hohmann transfer
a_hohmann = hohmann_params['transfer_orbit']['semi_major_axis_AU']
e_hohmann = hohmann_params['transfer_orbit']['eccentricity']
r_hohmann = a_hohmann * (1 - e_hohmann**2) / (1 + e_hohmann * np.cos(theta))
x_hohmann = r_hohmann * np.cos(theta)
y_hohmann = r_hohmann * np.sin(theta)
plt.plot(x_hohmann, y_hohmann, color='green', label='Hohmann Transfer')

# Plot a fast transfer (using the middle one from our calculations)
mid_idx = len(fast_transfer_results) // 2
fast_a = fast_transfer_results[mid_idx]['semi_major_axis']
fast_e = fast_transfer_results[mid_idx]['eccentricity']
r_fast = fast_a * (1 - fast_e**2) / (1 + fast_e * np.cos(theta))
x_fast = r_fast * np.cos(theta)
y_fast = r_fast * np.sin(theta)
plt.plot(x_fast, y_fast, color='red', linestyle='--', 
         label=f'Fast Transfer ({fast_transfer_results[mid_idx]["time_factor"]:.1f}x faster)')

# Plot a low-energy transfer (using the middle one from our calculations)
mid_idx = len(low_energy_results) // 2
low_e_a = low_energy_results[mid_idx]['semi_major_axis']
low_e_e = low_energy_results[mid_idx]['eccentricity']
r_low_e = low_e_a * (1 - low_e_e**2) / (1 + low_e_e * np.cos(theta))
x_low_e = r_low_e * np.cos(theta)
y_low_e = r_low_e * np.sin(theta)
plt.plot(x_low_e, y_low_e, color='purple', linestyle='--', 
         label=f'Low-Energy Transfer ({low_energy_results[mid_idx]["time_factor"]:.1f}x slower)')

# Plot a bi-elliptic transfer (first half)
mid_idx = len(bielliptic_results) // 2
r3 = bielliptic_results[mid_idx]['r3'] * AU
a1 = (r1 + r3) / 2 / AU
e1 = (r3 - r1) / (r3 + r1)
r_bi1 = a1 * (1 - e1**2) / (1 + e1 * np.cos(theta))
x_bi1 = r_bi1 * np.cos(theta)
y_bi1 = r_bi1 * np.sin(theta)

# Plot a bi-elliptic transfer (second half)
a2 = (r3 + r2) / 2 / AU
e2 = (r3 - r2) / (r3 + r2)
r_bi2 = a2 * (1 - e2**2) / (1 + e2 * np.cos(theta))
x_bi2 = r_bi2 * np.cos(theta)
y_bi2 = r_bi2 * np.sin(theta)

plt.plot(x_bi1, y_bi1, color='blue', linestyle='--', 
         label=f'Bi-elliptic Transfer (r3={bielliptic_results[mid_idx]["r3_factor"]:.1f}×Mars)')
plt.plot(x_bi2, y_bi2, color='blue', linestyle=':')

# Set labels and title
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title('Comparison of Different Transfer Orbits')
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.savefig('transfer_orbit_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

# 3. Create a 3D visualization of the optimal transfer
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the Sun
ax.scatter([0], [0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's orbit
theta = np.linspace(0, 2*np.pi, 1000)
r_earth = earth_a/AU * (1 - earth_e**2) / (1 + earth_e * np.cos(theta))
x_earth = r_earth * np.cos(theta)
y_earth = r_earth * np.sin(theta) * np.cos(earth_i)
z_earth = r_earth * np.sin(theta) * np.sin(earth_i)
ax.plot(x_earth, y_earth, z_earth, color='blue', label='Earth Orbit')

# Plot Mars's orbit
r_mars = mars_a/AU * (1 - mars_e**2) / (1 + mars_e * np.cos(theta))
x_mars = r_mars * np.cos(theta)
y_mars = r_mars * np.sin(theta) * np.cos(mars_i)
z_mars = r_mars * np.sin(theta) * np.sin(mars_i)
ax.plot(x_mars, y_mars, z_mars, color='red', label='Mars Orbit')

# Plot Earth and Mars positions
earth_pos_AU = orbital_params['earth']['position']
mars_pos_AU = orbital_params['mars']['position']
ax.scatter(earth_pos_AU[0], earth_pos_AU[1], earth_pos_AU[2], color='blue', s=100, label='Earth')
ax.scatter(mars_pos_AU[0], mars_pos_AU[1], mars_pos_AU[2], color='red', s=100, label='Mars')

# Plot Hohmann transfer (simplified to be in the ecliptic plane)
a_hohmann = hohmann_params['transfer_orbit']['semi_major_axis_AU']
e_hohmann = hohmann_params['transfer_orbit']['eccentricity']
r_hohmann = a_hohmann * (1 - e_hohmann**2) / (1 + e_hohmann * np.cos(theta))
x_hohmann = r_hohmann * np.cos(theta)
y_hohmann = r_hohmann * np.sin(theta)
z_hohmann = np.zeros_like(theta)  # Simplified to be in the ecliptic plane
ax.plot(x_hohmann, y_hohmann, z_hohmann, color='green', label='Hohmann Transfer')

# Set labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title('3D Visualization of Optimal Transfer Orbit')
ax.legend()
plt.savefig('optimal_transfer_3d.png', dpi=300, bbox_inches='tight')
plt.close()

# Save the results to a file
with open('mission_scenarios.txt', 'w') as f:
    f.write("Comparison of Different Mission Scenarios\n")
    f.write("=======================================\n\n")
    
    f.write("1. Hohmann Transfer (Baseline)\n")
    f.write("-----------------------------\n")
    f.write(f"Transfer time: {hohmann_transfer_time:.2f} days\n")
    f.write(f"Total delta-v: {hohmann_delta_v/1000:.2f} km/s\n")
    f.write(f"Semi-major axis: {hohmann_a/AU:.6f} AU\n")
    f.write(f"Eccentricity: {hohmann_e:.6f}\n\n")
    
    f.write("2. Bi-elliptic Transfers\n")
    f.write("----------------------\n")
    f.write("r3 factor | Transfer Time | Total Delta-v | Delta-v Difference\n")
    f.write("---------|---------------|--------------|------------------\n")
    for r in bielliptic_results:
        dv_diff = (r['total_delta_v'] - hohmann_delta_v) / 1000
        dv_diff_percent = (r['total_delta_v'] - hohmann_delta_v) / hohmann_delta_v * 100
        f.write(f"{r['r3_factor']:.2f}x     | {r['total_time']:.2f} days    | {r['total_delta_v']/1000:.2f} km/s     | {dv_diff:.2f} km/s ({dv_diff_percent:.2f}%)\n")
    f.write("\n")
    
    f.write("3. Fast Transfers\n")
    f.write("---------------\n")
    f.write("Time Factor | Transfer Time | Total Delta-v | Delta-v Difference\n")
    f.write("-----------|---------------|--------------|------------------\n")
    for r in fast_transfer_results:
        dv_diff = (r['total_delta_v'] - hohmann_delta_v) / 1000
        dv_diff_percent = (r['total_delta_v'] - hohmann_delta_v) / hohmann_delta_v * 100
        f.write(f"{r['time_factor']:.2f}x     | {r['transfer_time']:.2f} days    | {r['total_delta_v']/1000:.2f} km/s     | {dv_diff:.2f} km/s ({dv_diff_percent:.2f}%)\n")
    f.write("\n")
<response clipped><NOTE>To save on context only part of this file has been shown to you. You should retry this tool after you have searched inside the file with `grep -n` in order to find the line numbers of what you are looking for.</NOTE>