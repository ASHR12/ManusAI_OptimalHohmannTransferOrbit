"""
Determine Hohmann transfer parameters for Earth to Mars transfer.
This script calculates the parameters of the minimum-energy Hohmann transfer orbit
from Earth to Mars based on their current positions as of March 2025.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

# Load the orbital parameters calculated in the previous step
orbital_params = np.load('orbital_params.npy', allow_pickle=True).item()

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

# Constants
G = orbital_params['constants']['G']
M_sun = orbital_params['constants']['M_sun']
AU = orbital_params['constants']['AU']
day_to_sec = orbital_params['constants']['day_to_sec']

# Calculate Hohmann transfer parameters
# For a Hohmann transfer, we need the semi-major axis of the transfer orbit
r1 = np.linalg.norm(earth_pos)  # Earth's distance from Sun
r2 = np.linalg.norm(mars_pos)   # Mars's distance from Sun

# Semi-major axis of the transfer orbit
a_transfer = (r1 + r2) / 2

# Eccentricity of the transfer orbit
e_transfer = (r2 - r1) / (r2 + r1)

# Calculate the velocities at perihelion (Earth departure) and aphelion (Mars arrival)
# Velocity at perihelion (Earth departure)
v_perihelion = np.sqrt(G * M_sun * (2/r1 - 1/a_transfer))

# Velocity at aphelion (Mars arrival)
v_aphelion = np.sqrt(G * M_sun * (2/r2 - 1/a_transfer))

# Calculate Earth's orbital velocity at departure point
v_earth = np.linalg.norm(earth_vel)

# Calculate Mars's orbital velocity at arrival point
v_mars = np.linalg.norm(mars_vel)

# Calculate the delta-v required for departure (Earth to transfer orbit)
delta_v_departure = abs(v_perihelion - v_earth)

# Calculate the delta-v required for arrival (transfer orbit to Mars)
delta_v_arrival = abs(v_aphelion - v_mars)

# Total delta-v required for the mission
total_delta_v = delta_v_departure + delta_v_arrival

# Calculate the transfer time (half the period of the transfer orbit)
transfer_time_sec = np.pi * np.sqrt(a_transfer**3 / (G * M_sun))
transfer_time_days = transfer_time_sec / day_to_sec

# Calculate the phase angle required between Earth and Mars for optimal launch
# This is the angle that Mars needs to be ahead of Earth at launch
phase_angle_rad = np.pi * (1 - np.sqrt((r1 + r2)/(2 * r2)))
phase_angle_deg = np.degrees(phase_angle_rad)

# Calculate the current phase angle between Earth and Mars
# Get the angular positions of Earth and Mars
earth_angle = np.arctan2(earth_pos[1], earth_pos[0])
mars_angle = np.arctan2(mars_pos[1], mars_pos[0])

# Calculate the current phase angle
current_phase_angle_rad = (mars_angle - earth_angle) % (2 * np.pi)
current_phase_angle_deg = np.degrees(current_phase_angle_rad)

# Calculate the wait time until the next optimal launch window
# Earth's angular velocity (rad/day)
earth_angular_vel = 2 * np.pi / orbital_params['earth']['period_days']
# Mars's angular velocity (rad/day)
mars_angular_vel = 2 * np.pi / orbital_params['mars']['period_days']
# Relative angular velocity
relative_angular_vel = earth_angular_vel - mars_angular_vel

# Calculate the angle that needs to be covered
angle_to_cover_rad = (phase_angle_rad - current_phase_angle_rad) % (2 * np.pi)
# If the angle is very small, we might need to wait almost a full synodic period
if angle_to_cover_rad < 0.1:
    angle_to_cover_rad += 2 * np.pi

# Calculate the wait time
wait_time_days = angle_to_cover_rad / relative_angular_vel

# Calculate the synodic period (time between successive oppositions)
synodic_period_days = 2 * np.pi / relative_angular_vel

# Print the results
print("Hohmann Transfer Parameters:")
print(f"Transfer orbit semi-major axis: {a_transfer/AU:.6f} AU ({a_transfer/1000:.0f} km)")
print(f"Transfer orbit eccentricity: {e_transfer:.6f}")
print(f"Departure delta-v: {delta_v_departure/1000:.2f} km/s")
print(f"Arrival delta-v: {delta_v_arrival/1000:.2f} km/s")
print(f"Total delta-v: {total_delta_v/1000:.2f} km/s")
print(f"Transfer time: {transfer_time_days:.2f} days ({transfer_time_days/365.25:.2f} years)")
print(f"Optimal phase angle: {phase_angle_deg:.2f} degrees")
print(f"Current phase angle: {current_phase_angle_deg:.2f} degrees")
print(f"Synodic period: {synodic_period_days:.2f} days ({synodic_period_days/365.25:.2f} years)")
print(f"Wait time until next optimal launch window: {wait_time_days:.2f} days ({wait_time_days/365.25:.2f} years)")

# Save the results to a file
with open('hohmann_transfer_parameters.txt', 'w') as f:
    f.write("Hohmann Transfer Parameters:\n")
    f.write(f"Transfer orbit semi-major axis: {a_transfer/AU:.6f} AU ({a_transfer/1000:.0f} km)\n")
    f.write(f"Transfer orbit eccentricity: {e_transfer:.6f}\n")
    f.write(f"Departure delta-v: {delta_v_departure/1000:.2f} km/s\n")
    f.write(f"Arrival delta-v: {delta_v_arrival/1000:.2f} km/s\n")
    f.write(f"Total delta-v: {total_delta_v/1000:.2f} km/s\n")
    f.write(f"Transfer time: {transfer_time_days:.2f} days ({transfer_time_days/365.25:.2f} years)\n")
    f.write(f"Optimal phase angle: {phase_angle_deg:.2f} degrees\n")
    f.write(f"Current phase angle: {current_phase_angle_deg:.2f} degrees\n")
    f.write(f"Synodic period: {synodic_period_days:.2f} days ({synodic_period_days/365.25:.2f} years)\n")
    f.write(f"Wait time until next optimal launch window: {wait_time_days:.2f} days ({wait_time_days/365.25:.2f} years)\n")

# Create a visualization of the transfer orbit
def plot_orbit(a, e, color, label, n_points=1000):
    """Plot an orbit given its semi-major axis and eccentricity."""
    theta = np.linspace(0, 2*np.pi, n_points)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

# Create a 2D plot of the orbits in the ecliptic plane
plt.figure(figsize=(12, 12))

# Plot the Sun
plt.scatter([0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's orbit
earth_a_AU = earth_a / AU
x_earth, y_earth = plot_orbit(earth_a_AU, earth_e, 'blue', 'Earth Orbit')
plt.plot(x_earth, y_earth, color='blue', label='Earth Orbit')

# Plot Mars's orbit
mars_a_AU = mars_a / AU
x_mars, y_mars = plot_orbit(mars_a_AU, mars_e, 'red', 'Mars Orbit')
plt.plot(x_mars, y_mars, color='red', label='Mars Orbit')

# Plot the transfer orbit
transfer_a_AU = a_transfer / AU
x_transfer, y_transfer = plot_orbit(transfer_a_AU, e_transfer, 'green', 'Transfer Orbit')
plt.plot(x_transfer, y_transfer, color='green', label='Hohmann Transfer Orbit')

# Plot Earth and Mars at their current positions
earth_pos_AU = earth_pos / AU
mars_pos_AU = mars_pos / AU
plt.scatter(earth_pos_AU[0], earth_pos_AU[1], color='blue', s=100, label='Earth')
plt.scatter(mars_pos_AU[0], mars_pos_AU[1], color='red', s=100, label='Mars')

# Set labels and title
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title('Earth-Mars Hohmann Transfer Orbit')
plt.grid(True)
plt.axis('equal')
plt.legend()

# Save the figure
plt.savefig('hohmann_transfer_orbit.png', dpi=300, bbox_inches='tight')
plt.close()

# Create a 3D visualization
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the Sun
ax.scatter([0], [0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's orbit
theta = np.linspace(0, 2*np.pi, 1000)
r_earth = earth_a_AU * (1 - earth_e**2) / (1 + earth_e * np.cos(theta))
x_earth = r_earth * np.cos(theta)
y_earth = r_earth * np.sin(theta) * np.cos(earth_i)
z_earth = r_earth * np.sin(theta) * np.sin(earth_i)
ax.plot(x_earth, y_earth, z_earth, color='blue', label='Earth Orbit')

# Plot Mars's orbit
r_mars = mars_a_AU * (1 - mars_e**2) / (1 + mars_e * np.cos(theta))
x_mars = r_mars * np.cos(theta)
y_mars = r_mars * np.sin(theta) * np.cos(mars_i)
z_mars = r_mars * np.sin(theta) * np.sin(mars_i)
ax.plot(x_mars, y_mars, z_mars, color='red', label='Mars Orbit')

# Plot the transfer orbit (simplified to be in the ecliptic plane)
r_transfer = transfer_a_AU * (1 - e_transfer**2) / (1 + e_transfer * np.cos(theta))
x_transfer = r_transfer * np.cos(theta)
y_transfer = r_transfer * np.sin(theta)
z_transfer = np.zeros_like(theta)  # Simplified to be in the ecliptic plane
ax.plot(x_transfer, y_transfer, z_transfer, color='green', label='Hohmann Transfer Orbit')

# Plot Earth and Mars at their current positions
earth_pos_AU = earth_pos / AU
mars_pos_AU = mars_pos / AU
ax.scatter(earth_pos_AU[0], earth_pos_AU[1], earth_pos_AU[2], color='blue', s=100, label='Earth')
ax.scatter(mars_pos_AU[0], mars_pos_AU[1], mars_pos_AU[2], color='red', s=100, label='Mars')

# Set labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title('Earth-Mars Hohmann Transfer Orbit (3D)')
ax.legend()

# Save the figure
plt.savefig('hohmann_transfer_orbit_3d.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nHohmann transfer parameters calculation complete.")
print("Results saved to 'hohmann_transfer_parameters.txt'")
print("Visualizations saved to 'hohmann_transfer_orbit.png' and 'hohmann_transfer_orbit_3d.png'")

# Save the Hohmann transfer parameters for use in other scripts
hohmann_params = {
    'transfer_orbit': {
        'semi_major_axis': a_transfer,
        'semi_major_axis_AU': a_transfer / AU,
        'eccentricity': e_transfer
    },
    'delta_v': {
        'departure': delta_v_departure,
        'arrival': delta_v_arrival,
        'total': total_delta_v
    },
    'timing': {
        'transfer_time_days': transfer_time_days,
        'optimal_phase_angle_deg': phase_angle_deg,
        'current_phase_angle_deg': current_phase_angle_deg,
        'wait_time_days': wait_time_days,
        'synodic_period_days': synodic_period_days
    }
}

np.save('hohmann_params.npy', hohmann_params)
