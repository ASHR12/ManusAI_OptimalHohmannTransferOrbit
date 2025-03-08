"""
Calculate orbital parameters for Earth and Mars based on March 2025 positions.
This script uses the astropy library to retrieve accurate planetary ephemeris data.
"""

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Set the ephemeris to use
solar_system_ephemeris.set('builtin')

# Define the time for March 2025
time = Time('2025-03-08')  # Current date as specified in the problem

# Get the positions of Earth and Mars in the Solar System Barycentric frame
earth_pos = get_body_barycentric('earth', time)
mars_pos = get_body_barycentric('mars', time)

# Get positions at a slightly different time to calculate velocities
time_delta = 1 * u.hour
time_plus = time + time_delta
earth_pos_plus = get_body_barycentric('earth', time_plus)
mars_pos_plus = get_body_barycentric('mars', time_plus)

# Calculate velocities using finite difference
earth_vel = (earth_pos_plus - earth_pos) / time_delta.to(u.day)
mars_vel = (mars_pos_plus - mars_pos) / time_delta.to(u.day)

# Convert to Cartesian coordinates (in AU and AU/day)
earth_pos_cart = np.array([earth_pos.x.to(u.AU).value, earth_pos.y.to(u.AU).value, earth_pos.z.to(u.AU).value])
mars_pos_cart = np.array([mars_pos.x.to(u.AU).value, mars_pos.y.to(u.AU).value, mars_pos.z.to(u.AU).value])
earth_vel_cart = np.array([earth_vel.x.to(u.AU/u.day).value, earth_vel.y.to(u.AU/u.day).value, earth_vel.z.to(u.AU/u.day).value])
mars_vel_cart = np.array([mars_vel.x.to(u.AU/u.day).value, mars_vel.y.to(u.AU/u.day).value, mars_vel.z.to(u.AU/u.day).value])

# Constants
G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
M_sun = 1.989e30  # Mass of the Sun in kg
AU = 1.496e11    # 1 AU in meters
day_to_sec = 86400  # Conversion from days to seconds

# Convert positions to meters and velocities to m/s for calculations
earth_pos_m = earth_pos_cart * AU
mars_pos_m = mars_pos_cart * AU
earth_vel_ms = earth_vel_cart * AU / day_to_sec
mars_vel_ms = mars_vel_cart * AU / day_to_sec

# Calculate orbital parameters for Earth
earth_r = np.linalg.norm(earth_pos_m)  # Distance from Sun in meters
earth_v = np.linalg.norm(earth_vel_ms)  # Orbital velocity in m/s
earth_h = np.cross(earth_pos_m, earth_vel_ms)  # Specific angular momentum
earth_h_mag = np.linalg.norm(earth_h)
earth_e_vec = np.cross(earth_vel_ms, earth_h) / (G * M_sun) - earth_pos_m / earth_r
earth_e = np.linalg.norm(earth_e_vec)  # Eccentricity
earth_a = earth_h_mag**2 / (G * M_sun * (1 - earth_e**2))  # Semi-major axis
earth_p = earth_h_mag**2 / (G * M_sun)  # Semi-latus rectum
earth_i = np.arccos(earth_h[2] / earth_h_mag)  # Inclination

# Calculate orbital parameters for Mars
mars_r = np.linalg.norm(mars_pos_m)  # Distance from Sun in meters
mars_v = np.linalg.norm(mars_vel_ms)  # Orbital velocity in m/s
mars_h = np.cross(mars_pos_m, mars_vel_ms)  # Specific angular momentum
mars_h_mag = np.linalg.norm(mars_h)
mars_e_vec = np.cross(mars_vel_ms, mars_h) / (G * M_sun) - mars_pos_m / mars_r
mars_e = np.linalg.norm(mars_e_vec)  # Eccentricity
mars_a = mars_h_mag**2 / (G * M_sun * (1 - mars_e**2))  # Semi-major axis
mars_p = mars_h_mag**2 / (G * M_sun)  # Semi-latus rectum
mars_i = np.arccos(mars_h[2] / mars_h_mag)  # Inclination

# Calculate orbital periods
earth_period = 2 * np.pi * np.sqrt(earth_a**3 / (G * M_sun))  # in seconds
mars_period = 2 * np.pi * np.sqrt(mars_a**3 / (G * M_sun))  # in seconds

# Convert to more readable units
earth_a_AU = earth_a / AU
mars_a_AU = mars_a / AU
earth_period_days = earth_period / day_to_sec
mars_period_days = mars_period / day_to_sec

# Print the results
print(f"Date: {time.iso}")
print("\nEarth Orbital Parameters:")
print(f"Position (AU): [{earth_pos_cart[0]:.6f}, {earth_pos_cart[1]:.6f}, {earth_pos_cart[2]:.6f}]")
print(f"Velocity (AU/day): [{earth_vel_cart[0]:.6f}, {earth_vel_cart[1]:.6f}, {earth_vel_cart[2]:.6f}]")
print(f"Semi-major axis: {earth_a_AU:.6f} AU ({earth_a/1000:.0f} km)")
print(f"Eccentricity: {earth_e:.6f}")
print(f"Inclination: {np.degrees(earth_i):.6f} degrees")
print(f"Orbital period: {earth_period_days:.2f} days ({earth_period_days/365.25:.6f} years)")

print("\nMars Orbital Parameters:")
print(f"Position (AU): [{mars_pos_cart[0]:.6f}, {mars_pos_cart[1]:.6f}, {mars_pos_cart[2]:.6f}]")
print(f"Velocity (AU/day): [{mars_vel_cart[0]:.6f}, {mars_vel_cart[1]:.6f}, {mars_vel_cart[2]:.6f}]")
print(f"Semi-major axis: {mars_a_AU:.6f} AU ({mars_a/1000:.0f} km)")
print(f"Eccentricity: {mars_e:.6f}")
print(f"Inclination: {np.degrees(mars_i):.6f} degrees")
print(f"Orbital period: {mars_period_days:.2f} days ({mars_period_days/365.25:.6f} years)")

# Save the orbital parameters to a file
with open('orbital_parameters.txt', 'w') as f:
    f.write(f"Date: {time.iso}\n")
    f.write("\nEarth Orbital Parameters:\n")
    f.write(f"Position (AU): [{earth_pos_cart[0]:.6f}, {earth_pos_cart[1]:.6f}, {earth_pos_cart[2]:.6f}]\n")
    f.write(f"Velocity (AU/day): [{earth_vel_cart[0]:.6f}, {earth_vel_cart[1]:.6f}, {earth_vel_cart[2]:.6f}]\n")
    f.write(f"Semi-major axis: {earth_a_AU:.6f} AU ({earth_a/1000:.0f} km)\n")
    f.write(f"Eccentricity: {earth_e:.6f}\n")
    f.write(f"Inclination: {np.degrees(earth_i):.6f} degrees\n")
    f.write(f"Orbital period: {earth_period_days:.2f} days ({earth_period_days/365.25:.6f} years)\n")
    
    f.write("\nMars Orbital Parameters:\n")
    f.write(f"Position (AU): [{mars_pos_cart[0]:.6f}, {mars_pos_cart[1]:.6f}, {mars_pos_cart[2]:.6f}]\n")
    f.write(f"Velocity (AU/day): [{mars_vel_cart[0]:.6f}, {mars_vel_cart[1]:.6f}, {mars_vel_cart[2]:.6f}]\n")
    f.write(f"Semi-major axis: {mars_a_AU:.6f} AU ({mars_a/1000:.0f} km)\n")
    f.write(f"Eccentricity: {mars_e:.6f}\n")
    f.write(f"Inclination: {np.degrees(mars_i):.6f} degrees\n")
    f.write(f"Orbital period: {mars_period_days:.2f} days ({mars_period_days/365.25:.6f} years)\n")

# Create a visualization of the orbits
def plot_orbit(a, e, i, color, label):
    """Plot an orbit given its semi-major axis, eccentricity, and inclination."""
    theta = np.linspace(0, 2*np.pi, 1000)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta) * np.cos(i)
    z = r * np.sin(theta) * np.sin(i)
    return x, y, z

# Create a 3D plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the Sun
ax.scatter([0], [0], [0], color='yellow', s=100, label='Sun')

# Plot Earth's orbit
x_earth, y_earth, z_earth = plot_orbit(earth_a_AU, earth_e, earth_i, 'blue', 'Earth')
ax.plot(x_earth, y_earth, z_earth, color='blue', label='Earth Orbit')
ax.scatter(earth_pos_cart[0], earth_pos_cart[1], earth_pos_cart[2], color='blue', s=50, label='Earth')

# Plot Mars's orbit
x_mars, y_mars, z_mars = plot_orbit(mars_a_AU, mars_e, mars_i, 'red', 'Mars')
ax.plot(x_mars, y_mars, z_mars, color='red', label='Mars Orbit')
ax.scatter(mars_pos_cart[0], mars_pos_cart[1], mars_pos_cart[2], color='red', s=50, label='Mars')

# Set labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title(f'Earth and Mars Orbits - {time.iso}')
ax.legend()

# Set equal aspect ratio
max_range = np.array([x_mars.max()-x_mars.min(), y_mars.max()-y_mars.min(), z_mars.max()-z_mars.min()]).max() / 2.0
mid_x = (x_mars.max()+x_mars.min()) * 0.5
mid_y = (y_mars.max()+y_mars.min()) * 0.5
mid_z = (z_mars.max()+z_mars.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Save the figure
plt.savefig('planetary_orbits.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nOrbital parameters calculation complete. Results saved to 'orbital_parameters.txt'")
print("Visualization saved to 'planetary_orbits.png'")

# Store the calculated parameters as global variables for use in other scripts
orbital_params = {
    'earth': {
        'position': earth_pos_cart,
        'velocity': earth_vel_cart,
        'position_m': earth_pos_m,
        'velocity_ms': earth_vel_ms,
        'semi_major_axis': earth_a,
        'semi_major_axis_AU': earth_a_AU,
        'eccentricity': earth_e,
        'inclination': earth_i,
        'period_days': earth_period_days
    },
    'mars': {
        'position': mars_pos_cart,
        'velocity': mars_vel_cart,
        'position_m': mars_pos_m,
        'velocity_ms': mars_vel_ms,
        'semi_major_axis': mars_a,
        'semi_major_axis_AU': mars_a_AU,
        'eccentricity': mars_e,
        'inclination': mars_i,
        'period_days': mars_period_days
    },
    'constants': {
        'G': G,
        'M_sun': M_sun,
        'AU': AU,
        'day_to_sec': day_to_sec
    }
}

# Save the orbital parameters as a numpy array for later use
np.save('orbital_params.npy', orbital_params)
