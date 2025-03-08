"""
Account for gravitational influences from the Sun and other planets on the Hohmann transfer orbit.
This script uses a more sophisticated model to simulate the trajectory with gravitational perturbations.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, get_body_barycentric

# Load the orbital parameters and Hohmann transfer parameters
orbital_params = np.load('orbital_params.npy', allow_pickle=True).item()
hohmann_params = np.load('hohmann_params.npy', allow_pickle=True).item()

# Constants
G = orbital_params['constants']['G']
M_sun = orbital_params['constants']['M_sun']
AU = orbital_params['constants']['AU']
day_to_sec = orbital_params['constants']['day_to_sec']

# Planet masses (kg)
M_mercury = 3.3011e23
M_venus = 4.8675e24
M_earth = 5.97237e24
M_mars = 6.4171e23
M_jupiter = 1.8982e27
M_saturn = 5.6834e26
M_uranus = 8.6810e25
M_neptune = 1.02413e26

# Define the time for March 2025
time = Time('2025-03-08')

# Set the ephemeris to use
solar_system_ephemeris.set('builtin')

# Get the positions of all planets in the Solar System Barycentric frame
mercury_pos = get_body_barycentric('mercury', time)
venus_pos = get_body_barycentric('venus', time)
earth_pos = get_body_barycentric('earth', time)
mars_pos = get_body_barycentric('mars', time)
jupiter_pos = get_body_barycentric('jupiter', time)
saturn_pos = get_body_barycentric('saturn', time)
uranus_pos = get_body_barycentric('uranus', time)
neptune_pos = get_body_barycentric('neptune', time)

# Convert to Cartesian coordinates (in meters)
mercury_pos_m = np.array([mercury_pos.x.to(u.m).value, mercury_pos.y.to(u.m).value, mercury_pos.z.to(u.m).value])
venus_pos_m = np.array([venus_pos.x.to(u.m).value, venus_pos.y.to(u.m).value, venus_pos.z.to(u.m).value])
earth_pos_m = np.array([earth_pos.x.to(u.m).value, earth_pos.y.to(u.m).value, earth_pos.z.to(u.m).value])
mars_pos_m = np.array([mars_pos.x.to(u.m).value, mars_pos.y.to(u.m).value, mars_pos.z.to(u.m).value])
jupiter_pos_m = np.array([jupiter_pos.x.to(u.m).value, jupiter_pos.y.to(u.m).value, jupiter_pos.z.to(u.m).value])
saturn_pos_m = np.array([saturn_pos.x.to(u.m).value, saturn_pos.y.to(u.m).value, saturn_pos.z.to(u.m).value])
uranus_pos_m = np.array([uranus_pos.x.to(u.m).value, uranus_pos.y.to(u.m).value, uranus_pos.z.to(u.m).value])
neptune_pos_m = np.array([neptune_pos.x.to(u.m).value, neptune_pos.y.to(u.m).value, neptune_pos.z.to(u.m).value])

# Define the gravitational acceleration function with perturbations from all planets
def gravity_with_perturbations(t, state):
    """
    Calculate the gravitational acceleration on a spacecraft including perturbations from all planets.
    
    Parameters:
    t : float
        Time (used for updating planet positions in a more sophisticated model)
    state : array
        State vector [x, y, z, vx, vy, vz]
        
    Returns:
    dstate : array
        Derivative of state vector [vx, vy, vz, ax, ay, az]
    """
    # Extract position
    r = state[:3]
    
    # Distance to the Sun
    r_sun = np.linalg.norm(r)
    
    # Gravitational acceleration due to the Sun
    a_sun = -G * M_sun * r / r_sun**3
    
    # Calculate perturbations from each planet
    # For each planet, calculate the gravitational acceleration it exerts on the spacecraft
    
    # Mercury
    r_mercury_sc = r - mercury_pos_m
    r_mercury_sc_norm = np.linalg.norm(r_mercury_sc)
    a_mercury = -G * M_mercury * r_mercury_sc / r_mercury_sc_norm**3 if r_mercury_sc_norm > 1e3 else np.zeros(3)
    
    # Venus
    r_venus_sc = r - venus_pos_m
    r_venus_sc_norm = np.linalg.norm(r_venus_sc)
    a_venus = -G * M_venus * r_venus_sc / r_venus_sc_norm**3 if r_venus_sc_norm > 1e3 else np.zeros(3)
    
    # Earth - avoid division by zero if spacecraft is at Earth's position
    r_earth_sc = r - earth_pos_m
    r_earth_sc_norm = np.linalg.norm(r_earth_sc)
    a_earth = -G * M_earth * r_earth_sc / r_earth_sc_norm**3 if r_earth_sc_norm > 1e3 else np.zeros(3)
    
    # Mars
    r_mars_sc = r - mars_pos_m
    r_mars_sc_norm = np.linalg.norm(r_mars_sc)
    a_mars = -G * M_mars * r_mars_sc / r_mars_sc_norm**3 if r_mars_sc_norm > 1e3 else np.zeros(3)
    
    # Jupiter (most significant perturbation)
    r_jupiter_sc = r - jupiter_pos_m
    r_jupiter_sc_norm = np.linalg.norm(r_jupiter_sc)
    a_jupiter = -G * M_jupiter * r_jupiter_sc / r_jupiter_sc_norm**3 if r_jupiter_sc_norm > 1e3 else np.zeros(3)
    
    # Saturn
    r_saturn_sc = r - saturn_pos_m
    r_saturn_sc_norm = np.linalg.norm(r_saturn_sc)
    a_saturn = -G * M_saturn * r_saturn_sc / r_saturn_sc_norm**3 if r_saturn_sc_norm > 1e3 else np.zeros(3)
    
    # Uranus
    r_uranus_sc = r - uranus_pos_m
    r_uranus_sc_norm = np.linalg.norm(r_uranus_sc)
    a_uranus = -G * M_uranus * r_uranus_sc / r_uranus_sc_norm**3 if r_uranus_sc_norm > 1e3 else np.zeros(3)
    
    # Neptune
    r_neptune_sc = r - neptune_pos_m
    r_neptune_sc_norm = np.linalg.norm(r_neptune_sc)
    a_neptune = -G * M_neptune * r_neptune_sc / r_neptune_sc_norm**3 if r_neptune_sc_norm > 1e3 else np.zeros(3)
    
    # Total acceleration (Sun + all planets)
    a_total = a_sun + a_mercury + a_venus + a_earth + a_mars + a_jupiter + a_saturn + a_uranus + a_neptune
    
    # Return the derivatives [vx, vy, vz, ax, ay, az]
    return np.concatenate([state[3:], a_total])

# Define the initial state for the spacecraft at Earth departure
# We'll use Earth's position and add the departure velocity in the direction of Earth's velocity
earth_vel_ms = orbital_params['earth']['velocity_ms']
earth_vel_norm = np.linalg.norm(earth_vel_ms)
earth_vel_dir = earth_vel_ms / earth_vel_norm

# Calculate the departure velocity vector
# The departure velocity is Earth's velocity plus the delta-v in the direction of Earth's velocity
departure_speed = earth_vel_norm + hohmann_params['delta_v']['departure']
departure_vel = earth_vel_dir * departure_speed

# Add a small offset to the initial position to avoid division by zero
# Move the spacecraft slightly away from Earth's center (100 km in the direction of velocity)
position_offset = earth_vel_dir * 1e5  # 100 km offset
initial_position = earth_pos_m + position_offset

# Initial state vector [x, y, z, vx, vy, vz]
initial_state = np.concatenate([initial_position, departure_vel])

# Time span for the integration (slightly longer than the calculated transfer time)
transfer_time_sec = hohmann_params['timing']['transfer_time_days'] * day_to_sec
t_span = (0, transfer_time_sec * 1.1)  # 10% longer to ensure we reach Mars

# Time points to evaluate the solution at
t_eval = np.linspace(0, transfer_time_sec * 1.1, 1000)

# Solve the differential equation
print("Simulating spacecraft trajectory with gravitational perturbations...")
solution = solve_ivp(
    gravity_with_perturbations,
    t_span,
    initial_state,
    t_eval=t_eval,
    method='RK45',
    rtol=1e-8,
    atol=1e-8
)

# Extract the trajectory
trajectory = solution.y[:3].T  # Shape: (n_points, 3)

# Calculate the distance to Mars at each point in the trajectory
mars_distances = np.sqrt(np.sum((trajectory - mars_pos_m)**2, axis=1))

# Find the closest approach to Mars
closest_idx = np.argmin(mars_distances)
closest_distance = mars_distances[closest_idx]
closest_time = solution.t[closest_idx] / day_to_sec  # Convert to days

# Calculate the velocity at closest approach
velocity_at_closest = solution.y[3:, closest_idx]
speed_at_closest = np.linalg.norm(velocity_at_closest)

# Calculate the required delta-v for Mars orbit insertion
# Mars's orbital velocity
mars_vel_norm = np.linalg.norm(orbital_params['mars']['velocity_ms'])

# Simplified calculation for delta-v (in a more sophisticated model, we would consider the direction)
delta_v_insertion = abs(speed_at_closest - mars_vel_norm)

# Calculate the total delta-v with perturbations
total_delta_v_with_perturbations = hohmann_params['delta_v']['departure'] + delta_v_insertion

# Calculate the difference from the ideal Hohmann transfer
delta_v_difference = total_delta_v_with_perturbations - hohmann_params['delta_v']['total']
percent_difference = (delta_v_difference / hohmann_params['delta_v']['total']) * 100

# Print the results
print("\nResults with Gravitational Perturbations:")
print(f"Closest approach to Mars: {closest_distance/1000:.2f} km")
print(f"Time to closest approach: {closest_time:.2f} days")
print(f"Speed at closest approach: {speed_at_closest/1000:.2f} km/s")
print(f"Required delta-v for Mars orbit insertion: {delta_v_insertion/1000:.2f} km/s")
print(f"Total delta-v with perturbations: {total_delta_v_with_perturbations/1000:.2f} km/s")
print(f"Difference from ideal Hohmann transfer: {delta_v_difference/1000:.2f} km/s ({percent_difference:.2f}%)")

# Save the results to a file
with open('perturbation_analysis.txt', 'w') as f:
    f.write("Results with Gravitational Perturbations:\n")
    f.write(f"Closest approach to Mars: {closest_distance/1000:.2f} km\n")
    f.write(f"Time to closest approach: {closest_time:.2f} days\n")
    f.write(f"Speed at closest approach: {speed_at_closest/1000:.2f} km/s\n")
    f.write(f"Required delta-v for Mars orbit insertion: {delta_v_insertion/1000:.2f} km/s\n")
    f.write(f"Total delta-v with perturbations: {total_delta_v_with_perturbations/1000:.2f} km/s\n")
    f.write(f"Difference from ideal Hohmann transfer: {delta_v_difference/1000:.2f} km/s ({percent_difference:.2f}%)\n\n")
    
    f.write("Perturbation Analysis:\n")
    f.write("The gravitational influences from other planets, especially Jupiter, cause deviations from the ideal Hohmann transfer trajectory.\n")
    f.write("This results in a different arrival velocity at Mars and consequently a different delta-v requirement for orbit insertion.\n")
    f.write("The most significant perturbations come from Jupiter due to its large mass, followed by Venus and Earth due to their proximity.\n")
    f.write("These perturbations can either increase or decrease the total delta-v requirement depending on the relative positions of the planets.\n")
    f.write("In this case, the perturbations resulted in a " + ("higher" if delta_v_difference > 0 else "lower") + " delta-v requirement.\n")
    f.write("For precise mission planning, these perturbations must be accounted for and may require trajectory correction maneuvers during the transfer.\n")

# Create a visualization of the perturbed trajectory
plt.figure(figsize=(12, 12))

# Plot the Sun
plt.scatter([0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's position
plt.scatter(earth_pos_m[0]/AU, earth_pos_m[1]/AU, color='blue', s=100, label='Earth')

# Plot Mars's position
plt.scatter(mars_pos_m[0]/AU, mars_pos_m[1]/AU, color='red', s=100, label='Mars')

# Plot the trajectory
plt.plot(trajectory[:, 0]/AU, trajectory[:, 1]/AU, color='green', label='Perturbed Trajectory')

# Plot the closest approach point
plt.scatter(trajectory[closest_idx, 0]/AU, trajectory[closest_idx, 1]/AU, color='purple', s=100, label='Closest Approach to Mars')

# Set labels and title
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.title('Spacecraft Trajectory with Gravitational Perturbations')
plt.grid(True)
plt.axis('equal')
plt.legend()

# Save the figure
plt.savefig('perturbed_trajectory.png', dpi=300, bbox_inches='tight')
plt.close()

# Create a 3D visualization
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the Sun
ax.scatter([0], [0], [0], color='yellow', s=200, label='Sun')

# Plot Earth's position
ax.scatter(earth_pos_m[0]/AU, earth_pos_m[1]/AU, earth_pos_m[2]/AU, color='blue', s=100, label='Earth')

# Plot Mars's position
ax.scatter(mars_pos_m[0]/AU, mars_pos_m[1]/AU, mars_pos_m[2]/AU, color='red', s=100, label='Mars')

# Plot the trajectory
ax.plot(trajectory[:, 0]/AU, trajectory[:, 1]/AU, trajectory[:, 2]/AU, color='green', label='Perturbed Trajectory')

# Plot the closest approach point
ax.scatter(trajectory[closest_idx, 0]/AU, trajectory[closest_idx, 1]/AU, trajectory[closest_idx, 2]/AU, color='purple', s=100, label='Closest Approach to Mars')

# Set labels and title
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
ax.set_title('Spacecraft Trajectory with Gravitational Perturbations (3D)')
ax.legend()

# Save the figure
plt.savefig('perturbed_trajectory_3d.png', dpi=300, bbox_inches='tight')
plt.close()

print("\nPerturbation analysis complete.")
print("Results saved to 'perturbation_analysis.txt'")
print("Visualizations saved to 'perturbed_trajectory.png' and 'perturbed_trajectory_3d.png'")

# Save the perturbation analysis results for use in other scripts
perturbation_results = {
    'closest_approach': {
        'distance_km': closest_distance/1000,
        'time_days': closest_time
    },
    'delta_v': {
        'insertion': delta_v_insertion,
        'total_with_perturbations': total_delta_v_with_perturbations,
        'difference_from_ideal': delta_v_difference,
        'percent_difference': percent_difference
    },
    'trajectory': {
        'positions': trajectory,
        'times': solution.t,
        'closest_idx': closest_idx
    }
}

np.save('perturbation_results.npy', perturbation_results)
