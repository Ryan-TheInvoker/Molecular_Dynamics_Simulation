import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import os

##KINETIC ENERGY GRAPH##

# Initialize lists to store time and kinetic energy values
time = []
kinetic_energy = []

# Open the file and read the data
if os.path.exists("sim4_KE_list.txt"):
    with open("sim4_KE_list.txt", "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) == 2:  # Ensure there are two elements (time and kinetic energy)
                time.append(float(parts[0]))
                kinetic_energy.append(float(parts[1]))

    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(time, kinetic_energy, label='Kinetic Energy', color='r',  linestyle='-')
    plt.title('Kinetic Energy over Time')
    plt.xlabel('Time')
    plt.ylabel('Kinetic Energy')
    plt.grid(True)
    plt.legend()

    plt.savefig("sim4_Kinetic_Energy_Plot.png", format='png', dpi=300)
    plt.show()


##POTENTIAL ENERGY GRAPH##
time = []
potential_energy = []

# Open the file and read the data
if os.path.exists("sim4_PE_list.txt"):

    with open("sim4_PE_list.txt", "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) == 2:  # Ensure there are two elements (time and kinetic energy)
                time.append(float(parts[0]))
                potential_energy.append(float(parts[1]))

    # Create a plot
    plt.figure(figsize=(10, 6))
    plt.plot(time, potential_energy, label='Potential Energy', color='g',  linestyle='-')
    plt.title('Potential Energy over Time')
    plt.xlabel('Time')
    plt.ylabel('Potential Energy')
    plt.grid(True)
    plt.legend()

    plt.savefig("sim4_Potential_Energy_Plot.png", format='png', dpi=300)
    plt.show()





# CALCULATE TOTAL ENERGY ##

# Ensure kinetic and potential energy lists are the same length
if os.path.exists("sim4_PE_list.txt") and os.path.exists("sim2_KE_list.txt"):
    if len(kinetic_energy) != len(potential_energy):
        print("Error: Kinetic and Potential energy lists are of different lengths.")
        exit()

    # Calculate total energy as the sum of kinetic and potential energies
    total_energy = [ke + pe for ke, pe in zip(kinetic_energy, potential_energy)]

    ## PLOT TOTAL ENERGY ##

    plt.figure(figsize=(10, 6))
    plt.plot(time, total_energy, label='Total Energy', color='b', linestyle='-')
    plt.title('Total Energy over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Total Energy (J)')
    plt.grid(True)
    plt.legend()

    # Save the plot to a file and display it
    plt.savefig("sim4_Total_Energy_Plot.png", format='png', dpi=300)
    plt.show()



## Pressure vs Temperature Plot ##

# Constants
PARTICLE_NUMBER = 200
BOX_WIDTH = 40.0
BOX_HEIGHT = 40.0
KB = 1.38e-23  # Boltzmann constant in J/K
VOLUME = BOX_WIDTH * BOX_HEIGHT


def read_data(filename):
    temperatures = []
    pressures = []

    with open(filename, 'r') as file:
        for line in file:
            if line.strip():  # Ensure the line is not empty
                temp, pressure = map(float, line.split())
                temperatures.append(temp)
                pressures.append(pressure)
    
    return np.array(temperatures), np.array(pressures)

def plot_temperature_vs_pressure(temperatures, pressures):
    plt.figure(figsize=(10, 6))
    plt.plot(temperatures, pressures, 'o', label='P vs T')

        # Plot the ideal gas curve
    ideal_pressures = (PARTICLE_NUMBER * 1 * temperatures) / VOLUME #KB replaced with 1
    plt.plot(temperatures, ideal_pressures, 'r--', label='Ideal Gas P vs T')
    print("ideal pressure:",ideal_pressures)
    plt.xlabel('Temperature (T)')
    plt.ylabel('Pressure (P)')
    plt.title('Pressure vs Temperature')
    plt.legend()
    plt.grid(True)
    plt.savefig("sim4_TP_Plot.png", format='png', dpi=300)

    plt.show()

filename = 'sim4_TP_data.txt'
if os.path.exists('sim4_TP_data.txt'):
    temperatures, pressures = read_data(filename)
    plot_temperature_vs_pressure(temperatures, pressures)



##CALORIC CURVE##

def read_caloric_curve_data(filename):
    temperatures = []
    internal_energies = []

    with open(filename, 'r') as file:
        for line in file:
            if line.strip():  # Ensure the line is not empty
                temp, energy = map(float, line.split())
                temperatures.append(temp)
                internal_energies.append(energy)
    
    return np.array(temperatures), np.array(internal_energies)

def plot_caloric_curve(temperatures, internal_energies):
    plt.figure(figsize=(10, 6))
    plt.plot(temperatures, internal_energies, 'o-', label='Caloric Curve')
    
    plt.xlabel('Temperature (T)')
    plt.ylabel('Internal Energy (E)')
    plt.title('Caloric Curve')
    plt.legend()
    plt.grid(True)
    plt.savefig("caloric_curve.png", format='png', dpi=300)

    plt.show()

# Update the filename as per your C code's output file naming convention
filename = 'sim4_TE_list.txt'  # Replace 'your_filename' with the actual base filename used in your C code
if os.path.exists(filename):
    temperatures, internal_energies = read_caloric_curve_data(filename)
    plot_caloric_curve(temperatures, internal_energies)



##ANIMATION##

# Function to read data
def read_data(filename):
    with open(filename, 'r') as file:
        data = []
        for line in file:
            # Split the line into strings and filter out empty strings
            line_data = [x.strip() for x in line.split(',') if x.strip()]
            
            if line_data:  # Ensures we don't process empty lists
                data.append(np.array(line_data, dtype=float).reshape(-1, 2))
    return data


# Read the particle position data; here I have the file name attached manually
data = read_data('sim4_particle_positions.csv') 

# Create a figure and axis; note we MUST use same width and height as box
fig, ax = plt.subplots()
ax.set_xlim(0, 40)  
ax.set_ylim(0, 40)

# Initial plot objects, setting up as blue circles
lines, = ax.plot([], [], 'bo')

# Initialization function: plot the background of each frame
def init():
    lines.set_data([], [])
    return lines,

# Animation function: this is called sequentially
def animate(i):
    x = data[i][:, 0]  # X coordinates of all particles at time step i
    y = data[i][:, 1]  # Y coordinates of all particles at time step i
    lines.set_data(x, y)  # Update the data of the plot 
    return lines,

# Call the animator, pass the function and its parameters
ani = FuncAnimation(fig, animate, init_func=init, frames=len(data), interval=10, blit=True)



ani.save('particle_motion_2_ass2.mp4', writer='ffmpeg')

plt.show()




