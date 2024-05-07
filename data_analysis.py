import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


##KINETIC ENERGY GRAPH##

# Initialize lists to store time and kinetic energy values
time = []
kinetic_energy = []

# Open the file and read the data
with open("sim1_KE_list.txt", "r") as file:
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

plt.savefig("Kinetic_Energy_Plot.png", format='png', dpi=300)
plt.show()


##POTENTIAL ENERGY GRAPH##
time = []
potential_energy = []

# Open the file and read the data
with open("sim1_PE_list.txt", "r") as file:
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

plt.savefig("Potential_Energy_Plot.png", format='png', dpi=300)
plt.show()





# CALCULATE TOTAL ENERGY ##

# Ensure kinetic and potential energy lists are the same length
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
plt.savefig("Total_Energy_Plot.png", format='png', dpi=300)
plt.show()







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
data = read_data('sim1_particle_positions.csv') 

# Create a figure and axis; note we MUST use same width and height as box
fig, ax = plt.subplots()
ax.set_xlim(0, 10)  
ax.set_ylim(0, 10)

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



ani.save('particle_motion_1_ass2.mp4', writer='ffmpeg')

plt.show()




