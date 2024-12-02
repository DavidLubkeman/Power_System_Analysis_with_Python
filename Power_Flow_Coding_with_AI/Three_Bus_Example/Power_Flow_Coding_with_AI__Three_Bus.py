#!/usr/bin/env python
# coding: utf-8

# <center> Power System Analysis with Python 
# <br>
# Topic Area: Power System Coding with AI
# </center>

# <center>
# 
# # Power Flow Coding with AI - Three-Bus Example
# 
# <center>
#     
# ## Dr. David Lubkeman

# # Example 1: Request with Program Requirements
# Generate a Python program for solving a power system power flow for a three-bus system with a single source bus (SW) and all other buses being load buses (PQ).  
# 
# The three-bus system data is given in the Jupyter notebook markdown tables below:       
# # Three-Bus Power System Case
# # Bus Data
# |    | bus_name | bus_type |  voltage (kV) |  generation_power (MW) |  load_real_power (MW) |  load_reactive_power (MVAR) |
# |---:|:---------|:---------|--------------:|-----------------------:|----------------------:|----------------------------:|
# |  0 | Bus_1    | SW       |           230 |                     NA |                    50 |                       30.99 |
# |  1 | Bus_2    | PQ       |           230 |                      0 |                   170 |                      105.35 |
# |  2 | Bus_3    | PQ       |           230 |                      0 |                   200 |                      123.94 |
# 
# # Line Data
# |    | from_bus_name   | to_bus_name   |   r (Ohms)           |   x (Ohms)           |   capacitance (nF) |
# |---:|:----------------|:--------------|---------------------:|---------------------:|-------------------:|
# |  0 | Bus_1           | Bus_2         |              5.33232 |              26.6616 |            513.969 |
# |  1 | Bus_1           | Bus_3         |              3.93576 |              19.6788 |            388.611 |
# |  2 | Bus_2           | Bus_3         |              6.72888 |              33.6444 |            639.328 |
# 
# Program Requirements are as follow:
# 1. Utilize a gauss-seidel iterative solution technique,
# 2. Input bus quantities are in terms of line voltage and three-phase power, 
# 3. Algorithm should use a per-phase model in terms of phase voltage (line/sqrt(3)) and single-phase power (three-phase power/sqrt(3),
# 4. The accuracy for voltage should be within 0.01 percent of the actual solution,
# 5. Limit the number of iterations to 20,
# 6. Program power flow calculations should be in terms of volts, amperes and ohms,  
# 7. Do not use per unit system calculations.
# 8. Do not use matrix equation solvers.
#                                                                                                            
# Program Output is to include:
# 1. Bus line voltages with magnitudes in kilovolts and voltage angle in degrees, 
# 2. Three-phase real and reactive line power flow in MW and MVAR,
# 3. Total three-phase line real power loss in MW,
# 4. Three-phase real and reactive power delivered by the source voltage in MW and MVAR.
# 

# ## Example 1, Attempt 1

# In[3]:


# Example 1, Attempt 1
import numpy as np
import math

# Constants
MAX_ITER = 20
ACCURACY = 1e-4  # 0.01 percent accuracy
V_base_kV = 230  # Base voltage in kV

# Bus Data
bus_data = {
    'Bus_1': {'type': 'SW', 'voltage': 230 / math.sqrt(3), 'P_load': 50 / math.sqrt(3), 'Q_load': 30.99 / math.sqrt(3), 'generation': None},
    'Bus_2': {'type': 'PQ', 'voltage': 230 / math.sqrt(3), 'P_load': 170 / math.sqrt(3), 'Q_load': 105.35 / math.sqrt(3), 'generation': 0},
    'Bus_3': {'type': 'PQ', 'voltage': 230 / math.sqrt(3), 'P_load': 200 / math.sqrt(3), 'Q_load': 123.94 / math.sqrt(3), 'generation': 0}
}

# Line Data
line_data = [
    {'from': 'Bus_1', 'to': 'Bus_2', 'r': 5.33232, 'x': 26.6616, 'b': 2 * math.pi * 60 * 513.969e-9},
    {'from': 'Bus_1', 'to': 'Bus_3', 'r': 3.93576, 'x': 19.6788, 'b': 2 * math.pi * 60 * 388.611e-9},
    {'from': 'Bus_2', 'to': 'Bus_3', 'r': 6.72888, 'x': 33.6444, 'b': 2 * math.pi * 60 * 639.328e-9}
]

# Convert line parameters to impedance
for line in line_data:
    line['z'] = complex(line['r'], line['x'])
    line['y'] = 1 / line['z']  # Admittance

# Initialize voltages
V = {
    'Bus_1': complex(bus_data['Bus_1']['voltage'], 0),
    'Bus_2': complex(bus_data['Bus_2']['voltage'], 0),
    'Bus_3': complex(bus_data['Bus_3']['voltage'], 0)
}

# Gauss-Seidel Iteration
for iteration in range(MAX_ITER):
    V_prev = V.copy()
    
    # Update PQ bus voltages
    for bus_name in ['Bus_2', 'Bus_3']:
        P_load = bus_data[bus_name]['P_load']
        Q_load = bus_data[bus_name]['Q_load']
        
        # Initialize sum of voltages over impedances for all lines connected to this bus
        sum_V_over_Z = 0
        sum_Y = 0
        
        # Iterate over all lines connected to the bus
        for line in line_data:
            if line['from'] == bus_name:
                connected_bus = line['to']
            elif line['to'] == bus_name:
                connected_bus = line['from']
            else:
                continue
                
            sum_V_over_Z += V[connected_bus] / line['z']
            sum_Y += line['y']
        
        # Gauss-Seidel update formula
        V[bus_name] = (P_load - 1j * Q_load) / np.conj(V[bus_name]) + sum_V_over_Z / sum_Y

    # Check for convergence
    error = max(abs(V['Bus_2'] - V_prev['Bus_2']), abs(V['Bus_3'] - V_prev['Bus_3']))
    if error < ACCURACY:
        print(f"Convergence achieved in {iteration + 1} iterations.")
        break
else:
    print("Max iterations reached without convergence.")

# Calculate power flows and losses
total_loss_real = 0
total_loss_reactive = 0
for line in line_data:
    from_bus = line['from']
    to_bus = line['to']
    
    I_line = (V[from_bus] - V[to_bus]) / line['z']
    S_line = V[from_bus] * np.conj(I_line)
    
    P_line_real = S_line.real * math.sqrt(3)
    Q_line_reactive = S_line.imag * math.sqrt(3)
    
    total_loss_real += P_line_real
    total_loss_reactive += Q_line_reactive
    
    print(f"Power flow from {from_bus} to {to_bus}: P = {P_line_real:.2f} MW, Q = {Q_line_reactive:.2f} MVAR")

# Output results
for bus_name, voltage in V.items():
    voltage_magnitude_kV = abs(voltage) * math.sqrt(3)
    voltage_angle_deg = math.degrees(np.angle(voltage))
    print(f"Bus {bus_name}: Voltage magnitude = {voltage_magnitude_kV:.2f} kV, Voltage angle = {voltage_angle_deg:.2f} degrees")

# Output total power losses
print(f"Total three-phase real power loss = {total_loss_real:.2f} MW")
print(f"Total three-phase reactive power loss = {total_loss_reactive:.2f} MVAR")


# ## Example 1, Attempt 2

# In[4]:


import numpy as np

# Bus Data (Bus number, Bus type, Voltage (kV), Load real power (MW), Load reactive power (MVAR))
bus_data = {
    1: {'type': 'SW', 'V': 230.0, 'P_load': 50.0, 'Q_load': 30.99, 'angle': 0.0},  # Slack bus
    2: {'type': 'PQ', 'V': 230.0, 'P_load': 170.0, 'Q_load': 105.35, 'angle': 0.0},
    3: {'type': 'PQ', 'V': 230.0, 'P_load': 200.0, 'Q_load': 123.94, 'angle': 0.0},
}

# Line Data (from_bus, to_bus, r (Ohms), x (Ohms), capacitance (nF))
line_data = [
    (1, 2, 5.33232, 26.6616, 513.969),
    (1, 3, 3.93576, 19.6788, 388.611),
    (2, 3, 6.72888, 33.6444, 639.328),
]

# Gauss-Seidel Iterative Solver
def gauss_seidel(bus_data, line_data, max_iterations=20, tolerance=1e-4):
    # Initial bus voltages in terms of line-to-neutral
    V = {
        1: bus_data[1]['V'] / np.sqrt(3),
        2: bus_data[2]['V'] / np.sqrt(3),
        3: bus_data[3]['V'] / np.sqrt(3),
    }
    
    angles = {1: bus_data[1]['angle'], 2: bus_data[2]['angle'], 3: bus_data[3]['angle']}
    
    for iteration in range(max_iterations):
        V_old = V.copy()
        
        # Update voltage at each PQ bus
        for bus in bus_data:
            if bus_data[bus]['type'] == 'PQ':
                # Sum the admittance from the connected buses
                I = 0
                for line in line_data:
                    if line[0] == bus or line[1] == bus:
                        if line[0] == bus:
                            other_bus = line[1]
                        else:
                            other_bus = line[0]
                        
                        r = line[2]
                        x = line[3]
                        Z = r + 1j * x
                        Y = 1 / Z
                        
                        I += Y * (V_old[other_bus] - V_old[bus])
                
                # Load power for the bus
                P_load = bus_data[bus]['P_load'] / 3  # Single-phase load
                Q_load = bus_data[bus]['Q_load'] / 3  # Single-phase load
                
                S = (P_load - 1j * Q_load) / V[bus].conjugate()
                V[bus] = S / I
        
        # Check for convergence
        max_error = max(abs(V[bus] - V_old[bus]) for bus in bus_data if bus_data[bus]['type'] == 'PQ')
        if max_error < tolerance:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print("Did not converge within the maximum number of iterations.")
    
    return V, angles

# Calculate power flow for each line
def calculate_line_power_flow(V, angles, line_data):
    power_flow = []
    for line in line_data:
        from_bus = line[0]
        to_bus = line[1]
        r = line[2]
        x = line[3]
        
        Z = r + 1j * x
        I = (V[from_bus] - V[to_bus]) / Z
        
        S_from = V[from_bus] * I.conjugate() * 3
        S_to = V[to_bus] * I.conjugate() * 3
        
        power_flow.append({
            'from_bus': from_bus,
            'to_bus': to_bus,
            'S_from': S_from,
            'S_to': S_to,
        })
    
    return power_flow

# Run Gauss-Seidel Solver
V, angles = gauss_seidel(bus_data, line_data)

# Calculate Power Flows
power_flows = calculate_line_power_flow(V, angles, line_data)

# Output results
print("\nBus Voltages:")
for bus in V:
    magnitude = abs(V[bus]) * np.sqrt(3)
    angle = np.angle(V[bus], deg=True)
    print(f"Bus {bus}: Voltage = {magnitude:.2f} kV, Angle = {angle:.2f} degrees")

print("\nLine Power Flows:")
for flow in power_flows:
    print(f"From Bus {flow['from_bus']} to Bus {flow['to_bus']}:")
    print(f"  Real Power: {flow['S_from'].real:.2f} MW, Reactive Power: {flow['S_from'].imag:.2f} MVAR")

# Calculate Total Line Losses
total_loss_real = sum(flow['S_from'].real - flow['S_to'].real for flow in power_flows)
total_loss_reactive = sum(flow['S_from'].imag - flow['S_to'].imag for flow in power_flows)

print(f"\nTotal Line Losses: {total_loss_real:.2f} MW, {total_loss_reactive:.2f} MVAR")

# Calculate Power Delivered by the Source
S_source = sum(flow['S_from'] for flow in power_flows if flow['from_bus'] == 1)
print(f"\nPower Delivered by Source: {S_source.real:.2f} MW, {S_source.imag:.2f} MVAR")


# ## Example 2, Data Structure Options

# What are options for putting the following power system bus and line data tables into a Python data structure?  
# Ease of understanding the data structure is the primary consideration.  
# 
# The three-bus system data is given in the Jupyter notebook markdown tables below:       
# # Three-Bus Power System Case
# # Bus Data
# |    | bus_name | bus_type |  voltage (kV) |  generation_power (MW) |  load_real_power (MW) |  load_reactive_power (MVAR) |
# |---:|:---------|:---------|--------------:|-----------------------:|----------------------:|----------------------------:|
# |  0 | Bus_1    | SW       |           230 |                     NA |                    50 |                       30.99 |
# |  1 | Bus_2    | PQ       |           230 |                      0 |                   170 |                      105.35 |
# |  2 | Bus_3    | PQ       |           230 |                      0 |                   200 |                      123.94 |
# 
# # Line Data
# |    | from_bus_name   | to_bus_name   |   r (Ohms)           |   x (Ohms)           |   capacitance (nF) |
# |---:|:----------------|:--------------|---------------------:|---------------------:|-------------------:|
# |  0 | Bus_1           | Bus_2         |              5.33232 |              26.6616 |            513.969 |
# |  1 | Bus_1           | Bus_3         |              3.93576 |              19.6788 |            388.611 |
# |  2 | Bus_2           | Bus_3         |              6.72888 |              33.6444 |            639.328 |

# ## Example 2, Data Structure entry access

# In[5]:


bus_data = [
    {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
    {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
    {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
]
line_data = [
    {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
    {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
    {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
]
print(bus_data[0]["load_real_power_MW"])
print(line_data[2]["r_Ohms"])


# ## Example 2, Main(), input_data()

# In[ ]:


Write a Python program consisting of a Main function that calls function "input_data". The "input_data" function should return bus_data and line_data data structures, 
where these structures are defined by:
bus_data = [ 
{"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99}, 
{"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35}, 
{"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94} 
] 
line_data = [ 
{"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969}, 
{"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328}, 
{"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611} 
]


# In[6]:


# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Print Bus Data
    print("Bus Data:")
    for bus in bus_data:
        print(bus)
    
    # Print Line Data
    print("\nLine Data:")
    for line in line_data:
        print(line)

    # Example: Print the load real power for Bus_1
    for bus in bus_data:
        if bus["bus_name"] == "Bus_1":
            print(f'\nLoad Real Power for Bus_1: {bus["load_real_power_MW"]} MW')

    # Example: Print line data between Bus_1 and Bus_2
    for line in line_data:
        if line["from_bus_name"] == "Bus_1" and line["to_bus_name"] == "Bus_2":
            print(f'\nLine from Bus_1 to Bus_2 has resistance {line["r_Ohms"]} Ohms and reactance {line["x_Ohms"]} Ohms.')

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# ## Example 2, ybus_build()

# In[ ]:


Add a Python function to the following code called ybus_build which is to be called after bus_data and line_data are returned by input_data() in the main program.

ybus_build creates a bus admittance matrix called Ybus based on data contained in bus_data and line_data data structures.

The Ybus matrix should also utilize the line charging capacitance in line_data structure, which has units of nanoFarads.  The system frequency is set in the global parameter FREQUENCY which is initially set to 60.

After ybus_build is called, then also add code to print the Ybus matrix entries.

The code to be modified is as follows:

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[12]:


import numpy as np

# Define global system frequency
FREQUENCY = 60  # Hz

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Print the Ybus matrix
    print("Ybus Matrix:")
    for i in range(Ybus.shape[0]):
        for j in range(Ybus.shape[1]):
            print(f"Ybus[{i+1}][{j+1}] = {Ybus[i,j]:.4f}")
        print()  # New line after each row

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[ ]:


# Test ybus_build by creating entries manually
Z12=5.33232+26.6616j
Y12=377*513.969e-9j
Z23=6.7288+33.6444j
Y23=377*639.328e-9j
print(f"Y12= {-1/Z12:.4f}")
print(f"Y23= {-1/Z23:.4f}")
print(f"Y22' {1/Z12+1/Z23+Y12/2+Y23/2:.4f}")


# ## Example 2, gauss_seidel()

# In[ ]:


Add a Python function to the following code called gauss_seidel which is to be called after ybus_build()in the main program. 

gauss_seidel() returns a voltage matrix V based on data contained in bus_data, line_data and Ybus data structures.

The initial value for voltages should match the the voltage_kV field provided for bus_type=“SW” in the bus_data data structure.

Other requirements are as follows:
1. Utilize a gauss-seidel iterative solution technique,
2. Algorithm should use a per-phase model in terms of phase voltage (line/sqrt(3)) and single-phase power (three-phase power/sqrt(3),
3. The accuracy for voltage should be within 0.01 percent of the actual solution, with tolerance defined by global variable TOLERANCE,
4. Limit the number of iterations to 20, with iterations defined by MAX_ITERATIONS
5. Program power flow calculations should be in terms of volts, amperes and ohms,  
6. Do not use per unit system calculations.
7. Do not use matrix equation solvers.


The code to be modified is as follows:
import numpy as np

# Define global system frequency
FREQUENCY = 60  # Hz

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[8]:


# Ex. 2 Gauss_Seidel Function Original Code

import numpy as np

# Define global constants
FREQUENCY = 60  # Hz
TOLERANCE = 1e-6  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Gauss-Seidel power flow solver
def gauss_seidel(bus_data, Ybus):
    num_buses = len(bus_data)
    V = np.zeros(num_buses, dtype=complex)  # Voltage at each bus, complex numbers
    
    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V[i] = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
        else:
            V[i] = 1.0 + 0j  # Start with a flat voltage guess
    
    # Start Gauss-Seidel iterations
    for iteration in range(MAX_ITERATIONS):
        max_error = 0
        for i, bus in enumerate(bus_data):
            if bus["bus_type"] == "SW":
                continue  # Skip slack bus, its voltage is fixed
            
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / np.sqrt(3)  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / np.sqrt(3)  # Convert MVAR to VARs per-phase
            
            # Calculate the new voltage at bus i
            sum_YV = 0
            for j in range(num_buses):
                if i != j:
                    sum_YV += Ybus[i, j] * V[j]
            
            V_new = (P - 1j * Q) / np.conj(V[i]) - sum_YV
            V_new /= Ybus[i, i]
            
            # Calculate the error (voltage change)
            error = np.abs(V_new - V[i])
            max_error = max(max_error, error)
            
            # Update the voltage
            V[i] = V_new
        
        # Check convergence
        if max_error < TOLERANCE:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations.")
    
    return V

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Perform power flow calculation using Gauss-Seidel method
    V = gauss_seidel(bus_data, Ybus)
    
    # Print final bus voltages in kV per-phase
    print("Bus Voltages (in kV, per-phase):")
    for i, bus in enumerate(bus_data):
        print(f"{bus['bus_name']}: {np.abs(V[i]) * np.sqrt(3) / 1000:.4f} kV")

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[9]:


# Ex. 2 Gauss_Seidel Function Debugging

import numpy as np

# Define global constants
FREQUENCY = 60  # Hz

#TOLERANCE = 1e-6  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
# CHANGE 1, set tolerance back to original request
TOLERANCE = 1e-4  # Convergence tolerance (0.01% = 0.0001, using original specified value)

MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Gauss-Seidel power flow solver
def gauss_seidel(bus_data, Ybus):
    num_buses = len(bus_data)
    V = np.zeros(num_buses, dtype=complex)  # Voltage at each bus, complex numbers

    """    
    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V[i] = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
        else:
            V[i] = 1.0 + 0j  # Start with a flat voltage guess
    """
    # CHANGE 2, Initial voltage is swing bus voltage in units of Volts
    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V_Source = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
    for i, bus in enumerate(bus_data):
        V[i] = V_Source  # Start with a flat voltage guess based on source voltage     

    # Start Gauss-Seidel iterations
    for iteration in range(MAX_ITERATIONS):
        max_error = 0
        for i, bus in enumerate(bus_data):
            if bus["bus_type"] == "SW":
                continue  # Skip slack bus, its voltage is fixed
            
            """
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / np.sqrt(3)  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / np.sqrt(3)  # Convert MVAR to VARs per-phase
            """
            # CHANGE 3, single-phase is 1/3 of three-phase power
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / 3  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / 3  # Convert MVAR to VARs per-phase
                      
            # Calculate the new voltage at bus i
            sum_YV = 0
            for j in range(num_buses):
                if i != j:
                    sum_YV += Ybus[i, j] * V[j]
 
            """
            V_new = (P - 1j * Q) / np.conj(V[i]) - sum_YV
            """
            # CHANGE 4, sign on Sload incorrect
            V_new = -(P - 1j * Q) / np.conj(V[i]) - sum_YV   
            V_new /= Ybus[i, i]
      
            # Calculate the error (voltage change)
            """
            error = np.abs(V_new - V[i])            
            """
            # Change 5, compute a normalized error term
            error = np.abs((V_new - V[i])/V_new)
            max_error = max(max_error, error)
            
            # Update the voltage
            V[i] = V_new
        
        # Check convergence
        if max_error < TOLERANCE:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations.")
    
    return V

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Perform power flow calculation using Gauss-Seidel method
    V = gauss_seidel(bus_data, Ybus)
    
    """
    # Print final bus voltages in kV per-phase
    print("Bus Voltages (in kV, per-phase):")
    for i, bus in enumerate(bus_data):
        print(f"{bus['bus_name']}: {np.abs(V[i]) * np.sqrt(3) / 1000:.4f} kV")
    """
    # CHANGE 6, add degrees to voltage printout, adjust significant figures on magnitude
    print("Bus Voltages (in kV, per-phase) and Angles (in degrees):")
    for i, bus in enumerate(bus_data):
        voltage_magnitude_kV = np.abs(V[i]) * np.sqrt(3) / 1000  # Convert back to kV per-phase
        voltage_angle_deg = np.angle(V[i], deg=True)  # Get angle in degrees
        print(f"{bus['bus_name']}: {voltage_magnitude_kV:.1f} kV, {voltage_angle_deg:.2f} degrees")   

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# ## Example 2, output functions

# In[ ]:


Add the following Python functions  to the code provided below, which are to be called after gauss_seidel()in the main program: 
line_flow_caluculation()- returns line_flow data structure, which consists of sending and receiving end real and reactive power lineflows in MW and MVAR based on the voltage matrix V and line_data data structures
2. loss_calculation() – returns total_loss in MW, which is equal to the summation of real line powers stored in line_flow data structure
3. generator_injection_calculation()- returns generator_bus_power in MW and MVAR, which consists of real and reactive power injected by swing (SW) bus generator.
4. print_results() – prints the following:
 Bus voltages converted to magnitude in kV and angle in degrees
(b) Sending and receiving line flows in MW and MVAR
(c) Total line real power loss in MW
(d) Power injected by the swing bus (SW) generator in MW and MVAR

The code to be modified is provided below:
# Ex. 2 Debugging

import numpy as np

# Define global constants
FREQUENCY = 60  # Hz
TOLERANCE = 1e-4  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Gauss-Seidel power flow solver
def gauss_seidel(bus_data, Ybus):
    num_buses = len(bus_data)
    V = np.zeros(num_buses, dtype=complex)  # Voltage at each bus, complex numbers

    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V_Source = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
    for i, bus in enumerate(bus_data):
        V[i] = V_Source  # Start with a flat voltage guess based on source voltage     

    # Start Gauss-Seidel iterations
    for iteration in range(MAX_ITERATIONS):
        max_error = 0
        for i, bus in enumerate(bus_data):
            if bus["bus_type"] == "SW":
                continue  # Skip slack bus, its voltage is fixed
            
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / 3  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / 3  # Convert MVAR to VARs per-phase
                      
            # Calculate the new voltage at bus i
            sum_YV = 0
            for j in range(num_buses):
                if i != j:
                    sum_YV += Ybus[i, j] * V[j]
 
            V_new = -(P - 1j * Q) / np.conj(V[i]) - sum_YV   
            V_new /= Ybus[i, i]
        
            # Calculate the error (voltage change)
            error = np.abs((V_new - V[i])/V_new)
            max_error = max(max_error, error)
            
            # Update the voltage
            V[i] = V_new
        
        # Check convergence
        if max_error < TOLERANCE:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations.")
    
    return V

# Define the main function
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Perform power flow calculation using Gauss-Seidel method
    V = gauss_seidel(bus_data, Ybus)
    
    # Print final bus voltages in kV per-phase
    print("Bus Voltages (in kV, per-phase) and Angles (in degrees):")
    for i, bus in enumerate(bus_data):
        voltage_magnitude_kV = np.abs(V[i]) * np.sqrt(3) / 1000  # Convert back to kV per-phase
        voltage_angle_deg = np.angle(V[i], deg=True)  # Get angle in degrees
        print(f"{bus['bus_name']}: {voltage_magnitude_kV:.1f} kV, {voltage_angle_deg:.2f} degrees")   

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[10]:


# Ex. 2 - After Output Functions Added, Original Code

import numpy as np

# Define global constants
FREQUENCY = 60  # Hz
TOLERANCE = 1e-4  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Gauss-Seidel power flow solver
def gauss_seidel(bus_data, Ybus):
    num_buses = len(bus_data)
    V = np.zeros(num_buses, dtype=complex)  # Voltage at each bus, complex numbers

    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V_Source = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
    for i, bus in enumerate(bus_data):
        V[i] = V_Source  # Start with a flat voltage guess based on source voltage     

    # Start Gauss-Seidel iterations
    for iteration in range(MAX_ITERATIONS):
        max_error = 0
        for i, bus in enumerate(bus_data):
            if bus["bus_type"] == "SW":
                continue  # Skip slack bus, its voltage is fixed
            
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / 3  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / 3  # Convert MVAR to VARs per-phase
                      
            # Calculate the new voltage at bus i
            sum_YV = 0
            for j in range(num_buses):
                if i != j:
                    sum_YV += Ybus[i, j] * V[j]
 
            V_new = -(P - 1j * Q) / np.conj(V[i]) - sum_YV   
            V_new /= Ybus[i, i]
        
            # Calculate the error (voltage change)
            error = np.abs((V_new - V[i])/V_new)
            max_error = max(max_error, error)
            
            # Update the voltage
            V[i] = V_new
        
        # Check convergence
        if max_error < TOLERANCE:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations.")
    
    return V

# Function to calculate line flows
def line_flow_calculation(V, bus_data, line_data, Ybus):
    line_flow = []
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop over each line
    for line in line_data:
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Calculate the voltage difference across the line
        V_from = V[from_idx]
        V_to = V[to_idx]
        
        # Get the admittance of the line from the Ybus matrix
        Y_line = -Ybus[from_idx, to_idx]
        
        # Calculate the sending-end and receiving-end complex powers
        I_from_to = Y_line * (V_from - V_to)  # Current from bus 'from_idx' to 'to_idx'
        S_from_to = V_from * np.conj(I_from_to)  # Sending-end power

        I_to_from = Y_line * (V_to - V_from)  # Current from bus 'to_idx' to 'from_idx'
        S_to_from = V_to * np.conj(I_to_from)  # Receiving-end power
        
        # Store results in MW and MVAR
        line_flow.append({
            "from_bus_name": line["from_bus_name"],
            "to_bus_name": line["to_bus_name"],
            "sending_real_power_MW": S_from_to.real / 1e6,
            "sending_reactive_power_MVAR": S_from_to.imag / 1e6,
            "receiving_real_power_MW": S_to_from.real / 1e6,
            "receiving_reactive_power_MVAR": S_to_from.imag / 1e6
        })
    
    return line_flow

# Function to calculate total real power loss
def loss_calculation(line_flow):
    total_loss = 0
    for flow in line_flow:
        # Line losses are the difference between sending and receiving real power
        total_loss += flow["sending_real_power_MW"] - flow["receiving_real_power_MW"]
    return total_loss

# Function to calculate power injected by the swing bus generator
def generator_injection_calculation(V, Ybus, bus_data):
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            swing_bus_idx = i
            break

    # Calculate the injected power by the swing bus generator
    I_swing = np.dot(Ybus[swing_bus_idx, :], V)  # Current injection at the swing bus
    S_swing = V[swing_bus_idx] * np.conj(I_swing)  # Complex power injection

    generator_bus_power = {
        "real_power_MW": S_swing.real / 1e6,
        "reactive_power_MVAR": S_swing.imag / 1e6
    }
    
    return generator_bus_power

# Function to print the results
def print_results(bus_data, V, line_flow, total_loss, generator_bus_power):
    print("Bus Voltages (in kV, per-phase) and Angles (in degrees):")
    for i, bus in enumerate(bus_data):
        voltage_magnitude_kV = np.abs(V[i]) * np.sqrt(3) / 1000  # Convert back to kV per-phase
        voltage_angle_deg = np.angle(V[i], deg=True)  # Get angle in degrees
        print(f"{bus['bus_name']}: {voltage_magnitude_kV:.4f} kV, {voltage_angle_deg:.2f} degrees")
    
    print("\nLine Flows (MW and MVAR):")
    for flow in line_flow:
        print(f"{flow['from_bus_name']} -> {flow['to_bus_name']}:")
        print(f"  Sending End: {flow['sending_real_power_MW']:.4f} MW, {flow['sending_reactive_power_MVAR']:.4f} MVAR")
        print(f"  Receiving End: {flow['receiving_real_power_MW']:.4f} MW, {flow['receiving_reactive_power_MVAR']:.4f} MVAR")
    
    print(f"\nTotal Line Real Power Loss: {total_loss:.4f} MW")
    print(f"\nSwing Bus (SW) Generator Injection:")
    print(f"  Real Power: {generator_bus_power['real_power_MW']:.4f} MW")
    print(f"  Reactive Power: {generator_bus_power['reactive_power_MVAR']:.4f} MVAR")

# Updated main function to call the new functions
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Perform power flow calculation using Gauss-Seidel method
    V = gauss_seidel(bus_data, Ybus)
    
    # Calculate line flows
    line_flow = line_flow_calculation(V, bus_data, line_data, Ybus)
    
    # Calculate total line losses
    total_loss = loss_calculation(line_flow)
    
    # Calculate generator injection at the swing bus
    generator_bus_power = generator_injection_calculation(V, Ybus, bus_data)
    
    # Print the results
    print_results(bus_data, V, line_flow, total_loss, generator_bus_power)

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[11]:


# Ex. 2 Debugging Output Function Code

import numpy as np

# Define global constants
FREQUENCY = 60  # Hz
TOLERANCE = 1e-4  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    # Bus Data (list of dictionaries)
    bus_data = [
        {"bus_name": "Bus_1", "bus_type": "SW", "voltage_kV": 230, "generation_power_MW": None, "load_real_power_MW": 50, "load_reactive_power_MVAR": 30.99},
        {"bus_name": "Bus_2", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 170, "load_reactive_power_MVAR": 105.35},
        {"bus_name": "Bus_3", "bus_type": "PQ", "voltage_kV": 230, "generation_power_MW": 0, "load_real_power_MW": 200, "load_reactive_power_MVAR": 123.94}
    ]
    
    # Line Data (list of dictionaries)
    line_data = [
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_2", "r_Ohms": 5.33232, "x_Ohms": 26.6616, "capacitance_nF": 513.969},
        {"from_bus_name": "Bus_2", "to_bus_name": "Bus_3", "r_Ohms": 6.72888, "x_Ohms": 33.6444, "capacitance_nF": 639.328},
        {"from_bus_name": "Bus_1", "to_bus_name": "Bus_3", "r_Ohms": 3.93576, "x_Ohms": 19.6788, "capacitance_nF": 388.611}
    ]
    
    return bus_data, line_data

# Function to build the Ybus admittance matrix
def ybus_build(bus_data, line_data):
    # Number of buses
    num_buses = len(bus_data)
    
    # Initialize the Ybus matrix as a complex matrix (zero initially)
    Ybus = np.zeros((num_buses, num_buses), dtype=complex)
    
    # System frequency in radians per second
    omega = 2 * np.pi * FREQUENCY
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop through each line to calculate admittance and update Ybus
    for line in line_data:
        # Retrieve the bus indices from names
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Line impedance
        r = line["r_Ohms"]
        x = line["x_Ohms"]
        Z = complex(r, x)  # Impedance Z = R + jX
        Y = 1 / Z  # Admittance Y = 1 / Z
        
        # Line charging susceptance (B_shunt)
        capacitance = line["capacitance_nF"] * 1e-9  # Convert from nF to Farads
        B_shunt = omega * capacitance  # Shunt susceptance due to line charging
        
        # Update Ybus matrix for the line
        # Off-diagonal terms (mutual admittance between buses)
        Ybus[from_idx, to_idx] -= Y
        Ybus[to_idx, from_idx] -= Y
        
        # Diagonal terms (self-admittance at buses)
        Ybus[from_idx, from_idx] += Y + complex(0, B_shunt / 2)
        Ybus[to_idx, to_idx] += Y + complex(0, B_shunt / 2)
    
    return Ybus

# Gauss-Seidel power flow solver
def gauss_seidel(bus_data, Ybus):
    num_buses = len(bus_data)
    V = np.zeros(num_buses, dtype=complex)  # Voltage at each bus, complex numbers

    # Initialize voltage based on bus types
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            V_Source = (bus["voltage_kV"] * 1000) / np.sqrt(3)  # Convert kV to volts per-phase
    for i, bus in enumerate(bus_data):
        V[i] = V_Source  # Start with a flat voltage guess based on source voltage     

    # Start Gauss-Seidel iterations
    for iteration in range(MAX_ITERATIONS):
        max_error = 0
        for i, bus in enumerate(bus_data):
            if bus["bus_type"] == "SW":
                continue  # Skip slack bus, its voltage is fixed
            
            # Calculate power in terms of single-phase values
            P = bus["load_real_power_MW"] * 1e6 / 3  # Convert MW to watts per-phase
            Q = bus["load_reactive_power_MVAR"] * 1e6 / 3  # Convert MVAR to VARs per-phase
                      
            # Calculate the new voltage at bus i
            sum_YV = 0
            for j in range(num_buses):
                if i != j:
                    sum_YV += Ybus[i, j] * V[j]
 
            V_new = -(P - 1j * Q) / np.conj(V[i]) - sum_YV   
            V_new /= Ybus[i, i]
        
            # Calculate the error (voltage change)
            error = np.abs((V_new - V[i])/V_new)
            max_error = max(max_error, error)
            
            # Update the voltage
            V[i] = V_new
        
        # Check convergence
        if max_error < TOLERANCE:
            print(f"Converged in {iteration + 1} iterations.")
            break
    else:
        print(f"Did not converge after {MAX_ITERATIONS} iterations.")
    
    return V

# Function to calculate line flows
def line_flow_calculation(V, bus_data, line_data, Ybus):
    line_flow = []
    
    # Map bus names to indices
    bus_indices = {bus["bus_name"]: idx for idx, bus in enumerate(bus_data)}
    
    # Loop over each line
    for line in line_data:
        from_idx = bus_indices[line["from_bus_name"]]
        to_idx = bus_indices[line["to_bus_name"]]
        
        # Calculate the voltage difference across the line
        V_from = V[from_idx]
        V_to = V[to_idx]
        
        """
        # Get the admittance of the line from the Ybus matrix
        Y_line = -Ybus[from_idx, to_idx]
        """
        # CHANGE 1, get Z_line from line_data structure and also add  Y_cap term for line charging
        Z_line = complex(line["r_Ohms"] , line["x_Ohms"])
        Y_cap = complex(0.0, 2 * np.pi * FREQUENCY * line["capacitance_nF"] * 1e-9)
        
        # Calculate the sending-end and receiving-end complex powers

        """
        I_from_to = Y_line * (V_from - V_to)  # Current from bus 'from_idx' to 'to_idx'
        """
        # CHANGE 2, use Z_line and Y_cap to compute current
        I_from_to = (V_from - V_to)/Z_line + (Y_cap/2)*V_from 
        
        """
        S_from_to = V_from * np.conj(I_from_to)  # Sending-end power
        """
        # CHANGE 3, multiply by 3 to convert per-phase to three-phase power
        S_from_to = 3*V_from * np.conj(I_from_to)  # Sending-end power

        """
        I_to_from = Y_line * (V_to - V_from)  # Current from bus 'to_idx' to 'from_idx'
        """
        I_to_from =  (V_to - V_from)/Z_line + (Y_cap/2)*V_to
        
        """
        S_to_from = V_to * np.conj(I_to_from)  # Receiving-end power
        """
        S_to_from = 3*V_to * np.conj(I_to_from)  # Receiving-end power
        
        # Store results in MW and MVAR
        line_flow.append({
            "from_bus_name": line["from_bus_name"],
            "to_bus_name": line["to_bus_name"],
            "sending_real_power_MW": S_from_to.real / 1e6,
            "sending_reactive_power_MVAR": S_from_to.imag / 1e6,
            "receiving_real_power_MW": S_to_from.real / 1e6,
            "receiving_reactive_power_MVAR": S_to_from.imag / 1e6
        })
    
    return line_flow

# Function to calculate total real power loss
def loss_calculation(line_flow):
    total_loss = 0
    for flow in line_flow:
        # Line losses are the difference between sending and receiving real power
        """
        total_loss += flow["sending_real_power_MW"] - flow["receiving_real_power_MW"]
        """
        # CHANGE 4, need too add flows, since they are both directed inwards
        total_loss += flow["sending_real_power_MW"] + flow["receiving_real_power_MW"]        
    return total_loss

# Function to calculate power injected by the swing bus generator
def generator_injection_calculation(V, Ybus, bus_data):
    for i, bus in enumerate(bus_data):
        if bus["bus_type"] == "SW":
            swing_bus_idx = i

            # CHANGE 5, Create a term to represent swing bus load
            swing_bus_load = complex(bus["load_real_power_MW"]*1e6,bus["load_reactive_power_MVAR"]*1e6)
            break

 
    # Calculate the injected power by the swing bus generator
    I_swing = np.dot(Ybus[swing_bus_idx, :], V)  # Current injection at the swing bus
    # CHANGE 6, Need to multiply by 3 and add the load at the swing bus      
    S_swing = 3*V[swing_bus_idx] * np.conj(I_swing)+swing_bus_load  # Complex power injection

    generator_bus_power = {
        "real_power_MW": S_swing.real / 1e6,
        "reactive_power_MVAR": S_swing.imag / 1e6
    }
    
    return generator_bus_power

# Function to print the results
def print_results(bus_data, V, line_flow, total_loss, generator_bus_power):
    print("Bus Voltages (in kV, per-phase) and Angles (in degrees):")
    for i, bus in enumerate(bus_data):
        voltage_magnitude_kV = np.abs(V[i]) * np.sqrt(3) / 1000  # Convert back to kV per-phase
        voltage_angle_deg = np.angle(V[i], deg=True)  # Get angle in degrees
        print(f"{bus['bus_name']}: {voltage_magnitude_kV:.4f} kV, {voltage_angle_deg:.2f} degrees")
    
    print("\nLine Flows (MW and MVAR):")
    for flow in line_flow:
        print(f"{flow['from_bus_name']} -> {flow['to_bus_name']}:")
        print(f"  Sending End: {flow['sending_real_power_MW']:.4f} MW, {flow['sending_reactive_power_MVAR']:.4f} MVAR")
        print(f"  Receiving End: {flow['receiving_real_power_MW']:.4f} MW, {flow['receiving_reactive_power_MVAR']:.4f} MVAR")
    
    print(f"\nTotal Line Real Power Loss: {total_loss:.4f} MW")
    print(f"\nSwing Bus (SW) Generator Injection:")
    print(f"  Real Power: {generator_bus_power['real_power_MW']:.4f} MW")
    print(f"  Reactive Power: {generator_bus_power['reactive_power_MVAR']:.4f} MVAR")

# Updated main function to call the new functions
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()
    
    # Build the Ybus matrix
    Ybus = ybus_build(bus_data, line_data)
    
    # Perform power flow calculation using Gauss-Seidel method
    V = gauss_seidel(bus_data, Ybus)
    
    # Calculate line flows
    line_flow = line_flow_calculation(V, bus_data, line_data, Ybus)
    
    # Calculate total line losses
    total_loss = loss_calculation(line_flow)
    
    # Calculate generator injection at the swing bus
    generator_bus_power = generator_injection_calculation(V, Ybus, bus_data)
    
    # Print the results
    print_results(bus_data, V, line_flow, total_loss, generator_bus_power)

# Call the main function if the script is executed directly
if __name__ == "__main__":
    main()


# In[ ]:




