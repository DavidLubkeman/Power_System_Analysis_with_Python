# Original program code for four-bus power flow video with generator_injection_calculation and print_results modified
# This is for a multi-source system with constant power loads
# 

import numpy as np
import csv

# Define global constants
FREQUENCY = 60  # Hz
TOLERANCE = 1e-4  # Convergence tolerance (0.01% = 0.0001, but using a smaller value for accuracy)
MAX_ITERATIONS = 20  # Maximum number of iterations

# Define the input_data function
def input_data():
    bus_data = []
    line_data = []
    
    with open("combined_data.txt", newline='') as csvfile:
              
        reader = csv.DictReader(csvfile)
        is_line_data = False
        
        for row in reader:
                     
            if row["section"] == "BUS_DATA":
                # Parse bus data
                bus_data.append({
                    "bus_name": row["bus_name"],
                    "bus_type": row["bus_type"],
                    "voltage_kV": float(row["voltage_kV"]),
                    "generation_power_MW": float(row["generation_power_MW"]) if row["generation_power_MW"] else None,
                    "load_real_power_MW": float(row["load_real_power_MW"]),
                    "load_reactive_power_MVAR": float(row["load_reactive_power_MVAR"]),
                })
            elif row["section"] == "LINE_DATA":
                # Parse line data
                line_data.append({
                    "from_bus_name": row["from_bus_name"],
                    "to_bus_name": row["to_bus_name"],
                    "r_Ohms": float(row["r_Ohms"]),
                    "x_Ohms": float(row["x_Ohms"]),
                    "capacitance_nF": float(row["capacitance_nF"]),
                })
            else:
                print("\n Invalid Data Entry")
    
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

            if bus["bus_type"] == "PV":
                # Fixed voltage magnitude for PV bus
                V_mag = (bus["voltage_kV"] * 1000) / np.sqrt(3)
                P_gen = bus["generation_power_MW"] * 1e6 / 3  # Convert MW to watts per-phase
                P_load = bus["load_real_power_MW"] * 1e6 / 3
                Q_load = bus["load_reactive_power_MVAR"] * 1e6 / 3
                P_net = P_gen - P_load

                # Compute reactive power
                sum_YV = 0
                for j in range(num_buses):
                    if i != j:
                        sum_YV += Ybus[i, j] * V[j]

                """
                Q_net = np.imag(np.conj(V[i]) * (Ybus[i, i] * V[i] + sum_YV))
                """
                # CHANGE Qnet=imag{V_i*cong(I_i)}
                Q_net = np.imag(V[i] * np.conj((Ybus[i, i] * V[i] + sum_YV)))

                # Update voltage phase
                V_new = (P_net - 1j * Q_net) / np.conj(V[i]) - sum_YV
                V_new /= Ybus[i, i]
                V_new = V_new / np.abs(V_new) * V_mag  # Maintain fixed magnitude

            else:  # PQ bus
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
            error = np.abs((V_new - V[i]) / V_new)
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
        
        # Get Z_line from line_data structure and also add  Y_cap term for line charging
        Z_line = complex(line["r_Ohms"] , line["x_Ohms"])
        Y_cap = complex(0.0, 2 * np.pi * FREQUENCY * line["capacitance_nF"] * 1e-9)
        
        # Calculate the sending-end and receiving-end complex powers

        # Use Z_line and Y_cap to compute current
        I_from_to = (V_from - V_to)/Z_line + (Y_cap/2)*V_from 
        
        # Multiply by 3 to convert per-phase to three-phase power
        S_from_to = 3*V_from * np.conj(I_from_to)  # Sending-end power

        I_to_from =  (V_to - V_from)/Z_line + (Y_cap/2)*V_to
        
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
        # Need too add flows, since they are both directed inwards
        total_loss += flow["sending_real_power_MW"] + flow["receiving_real_power_MW"]        
    return total_loss

# Function to calculate power injected by all generators
def generator_injection_calculation(V, Ybus, bus_data):
    generator_injections = []

    for i, bus in enumerate(bus_data):
        if bus["bus_type"] in ["SW", "PV"]:
            bus_type = bus["bus_type"]

            # Load at the generator bus (real and reactive)
            load_real_power = bus["load_real_power_MW"] * 1e6  # Convert MW to watts
            load_reactive_power = bus["load_reactive_power_MVAR"] * 1e6  # Convert MVAR to vars
            bus_load = complex(load_real_power, load_reactive_power)

            # Current injection at the bus
            I_injected = np.dot(Ybus[i, :], V)

            # Complex power injection (including bus load)
            S_injected = 3 * V[i] * np.conj(I_injected) + bus_load  # Multiply by 3 for three-phase

            generator_injections.append({
                "bus_name": bus["bus_name"],
                "bus_type": bus_type,
                "real_power_MW": S_injected.real / 1e6,  # Convert to MW
                "reactive_power_MVAR": S_injected.imag / 1e6  # Convert to MVAR
            })

    return generator_injections

# Updated print_results function to include generator injections
def print_results(bus_data, V, line_flow, total_loss, generator_injections):
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
    
    print("\nGenerator Injections:")
    for gen in generator_injections:
        print(f"{gen['bus_name']} ({gen['bus_type']}):")
        print(f"  Real Power: {gen['real_power_MW']:.4f} MW")
        print(f"  Reactive Power: {gen['reactive_power_MVAR']:.4f} MVAR")


# Updated main function to call the new functions
def main():
    # Call the input_data function to retrieve bus and line data
    bus_data, line_data = input_data()

    """
    # Print Bus Data
    print("Bus Data:")
    for bus in bus_data:
        print(bus)
    
    # Print Line Data
    print("\nLine Data:")
    for line in line_data:
        print(line)
    """
    
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


