#!/usr/bin/env python
# coding: utf-8

# <center> Power System Analysis with Python 
# <br>
# Topic Area: Power System Coding with AI
# </center>

# <center>
# 
# # Power Flow Coding with AI - Two-Bus Example
# 
# <center>
#     
# ## Dr. David Lubkeman

# # Example 1: Request with Minimal Requirements
# 
# Generate a Python program for analyzing a simple two-bus three-phase power system with buses Bus_1 and Bus_2.  Bus_1 connects to an ideal generation source while Bus_2 connects to a constant power load.  Bus_1 and Bus_2 are interconnected by a transmission line. 
# 
# The system parameters are as follows:
# 1. Bus_1 Generation source line voltage V_source: 230 kV
# 2. Bus_1 to Bus_2 Transmission line resistance R: 5.33232 Ohms
# 3. Bus_1 to Bus_2 Transmission line reactance X: 26.6616 Ohms
# 4. Bus_2 load real power consumption P_load: 170.0 MW
# 5. Bus_2 load reactive power consumption Q_load: 105.35 MVAR
# 
# Program Requirements:
# 
# 1.The accuracy for calculated voltage at Bus_2 should be within 0.01 percent of the actual solution.
# 
# Program Output:
# 1. Load bus line voltage in kilovolts, 
# 2. Total three-phase line real power loss in MW,
# 3. Three-phase real and reactive power delivered by the voltage source in MW and MVAR.
# 

# ## Example 1, Try #1

# In[37]:


# Example 1, Try #1

import numpy as np

# System parameters
V_source = 230e3  # Source line voltage in volts (230 kV)
R = 5.33232  # Transmission line resistance in ohms
X = 26.6616  # Transmission line reactance in ohms
P_load = 170e6  # Load real power consumption in watts (170 MW)
Q_load = 105.35e6  # Load reactive power consumption in vars (105.35 MVAR)

# Compute the impedance of the transmission line
Z_line = complex(R, X)  # Transmission line impedance

# Compute the load apparent power
S_load = complex(P_load, Q_load)  # Load apparent power in volt-amperes

# Calculate the current flowing through the transmission line
V_source_ph = V_source / np.sqrt(3)  # Phase voltage of the source
I_line = np.conj(S_load) / (V_source_ph)  # Line current

# Determine the voltage at Bus_2
V_drop = I_line * Z_line  # Voltage drop across the transmission line
V_bus2_ph = V_source_ph - V_drop  # Phase voltage at Bus_2
V_bus2 = np.abs(V_bus2_ph) * np.sqrt(3)  # Line voltage at Bus_2

# Calculate the power losses in the transmission line
P_loss = 3 * (I_line.real**2 * R + I_line.imag**2 * X)  # Three-phase real power loss

# Compute the power delivered by the source
S_source = 3 * V_source_ph * np.conj(I_line)  # Source apparent power
P_source = S_source.real  # Source real power
Q_source = S_source.imag  # Source reactive power

# Display the results
print(f"Load bus line voltage (Bus_2) in kilovolts: {V_bus2 / 1e3:.3f} kV")
print(f"Total three-phase line real power loss: {P_loss / 1e6:.3f} MW")
print(f"Three-phase real power delivered by the voltage source: {P_source / 1e6:.3f} MW")
print(f"Three-phase reactive power delivered by the voltage source: {Q_source / 1e6:.3f} MVAR")


# ## Example 1, Try #2

# In[39]:


# Example 1, Try #2

import numpy as np

# System parameters
V_source_kV = 230.0  # Source voltage in kV
R = 5.33232  # Line resistance in Ohms
X = 26.6616  # Line reactance in Ohms
P_load_MW = 170.0  # Load real power in MW
Q_load_MVAR = 105.35  # Load reactive power in MVAR

# Convert kV to V
V_source_V = V_source_kV * 1e3

# Calculate load impedance in per-unit
S_load = P_load_MW + 1j * Q_load_MVAR  # Load complex power in MVA
Z_line = R + 1j * X  # Line impedance in Ohms

# Initial guess for V2 (load bus voltage in per-unit)
V2_pu = V_source_V / np.sqrt(3)  # Assume load voltage is close to source voltage

# Iteratively solve for load bus voltage
tolerance = 1e-4  # Voltage accuracy tolerance
max_iterations = 1000  # Maximum number of iterations

for _ in range(max_iterations):
    I_line = np.conj(S_load / V2_pu)  # Line current in per-unit
    V2_new = V_source_V - I_line * Z_line  # Calculate new load voltage in per-unit
    
    if np.abs(V2_new - V2_pu) < tolerance:
        break
    V2_pu = V2_new

# Calculate final load bus voltage in kV
V2_kV = np.abs(V2_pu) / 1e3

# Calculate total three-phase line real power loss
I_line = np.conj(S_load / V2_pu)
P_loss_total = 3 * (I_line * np.conj(I_line)).real * R

# Calculate real and reactive power delivered by the source
S_source = 3 * V_source_V * np.conj(I_line) / 1e6  # Source power in MVA
P_source_MW = S_source.real  # Real power in MW
Q_source_MVAR = S_source.imag  # Reactive power in MVAR

# Output results
print(f"Load bus line voltage: {V2_kV:.2f} kV")
print(f"Total three-phase line real power loss: {P_loss_total:.2f} MW")
print(f"Real power delivered by the source: {P_source_MW:.2f} MW")
print(f"Reactive power delivered by the source: {Q_source_MVAR:.2f} MVAR")



# ## Example 1, Try #2 - Debugged

# In[41]:


# Example 1, Try #2 Debugging

import numpy as np

# System parameters
V_source_kV = 230.0  # Source voltage in kV
R = 5.33232  # Line resistance in Ohms
X = 26.6616  # Line reactance in Ohms
P_load_MW = 170.0  # Load real power in MW
Q_load_MVAR = 105.35  # Load reactive power in MVAR

# Convert kV to V
# *** Change_1 divide by sqrt(3) to get per-phase voltage
V_source_V = V_source_kV * 1e3/np.sqrt(3)


# *** Change_2, this comment not correct: Calculate load impedance in per-unit
# *** Change_3 S_load is per-phase complex power in Watts, VARs
S_load = ((P_load_MW + 1j * Q_load_MVAR)/3)*1e6  # Load complex power in MVA
Z_line = R + 1j * X  # Line impedance in Ohms

# Initial guess for V2 
# *** Change_4  set V2_pu to V2_old , modify comment to: Initial guess for V2
V2_old = V_source_V / np.sqrt(3)  # Assume load voltage is close to source voltage

# Iteratively solve for load bus voltage
tolerance = 1e-4  # Voltage accuracy tolerance
max_iterations = 1000  # Maximum number of iterations

for _ in range(max_iterations):
    # *** Change_5 V2_pu->V2_old
    I_line = np.conj(S_load / V2_old)  # Line current in per-unit
    V2_new = V_source_V - I_line * Z_line  # Calculate new load voltage in per-unit
 
    # *** Change_6 V2_pu -> V2_old, normalize error by V2_new, update V2_old
    if np.abs((V2_new - V2_old)/V2_new) < tolerance:
        break
    V2_old = V2_new

# Calculate final load bus voltage in kV
# *** Change_7 multiply by sqrt(3) for line value, V2_pu->V2_old  
V2_kV = np.sqrt(3)*np.abs(V2_old) / 1e3

# Calculate total three-phase line real power loss
# *** Change_8 V2_pu->V2_old
I_line = np.conj(S_load / V2_old)

# *** Change_9 Convert loss to MW by diving by 1e6
P_loss_total = 3 * (I_line * np.conj(I_line)).real * R/1e6

# Calculate real and reactive power delivered by the source
S_source = 3 * V_source_V * np.conj(I_line) / 1e6  # Source power in MVA
P_source_MW = S_source.real  # Real power in MW
Q_source_MVAR = S_source.imag  # Reactive power in MVAR

# Output results
print(f"Load bus line voltage: {V2_kV:.2f} kV")
print(f"Total three-phase line real power loss: {P_loss_total:.2f} MW")
print(f"Real power delivered by the source: {P_source_MW:.2f} MW")
print(f"Reactive power delivered by the source: {Q_source_MVAR:.2f} MVAR")



# # Example 2: Request with Detailed Requirements
# 
# Generate a Python program for analyzing a simple two-bus three-phase power system with buses Bus_1 and Bus_2.  Bus_1 connects to an ideal generation source while Bus_2 connects to a constant power load.  Bus_1 and Bus_2 are interconnected by a transmission line.
# 
# The system parameters are as follows:
# 1. Bus_1 Generation source line voltage V_source: 230 kV
# 2. Bus_1 to Bus_2 Transmission line resistance R: 5.33232 Ohms
# 3. Bus_1 to Bus_2 Transmission line reactance X: 26.6616 Ohms
# 4. Bus_2 load real power consumption P_load: 170.0 MW
# 5. Bus_2 load reactive power consumption Q_load: 105.35 MVAR
# 
# Program Requirements:
# 1. Utilize a Gauss-Seidel iterative solution technique,
# 2. Input voltage and power quantities are in terms of line voltage in kV and three-phase power in MW and MVAR, 
# 3. Algorithm should use a per-phase model in terms of phase voltage (line/sqrt(3)) and single-phase power (three-phase power/sqrt(3),
# 4. The accuracy for voltage should be within 0.01 percent of the actual solution,
# 5. Limit the number of iterations to 10,
# 6. Program power flow calculations should be in terms of volts, amperes and ohms,  
# 7. Do not use per unit system calculations.
# 8. Do not use matrix equation solvers.
# 9. Do not use a Ybus matrix.
# 
# Program Output:
# 1. Load bus line voltage in kilovolts, 
# 2. Total three-phase line real power loss in MW,
# 3. Three-phase real and reactive power delivered by the voltage source in MW and MVAR.
# 

# ## Example 2, Initial Code

# In[44]:


# Example 2, Initial Code

import cmath

# Constants
V_source_kV = 230  # Line voltage in kV
R = 5.33232  # Ohms
X = 26.6616  # Ohms
P_load_MW = 170.0  # MW
Q_load_MVAR = 105.35  # MVAR
max_iterations = 10
tolerance = 0.01 / 100  # 0.01 percent tolerance

# Convert three-phase quantities to per-phase quantities
V_source = V_source_kV * 1000 / cmath.sqrt(3)  # Phase voltage in volts
S_load = (P_load_MW * 1e6) / 3 + 1j * (Q_load_MVAR * 1e6) / 3  # Single-phase power in VA
Z_line = complex(R, X)  # Impedance of the transmission line

# Initial guesses
V2 = V_source
I2 = S_load.conjugate() / V2

# Gauss-Seidel iterations
for iteration in range(max_iterations):
    V2_new = V_source - I2 * Z_line
    I2_new = S_load.conjugate() / V2_new
    
    if abs((V2_new - V2) / V2_new) < tolerance:
        break
    
    V2 = V2_new
    I2 = I2_new

# Calculate final values
V2_kV = abs(V2) / 1000 * cmath.sqrt(3)  # Line voltage in kV
S_source = 3 * V_source * I2.conjugate()  # Three-phase power from the source in VA
P_loss = 3 * (I2 * I2.conjugate()).real * R  # Real power loss in W

# Convert to desired units
P_loss_MW = P_loss / 1e6  # Real power loss in MW
P_source_MW = S_source.real / 1e6  # Real power from source in MW
Q_source_MVAR = S_source.imag / 1e6  # Reactive power from source in MVAR

# Output results
print(f"Load bus line voltage: {V2_kV:.3f} kV")
print(f"Total three-phase line real power loss: {P_loss_MW:.3f} MW")
print(f"Three-phase real power delivered by the voltage source: {P_source_MW:.3f} MW")
print(f"Three-phase reactive power delivered by the voltage source: {Q_source_MVAR:.3f} MVAR")



# ## Example 2, Code Debugged

# In[46]:


# Example 2, Debugging

import cmath

# Constants
V_source_kV = 230  # Line voltage in kV
R = 5.33232  # Ohms
X = 26.6616  # Ohms
P_load_MW = 170.0  # MW
Q_load_MVAR = 105.35  # MVAR
max_iterations = 10
tolerance = 0.01 / 100  # 0.01 percent tolerance

# Convert three-phase quantities to per-phase quantities
V_source = V_source_kV * 1000 / cmath.sqrt(3)  # Phase voltage in volts
S_load = (P_load_MW * 1e6) / 3 + 1j * (Q_load_MVAR * 1e6) / 3  # Single-phase power in VA
Z_line = complex(R, X)  # Impedance of the transmission line

# Initial guesses
V2 = V_source
I2 = S_load.conjugate() / V2

# Gauss-Seidel iterations
for iteration in range(max_iterations):
    V2_new = V_source - I2 * Z_line
    
    # Change_1 Take complex conjugate of V2_new
    I2_new = S_load.conjugate() / V2_new.conjugate()
    
    if abs((V2_new - V2) / V2_new) < tolerance:
        break
    
    V2 = V2_new
    I2 = I2_new

# Calculate final values
V2_kV = abs(V2) / 1000 * cmath.sqrt(3)  # Line voltage in kV
S_source = 3 * V_source * I2.conjugate()  # Three-phase power from the source in VA
P_loss = 3 * (I2 * I2.conjugate()).real * R  # Real power loss in W

# Convert to desired units
P_loss_MW = P_loss / 1e6  # Real power loss in MW
P_source_MW = S_source.real / 1e6  # Real power from source in MW
Q_source_MVAR = S_source.imag / 1e6  # Reactive power from source in MVAR

# Output results
print(f"Load bus line voltage: {V2_kV:.3f} kV")
print(f"Total three-phase line real power loss: {P_loss_MW:.3f} MW")
print(f"Three-phase real power delivered by the voltage source: {P_source_MW:.3f} MW")
print(f"Three-phase reactive power delivered by the voltage source: {Q_source_MVAR:.3f} MVAR")



# # Example 3: Request with Detailed Flowchart
# 
# Generate a Python program for analyzing a simple two-bus three-phase power system with buses Bus_1 and Bus_2.  Bus_1 connects to an ideal generation source while Bus_2 connects to a constant power load.  Bus_1 and Bus_2 are interconnected by a transmission line.
# 
# The system parameters are as follows:
# 1. Bus_1 Generation source line voltage V_source: 230 kV
# 2. Bus_1 to Bus_2 Transmission line resistance R: 5.33232 Ohms
# 3. Bus_1 to Bus_2 Transmission line reactance X: 26.6616 Ohms
# 4. Bus_2 load real power consumption P_load: 170.0 MW
# 5. Bus_2 load reactive power consumption Q_load: 105.35 MVAR
# 
# Program Flowchart using the input data given above:
# 1. Input V_source, R, X, P_load, Q_load
# 2. Set
# 	Z = R+jX, 
# 	S_load=P_load+jQ_load, 
# 	V_1 = V_source/sqrt(3) 
# 	error_tolerance = 1e-4
# 3. Set intial voltage estimate V_2__est = V_1
# 4. I_load = conj(S_load)/3)/conj(V_2__est)
# 5. I_12 = I_load
# 6. V_2__update = V_1 – Z*I_12
# 7. error = abs((V_2__update - V_2__est)/V_2__update)
# 8. V_2__est = V_2__update
# 9. if error > error_tolerance go to Step 4
# 10. Output Load Voltage in kV =  sqrt(3) * abs(V_2__est)/1e3
# 11. Output Line Loss in MW = 3 * (abs(I_12)**2)*R/1e6
# 12. Output Source Power in MW, MVAr = 3 * V_1 * conj(I_12)/1e6

# ## Example 3, Initial Code

# In[49]:


# Example 3, Initial Code

import cmath

# Input parameters
V_source = 230e3  # Volts
R = 5.33232  # Ohms
X = 26.6616  # Ohms
P_load = 170e6  # Watts
Q_load = 105.35e6  # VARs

# Derived parameters
Z = complex(R, X)  # Impedance
S_load = complex(P_load, Q_load)  # Load power
V_1 = V_source / cmath.sqrt(3)  # Phase voltage at Bus_1
error_tolerance = 1e-4  # Error tolerance

# Initial voltage estimate
V_2_est = V_1

# Iterative process
while True:
    # Load current
    I_load = cmath.conj(S_load / 3) / cmath.conj(V_2_est)
    I_12 = I_load

    # Updated voltage at Bus_2
    V_2_update = V_1 - Z * I_12

    # Error calculation
    error = abs((V_2_update - V_2_est) / V_2_update)

    # Update voltage estimate
    V_2_est = V_2_update

    # Check for convergence
    if error <= error_tolerance:
        break

# Output results
load_voltage_kV = cmath.sqrt(3) * abs(V_2_est) / 1e3  # kV
line_loss_MW = 3 * (abs(I_12) ** 2) * R / 1e6  # MW
source_power_MVA = 3 * V_1 * cmath.conj(I_12) / 1e6  # MVA
source_power_MW = source_power_MVA.real  # MW
source_power_MVAR = source_power_MVA.imag  # MVAr

# Display results
print(f"Load Voltage: {load_voltage_kV:.4f} kV")
print(f"Line Loss: {line_loss_MW:.4f} MW")
print(f"Source Power: {source_power_MW:.4f} MW, {source_power_MVAR:.4f} MVAr")


# ## Example 3, Try 2 with function call error fixed

# In[52]:


# Example 3, Try #2 with function call error fixed

import cmath

# Input parameters
V_source = 230e3  # Volts
R = 5.33232  # Ohms
X = 26.6616  # Ohms
P_load = 170e6  # Watts
Q_load = 105.35e6  # VARs

# Derived parameters
Z = complex(R, X)  # Impedance
S_load = complex(P_load, Q_load)  # Load power
V_1 = V_source / cmath.sqrt(3)  # Phase voltage at Bus_1
error_tolerance = 1e-4  # Error tolerance

# Initial voltage estimate
V_2_est = V_1

# Iterative process
while True:
    # Load current
    I_load = S_load.conjugate() / 3 / V_2_est.conjugate()
    I_12 = I_load

    # Updated voltage at Bus_2
    V_2_update = V_1 - Z * I_12

    # Error calculation
    error = abs((V_2_update - V_2_est) / V_2_update)

    # Update voltage estimate
    V_2_est = V_2_update

    # Check for convergence
    if error <= error_tolerance:
        break

# Output results=
load_voltage_kV = cmath.sqrt(3) * abs(V_2_est) / 1e3  # kV
line_loss_MW = 3 * (abs(I_12) ** 2) * R / 1e6  # MW
source_power_MVA = 3 * V_1 * I_12.conjugate() / 1e6  # MVA
source_power_MW = source_power_MVA.real  # MW
source_power_MVAR = source_power_MVA.imag  # MVAr

# Display results
print(f"Load Voltage: {load_voltage_kV:.4f} kV")
print(f"Line Loss: {line_loss_MW:.4f} MW")
print(f"Source Power: {source_power_MW:.4f} MW, {source_power_MVAR:.4f} MVAr")


# ## Request AI to refactor (cleanup) code

# In[54]:


# Request AI to refactor (cleanup) code 

import cmath

# Input parameters
V_source = 230e3  # Line voltage at Bus_1 in Volts
R = 5.33232  # Transmission line resistance in Ohms
X = 26.6616  # Transmission line reactance in Ohms
P_load = 170e6  # Real power consumption at Bus_2 in Watts
Q_load = 105.35e6  # Reactive power consumption at Bus_2 in VARs

# Derived parameters
Z = complex(R, X)  # Impedance of the transmission line
S_load = complex(P_load, Q_load)  # Complex power of the load
V_1 = V_source / cmath.sqrt(3)  # Phase voltage at Bus_1
error_tolerance = 1e-4  # Error tolerance for iterative process

# Initial voltage estimate at Bus_2
V_2_est = V_1

# Iterative process to find the voltage at Bus_2
while True:
    # Calculate the load current
    I_load = S_load.conjugate() / 3 / V_2_est.conjugate()
    I_12 = I_load

    # Update the voltage at Bus_2
    V_2_update = V_1 - Z * I_12

    # Calculate the error
    error = abs((V_2_update - V_2_est) / V_2_update)

    # Update the voltage estimate for the next iteration
    V_2_est = V_2_update

    # Check for convergence
    if error <= error_tolerance:
        break

# Calculate the output results
load_voltage_kV = cmath.sqrt(3) * abs(V_2_est) / 1e3  # Load voltage in kV
line_loss_MW = 3 * (abs(I_12) ** 2) * R / 1e6  # Line loss in MW
source_power_MVA = 3 * V_1 * I_12.conjugate() / 1e6  # Source power in MVA
source_power_MW = source_power_MVA.real  # Source real power in MW
source_power_MVAR = source_power_MVA.imag  # Source reactive power in MVAr

# Display the results
print(f"Load Voltage: {load_voltage_kV:.4f} kV")
print(f"Line Loss: {line_loss_MW:.4f} MW")
print(f"Source Power: {source_power_MW:.4f} MW, {source_power_MVAR:.4f} MVAr")


# In[ ]:




