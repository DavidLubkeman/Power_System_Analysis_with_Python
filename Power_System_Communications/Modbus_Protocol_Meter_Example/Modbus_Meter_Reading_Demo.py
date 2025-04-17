import time
from pymodbus.client import ModbusTcpClient
import tkinter as tk
from tkinter import ttk

# Global variables
SERVER_IP = '192.168.0.40'
CLIENT_IP = '192.168.0.11'
GATEWAY = '192.168.0.1'
PORT = 502
START_ADDRESS = 0
SAMPLING_RATES_IN_SEC = [0.1, 0.2, 0.5, 1, 2, 5]

REGISTER_POINTS = [
    {"name": "Voltage", "scaling_factor": 0.01, "units": "V", "value": 0},
    {"name": "Current", "scaling_factor": 0.001, "units": "A", "value": 0},   
    {"name": "Phase Angle Volts", "scaling_factor": 0.001, "units": "Degrees", "value": 0},
    {"name": "Phase Angle Amps", "scaling_factor": 0.01, "units": "Degrees", "value": 0},
    {"name": "Frequency", "scaling_factor": 0.01, "units": "Hz", "value": 0},
    {"name": "Watts", "scaling_factor": 0.1, "units": "W", "value": 0},
    {"name": "VARs", "scaling_factor": 0.1, "units": "VAR", "value": 0},    
    {"name": "VAs", "scaling_factor": 0.1, "units": "VA", "value": 0},
    {"name": "Power Factor", "scaling_factor": 0.001, "units": "PF", "value": 0},
    {"name": "Kilowatt-hours", "scaling_factor": 1.0, "units": "kWh", "value": 0},
    {"name": "Energy Reset", "scaling_factor": 1, "units": "", "value": 0},
]

NUM_POINTS = len(REGISTER_POINTS)  # Number of points to read

def read_modbus_registers():
    client = ModbusTcpClient(SERVER_IP, port=PORT)
    if not client.connect():
        print("Unable to connect to the Modbus server")
        return

    try:                    
        result = client.read_holding_registers(START_ADDRESS, count=NUM_POINTS)
        if not result.isError():
            timestamp = time.strftime("%H:%M:%S", time.localtime()) + f".{int(time.time() * 100) % 100:02d}"
            timestamp_label.config(text=f"Timestamp: {timestamp}")
            for i, value in enumerate(result.registers):
                # Convert 16-bit unsigned integer to signed integer
                if value > 32767:   
                    value -= 65536
                scaled_value = value * REGISTER_POINTS[i]["scaling_factor"]
                REGISTER_POINTS[i]["value"] = scaled_value
                if i < 10:
                    name_labels[i].config(text=f"{REGISTER_POINTS[i]['name']}", anchor='e', fg='dark blue')
                    value_labels[i].config(text=f"{scaled_value:.3f}", anchor='e')
                    unit_labels[i].config(text=f"{REGISTER_POINTS[i]['units']}", anchor='w', fg='brown')
                else:
                    reset_button.config(text=f"{REGISTER_POINTS[i]['name']}: {int(scaled_value)} {REGISTER_POINTS[i]['units']}")
        else:
            print("Error reading registers")
    except Exception as e:
        print(f"Exception: {e}")
    finally:
        client.close()

def send_reset_command():
    client = ModbusTcpClient(SERVER_IP, port=PORT)
    if not client.connect():
        print("Unable to connect to the Modbus server")
        return

    try:
        client.write_register(10, 1)
        print("Sent reset command to Modbus server")
    except Exception as e:
        print(f"Exception: {e}")
    finally:
        client.close()

def update_gui():
    read_modbus_registers()
    root.after(int(float(sampling_rate.get()) * 1000), update_gui)

root = tk.Tk()
root.title("Bitronics Meter")

# Create frames for layout
frame_upper_left = tk.Frame(root)
frame_upper_right = tk.Frame(root)
frame_lower_left = tk.Frame(root)
frame_lower_right = tk.Frame(root)
frame_upper_left.grid(row=0, column=0, padx=13, pady=13)
frame_upper_right.grid(row=0, column=1, padx=13, pady=13)
frame_lower_left.grid(row=1, column=0, padx=13, pady=13)
frame_lower_right.grid(row=1, column=1, padx=13, pady=13)

# Create labels for each register point
name_labels = []
value_labels = []
unit_labels = []

for i in range(5):
    name_label = tk.Label(frame_upper_left, text=f"{REGISTER_POINTS[i]['name']}", anchor='e', fg='dark blue', font=("Helvetica", 16))
    name_label.grid(row=i, column=0, sticky='e')
    name_labels.append(name_label)
    
    value_label = tk.Label(frame_upper_left, text=f"{REGISTER_POINTS[i]['value']:.3f}", anchor='e', relief='sunken', width=13, fg='red', bg='white', font=("Helvetica", 16))
    value_label.grid(row=i, column=1, sticky='e')
    value_labels.append(value_label)
    
    unit_label = tk.Label(frame_upper_left, text=f"{REGISTER_POINTS[i]['units']}", anchor='w', fg='brown', font=("Helvetica", 16))
    unit_label.grid(row=i, column=2, sticky='w')
    unit_labels.append(unit_label)

for i in range(5, 10):
    name_label = tk.Label(frame_upper_right, text=f"{REGISTER_POINTS[i]['name']}", anchor='e', fg='dark blue', font=("Helvetica", 16))
    name_label.grid(row=i-5, column=0, sticky='e')
    name_labels.append(name_label)
    
    value_label = tk.Label(frame_upper_right, text=f"{REGISTER_POINTS[i]['value']:.3f}", anchor='e', relief='sunken', width=13, fg='red', bg='white', font=("Helvetica", 16))
    value_label.grid(row=i-5, column=1, sticky='e')
    value_labels.append(value_label)
    
    unit_label = tk.Label(frame_upper_right, text=f"{REGISTER_POINTS[i]['units']}", anchor='w', fg='brown', font=("Helvetica", 16))
    unit_label.grid(row=i-5, column=2, sticky='w')
    unit_labels.append(unit_label)

reset_button = tk.Button(frame_lower_right, text=f"{REGISTER_POINTS[10]['name']}: {int(REGISTER_POINTS[10]['value'])} {REGISTER_POINTS[10]['units']}", command=send_reset_command, font=("Helvetica", 16))
reset_button.pack(anchor='w')

# Create timestamp label
timestamp_label = tk.Label(frame_lower_left, text="Timestamp: ", font=("Helvetica", 16))
timestamp_label.pack(anchor='w')

# Create sampling rate dropdown
sampling_rate_frame = tk.Frame(frame_lower_left)
sampling_rate_frame.pack(anchor='w')
sampling_rate_label = tk.Label(sampling_rate_frame, text="Sampling Rate:", font=("Helvetica", 16))
sampling_rate_label.pack(side='left')
sampling_rate = tk.StringVar(value=SAMPLING_RATES_IN_SEC[3])
sampling_rate_dropdown = ttk.Combobox(sampling_rate_frame, textvariable=sampling_rate, values=SAMPLING_RATES_IN_SEC, font=("Helvetica", 16), width=5)
sampling_rate_dropdown.pack(side='left')
seconds_label = tk.Label(sampling_rate_frame, text="Seconds", font=("Helvetica", 16))
seconds_label.pack(side='left')

# Start the GUI update loop
root.after(int(float(sampling_rate.get()) * 1000), update_gui)
root.mainloop()