import tkinter as tk
from tkinter import ttk
from datetime import datetime
from pymodbus.client import ModbusTcpClient
import threading
import time

# Channel data structure
class Channel:
    def __init__(self, name, scale, unit):
        self.name = name
        self.scale = scale
        self.unit = unit

channel = [
    Channel("VA_MAG", 0.01, "V"),
    Channel("FREQ", 0.001, "Hz"),
    Channel("IA_MAG", 0.01, "mA"),
    Channel("IA_ANG", 0.01, "Degrees"),
    Channel("IARMS", 0.01, "mA"),
    Channel("PA", 0.01, "W"),
    Channel("QA", 0.01, "VAR"),
    Channel("SA", 0.01, "VA"),
    Channel("PFA", 0.001, ""),
    Channel("MWH3P", 0.00001, "kWh")
]

def combine_high_low(high, low):
    signed_high = high if high < 32768 else high - 65536
    combined = (signed_high << 16) | low
    return combined

class MeterGUI:
    def __init__(self, master):
        self.master = master
        master.title("Relay Example")
        master.geometry("400x600")
        master.configure(bg="#4472C4")
        master.resizable(False, False)

        # Display Panel
        self.display_frame = tk.Frame(master, bg="white", bd=2, relief="solid")
        self.display_frame.place(x=20, y=20, width=360, height=340)
        self.display_label = tk.Label(self.display_frame, text="Meter Points", font=("Arial", 24), bg="white")
        self.display_label.pack(anchor="n", pady=10)
        # Only one timestamp label, just below the title
        self.timestamp_label = tk.Label(self.display_frame, text="", font=("Arial", 14), bg="white", anchor="w")
        self.timestamp_label.pack(anchor="n", pady=(0, 10))
        self.data_text = tk.Text(self.display_frame, font=("Courier New", 14), bg="white", bd=0, height=12, width=38)
        self.data_text.pack(pady=(0,0))
        self.data_text.tag_configure("indent", lmargin1=36, lmargin2=36)

        # Sampling Rate Frame
        self.sampling_frame = tk.Frame(master, bg="white", bd=2, relief="solid")
        self.sampling_frame.place(x=20, y=380, width=220, height=120)
        tk.Label(self.sampling_frame, text="Sampling Rate:", font=("Arial", 16), bg="white").pack(anchor="nw", padx=10, pady=(10,0))
        self.sampling_var = tk.StringVar(value="1")
        self.sampling_menu = ttk.Combobox(self.sampling_frame, textvariable=self.sampling_var, state="readonly", font=("Arial", 14))
        self.sampling_menu['values'] = ("1", "2", "5", "10")
        self.sampling_menu.pack(anchor="nw", padx=10, pady=5)
        tk.Label(self.sampling_frame, text="Seconds", font=("Arial", 14), bg="white").pack(anchor="nw", padx=10)

        # Control Buttons
        self.close_btn = tk.Button(master, text="CLOSE", font=("Arial", 16), width=8, height=2, bg="white", relief="groove", command=self.close_action)
        self.close_btn.place(x=260, y=380, width=120, height=55)
        self.trip_btn = tk.Button(master, text="TRIP", font=("Arial", 16), width=8, height=2, bg="white", relief="groove", command=self.trip_action)
        self.trip_btn.place(x=260, y=445, width=120, height=55)

        # Coil Toggle Button
        self.coil_toggle_btn = tk.Button(master, text="COIL TOGGLE", font=("Arial", 16), width=16, height=2, bg="white", relief="groove", command=self.toggle_coil)
        self.coil_toggle_btn.place(x=20, y=530, width=360, height=55)
        self.update_coil_button_color()  # Set initial color

        # Modbus client
        self.client = ModbusTcpClient('192.168.0.30')
        self.running = True
        self.update_thread = threading.Thread(target=self.update_data, daemon=True)
        self.update_thread.start()

        # Handle closing
        master.protocol("WM_DELETE_WINDOW", self.on_close)

    def update_data(self):
        while self.running:
            try:
                if not self.client.connect():
                    self.display_error("Modbus connection failed.")
                    time.sleep(2)
                    continue

                # Read analog input registers
                rr = self.client.read_input_registers(address=540, count=20, slave=1)
                if not rr or not hasattr(rr, 'registers'):
                    self.display_error("Read failed.")
                    time.sleep(2)
                    continue

                registers = rr.registers
                values = []
                label_width = max(len(ch.name) for ch in channel) + 2  # +2 for spacing
                for i in range(0, 20, 2):
                    high = registers[i]
                    low = registers[i+1]
                    raw = combine_high_low(high, low)
                    ch = channel[i//2]
                    scaled = raw * ch.scale
                    values.append(f"{ch.name:<{label_width}} {scaled:>10.4f} {ch.unit}")

                # Read discrete input register 8 (bit 8, zero-based)
                dr = self.client.read_discrete_inputs(address=8, count=1, slave=1)
                if dr and hasattr(dr, 'bits'):
                    breaker_contact = dr.bits[0]
                else:
                    breaker_contact = None

                # Update button colors based on 52A_breaker_contact
                if breaker_contact is True:
                    self.close_btn.config(bg="red")
                    self.trip_btn.config(bg="white")
                elif breaker_contact is False:
                    self.close_btn.config(bg="white")
                    self.trip_btn.config(bg="red")
                else:
                    # If read failed, set both to default
                    self.close_btn.config(bg="white")
                    self.trip_btn.config(bg="white")

                # Update display
                now = datetime.now().strftime("%H:%M:%S")  # Only show time
                self.timestamp_label.config(text=f"Timestamp: {now}")
                self.data_text.delete(1.0, tk.END)
                self.data_text.insert(tk.END, "\n".join(values), "indent")

            except Exception as e:
                self.display_error(str(e))

            # Wait for next interval
            try:
                interval = int(self.sampling_var.get())
            except Exception:
                interval = 1
            time.sleep(interval)
            self.update_coil_button_color()

    def display_error(self, msg):
        self.data_text.delete(1.0, tk.END)
        self.data_text.insert(tk.END, f"Error: {msg}")

    def close_action(self):
        # Implement relay CLOSE action here
        tk.messagebox.showinfo("Relay", "CLOSE command sent.")

    def trip_action(self):
        # Implement relay TRIP action here
        tk.messagebox.showinfo("Relay", "TRIP command sent.")

    def update_coil_button_color(self):
        """Read coil 0 and update the button color."""
        try:
            rr = self.client.read_coils(address=0, count=1, slave=1)
            if rr and hasattr(rr, 'bits'):
                coil_status = rr.bits[0]
                if coil_status:
                    self.coil_toggle_btn.config(bg="green")
                else:
                    self.coil_toggle_btn.config(bg="white")
            else:
                self.coil_toggle_btn.config(bg="white")
        except Exception:
            self.coil_toggle_btn.config(bg="white")

    def toggle_coil(self):
        """Toggle coil 0: read, then write the opposite value."""
        try:
            rr = self.client.read_coils(address=0, count=1, slave=1)
            if rr and hasattr(rr, 'bits'):
                current = rr.bits[0]
                new_value = not current
                self.client.write_coil(address=0, value=new_value, slave=1)
                # Update button color after write
                self.coil_toggle_btn.config(bg="green" if new_value else "white")
            else:
                tk.messagebox.showerror("Coil Toggle", "Failed to read coil status.")
        except Exception as e:
            tk.messagebox.showerror("Coil Toggle", f"Error: {e}")

    def on_close(self):
        self.running = False
        self.client.close()
        self.master.destroy()

if __name__ == "__main__":
    root = tk.Tk()
    app = MeterGUI(root)
    root.mainloop()