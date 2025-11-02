#!/usr/bin/env python3
"""
tm_sliding_tk.py

Tkinter app that:
 - Loads a FASTA file (first sequence found)
 - Scans the sequence with a sliding window of size 9 (fixed)
 - Computes Tm for each 9-mer by:
     * Wallace/simple: Tm = 4*(G+C) + 2*(A+T)
     * Salt-adjusted:  Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
 - Plots both Tm traces along the sequence (embedded matplotlib)
 - Allows saving the chart and setting [Na+] in mM
 - Hovering the plot shows the nearest point values in the status bar
"""

import math
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

WINDOW_SIZE = 9  # fixed sliding window size as requested


def parse_fasta_first_sequence(path):
    """Return the first sequence (concatenated, no headers) from a FASTA file."""
    seq_lines = []
    with open(path, "r") as f:
        found_header = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found_header:
                    # we already got header and entries; stop at next header
                    break
                found_header = True
                continue
            # accept sequence lines even if no header (robustness)
            seq_lines.append(line)
    if not seq_lines:
        raise ValueError("No sequence lines found in FASTA.")
    seq = "".join(seq_lines).upper()
    # Basic validation: allow A,T,G,C,N only
    for ch in seq:
        if ch not in "ATGCN":
            raise ValueError(f"Invalid base '{ch}' in sequence (only A,T,G,C,N allowed).")
    return seq


def tm_wallace(seq9):
    """Simple Wallace rule for a short oligo (works for any length algebraically)."""
    g = seq9.count("G")
    c = seq9.count("C")
    a = seq9.count("A")
    t = seq9.count("T")
    return 4 * (g + c) + 2 * (a + t)


def tm_salt_adjusted(seq9, na_mM):
    """Salt-adjusted formula. na_mM is in millimolar; converted to M inside."""
    length = len(seq9)
    if length == 0:
        return float("nan")
    na_M = float(na_mM) / 1000.0
    if na_M <= 0:
        raise ValueError("[Na+] must be > 0.")
    g = seq9.count("G")
    c = seq9.count("C")
    percent_gc = (g + c) / length * 100.0
    return 81.5 + 16.6 * math.log10(na_M) + 0.41 * percent_gc - 600.0 / length


class TmSlidingApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-window Tm scanner — nine positions of fate")
        self.geometry("950x650")

        self.seq = ""            # loaded sequence
        self.x_positions = None  # x values for plot (1-based center positions)
        self.tm_simple = None
        self.tm_salt = None

        self._build_ui()

    def _build_ui(self):
        # Top frame: Controls
        control = ttk.Frame(self, padding=(8, 6))
        control.pack(side="top", fill="x")

        ttk.Label(control, text="FASTA file:").grid(row=0, column=0, sticky="w")
        self.fasta_label = ttk.Label(control, text="(no file)", width=50)
        self.fasta_label.grid(row=0, column=1, sticky="w", padx=(4, 10))

        btn_load = ttk.Button(control, text="Load FASTA...", command=self.load_fasta)
        btn_load.grid(row=0, column=2, padx=4)

        ttk.Label(control, text="[Na+] (mM):").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.na_entry = ttk.Entry(control, width=10)
        self.na_entry.grid(row=1, column=1, sticky="w", padx=(4, 0), pady=(6, 0))
        self.na_entry.insert(0, "50.0")

        btn_plot = ttk.Button(control, text="Plot sliding Tm (window=9)", command=self.plot_tm)
        btn_plot.grid(row=1, column=2, padx=4, pady=(6, 0))

        btn_save = ttk.Button(control, text="Save chart...", command=self.save_plot)
        btn_save.grid(row=1, column=3, padx=4, pady=(6, 0))

        # Info frame
        info = ttk.Frame(self, padding=(8, 6))
        info.pack(side="top", fill="x")
        self.info_label = ttk.Label(info, text="Load a FASTA file to begin.", anchor="w")
        self.info_label.pack(fill="x")

        # Matplotlib figure area
        fig_frame = ttk.Frame(self)
        fig_frame.pack(fill="both", expand=True, padx=8, pady=8)

        self.fig = Figure(figsize=(8.5, 5.5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_title("Melting temperature along the sequence (sliding 9-mer)")
        self.ax.set_xlabel("Position (center of 9-mer, 1-based)")
        self.ax.set_ylabel("Tm (°C)")
        self.ax.grid(True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)

        # Status bar at bottom
        status = ttk.Frame(self, relief="sunken")
        status.pack(side="bottom", fill="x")
        self.status_var = tk.StringVar(value="Idle.")
        self.status_label = ttk.Label(status, textvariable=self.status_var, anchor="w")
        self.status_label.pack(fill="x")

        # Connect mouse move on canvas to show nearest point
        self.canvas.mpl_connect("motion_notify_event", self._on_motion_event)

    def load_fasta(self):
        path = filedialog.askopenfilename(title="Open FASTA file", filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("All files", "*.*")])
        if not path:
            return
        try:
            seq = parse_fasta_first_sequence(path)
        except Exception as e:
            messagebox.showerror("Error reading FASTA", str(e))
            return
        self.seq = seq
        self.fasta_label.config(text=path.split("/")[-1])
        self.info_label.config(text=f"Loaded sequence length = {len(seq)}. Window size fixed = {WINDOW_SIZE}.")
        self.status_var.set(f"Loaded sequence ({len(seq)} nt). Ready to plot.")
        # Clear previous plot
        self.ax.clear()
        self.ax.set_title("Melting temperature along the sequence (sliding 9-mer)")
        self.ax.set_xlabel("Position (center of 9-mer, 1-based)")
        self.ax.set_ylabel("Tm (°C)")
        self.ax.grid(True)
        self.canvas.draw_idle()

    def _compute_sliding_tm(self, na_mM):
        if not self.seq:
            raise RuntimeError("No sequence loaded.")
        seq = self.seq
        L = len(seq)
        if L < WINDOW_SIZE:
            raise RuntimeError(f"Sequence length {L} is shorter than window size {WINDOW_SIZE}.")
        # positions: use center position (1-based)
        centers = []
        simple_vals = []
        salt_vals = []
        for i in range(0, L - WINDOW_SIZE + 1):
            w = seq[i:i + WINDOW_SIZE]
            center_pos = i + (WINDOW_SIZE // 2) + 1  # 1-based
            centers.append(center_pos)
            simple_vals.append(tm_wallace(w))
            salt_vals.append(tm_salt_adjusted(w, na_mM))
        return np.array(centers), np.array(simple_vals), np.array(salt_vals)

    def plot_tm(self):
        if not self.seq:
            messagebox.showwarning("No sequence", "Please load a FASTA file first.")
            return
        # read [Na+]
        try:
            na_mM = float(self.na_entry.get())
        except ValueError:
            messagebox.showerror("Invalid [Na+]", "Please enter a valid numeric Na+ concentration (mM).")
            return
        if na_mM <= 0:
            messagebox.showerror("Invalid [Na+]", "Please supply a positive [Na+] in mM.")
            return

        try:
            x, svals, tvals = self._compute_sliding_tm(na_mM)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        self.x_positions = x
        self.tm_simple = svals
        self.tm_salt = tvals

        # clear and plot
        self.ax.clear()
        self.ax.set_title(f"Melting temperature along the sequence (sliding {WINDOW_SIZE}-mer)")
        self.ax.set_xlabel("Position (center of 9-mer, 1-based)")
        self.ax.set_ylabel("Tm (°C)")
        self.ax.grid(True)

        # Plot both traces. Keep visuals clean and legible.
        line1, = self.ax.plot(x, svals, label="Wallace (simple)", linewidth=1.6)
        line2, = self.ax.plot(x, tvals, label=f"Salt-adjusted ([Na+] = {na_mM} mM)", linewidth=1.6, linestyle='--')
        self.ax.legend(loc="best")
        # Provide vertical marker at sequence midpoint for visual interest
        mid = len(self.seq) // 2
        self.ax.axvline(mid + 1, alpha=0.12)

        # Rescale nicely
        ymin = min(np.min(svals), np.min(tvals))
        ymax = max(np.max(svals), np.max(tvals))
        yrange = ymax - ymin
        if yrange == 0:
            ymin -= 1
            ymax += 1
        else:
            ymin -= 0.08 * yrange
            ymax += 0.08 * yrange
        self.ax.set_ylim(ymin, ymax)

        self.canvas.draw_idle()
        self.status_var.set(f"Plotted {len(x)} windows (window size={WINDOW_SIZE}). Hover to inspect points.")

    def save_plot(self):
        if self.x_positions is None:
            messagebox.showwarning("Nothing to save", "No plot available. Generate a plot first.")
            return
        path = filedialog.asksaveasfilename(title="Save chart", defaultextension=".png",
                                            filetypes=[("PNG image", "*.png"), ("PDF", "*.pdf"), ("All files", "*.*")])
        if not path:
            return
        try:
            self.fig.savefig(path, dpi=300, bbox_inches="tight")
            self.status_var.set(f"Saved chart to {path}")
        except Exception as e:
            messagebox.showerror("Save error", str(e))

    def _on_motion_event(self, event):
        """Show nearest plotted point values when hovering over the axes."""
        if event.inaxes is not self.ax:
            return
        if self.x_positions is None:
            return
        # Find nearest x index
        try:
            xmouse = event.xdata
        except Exception:
            return
        # Use numpy to find nearest
        idx = np.abs(self.x_positions - xmouse).argmin()
        xpos = int(self.x_positions[idx])
        sval = float(self.tm_simple[idx])
        tval = float(self.tm_salt[idx])
        # Display concise info in status bar
        self.status_var.set(f"Pos {xpos}: Wallace Tm = {sval:.2f} °C | Salt-adjusted Tm = {tval:.2f} °C")


if __name__ == "__main__":
    app = TmSlidingApp()
    app.mainloop()
