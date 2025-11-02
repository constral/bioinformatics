#!/usr/bin/env python3
"""
tm_sliding_threshold_tk.py

Tkinter app that:
 - Loads a FASTA file (first sequence found)
 - Scans the sequence with a sliding window of size 9 (fixed)
 - Computes Tm for each 9-mer by:
     * Wallace/simple: Tm = 4*(G+C) + 2*(A+T)
     * Salt-adjusted:  Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length
 - Plots both Tm traces along the sequence (embedded matplotlib)
 - Shows min / max values for both traces
 - Allows the user to set a threshold; displays regions above threshold on a second chart
 - The second chart shows two horizontal bars (one per signal). Bars are filled only
   where the signal is above the threshold. Empty elsewhere.
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
    seq_lines = []
    with open(path, "r") as f:
        found_header = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if found_header:
                    break
                found_header = True
                continue
            seq_lines.append(line)
    if not seq_lines:
        raise ValueError("No sequence lines found in FASTA.")
    seq = "".join(seq_lines).upper()
    for ch in seq:
        if ch not in "ATGCN":
            raise ValueError(f"Invalid base '{ch}' in sequence (only A,T,G,C,N allowed).")
    return seq


def tm_wallace(seq9):
    g = seq9.count("G")
    c = seq9.count("C")
    a = seq9.count("A")
    t = seq9.count("T")
    return 4 * (g + c) + 2 * (a + t)


def tm_salt_adjusted(seq9, na_mM):
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


def find_contiguous_regions(mask):
    """
    Given a boolean mask (list/ndarray) of length N, return list of (start_idx, end_idx)
    inclusive for runs where mask is True. Indices are window start indices (0-based).
    """
    regions = []
    N = len(mask)
    i = 0
    while i < N:
        if mask[i]:
            start = i
            j = i + 1
            while j < N and mask[j]:
                j += 1
            end = j - 1
            regions.append((start, end))
            i = j
        else:
            i += 1
    return regions


class TmSlidingThresholdApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Sliding-window Tm scanner — threshold regions (window=9)")
        self.geometry("1000x700")

        self.seq = ""
        self.window_starts = None  # list of window start indices (0-based)
        self.centers = None
        self.tm_simple = None
        self.tm_salt = None

        self._build_ui()

    def _build_ui(self):
        # Controls frame
        control = ttk.Frame(self, padding=(8, 6))
        control.pack(side="top", fill="x")

        ttk.Label(control, text="FASTA file:").grid(row=0, column=0, sticky="w")
        self.fasta_label = ttk.Label(control, text="(no file)", width=40)
        self.fasta_label.grid(row=0, column=1, sticky="w", padx=(4, 10))

        btn_load = ttk.Button(control, text="Load FASTA...", command=self.load_fasta)
        btn_load.grid(row=0, column=2, padx=4)

        ttk.Label(control, text="[Na+] (mM):").grid(row=1, column=0, sticky="w", pady=(6, 0))
        self.na_entry = ttk.Entry(control, width=12)
        self.na_entry.grid(row=1, column=1, sticky="w", padx=(4, 0), pady=(6, 0))
        self.na_entry.insert(0, "50.0")

        btn_plot = ttk.Button(control, text="Compute & Plot", command=self.plot_tm)
        btn_plot.grid(row=1, column=2, padx=4, pady=(6, 0))

        ttk.Label(control, text="Threshold (°C):").grid(row=2, column=0, sticky="w", pady=(6, 0))
        self.threshold_entry = ttk.Entry(control, width=12)
        self.threshold_entry.grid(row=2, column=1, sticky="w", padx=(4, 0), pady=(6, 0))
        self.threshold_entry.insert(0, "50.0")

        btn_apply_thresh = ttk.Button(control, text="Apply Threshold", command=self.apply_threshold)
        btn_apply_thresh.grid(row=2, column=2, padx=4, pady=(6, 0))

        btn_save = ttk.Button(control, text="Save chart...", command=self.save_plot)
        btn_save.grid(row=1, column=3, padx=6)

        # Info / min-max
        info = ttk.Frame(self, padding=(8, 6))
        info.pack(side="top", fill="x")
        self.info_label = ttk.Label(info, text="Load a FASTA file to begin.", anchor="w")
        self.info_label.pack(fill="x")
        summary = ttk.Frame(self)
        summary.pack(side="top", fill="x", padx=8)
        self.simple_minmax = ttk.Label(summary, text="Wallace: min=--  max=--")
        self.simple_minmax.pack(side="left", padx=(0,20))
        self.salt_minmax = ttk.Label(summary, text="Salt-adjusted: min=--  max=--")
        self.salt_minmax.pack(side="left")

        # Figure (two subplots stacked)
        fig_frame = ttk.Frame(self)
        fig_frame.pack(fill="both", expand=True, padx=8, pady=8)

        self.fig = Figure(figsize=(9.5, 6.5), dpi=100)
        # top: Tm traces; bottom: threshold bars
        self.ax_top = self.fig.add_subplot(211)
        self.ax_bottom = self.fig.add_subplot(212, sharex=self.ax_top)

        self.ax_top.set_title("Melting temperature along the sequence (sliding 9-mer)")
        self.ax_top.set_ylabel("Tm (°C)")
        self.ax_top.grid(True)

        self.ax_bottom.set_title("Regions above threshold (horizontal bars)")
        self.ax_bottom.set_xlabel("Position (1-based, sequence coordinates)")
        # We'll configure y ticks for bottom to label the two signals
        self.ax_bottom.set_yticks([15, 35])
        self.ax_bottom.set_yticklabels(["Wallace", "Salt-adjusted"])
        self.ax_bottom.set_ylim(0, 50)
        self.ax_bottom.set_xlim(0, 1)  # will be updated
        self.ax_bottom.grid(False)

        self.canvas = FigureCanvasTkAgg(self.fig, master=fig_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)

        # Status bar
        status = ttk.Frame(self, relief="sunken")
        status.pack(side="bottom", fill="x")
        self.status_var = tk.StringVar(value="Idle.")
        self.status_label = ttk.Label(status, textvariable=self.status_var, anchor="w")
        self.status_label.pack(fill="x")

        # Hover
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
        self.status_var.set(f"Loaded sequence ({len(seq)} nt). Ready to compute.")
        # clear plots
        self._clear_plots()

    def _clear_plots(self):
        self.ax_top.clear()
        self.ax_bottom.clear()
        self.ax_top.set_title("Melting temperature along the sequence (sliding 9-mer)")
        self.ax_top.set_ylabel("Tm (°C)")
        self.ax_top.grid(True)
        self.ax_bottom.set_title("Regions above threshold (horizontal bars)")
        self.ax_bottom.set_yticks([15, 35])
        self.ax_bottom.set_yticklabels(["Wallace", "Salt-adjusted"])
        self.ax_bottom.set_ylim(0, 50)
        self.canvas.draw_idle()
        self.x_positions = None
        self.window_starts = None
        self.tm_simple = None
        self.tm_salt = None
        self.simple_minmax.config(text="Wallace: min=--  max=--")
        self.salt_minmax.config(text="Salt-adjusted: min=--  max=--")

    def _compute_sliding_tm(self, na_mM):
        if not self.seq:
            raise RuntimeError("No sequence loaded.")
        seq = self.seq
        L = len(seq)
        if L < WINDOW_SIZE:
            raise RuntimeError(f"Sequence length {L} is shorter than window size {WINDOW_SIZE}.")
        starts = []   # window start indices (0-based)
        centers = []  # center positions (1-based)
        simple_vals = []
        salt_vals = []
        for i in range(0, L - WINDOW_SIZE + 1):
            w = seq[i:i + WINDOW_SIZE]
            starts.append(i)
            centers.append(i + (WINDOW_SIZE // 2) + 1)  # 1-based center
            simple_vals.append(tm_wallace(w))
            salt_vals.append(tm_salt_adjusted(w, na_mM))
        return np.array(starts), np.array(centers), np.array(simple_vals), np.array(salt_vals)

    def plot_tm(self):
        if not self.seq:
            messagebox.showwarning("No sequence", "Please load a FASTA file first.")
            return
        try:
            na_mM = float(self.na_entry.get())
        except ValueError:
            messagebox.showerror("Invalid [Na+]", "Please enter a valid numeric Na+ concentration (mM).")
            return
        if na_mM <= 0:
            messagebox.showerror("Invalid [Na+]", "Please supply a positive [Na+] in mM.")
            return

        try:
            starts, centers, svals, tvals = self._compute_sliding_tm(na_mM)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        self.window_starts = starts
        self.centers = centers
        self.tm_simple = svals
        self.tm_salt = tvals

        # Top plot: Tm traces
        self.ax_top.clear()
        self.ax_top.plot(centers, svals, label="Wallace (simple)", linewidth=1.6, zorder=2)
        self.ax_top.plot(centers, tvals, label=f"Salt-adjusted ([Na+] = {na_mM} mM)", linewidth=1.6, linestyle="--", zorder=2)
        self.ax_top.set_title(f"Melting temperature along the sequence (sliding {WINDOW_SIZE}-mer)")
        self.ax_top.set_ylabel("Tm (°C)")
        self.ax_top.set_xlim(1, len(self.seq))
        self.ax_top.grid(True)
        self.ax_top.legend(loc="best")

        # Min / max widgets
        smin, smax = float(np.min(svals)), float(np.max(svals))
        tmin, tmax = float(np.min(tvals)), float(np.max(tvals))
        self.simple_minmax.config(text=f"Wallace: min={smin:.2f}  max={smax:.2f}")
        self.salt_minmax.config(text=f"Salt-adjusted: min={tmin:.2f}  max={tmax:.2f}")

        # Bottom plot: empty for now; user must apply threshold to see bars
        self.ax_bottom.clear()
        self.ax_bottom.set_title("Regions above threshold (horizontal bars)")
        self.ax_bottom.set_yticks([15, 35])
        self.ax_bottom.set_yticklabels(["Wallace", "Salt-adjusted"])
        self.ax_bottom.set_ylim(0, 50)
        self.ax_bottom.set_xlim(1, len(self.seq))
        self.ax_bottom.set_xlabel("Position (1-based, sequence coordinates)")

        self.canvas.draw_idle()
        self.status_var.set(f"Computed {len(centers)} windows. Min/Max updated. Set threshold and click 'Apply Threshold'.")

    def apply_threshold(self):
        if self.tm_simple is None:
            messagebox.showwarning("No data", "Compute Tm traces first by clicking 'Compute & Plot'.")
            return
        try:
            threshold = float(self.threshold_entry.get())
        except ValueError:
            messagebox.showerror("Invalid threshold", "Please enter a valid numeric threshold (°C).")
            return

        # Build boolean masks based on window starts (window count = N)
        N = len(self.window_starts)
        mask_simple = np.array([v > threshold for v in self.tm_simple], dtype=bool)
        mask_salt   = np.array([v > threshold for v in self.tm_salt], dtype=bool)

        # convert masks into contiguous regions of window start indices
        regions_simple = find_contiguous_regions(mask_simple)
        regions_salt   = find_contiguous_regions(mask_salt)

        # Convert region start/end (window indices) -> sequence coordinates (1-based)
        # For window starting at s (0-based), window covers positions [s+1 .. s+WINDOW_SIZE]
        bars_simple = []
        for (s, e) in regions_simple:
            xstart = s + 1
            xend = e + WINDOW_SIZE
            width = xend - xstart + 1
            bars_simple.append((xstart, width))

        bars_salt = []
        for (s, e) in regions_salt:
            xstart = s + 1
            xend = e + WINDOW_SIZE
            width = xend - xstart + 1
            bars_salt.append((xstart, width))

        # Redraw bottom axes: two horizontal bars (y positions chosen to separate visually)
        self.ax_bottom.clear()
        self.ax_bottom.set_title(f"Regions above threshold = {threshold:.2f} °C")
        self.ax_bottom.set_xlim(1, len(self.seq))
        self.ax_bottom.set_ylim(0, 50)
        self.ax_bottom.set_yticks([15, 35])
        self.ax_bottom.set_yticklabels(["Wallace", "Salt-adjusted"])
        self.ax_bottom.set_xlabel("Position (1-based, sequence coordinates)")

        # For readability, draw background "empty" bars outline (light)
        self.ax_bottom.broken_barh([(1, len(self.seq))], (10, 8), facecolors='none', edgecolors='lightgray', linewidth=0.6)
        self.ax_bottom.broken_barh([(1, len(self.seq))], (30, 8), facecolors='none', edgecolors='lightgray', linewidth=0.6)

        # Now draw filled regions for each signal where above threshold
        if bars_simple:
            self.ax_bottom.broken_barh(bars_simple, (10, 8), facecolors=("tab:blue",), edgecolors=("tab:blue",), alpha=0.85)
        if bars_salt:
            self.ax_bottom.broken_barh(bars_salt, (30, 8), facecolors=("tab:orange",), edgecolors=("tab:orange",), alpha=0.85)

        # Annotate counts
        count_simple = len(bars_simple)
        count_salt = len(bars_salt)
        total_windows = len(self.window_starts)
        self.canvas.draw_idle()
        self.status_var.set(f"Threshold applied ({threshold:.2f} °C). Wallace regions: {len(regions_simple)}  Salt regions: {len(regions_salt)}  (of {total_windows} windows)")

    def save_plot(self):
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
        if event.inaxes is not self.ax_top:
            return
        if self.centers is None:
            return
        try:
            xmouse = event.xdata
        except Exception:
            return
        idx = np.abs(self.centers - xmouse).argmin()
        xpos = int(self.centers[idx])
        sval = float(self.tm_simple[idx])
        tval = float(self.tm_salt[idx])
        self.status_var.set(f"Pos {xpos}: Wallace Tm = {sval:.2f} °C | Salt-adjusted Tm = {tval:.2f} °C")


if __name__ == "__main__":
    app = TmSlidingThresholdApp()
    app.mainloop()
