"""
DNA Sliding-Window Frequency Analyzer
Single-file Python application with a GUI (Tkinter) that allows the user to select a FASTA file,
compute relative frequencies of symbols in a sliding window (default size 30) across the sequence,
and plot the per-symbol percentage as lines on a chart (matplotlib) embedded in the GUI.

Features:
- Choose FASTA file using a file dialog
- Optionally change sliding window size
- Automatically detect sequence alphabet (e.g. A,C,G,T,N or others) and plot one line per symbol
- Show x-axis in window centers and y-axis as percent frequency (0-100)
- Save plotted chart as PNG

Requirements:
- Python 3.8+
- matplotlib
- numpy

Run:
python dna_sliding_window_analyzer.py

"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
from collections import Counter, OrderedDict
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


WINDOW_DEFAULT = 30


def parse_fasta(filepath):
    """Return concatenated sequence (uppercase) from the first entry in FASTA file."""
    seq_parts = []
    try:
        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # skip header lines
                    continue
                # keep letters (A-Z), convert to uppercase
                seq_parts.append(line.upper())
    except Exception as e:
        raise
    return "".join(seq_parts)


def sliding_window_frequencies(seq, window_size):
    """
    Compute relative frequencies (percent) of each symbol in seq for each sliding window.
    Returns (alphabet, x_positions, freq_matrix)
    - alphabet: list of symbols (sorted)
    - x_positions: list of window center indices (0-based)
    - freq_matrix: dict symbol -> list of percentages (same length as x_positions)
    """
    n = len(seq)
    if n == 0:
        return [], [], {}
    if window_size <= 0:
        raise ValueError("window_size must be > 0")

    # Determine alphabet from sequence (unique symbols in the order of appearance)
    seen = OrderedDict()
    for ch in seq:
        if ch not in seen:
            seen[ch] = True
    alphabet = list(seen.keys())

    # prepare container for frequencies
    freq = {sym: [] for sym in alphabet}
    x_pos = []

    # if sequence shorter than window -> compute a single window covering full sequence
    if n < window_size:
        win = seq
        counts = Counter(win)
        for sym in alphabet:
            pct = (counts.get(sym, 0) / len(win)) * 100
            freq[sym].append(pct)
        x_pos.append((0 + n - 1) / 2.0)
        return alphabet, x_pos, freq

    # sliding windows: start from 0 to n - window_size
    for start in range(0, n - window_size + 1):
        win = seq[start:start + window_size]
        counts = Counter(win)
        for sym in alphabet:
            pct = (counts.get(sym, 0) / window_size) * 100
            freq[sym].append(pct)
        center = start + (window_size - 1) / 2.0
        x_pos.append(center)

    return alphabet, x_pos, freq


class DNAAnalyzerApp:
    def __init__(self, master):
        self.master = master
        master.title("Sliding-window FASTA frequency analyzer")
        master.geometry("900x650")

        # Top frame: controls
        ctrl = ttk.Frame(master, padding=(10, 10))
        ctrl.pack(side=tk.TOP, fill=tk.X)

        ttk.Label(ctrl, text="FASTA file:").grid(row=0, column=0, sticky=tk.W)
        self.file_var = tk.StringVar()
        self.file_entry = ttk.Entry(ctrl, textvariable=self.file_var, width=60)
        self.file_entry.grid(row=0, column=1, padx=6)
        ttk.Button(ctrl, text="Browse...", command=self.browse_file).grid(row=0, column=2)

        # Window size control (slider)
        ttk.Label(ctrl, text="Window size:").grid(row=1, column=0, sticky=tk.W, pady=(8,0))
        self.window_var = tk.IntVar(value=WINDOW_DEFAULT)

        # Horizontal slider (Scale). Range can be adjusted as needed.
        # Using tk.Scale for integer stepping; showvalue=False because we display value in a separate label.
        self.win_scale = tk.Scale(ctrl,
                                  from_=1, to=1000,
                                  orient=tk.HORIZONTAL,
                                  variable=self.window_var,
                                  showvalue=False,
                                  length=300,
                                  resolution=1)
        self.win_scale.grid(row=1, column=1, sticky=tk.W, pady=(8,0))

        # Label that displays the current slider value (keeps visual parity with previous Spinbox)
        self.win_value_label = ttk.Label(ctrl, text=str(self.window_var.get()), width=6, anchor=tk.CENTER)
        self.win_value_label.grid(row=1, column=2, sticky=tk.W, padx=(6,0), pady=(8,0))

        # Update the small label when the slider moves
        def _on_scale_change(var, idx, mode):
            self.win_value_label.config(text=str(self.window_var.get()))
        self.window_var.trace_add('write', _on_scale_change)


        self.run_btn = ttk.Button(ctrl, text="Analyze & Plot", command=self.analyze_and_plot)
        self.run_btn.grid(row=1, column=2, pady=(8,0))

        ttk.Separator(master, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=8)

        # Main frame: left = plot, right = legend/table
        main = ttk.Frame(master, padding=(10, 0))
        main.pack(fill=tk.BOTH, expand=True)

        # Plot area
        plot_frame = ttk.Frame(main)
        plot_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.fig = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_xlabel("Sequence position (window center)")
        self.ax.set_ylabel("Percentage (%)")
        self.ax.grid(True, linestyle=':', linewidth=0.5)

        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        bottom_plot_controls = ttk.Frame(plot_frame)
        bottom_plot_controls.pack(fill=tk.X)
        ttk.Button(bottom_plot_controls, text="Save plot as PNG", command=self.save_plot).pack(side=tk.LEFT)

        # Right side: a Treeview showing the alphabet and last computed mean frequency
        right = ttk.Frame(main, width=220)
        right.pack(side=tk.RIGHT, fill=tk.Y)

        ttk.Label(right, text="Symbols (alphabet)").pack(anchor=tk.NW)
        columns = ("Symbol", "Mean %")
        self.tree = ttk.Treeview(right, columns=columns, show='headings', height=18)
        self.tree.heading('Symbol', text='Symbol')
        self.tree.heading('Mean %', text='Mean %')
        self.tree.column('Symbol', width=60, anchor=tk.CENTER)
        self.tree.column('Mean %', width=100, anchor=tk.E)
        self.tree.pack(fill=tk.Y, expand=True)

        # status bar
        self.status = tk.StringVar(value="Ready")
        ttk.Label(master, textvariable=self.status, relief=tk.SUNKEN, anchor=tk.W).pack(side=tk.BOTTOM, fill=tk.X)

        # store last results
        self.last_alphabet = []
        self.last_x = []
        self.last_freq = {}

    def browse_file(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.fna *.ffn *.faa"), ("All files","*.*")])
        if path:
            self.file_var.set(path)

    def analyze_and_plot(self):
        filepath = self.file_var.get().strip()
        if not filepath or not os.path.exists(filepath):
            messagebox.showerror("File missing", "Please choose a valid FASTA file.")
            return

        try:
            seq = parse_fasta(filepath)
        except Exception as e:
            messagebox.showerror("Read error", f"Failed to parse FASTA file:\n{e}")
            return

        if len(seq) == 0:
            messagebox.showwarning("Empty sequence", "The FASTA file contains no sequence data after headers.")
            return

        window_size = int(self.window_var.get())
        try:
            alphabet, x_pos, freq = sliding_window_frequencies(seq, window_size)
        except Exception as e:
            messagebox.showerror("Error", str(e))
            return

        if not alphabet:
            messagebox.showwarning("No symbols", "No sequence symbols were detected.")
            return

        self.last_alphabet = alphabet
        self.last_x = x_pos
        self.last_freq = freq

        # Clear plot and draw new lines
        self.ax.clear()
        self.ax.set_xlabel("Sequence position (window center)")
        self.ax.set_ylabel("Percentage (%)")
        self.ax.set_ylim(0, 100)
        self.ax.grid(True, linestyle=':', linewidth=0.5)

        for sym in alphabet:
            y = freq[sym]
            # plot as a line; use marker for clarity
            self.ax.plot(x_pos, y, label=sym, linewidth=1.5, marker='o', markersize=3)

        self.ax.legend(loc='upper right', title='Symbols')
        self.fig.tight_layout()
        self.canvas.draw()

        # update treeview with mean percentage per symbol
        for r in self.tree.get_children():
            self.tree.delete(r)
        for sym in alphabet:
            mean_pct = np.mean(freq[sym]) if len(freq[sym]) > 0 else 0.0
            self.tree.insert('', tk.END, values=(sym, f"{mean_pct:.2f}"))

        self.status.set(f"Analyzed {len(seq)} bases; window={window_size}; windows={len(x_pos)}")

    def save_plot(self):
        if not self.last_x:
            messagebox.showinfo("No plot", "Please run an analysis and plot first.")
            return
        fpath = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG image','*.png')])
        if fpath:
            try:
                self.fig.savefig(fpath)
                messagebox.showinfo("Saved", f"Plot saved to {fpath}")
            except Exception as e:
                messagebox.showerror("Save failed", str(e))


if __name__ == '__main__':
    root = tk.Tk()
    app = DNAAnalyzerApp(root)
    root.mainloop()
