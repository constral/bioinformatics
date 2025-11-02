#!/usr/bin/env python3
"""
fasta_analyzer.py

tkinter GUI to read a fasta file, analyze letters, and show results.
- reads a fasta (default fallback: ./fasta.fa)
- analyzes unique letters, counts and percentages, gc% etc.
- runs analysis in a background thread and updates a progress bar
"""

import tkinter as tk
from tkinter import filedialog, ttk
import os
import threading
import time
from collections import Counter
from typing import Tuple, Callable, Optional

DEFAULT_FASTA = "fasta.fa"   # fallback file if user cancels dialog

# ---------- core utilities ----------

def read_fasta(filepath: str) -> Tuple[str, str]:
    """Read a single-record fasta file and return (header, sequence).
    header and sequence are empty strings on error (caller handles message)."""
    header = ""
    seq_lines = []
    try:
        with open(filepath, "r") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    # take the first header we encounter
                    if not header:
                        header = line[1:].strip()
                    # continue reading sequence lines after header
                elif header:
                    seq_lines.append(line)
        sequence = "".join(seq_lines).upper()
        return header, sequence
    except Exception:
        return "", ""

def gc_percent(sequence: str) -> float:
    """Return percentage of G + C over total bases (0-100)."""
    if not sequence:
        return 0.0
    g = sequence.count("G")
    c = sequence.count("C")
    return 100.0 * (g + c) / len(sequence)

# ---------- analysis routine (runs in worker thread) ----------

def analyze_sequence(sequence: str, progress_cb: Callable[[int], None]) -> str:
    """Analyze the sequence and call progress_cb(percent) during the run.
    Returns a multiline textual result for display.
    """
    # immediate checks
    progress_cb(0)
    if not sequence:
        progress_cb(100)
        return "No sequence data found."

    # step 1: clean and collect letters
    progress_cb(5)
    cleaned = "".join(ch for ch in sequence if ch.isalpha()).upper()
    total = len(cleaned)
    if total == 0:
        progress_cb(100)
        return "Sequence found, but it contains no alphabetic characters."

    # update progress
    progress_cb(20)
    time.sleep(0.05)  # tiny pause so the UI can reflect progress

    # step 2: unique alphabet (sorted)
    unique_letters = sorted(set(cleaned))
    alphabet_str = "".join(unique_letters)
    progress_cb(40)
    time.sleep(0.05)

    # step 3: counts and percentages
    counts = Counter(cleaned)
    progress_cb(65)
    time.sleep(0.05)

    # step 4: gc%
    # note: if sequence contains non-ACGT letters they are ignored in gc calculation
    seq_for_gc = "".join(ch for ch in cleaned if ch in "ACGT")
    gc = gc_percent(seq_for_gc) if seq_for_gc else gc_percent(cleaned)
    progress_cb(80)
    time.sleep(0.05)

    # format results
    lines = []
    lines.append("--- FASTA Sequence Analysis ---")
    lines.append(f"Total letters analyzed: {total}")
    lines.append("")
    lines.append("[1] unique letters (alphabet)")
    lines.append(alphabet_str)
    lines.append("")
    lines.append("[2] letter frequency (percentage of total letters)")
    for letter in sorted(counts.keys()):
        pct = (counts[letter] * 100.0) / total
        lines.append(f"  {letter} : {pct:.2f}%")
    lines.append("")
    lines.append(f"[3] gc% (computed over A/C/G/T if present): {gc:.2f}%")

    progress_cb(100)
    return "\n".join(lines)

# ---------- GUI class ----------

class FastaAnalyzerGUI:
    def __init__(self, master: tk.Tk):
        self.master = master
        master.title("FASTA File Analyzer")
        master.geometry("700x540")

        # top button
        self.select_btn = tk.Button(master, text="Choose FASTA File (or use default)", command=self.on_select,
                                    font=("Arial", 11), bg="#2E86AB", fg="white", padx=8, pady=6)
        self.select_btn.pack(padx=12, pady=12, fill="x")

        # progress bar (hidden until used)
        self.progress = ttk.Progressbar(master, orient="horizontal", mode="determinate", length=640)
        self.progress.pack(padx=12, pady=(0,8))
        self.progress.pack_forget()

        # results label + text
        self.label = tk.Label(master, text="Analysis Results:", font=("Arial", 10, "bold"))
        self.label.pack(anchor="w", padx=12)

        self.results = tk.Text(master, wrap="word", height=24, font=("Courier", 10), bg="#fafafa", relief="sunken")
        self.results.pack(padx=12, pady=(4,12), fill="both", expand=True)
        self.results.insert("1.0", "press the button above to select a FASTA file (or click cancel to use fasta.fa).")

        # thread control
        self._worker_thread: Optional[threading.Thread] = None
        self._result_text: Optional[str] = None

    def on_select(self):
        """Handle user clicking the select button. Starts background analysis."""
        # ask file
        path = filedialog.askopenfilename(
            title="Select FASTA sequence file",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta *.fna *.fa"), ("All files", "*.*")]
        )
        # if user cancels, try default file in cwd
        if not path:
            if os.path.exists(DEFAULT_FASTA):
                path = DEFAULT_FASTA
            else:
                self._set_result_text("no file selected and default fasta.fa not found in current directory.")
                return

        # read fasta
        header, seq = read_fasta(path)
        if not seq:
            self._set_result_text(f"could not parse sequence from {os.path.basename(path)}. check file format.")
            return

        # display header and basic info
        self.results.delete("1.0", tk.END)
        self.results.insert(tk.END, f"selected file: {os.path.basename(path)}\n")
        self.results.insert(tk.END, f"header: {header}\n")
        self.results.insert(tk.END, f"raw buffer length: {len(seq)}\n\n")

        # show progress bar and disable the button while working
        self.progress.pack(pady=(0, 8))
        self.progress["value"] = 0
        self.select_btn.config(state="disabled")

        # start worker thread
        self._result_text = None
        self._worker_thread = threading.Thread(target=self._worker_run, args=(seq,), daemon=True)
        self._worker_thread.start()
        # poll for completion
        self.master.after(100, self._poll_worker)

    def _worker_run(self, seq: str):
        """Background worker: runs analysis and stores result in self._result_text."""
        # progress callback updates progress bar in the main thread via after
        def progress_cb(pct: int):
            # clamp and schedule in main thread
            pct = max(0, min(100, int(pct)))
            self.master.after(0, lambda: self.progress.configure(value=pct))

        result = analyze_sequence(seq, progress_cb)
        self._result_text = result

    def _poll_worker(self):
        """Poll for the worker completion and finalize UI when done."""
        if self._worker_thread and self._worker_thread.is_alive():
            # still working; reschedule check
            self.master.after(120, self._poll_worker)
            return

        # worker finished
        self.select_btn.config(state="normal")
        self.progress.pack_forget()
        if self._result_text is None:
            self._set_result_text("analysis finished but no result was produced.")
        else:
            # append the analysis text
            self.results.insert(tk.END, self._result_text)
            self.results.see("end")

    def _set_result_text(self, text: str):
        """Helper to replace contents of results text widget."""
        self.results.delete("1.0", tk.END)
        self.results.insert(tk.END, text)

# ---------- main ----------

def main():
    root = tk.Tk()
    app = FastaAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
