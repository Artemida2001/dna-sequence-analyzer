#!/usr/bin/env python3
"""
DNA Sequence Analyzer - A GUI tool for DNA analysis.
Features: FASTA/NCBI loading, motif search, ORF detection, alignment, heatmap, CSV/PDF export.
"""

import re
import tempfile
import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, scrolledtext, ttk, Toplevel
from typing import Dict, List, Optional, Tuple

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from fpdf import FPDF
from fpdf.enums import XPos, YPos
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


__author__ = "Artemida Chadzinikolau"
__version__ = "1.1"

# App configuration
Entrez.email = "your_email@example.com"  # Required by NCBI
DEFAULT_SEGMENT_SIZE = 100
MIN_ORF_AA_LENGTH = 10
LONG_SEQUENCE_THRESHOLD = 10_000
FASTA_LINE_WIDTH = 70

# GUI configuration
STYLES = {
    "TLabel": {"font": ("Helvetica", 10)},
    "TButton": {"font": ("Helvetica", 10, "bold")},
    "TLabelFrame.Label": {"font": ("Helvetica", 12, "bold")},
    "Treeview.Heading": {"font": ("Helvetica", 10, "bold")},
    "Export.TButton": {"background": "#DCEDC8", "foreground": "#33691E"},
    "PDF.TButton": {"background": "#BBDEFB", "foreground": "#0D47A1"},
    "Heatmap.TButton": {"background": "#FFC107"},
}
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8,             # base text size
    "axes.titlesize": 8,        # chart title
    "axes.labelsize": 8,        # X/Y labels
    "xtick.labelsize": 8,       # tick labels
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 8,
    "figure.dpi": 100,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",    # prevents cutting off labels in PDF
    "savefig.pad_inches": 0.2,
})

class App:
    """Main application class for DNA sequence analysis."""

    def __init__(self, master: tk.Tk) -> None:
        self.master = master
        master.title(f"🧬 DNA Sequence Analyzer v{__version__} - {__author__}")
        master.geometry("1200x700")

        # State
        self.sequences: Dict[str, Seq] = {}
        self.current_seq_id: Optional[str] = None
        self.analysis_results: Dict[str, Dict] = {}
        self.last_detected_orfs: List[Dict] = []
        self.segment_size: int = DEFAULT_SEGMENT_SIZE

        # UI Setup
        self._configure_styles()
        self._build_layout()
        self._setup_tabs()

    # === UI Setup ===
    def _configure_styles(self) -> None:
        style = ttk.Style()
        style.theme_use("clam")
        for name, config in STYLES.items():
            style.configure(name, **config)

    def _build_layout(self) -> None:
        self.paned_window = ttk.PanedWindow(self.master, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # Left Panel
        self.left_panel = ttk.Frame(self.paned_window, width=480)
        self.paned_window.add(self.left_panel, weight=5)
        self._setup_left_panel()

        # Right Panel
        self.right_panel = ttk.Frame(self.paned_window)
        self.paned_window.add(self.right_panel, weight=8)

    def _setup_left_panel(self) -> None:
        # Load Frame
        load_frame = ttk.LabelFrame(self.left_panel, text="📥 Load Sequences")
        load_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(load_frame, text="From FASTA/TXT file", command=self.load_from_file).pack(
            fill=tk.X, padx=5, pady=5
        )

        ncbi_frame = ttk.Frame(load_frame)
        ncbi_frame.pack(fill=tk.X, padx=5, pady=5)
        self.ncbi_entry = ttk.Entry(ncbi_frame, width=20)
        self.ncbi_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))
        self.ncbi_entry.insert(0, "NR_102783, MT012657")
        ttk.Button(ncbi_frame, text="Fetch from NCBI", command=self.load_from_ncbi).pack(side=tk.LEFT)

        self.status_label = ttk.Label(load_frame, text="Status: Ready", foreground="blue")
        self.status_label.pack(fill=tk.X, padx=5, pady=5)

        # Stats Table
        stats_frame = ttk.LabelFrame(self.left_panel, text="Loaded Sequences")
        stats_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        columns = ("ID", "Length (bp)", "GC Content (%)")
        self.stats_tree = ttk.Treeview(stats_frame, columns=columns, show="headings", height=12)
        for col in columns:
            self.stats_tree.heading(col, text=col)
        self.stats_tree.column("ID", width=80, anchor=tk.W)
        self.stats_tree.column("Length (bp)", width=80, anchor=tk.E)
        self.stats_tree.column("GC Content (%)", width=90, anchor=tk.E)

        scrollbar = ttk.Scrollbar(stats_frame, orient=tk.VERTICAL, command=self.stats_tree.yview)
        self.stats_tree.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.stats_tree.pack(fill=tk.BOTH, expand=True)
        self.stats_tree.bind("<<TreeviewSelect>>", self._on_sequence_select)

    def _setup_tabs(self) -> None:
        self.notebook = ttk.Notebook(self.right_panel)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.tab_analysis = ttk.Frame(self.notebook)
        self.tab_orf = ttk.Frame(self.notebook)
        self.tab_alignment = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_analysis, text="🔎 Motif Analysis")
        self.notebook.add(self.tab_orf, text="▶️ ORF Detection")
        self.notebook.add(self.tab_alignment, text="🧬 Alignment")

        self._setup_analysis_tab()
        self._setup_orf_tab()
        self._setup_alignment_tab()

    def _setup_analysis_tab(self) -> None:
        top_frame = ttk.Frame(self.tab_analysis)
        top_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Label(top_frame, text="Motifs:").pack(side=tk.LEFT, padx=5)
        self.motif_entry = ttk.Entry(top_frame, width=30)
        self.motif_entry.insert(0, "ATG,CATA,CGCG")
        self.motif_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        ttk.Label(top_frame, text="Segment (bp):").pack(side=tk.LEFT, padx=(15, 5))
        self.segment_entry = ttk.Entry(top_frame, width=6)
        self.segment_entry.insert(0, str(self.segment_size))
        self.segment_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(top_frame, text="Analyze", command=self.analyze_sequence).pack(side=tk.LEFT, padx=10)
        ttk.Button(
            top_frame,
            text="Generate Heatmap",
            command=self.show_comparison_heatmap,
            style="Heatmap.TButton",
        ).pack(side=tk.LEFT, padx=10)

        # Results Paned Window
        results_paned = ttk.PanedWindow(self.tab_analysis, orient=tk.VERTICAL)
        results_paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        # Plot Canvas
        self.canvas_frame = ttk.Frame(results_paned)
        results_paned.add(self.canvas_frame, weight=2)
        self.fig, self.ax = plt.subplots(figsize=(8, 2))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._update_plot(message="Select a sequence and motifs, then click 'Analyze'.")

        # Bottom panel: Export + Text
        bottom_frame = ttk.Frame(results_paned)
        results_paned.add(bottom_frame, weight=5)

        export_frame = ttk.Frame(bottom_frame)
        export_frame.pack(fill=tk.X, pady=5)
        ttk.Button(export_frame, text="Export to CSV", command=self.export_csv,
            style="Export.TButton").pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(
            export_frame, text="Export PDF Report", command=self.export_pdf_report, style="PDF.TButton"
        ).pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        text_frame = ttk.Frame(bottom_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)
        self.results_text = scrolledtext.ScrolledText(text_frame, height=8, font=("Courier", 9),
            wrap=tk.NONE)
        self.results_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 5), pady=5)
        self.seq_text = scrolledtext.ScrolledText(text_frame, height=8, font=("Courier", 10),
            wrap=tk.WORD)
        self.seq_text.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(5, 0), pady=5)
        self.seq_text.insert(tk.INSERT, "Sequence with highlighted motifs will appear here...")
        self.seq_text.config(state=tk.DISABLED)

    def _setup_orf_tab(self) -> None:
        control_frame = ttk.Frame(self.tab_orf, padding=10)
        control_frame.pack(fill=tk.X)

        ttk.Button(
            control_frame,
            text=f"Detect ORFs (min {MIN_ORF_AA_LENGTH} aa)",
            command=self.detect_orfs,
        ).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))
        ttk.Button(
            control_frame,
            text="Export ORFs (FASTA)",
            command=self.export_orfs_to_fasta,
            style="Export.TButton",
        ).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 0))

        self.orf_text = scrolledtext.ScrolledText(
            self.tab_orf, wrap=tk.WORD, font=("Courier", 10), height=20
        )
        self.orf_text.tag_config("AA_SEQ", foreground="#00008B", font=("Courier", 10, "bold"))
        self.orf_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def _setup_alignment_tab(self) -> None:
        control_frame = ttk.Frame(self.tab_alignment, padding=10)
        control_frame.pack(fill=tk.X)

        ttk.Label(control_frame, text="Seq 1:").pack(side=tk.LEFT, padx=(0, 5))
        self.align_seq1_selector = ttk.Combobox(control_frame, state="readonly", width=25)
        self.align_seq1_selector.pack(side=tk.LEFT, padx=5)

        ttk.Label(control_frame, text="Seq 2:").pack(side=tk.LEFT, padx=(20, 5))
        self.align_seq2_selector = ttk.Combobox(control_frame, state="readonly", width=25)
        self.align_seq2_selector.pack(side=tk.LEFT, padx=5)

        ttk.Button(control_frame, text="Run Alignment", command=self.perform_alignment,
            style="Align.TButton").pack(side=tk.LEFT, padx=20)

        self.align_results_text = scrolledtext.ScrolledText(
            self.tab_alignment, wrap=tk.NONE, font=("Courier", 10)
        )
        self.align_results_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.align_results_text.insert(
            tk.INSERT,
            "Select two sequences and run alignment...\n\n"
            f"WARNING: Alignment of sequences >{LONG_SEQUENCE_THRESHOLD:,} bp may be very slow.",
        )

    # Sequence load
    def load_from_file(self) -> None:
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA/TXT", "*.fasta *.fa *.txt"), ("All files", "*.*")]
        )
        if not filepath:
            return

        path = Path(filepath)
        try:
            new_seqs = {record.id: record.seq for record in SeqIO.parse(path, "fasta")}
            if not new_seqs:
                raw = path.read_text().replace(" ", "").replace("\n", "").upper()
                if raw:
                    new_seqs[f"TXT_{len(self.sequences) + 1}"] = Seq(raw)

            if not new_seqs:
                raise ValueError("Empty or invalid file.")

            self.sequences.update(new_seqs)
            self._refresh_stats_table()
            self._update_status(f"Loaded {len(new_seqs)} sequence(s) from file.", "green")

        except Exception as e:
            messagebox.showerror("File Error", f"Failed to load file:\n{e}")
            self._update_status("ERROR", "red")

    def load_from_ncbi(self) -> None:
        raw = self.ncbi_entry.get().strip()
        if not raw:
            messagebox.showwarning(
                "Input Required", "Enter NCBI accession(s), e.g., NR_102783,MT012657")
            return

        accessions = [a.strip() for a in re.split(r"[,\s]+", raw) if a.strip()]
        if not accessions:
            return

        try:
            self._update_status("Connecting to NCBI...", "orange")
            self.master.update_idletasks()
            ids = ",".join(accessions)
            with Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text") as handle:
                new_seqs = {rec.id: rec.seq for rec in SeqIO.parse(handle, "fasta")}

            if not new_seqs:
                raise ValueError("No sequences returned.")

            self.sequences.update(new_seqs)
            self._refresh_stats_table()
            self._update_status(f"Loaded {len(new_seqs)} sequence(s) from NCBI.", "green")

        except Exception as e:
            messagebox.showerror("NCBI Error", f"Failed to fetch:\n{e}")
            self._update_status("NCBI ERROR", "red")

    # === Stats & Selection ===
    def _calculate_stats(self, seq: Seq) -> Dict[str, float]:
        s = str(seq).upper()
        length = len(s)
        counts = {b: s.count(b) for b in "GCAT"}
        gc = (counts["G"] + counts["C"]) / sum(counts.values()) * 100 if sum(counts.values()) else 0.0
        return {"length": length, "gc": gc}

    def _refresh_stats_table(self) -> None:
        for item in self.stats_tree.get_children():
            self.stats_tree.delete(item)

        seq_ids = []
        for seq_id, seq in self.sequences.items():
            stats = self._calculate_stats(seq)
            self.stats_tree.insert(
                "", "end", values=(seq_id, f"{stats['length']:,}", f"{stats['gc']:.1f}")
            )
            seq_ids.append(seq_id)

        self.align_seq1_selector["values"] = seq_ids
        self.align_seq2_selector["values"] = seq_ids
        if seq_ids:
            self.align_seq1_selector.set(seq_ids[0])
            self.align_seq2_selector.set(seq_ids[-1] if len(seq_ids) > 1 else seq_ids[0])

    def _on_sequence_select(self, event: tk.Event) -> None:
        selection = self.stats_tree.selection()
        if selection:
            self.current_seq_id = self.stats_tree.item(selection[0])["values"][0]
            self._update_status(f"Active: {self.current_seq_id}", "blue")

    # Motif analysis
    def analyze_sequence(self) -> None:
        if not self.current_seq_id:
            messagebox.showwarning("No Sequence", "Select a sequence first.")
            return

        seq = self.sequences[self.current_seq_id]
        seq_str = str(seq).upper()
        motifs = [m.strip().upper() for m in self.motif_entry.get().split(",") if m.strip()]

        try:
            self.segment_size = int(self.segment_entry.get().strip())
            if self.segment_size <= 0:
                raise ValueError
        except ValueError:
            messagebox.showwarning("Invalid Segment", "Segment size must be a positive integer.")
            return

        if not motifs:
            messagebox.showwarning("No Motifs", "Enter at least one motif.")
            return

        self.analysis_results = {}
        num_segments = (len(seq_str) + self.segment_size - 1) // self.segment_size

        for motif in motifs:
            try:
                positions = [m.start() for m in re.finditer(re.escape(motif), seq_str)]
            except re.error as e:
                messagebox.showerror("Regex Error", f"Invalid motif '{motif}': {e}")
                return

            segment_counts = [0] * num_segments
            for pos in positions:
                idx = pos // self.segment_size
                if idx < num_segments:
                    segment_counts[idx] += 1

            self.analysis_results[motif] = {
                "count": len(positions),
                "frequency": len(positions) / len(seq_str),
                "positions": positions,
                "segment_counts": segment_counts,
            }

        self._display_motif_results(seq_str)
        self._update_plot(seq_len=len(seq_str))

    def _display_motif_results(self, seq_str: str) -> None:
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete("1.0", tk.END)
        self.results_text.insert(tk.END, f"ID: {self.current_seq_id}, Length: {len(seq_str):,} bp\n")
        total = sum(d["count"] for d in self.analysis_results.values())
        self.results_text.insert(tk.END, f"Total Motifs: {total:,}\n" + "=" * 60 + "\n")
        self.results_text.insert(tk.END, f"{'Motif':<10} | {'Count':<12} | {'Freq (%)':<10} | Max/Seg\n")
        self.results_text.insert(tk.END, "-" * 60 + "\n")

        for motif, data in self.analysis_results.items():
            max_seg = max(data["segment_counts"]) if data["segment_counts"] else 0
            self.results_text.insert(
                tk.END,
                f"{motif:<10} | {data['count']:<12,} | {data['frequency']*100:>8.4f}% | {max_seg}\n",
            )
        self.results_text.config(state=tk.DISABLED)

        self.seq_text.config(state=tk.NORMAL)
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert(tk.END, seq_str)

        colors = ["red", "blue", "green", "orange", "purple", "brown"]
        for i, motif in enumerate(self.analysis_results):
            color = colors[i % len(colors)]
            self.seq_text.tag_config(motif, foreground=color, font=("Courier", 10, "bold"))
            for pos in self.analysis_results[motif]["positions"]:
                start = f"1.0 + {pos} chars"
                end = f"1.0 + {pos + len(motif)} chars"
                self.seq_text.tag_add(motif, start, end)
        self.seq_text.config(state=tk.DISABLED)

    def _update_plot(self, seq_len: Optional[int] = None, message: Optional[str] = None) -> None:
        self.ax.clear()
        if message:
            self.ax.text(0.5, 0.5, message, ha="center", va="center", transform=self.ax.transAxes)
        elif not self.analysis_results or not seq_len:
            self.ax.text(0.5, 0.5, "No results.", ha="center", va="center", transform=self.ax.transAxes)
        else:
            colors = plt.colormaps["tab10"]
            max_count = 0
            scatter_bottom = 0
            for i, (motif, data) in enumerate(self.analysis_results.items()):
                color = colors(i)
                centers = [
                    j * self.segment_size + self.segment_size / 2
                    for j in range(len(data["segment_counts"]))
                ]
                self.ax.bar(
                    centers,
                    data["segment_counts"],
                    width=self.segment_size * 0.9,
                    color=color,
                    alpha=0.6,
                    label=f"{motif} (Seg)",
                )
                y = -(i + 1)
                self.ax.scatter(
                    data["positions"], [y] * len(data["positions"]), color=color, marker="|", s=50
                )
                max_count = max(max_count, max(data["segment_counts"]))
                scatter_bottom = min(scatter_bottom, y)

            self.ax.set_title(f"Motif Distribution in {self.current_seq_id}")
            self.ax.set_xlabel("Position (bp)")
            self.ax.set_ylabel("Count in Segment")
            self.ax.set_xlim(0, seq_len)
            self.ax.set_ylim(scatter_bottom - 1, max_count * 1.2)
            self.ax.axhline(0, color="black", linewidth=0.5)
            self.ax.legend(fontsize="small", loc="upper right")
        self.fig.tight_layout(pad=3.0)
        self.canvas.draw()

    # Heatmap
    def show_comparison_heatmap(self) -> None:
        motifs = [m.strip().upper() for m in self.motif_entry.get().split(",") if m.strip()]
        if len(self.sequences) < 2:
            messagebox.showwarning("Not Enough Data", "Load at least 2 sequences.")
            return
        if not motifs:
            messagebox.showwarning("No Motifs", "Enter motifs.")
            return

        # Build DataFrame: motifs = rows (index), sequences = columns
        data = {}
        for seq_id, seq in self.sequences.items():
            counts = self._count_motifs_in_seq(str(seq).upper(), motifs)
            for motif, count in zip(motifs, counts):
                data.setdefault(motif, {})[seq_id] = count

        df = pd.DataFrame(data).T.fillna(0).astype(int)

        # Plot
        fig = plt.figure(figsize=(11, max(3, len(motifs) * 0.7 + 1.5)))  # Auto-height
        ax = fig.add_subplot(111)
        sns.heatmap(df, annot=True, fmt="d", cmap="viridis", ax=ax,
                    linewidths=0.5, cbar_kws={"shrink": 0.8})

        ax.set_title("Motif Occurrence Heatmap", pad=20)
        ax.set_xlabel("Sequence", labelpad=12)
        ax.set_ylabel("Motif", labelpad=12)

        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)

        fig.tight_layout(pad=3.0)

        # Show in window
        win = Toplevel(self.master)
        win.title("Comparison Heatmap")
        win.geometry("900x650")
        canvas = FigureCanvasTkAgg(fig, win)
        canvas.draw()
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _count_motifs_in_seq(self, seq: str, motifs: List[str]) -> List[int]:
        return [len(re.findall(re.escape(m), seq)) for m in motifs]

    # Alignment
    def perform_alignment(self) -> None:
        id1, id2 = self.align_seq1_selector.get(), self.align_seq2_selector.get()
        if not id1 or not id2 or id1 == id2:
            messagebox.showwarning("Invalid Selection", "Select two different sequences.")
            return

        seq1, seq2 = str(self.sequences[id1]), str(self.sequences[id2])
        if max(len(seq1), len(seq2)) > LONG_SEQUENCE_THRESHOLD:
            if not messagebox.askyesno("Slow Operation", "Long sequences. Continue?"):
                return

        self.align_results_text.config(state=tk.NORMAL)
        self.align_results_text.delete("1.0", tk.END)
        self.align_results_text.insert(tk.END, "Aligning...")
        self.master.update_idletasks()

        try:
            aligner = PairwiseAligner()
            aligner.mode = "fogsaa"
            alignments = aligner.align(seq1, seq2)
            self.align_results_text.delete("1.0", tk.END)
            if alignments:
                algn = alignments[0]
                output = f"--- GLOBAL ALIGNMENT ---\nScore: {algn.score}\n\n{algn}"
            else:
                output = "No alignment found."
            self.align_results_text.insert(tk.END, output)
        except Exception as e:
            self.align_results_text.insert(tk.END, f"Error: {e}")
        finally:
            self.align_results_text.config(state=tk.DISABLED)

    # ORF Detection
    def detect_orfs(self) -> None:
        if not self.current_seq_id:
            messagebox.showwarning("No Sequence", "Select a sequence.")
            return

        seq = self.sequences[self.current_seq_id]
        self.last_detected_orfs = []
        self.orf_text.config(state=tk.NORMAL)
        self.orf_text.delete("1.0", tk.END)
        self.orf_text.insert(tk.END, f"Detecting ORFs (min {MIN_ORF_AA_LENGTH} aa)...\n\n")

        orf_id = 1
        for strand, nuc in [(1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                text_line = f"\n--- Frame {frame+1}, Strand {'+' if strand == 1 else '-'} ---\n"
                self.orf_text.insert(tk.END, text_line)
                prot = str(nuc[frame:].translate())
                for match in re.finditer(r"M([^*]*)\*", prot):
                    aa_seq = match.group(0)
                    if len(aa_seq) - 1 < MIN_ORF_AA_LENGTH:
                        continue

                    dna_start = frame + match.start() * 3
                    dna_end = frame + match.end() * 3
                    if strand == -1:
                        start, end = len(seq) - dna_end, len(seq) - dna_start
                    else:
                        start, end = dna_start, dna_end

                    orf_name = f"ORF_{orf_id}_S{strand}_F{frame+1}"
                    pos = f"{start+1}-{end}"
                    text_line = f"{orf_name} [{pos}] ({end-start} bp, {len(aa_seq)-1} aa)\n"
                    self.orf_text.insert(tk.END, text_line)
                    self.orf_text.insert(tk.END, f"{aa_seq}\n\n", "AA_SEQ")

                    self.last_detected_orfs.append({
                        "id": orf_name,
                        "pos": pos,
                        "seq": aa_seq,
                        "parent": self.current_seq_id,
                    })
                    orf_id += 1
        self.orf_text.config(state=tk.DISABLED)

    def export_orfs_to_fasta(self) -> None:
        if not self.last_detected_orfs:
            messagebox.showwarning("No ORFs", "Run ORF detection first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA", "*.fasta *.faa")],
            initialfile=f"orfs_{self.current_seq_id}.fasta",
        )
        if not filepath:
            return

        try:
            with open(filepath, "w") as f:
                for orf in self.last_detected_orfs:
                    f.write(f">{orf['id']} [pos={orf['pos']}] [src={orf['parent']}]\n")
                    seq = orf["seq"]
                    for i in range(0, len(seq), FASTA_LINE_WIDTH):
                        f.write(seq[i:i + FASTA_LINE_WIDTH] + "\n")
            messagebox.showinfo("Success", f"Exported {len(self.last_detected_orfs)} ORFs to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Failed to save:\n{e}")

    def export_csv(self) -> None:
        if not self.analysis_results:
            messagebox.showwarning("No Results", "Run analysis first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV", "*.csv")],
            initialfile=f"motifs_{self.current_seq_id}.csv"
        )
        if not filepath:
            return

        rows = []
        stats = self._calculate_stats(self.sequences[self.current_seq_id])
        for motif, data in self.analysis_results.items():
            rows.append({
                "Sequence ID": self.current_seq_id,
                "Length": stats["length"],
                "GC%": f"{stats['gc']:.2f}",
                "Motif": motif,
                "Count": data["count"],
                "Frequency": f"{data['frequency']:.6f}",
                "Positions": ";".join(map(str, data["positions"])),
            })

        pd.DataFrame(rows).to_csv(filepath, index=False)
        messagebox.showinfo("Success", f"Exported to:\n{filepath}")

    def export_pdf_report(self) -> None:
        if not self.analysis_results:
            messagebox.showwarning("No Results", "Run analysis first.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".pdf", filetypes=[("PDF", "*.pdf")],
            initialfile=f"report_{self.current_seq_id}.pdf"
        )
        if not filepath:
            return

        try:
            pdf = FPDF()
            pdf.add_page()
            pdf.set_font("Helvetica", size=12)
            pdf.cell(0, 10, "DNA SEQUENCE ANALYSIS REPORT", new_x=XPos.LMARGIN, new_y=YPos.NEXT, align="C")
            pdf.ln(8)

            seq = self.sequences[self.current_seq_id]
            stats = self._calculate_stats(seq)
            pdf.cell(0, 8, f"Sequence: {self.current_seq_id}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.cell(0, 8, f"Length: {stats['length']:,} bp", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.cell(0, 8, f"GC Content: {stats['gc']:.2f}%", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.ln(8)

            # Table
            pdf.cell(60, 8, "Motif", border=1)
            pdf.cell(40, 8, "Count", border=1)
            pdf.cell(50, 8, "Freq (%)", border=1)
            pdf.ln(8)
            for motif, data in self.analysis_results.items():
                pdf.cell(60, 8, motif, border=1)
                pdf.cell(40, 8, f"{data['count']:,}", border=1)
                pdf.cell(50, 8, f"{data['frequency']*100:.4f}", border=1)
                pdf.ln(8)

            # Plot
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                self.fig.savefig(tmp.name, format="png", dpi=150)
                pdf.image(tmp.name, w=pdf.w * 0.8)
                Path(tmp.name).unlink()

            # Heatmap
            fig_hm = self._generate_heatmap_figure()
            if fig_hm:
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                    fig_hm.savefig(tmp.name, format="png", dpi=150)
                    pdf.add_page()
                    pdf.image(tmp.name, w=pdf.w * 0.8)
                    Path(tmp.name).unlink()
                plt.close(fig_hm)

            pdf.output(filepath)
            messagebox.showinfo("Success", f"Report saved to:\n{filepath}")
        except Exception as e:
            messagebox.showerror("PDF Error", f"Failed to generate PDF:\n{e}")

    def _generate_heatmap_figure(self):
        motifs = [m.strip().upper() for m in self.motif_entry.get().split(",") if m.strip()]
        if len(self.sequences) < 2 or not motifs:
            return None

        data = {}
        for seq_id, seq in self.sequences.items():
            counts = self._count_motifs_in_seq(str(seq).upper(), motifs)
            for motif, count in zip(motifs, counts):
                data.setdefault(motif, {})[seq_id] = count

        df = pd.DataFrame(data).T.fillna(0).astype(int)

        fig = plt.figure(figsize=(11, max(3, len(motifs) * 0.7 + 1.5)))
        ax = fig.add_subplot(111)
        sns.heatmap(df, annot=True, fmt="d", cmap="viridis", ax=ax,
                    linewidths=0.5, cbar_kws={"shrink": 0.8})

        ax.set_title("Motif Heatmap", pad=15)
        plt.xticks(rotation=45, ha="right")
        plt.yticks(rotation=0)
        fig.tight_layout(pad=2.5)
        return fig

    # Helpers
    def _update_status(self, text: str, color: str) -> None:
        self.status_label.config(text=f"Status: {text}", foreground=color)


# Run app
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()