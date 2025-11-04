import io
import re
import sys
import tempfile
import tkinter as tk
from tkinter import (
    filedialog, scrolledtext, messagebox, ttk, Toplevel
)

from Bio import SeqIO, Entrez
from Bio import Align
from Bio.Seq import Seq

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import seaborn as sns

from fpdf import FPDF
from fpdf.enums import XPos, YPos

# Ustawienie e-maila dla Entrez (wymagane przez NCBI)
__version__ = '1.0'
Entrez.email = "your_email@example.com"


class DNAAnalyzerAppPro:
    def __init__(self, master):
        self.master = master
        master.title(f"🧬 DNA Sequence Analyzer v{__version__} - Artemida Chadzinikolau")
        master.geometry("1200x700")

        # --- Style ---
        self.style = ttk.Style()
        self.style.theme_use('clam')
        self.style.configure("TLabel", font=('Helvetica', 10))
        self.style.configure("TButton", font=('Helvetica', 10, 'bold'))
        self.style.configure("TLabelFrame.Label", font=('Helvetica', 12, 'bold'))
        self.style.configure("Treeview.Heading", font=('Helvetica', 10, 'bold'))
        self.style.configure("Export.TButton", background="#DCEDC8", foreground="#33691E")  # Zielony
        self.style.configure("PDF.TButton", background="#BBDEFB", foreground="#0D47A1")  # Niebieski

        # --- Stan Aplikacji ---
        self.sequences = {}  # {'ID': SeqObject}
        self.current_seq_id = None
        self.analysis_results = {}
        self.last_detected_orfs = []  # Przechowuje ORF do eksportu
        self.segment_size = 100

        # --- Główny Układ ---
        self.paned_window = ttk.PanedWindow(master, orient=tk.HORIZONTAL)
        self.paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        # --- Panel Lewy ---
        self.left_panel = ttk.Frame(self.paned_window, width=480)
        self.paned_window.add(self.left_panel, weight=2)
        self._setup_left_panel()

        # --- Panel Prawy ---
        self.right_panel = ttk.Frame(self.paned_window)
        self.paned_window.add(self.right_panel, weight=4)
        self._setup_right_panel()

    # -----------------------------------------------------------------
    # Inicjalizacja Paneli UI
    # -----------------------------------------------------------------

    def _setup_left_panel(self):
        """Tworzy widżety dla lewego panelu (Wczytywanie, Statystyki)."""
        load_frame = ttk.LabelFrame(self.left_panel, text="📥 Wczytaj Sekwencje")
        load_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Button(load_frame, text="Z Pliku FASTA/TXT", command=self.load_sequence_from_file).pack(fill=tk.X, padx=5,
                                                                                                    pady=5)

        ncbi_frame = ttk.Frame(load_frame)
        ncbi_frame.pack(fill=tk.X, padx=5, pady=5)
        self.ncbi_entry = ttk.Entry(ncbi_frame, width=20)
        self.ncbi_entry.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 5))
        self.ncbi_entry.insert(0, "Accession (np. NR_102783)")
        ttk.Button(ncbi_frame, text="Pobierz z NCBI", command=self.load_sequence_from_ncbi).pack(side=tk.LEFT)

        self.status_label = ttk.Label(load_frame, text="Status: Oczekuje...", foreground="blue")
        self.status_label.pack(fill=tk.X, padx=5, pady=5)

        stats_frame = ttk.LabelFrame(self.left_panel, text="📊 Załadowane Sekwencje")
        stats_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        cols = ('ID', 'Długość (bp)', 'GC Content')
        self.stats_tree = ttk.Treeview(stats_frame, columns=cols, show='headings')
        for col in cols:
            self.stats_tree.heading(col, text=col)

        self.stats_tree.column('ID', anchor=tk.W, width=80)
        self.stats_tree.column('Długość (bp)', anchor=tk.E, width=80)
        self.stats_tree.column('GC Content', anchor=tk.E, width=80)

        scrollbar = ttk.Scrollbar(stats_frame, orient=tk.VERTICAL, command=self.stats_tree.yview)
        self.stats_tree.configure(yscroll=scrollbar.set)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.stats_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.stats_tree.bind('<<TreeviewSelect>>', self.on_sequence_select)

    def _setup_right_panel(self):
        """Tworzy notatnik (zakładki) dla prawego panelu."""
        self.notebook = ttk.Notebook(self.right_panel)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.tab_analysis = ttk.Frame(self.notebook)
        self.tab_orf = ttk.Frame(self.notebook)
        self.tab_alignment = ttk.Frame(self.notebook)

        self.notebook.add(self.tab_analysis, text="🔎 Analiza Motywów")
        self.notebook.add(self.tab_orf, text="▶️ Wykrywanie ORF")
        self.notebook.add(self.tab_alignment, text="🧬 Alignment")

        self._setup_analysis_tab()
        self._setup_orf_tab()
        self._setup_alignment_tab()

    def _setup_analysis_tab(self):
        """Konfiguracja zakładki Analiza Motywów."""
        top_frame = ttk.Frame(self.tab_analysis)
        top_frame.pack(fill=tk.X, padx=10, pady=10)

        ttk.Label(top_frame, text="Motywy:").pack(side=tk.LEFT, padx=5)
        self.motif_entry = ttk.Entry(top_frame, width=30)
        self.motif_entry.insert(0, "ATG,TATA,CGCG")
        self.motif_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        ttk.Label(top_frame, text="Segment (bp):").pack(side=tk.LEFT, padx=(10, 5))
        self.segment_entry = ttk.Entry(top_frame, width=6)
        self.segment_entry.insert(0, str(self.segment_size))
        self.segment_entry.pack(side=tk.LEFT, padx=5)

        ttk.Button(top_frame, text="Analizuj", command=self.analyze_sequence).pack(side=tk.LEFT, padx=10)

        self.style.configure("Heatmap.TButton", background="#FFC107")
        ttk.Button(top_frame, text="Generuj Heatmapę",
                   command=self.generate_comparison_heatmap,
                   style="Heatmap.TButton").pack(side=tk.LEFT, padx=10)

        # --- Panel Wyników ---
        results_paned = ttk.PanedWindow(self.tab_analysis, orient=tk.VERTICAL)
        results_paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        self.canvas_frame = ttk.Frame(results_paned)
        results_paned.add(self.canvas_frame, weight=3)
        self.fig, self.ax = plt.subplots(figsize=(8, 3))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)
        self._update_plot_canvas(None, message="Wybierz sekwencję i motywy, a następnie kliknij 'Analizuj'.")

        bottom_frame = ttk.Frame(results_paned)
        results_paned.add(bottom_frame, weight=3)

        # --- Przyciski Eksportu (NOWA LOKALIZACJA) ---
        export_frame = ttk.Frame(bottom_frame)
        export_frame.pack(fill=tk.X, padx=0, pady=5)

        ttk.Button(export_frame, text="Eksportuj do CSV", command=self.export_csv,
                   style="Export.TButton").pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        ttk.Button(export_frame, text="Eksportuj Raport PDF (z Wykresami)", command=self.export_pdf_report,
                   style="PDF.TButton").pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)

        # --- Pola tekstowe wyników ---
        text_frame = ttk.Frame(bottom_frame)
        text_frame.pack(fill=tk.BOTH, expand=True)

        self.results_text = scrolledtext.ScrolledText(text_frame, height=8, font=('Courier', 9), wrap=tk.NONE)
        self.results_text.pack(fill=tk.BOTH, expand=True, padx=(0, 5), pady=5, side=tk.LEFT)

        self.seq_text = scrolledtext.ScrolledText(text_frame, wrap=tk.WORD, height=8, font=('Courier', 10))
        self.seq_text.pack(fill=tk.BOTH, expand=True, padx=(5, 0), pady=5, side=tk.RIGHT)
        self.seq_text.insert(tk.INSERT, "Tutaj pojawi się sekwencja z zaznaczonymi motywami...")
        self.seq_text.config(state=tk.DISABLED)

    def _setup_orf_tab(self):
        """Konfiguracja zakładki Wykrywanie ORF."""
        control_frame = ttk.Frame(self.tab_orf, padding=10)
        control_frame.pack(fill=tk.X)

        ttk.Button(control_frame, text="Wykryj ORF (min. 10 aa) dla aktywnej sekwencji",
                   command=self.detect_orfs).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 10))

        # NOWY PRZYCISK: Eksport ORF
        ttk.Button(control_frame, text="Eksportuj Sekwencje ORF (FASTA)",
                   command=self.export_orfs_to_fasta,
                   style="Export.TButton").pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(10, 0))

        self.orf_text = scrolledtext.ScrolledText(self.tab_orf, wrap=tk.WORD, font=('Courier', 10))
        self.orf_text.tag_config("AA_SEQ", font=('Courier', 10), foreground="#00008B")
        self.orf_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    def _setup_alignment_tab(self):
        """Konfiguracja zakładki Alignmentu."""
        control_frame = ttk.Frame(self.tab_alignment, padding=10)
        control_frame.pack(fill=tk.X)

        ttk.Label(control_frame, text="Sekwencja 1:").pack(side=tk.LEFT, padx=(0, 5))
        self.align_seq1_selector = ttk.Combobox(control_frame, state="readonly", width=25)
        self.align_seq1_selector.pack(side=tk.LEFT, padx=5)

        ttk.Label(control_frame, text="Sekwencja 2:").pack(side=tk.LEFT, padx=(20, 5))
        self.align_seq2_selector = ttk.Combobox(control_frame, state="readonly", width=25)
        self.align_seq2_selector.pack(side=tk.LEFT, padx=5)

        ttk.Button(control_frame, text="Uruchom Alignment", command=self.perform_alignment).pack(side=tk.LEFT, padx=20)

        self.align_results_text = scrolledtext.ScrolledText(
            self.tab_alignment, wrap=tk.NONE, font=('Courier', 10)
        )
        self.align_results_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        self.align_results_text.insert(tk.INSERT,
                                       "Wybierz dwie sekwencje i uruchom alignment...\n\nOSTRZEŻENIE: Alignment bardzo długich sekwencji (>10 000 bp) może być BARDZO powolny.")

    # -----------------------------------------------------------------
    # Metody Wczytywania i Zarządzania Danymi
    # -----------------------------------------------------------------

    def load_sequence_from_file(self):
        """Wczytuje sekwencje z pliku FASTA lub TXT."""
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA/TXT files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
        )
        if not filepath: return
        try:
            new_sequences = {}
            for record in SeqIO.parse(filepath, "fasta"):
                new_sequences[record.id] = record.seq

            if not new_sequences:
                with open(filepath, 'r') as f:
                    raw_seq = "".join(f.read().split()).upper().replace(' ', '')
                    if raw_seq: new_sequences[f"Plik_TXT_{len(self.sequences) + 1}"] = Seq(raw_seq)

            if not new_sequences: raise ValueError("Plik jest pusty lub nie zawiera prawidłowej sekwencji.")

            self.sequences.update(new_sequences)
            self.refresh_global_stats_table()
            self.status_label.config(text=f"Status: Wczytano {len(new_sequences)} sekwencj(i).", foreground="green")

        except Exception as e:
            messagebox.showerror("Błąd Wczytywania", f"Nie udało się wczytać pliku: {e}")
            self.status_label.config(text="Status: BŁĄD", foreground="red")

    def load_sequence_from_ncbi(self):
        """Pobiera sekwencję z bazy NCBI."""
        accession = self.ncbi_entry.get().strip()
        if not accession or accession == "Accession (np. NR_102783)":
            messagebox.showwarning("Brak Danych", "Wprowadź numer dostępu NCBI.")
            return
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()

            self.sequences[record.id] = record.seq
            self.refresh_global_stats_table()
            self.status_label.config(text=f"Status: Wczytano NCBI {record.id}.", foreground="blue")

        except Exception as e:
            messagebox.showerror("Błąd NCBI", f"Nie udało się pobrać sekwencji dla {accession}: {e}")
            self.status_label.config(text="Status: BŁĄD NCBI", foreground="red")

    def calculate_sequence_stats(self, seq_str):
        """Manualne obliczanie statystyk (GC, N)."""
        seq_str = seq_str.upper()
        length = len(seq_str)
        g_c_a_t = (seq_str.count('G'), seq_str.count('C'), seq_str.count('A'), seq_str.count('T'))

        gc_denominator = sum(g_c_a_t)
        gc_percent = ((g_c_a_t[0] + g_c_a_t[1]) / gc_denominator) * 100 if gc_denominator > 0 else 0.0

        return {'length': length, 'gc': gc_percent}

    def refresh_global_stats_table(self):
        """Odświeża tabelę statystyk ORAZ listy w zakładce Alignmentu."""
        self.stats_tree.delete(*self.stats_tree.get_children())
        seq_ids = []

        for seq_id, seq_obj in self.sequences.items():
            stats = self.calculate_sequence_stats(str(seq_obj))
            self.stats_tree.insert(
                '', 'end',
                values=(
                    seq_id, f"{stats['length']:,d}", f"{stats['gc']:.2f}%"
                )
            )
            seq_ids.append(seq_id)

        self.align_seq1_selector['values'] = seq_ids
        self.align_seq2_selector['values'] = seq_ids
        if seq_ids:
            self.align_seq1_selector.set(seq_ids[0])
            if len(seq_ids) > 1:
                self.align_seq2_selector.set(seq_ids[1])
            else:
                self.align_seq2_selector.set(seq_ids[0])

    def on_sequence_select(self, event):
        """Aktualizuje `current_seq_id` po kliknięciu w tabelę."""
        selected_item = self.stats_tree.focus()
        if selected_item:
            self.current_seq_id = self.stats_tree.item(selected_item)['values'][0]
            self.status_label.config(text=f"Aktywna sekwencja: {self.current_seq_id}", foreground="blue")

    # -----------------------------------------------------------------
    # Metody Analityczne (Analiza Motywów)
    # -----------------------------------------------------------------

    def analyze_sequence(self):
        """Główna funkcja analityczna: wyszukiwanie motywów i wizualizacja."""
        if not self.current_seq_id or self.current_seq_id not in self.sequences:
            messagebox.showwarning("Brak Sekwencji", "Najpierw wybierz aktywną sekwencję z tabeli.")
            return

        seq = self.sequences[self.current_seq_id]
        seq_str = str(seq).upper()

        motifs_raw = self.motif_entry.get().upper().replace(' ', '').split(',')
        motifs = [m for m in motifs_raw if m]

        try:
            self.segment_size = int(self.segment_entry.get())
            if self.segment_size <= 0: raise ValueError
        except ValueError:
            messagebox.showwarning("Błąd Segmentu", "Rozmiar segmentu musi być liczbą dodatnią.")
            return
        if not motifs:
            messagebox.showwarning("Brak Motywów", "Wprowadź motywy do analizy.")
            return

        self.analysis_results = {}
        num_segments = (len(seq_str) + self.segment_size - 1) // self.segment_size

        for motif in motifs:
            try:
                positions = [m.start() for m in re.finditer(f'(?={re.escape(motif)})', seq_str)]
            except re.error as e:
                messagebox.showerror("Błąd Wyrażenia Regularnego", f"Błąd w motywie '{motif}': {e}")
                return

            if positions:
                segment_counts = [0] * num_segments
                for pos in positions:
                    segment_index = pos // self.segment_size
                    if segment_index < num_segments: segment_counts[segment_index] += 1

                self.analysis_results[motif] = {
                    "liczba": len(positions),
                    "czestotliwosc": len(positions) / len(seq_str),
                    "pozycje": positions,
                    "segment_counts": segment_counts
                }

        self._display_analysis_results(seq_str, self.analysis_results)
        self._update_plot_canvas(self.analysis_results, seq_len=len(seq_str))

    def _display_analysis_results(self, seq_str, results):
        """Wyświetla wyniki tekstowe (tabelaryczne) i pokolorowaną sekwencję."""
        self.results_text.config(state=tk.NORMAL)
        self.results_text.delete(1.0, tk.END)
        header = f"ID: {self.current_seq_id}, Długość: {len(seq_str):,d} bp\n"
        self.results_text.insert(tk.INSERT, header)

        total_motifs = sum(data['liczba'] for data in results.values())
        self.results_text.insert(tk.INSERT, f"Suma Wszystkich Motywów: {total_motifs:,d}\n" + "=" * 60 + "\n")

        # POPRAWKA 3: Zmieniony nagłówek dla jasności
        self.results_text.insert(tk.INSERT,
                                 f"{'Motyw':<10} | {'Suma (Całk.)':<12} | {'Częst. (%)':<10} | Max w Segmencie\n")
        self.results_text.insert(tk.INSERT, "-" * 60 + "\n")

        for motif, data in results.items():
            max_seg = max(data.get('segment_counts', [0]))
            self.results_text.insert(tk.INSERT,
                                     f"{motif:<10} | {data['liczba']:<12,d} | {data['czestotliwosc'] * 100:<9.4f}% | {max_seg}\n"
                                     )
        self.results_text.config(state=tk.DISABLED)

        self.seq_text.config(state=tk.NORMAL)
        self.seq_text.delete(1.0, tk.END)
        self.seq_text.insert(tk.INSERT, seq_str)

        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown']
        color_map = {motif: colors[i % len(colors)] for i, motif in enumerate(results.keys())}

        for motif, color in color_map.items():
            self.seq_text.tag_config(motif, foreground=color, font=('Courier', 10, 'bold'))
        for motif, data in results.items():
            for pos in data['pozycje']:
                self.seq_text.tag_add(motif, f"1.0+{pos}c", f"1.0+{pos + len(motif)}c")

        self.seq_text.config(state=tk.DISABLED)

    def _update_plot_canvas(self, results, seq_len=None, message=None):
        """Wykres hybrydowy (słupki + punkty)."""
        self.ax.clear()

        if message:
            self.ax.text(0.5, 0.5, message, ha='center', va='center', transform=self.ax.transAxes, wrap=True)
        elif not results or not seq_len:
            self.ax.text(0.5, 0.5, 'Brak wyników do wizualizacji.', ha='center', va='center',
                         transform=self.ax.transAxes)
        else:
            colors = plt.colormaps['tab10']
            all_counts = [0]
            y_scatter_pos = 0

            for i, (motif, data) in enumerate(results.items()):
                color = colors(i)
                segment_centers = [j * self.segment_size + self.segment_size / 2 for j in
                                   range(len(data['segment_counts']))]
                self.ax.bar(segment_centers, data['segment_counts'], width=self.segment_size * 0.9, color=color,
                            alpha=0.6, label=f'{motif} (Segmenty)')

                y_val = - (i + 1)
                self.ax.scatter(data['pozycje'], [y_val] * len(data['pozycje']), color=color, marker='|', s=50,
                                label=f'{motif} (Pozycje)')

                all_counts.extend(data['segment_counts'])
                y_scatter_pos = y_val

            self.ax.set_title(f'Rozmieszczenie Motywów w Sekwencji {self.current_seq_id}')
            self.ax.set_xlabel('Pozycja w Sekwencji (bp)')
            self.ax.set_ylabel('Liczba Wystąpień w Segmencie')
            self.ax.set_xlim(0, seq_len)
            self.ax.set_ylim(bottom=y_scatter_pos - 1, top=max(all_counts) * 1.1 + 1)
            self.ax.axhline(0, color='black', linewidth=0.5)
            ticks = self.ax.get_yticks()
            self.ax.set_yticks([tick for tick in ticks if tick >= 0])
            self.ax.legend(loc='upper right', fontsize='small')

        self.fig.tight_layout()
        self.canvas.draw()

    # -----------------------------------------------------------------
    # Metoda Porównania (Heatmapa)
    # -----------------------------------------------------------------

    def generate_comparison_heatmap(self, return_fig=False):
        """Generuje heatmapę porównawczą. Opcjonalnie zwraca obiekt figury (dla PDF)."""
        motifs_raw = self.motif_entry.get().upper().replace(' ', '').split(',')
        motifs = [m for m in motifs_raw if m]

        if not motifs:
            messagebox.showwarning("Brak Motywów", "Wprowadź motywy do analizy, aby wygenerować heatmapę.")
            return None
        if len(self.sequences) < 2:
            messagebox.showwarning("Niewystarczające Dane",
                                   "Wczytaj co najmniej dwie sekwencje, aby wygenerować heatmapę.")
            return None

        data = []
        seq_ids = []
        for seq_id, seq_obj in self.sequences.items():
            counts = self._run_motif_search(str(seq_obj).upper(), motifs)
            data.append(counts)
            seq_ids.append(seq_id)

        try:
            df = pd.DataFrame(data, index=seq_ids).T  # Transpozycja
        except Exception as e:
            messagebox.showerror("Błąd Tworzenia Danych", f"Nie udało się przetworzyć danych do heatmapy: {e}")
            return None

        fig_hm = plt.figure(figsize=(10, 8))
        ax_hm = fig_hm.add_subplot(111)
        sns.heatmap(df, ax=ax_hm, annot=True, fmt='d', cmap='viridis', linewidths=.5)
        ax_hm.set_title('Heatmapa Porównawcza Wystąpień Motywów')
        ax_hm.set_xlabel('Sekwencja')
        ax_hm.set_ylabel('Motyw')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        fig_hm.tight_layout()

        if return_fig:
            return fig_hm  # Zwraca figurę dla PDF

        # Jeśli nie dla PDF, pokaż w nowym oknie
        hm_window = Toplevel(self.master)
        hm_window.title("Heatmapa Porównawcza Motywów")
        hm_window.geometry("800x600")

        canvas_hm = FigureCanvasTkAgg(fig_hm, master=hm_window)
        canvas_hm.draw()
        toolbar = NavigationToolbar2Tk(canvas_hm, hm_window)
        toolbar.update()
        canvas_hm.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    def _run_motif_search(self, seq_str, motifs):
        """Pomocnicza funkcja do liczenia motywów."""
        results = {}
        for motif in motifs:
            try:
                count = len(re.findall(f'(?={re.escape(motif)})', seq_str))
                results[motif] = count
            except re.error:
                results[motif] = 0
        return results

    # -----------------------------------------------------------------
    # Metody Zakładki Alignmentu
    # -----------------------------------------------------------------

    def perform_alignment(self):
        """Uruchamia alignment globalny i lokalny dla dwóch wybranych sekwencji."""
        id1 = self.align_seq1_selector.get()
        id2 = self.align_seq2_selector.get()

        if not id1 or not id2:
            messagebox.showwarning("Brak Sekwencji", "Wybierz dwie sekwencje do alignmentu.")
            return
        if id1 == id2:
            messagebox.showwarning("Błąd", "Wybierz dwie *różne* sekwencje.")
            return

        try:
            seq1 = str(self.sequences[id1])
            seq2 = str(self.sequences[id2])
        except KeyError:
            messagebox.showerror("Błąd", "Nie można odnaleźć wybranych sekwencji.")
            return

        if len(seq1) > 10000 or len(seq2) > 10000:
            if not messagebox.askyesno("Ostrzeżenie o Wydajności",
                                       "Jedna z sekwencji jest bardzo długa (>10 000 bp). "
                                       "Alignment może zająć BARDZO dużo czasu. Kontynuować?"):
                return

        self.align_results_text.config(state=tk.NORMAL)
        self.align_results_text.delete(1.0, tk.END)
        self.align_results_text.insert(tk.INSERT, "Przetwarzanie... To może potrwać chwilę...")
        self.master.update_idletasks()

        try:

            global_aligner = Align.PairwiseAligner()
            global_aligner.mode = 'global'
            global_aligns = global_aligner.align(seq1, seq2)

            local_aligner = Align.PairwiseAligner()
            local_aligner.mode = 'local'
            local_aligns = local_aligner.align(seq1, seq2)

            self.align_results_text.delete(1.0, tk.END)

            if global_aligns:
                g_align = global_aligns[0]
                output_global = (
                    f"--- ALIGNMENT GLOBALNY (globalxx) ---\n"
                    f"Score: {g_align.score}\n\n"
                    f"{g_align}"
                )
            else:
                output_global = "--- ALIGNMENT GLOBALNY ---\nNie znaleziono alignmentu.\n"

            if local_aligns:
                l_align = local_aligns[0]
                output_local = (
                    f"--- ALIGNMENT LOKALNY (localxx) ---\n"
                    f"Score: {l_align.score}\n\n"
                    f"{l_align}"
                )
            else:
                output_local = "--- ALIGNMENT LOKALNY ---\nNie znaleziono alignmentu.\n"

            self.align_results_text.insert(tk.INSERT, output_global + "\n" + "=" * 70 + "\n\n" + output_local)

        except Exception as e:
            self.align_results_text.delete(1.0, tk.END)
            self.align_results_text.insert(tk.INSERT, f"Wystąpił błąd podczas alignmentu: {e}")

        self.align_results_text.config(state=tk.DISABLED)

    # -----------------------------------------------------------------
    # Metody Zakładki ORF (Wykrywanie i Eksport)
    # -----------------------------------------------------------------

    def detect_orfs(self):
        """Wykrywa i wyświetla ORF oraz zapisuje je do eksportu."""
        if not self.current_seq_id or self.current_seq_id not in self.sequences:
            messagebox.showwarning("Brak Sekwencji", "Najpierw wybierz aktywną sekwencję z tabeli.")
            return

        seq = self.sequences[self.current_seq_id]
        self.last_detected_orfs = []  # Resetowanie listy

        self.orf_text.config(state=tk.NORMAL)
        self.orf_text.delete(1.0, tk.END)

        min_orf_len_aa = 10  # Minimalna długość

        self.orf_text.insert(tk.INSERT, f"Wykrywanie ORF dla {self.current_seq_id} (min. {min_orf_len_aa} aa)...\n\n")

        orf_counter = 1
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                self.orf_text.insert(tk.INSERT, f"\n--- Ramka Odczytu {frame + 1}, Nić {strand} ---\n")
                trans_seq = str(nuc[frame:].translate())

                for match in re.finditer(r'M([^*]*)\*', trans_seq):
                    orf_aa_seq = match.group(0)
                    orf_aa_len_protein = len(orf_aa_seq) - 1  # Bez stop

                    if orf_aa_len_protein >= min_orf_len_aa:
                        orf_start_pos_dna = frame + match.start() * 3
                        orf_end_pos_dna = frame + match.end() * 3

                        if strand == -1:
                            start_final = len(seq) - orf_end_pos_dna
                            end_final = len(seq) - orf_start_pos_dna
                        else:
                            start_final = orf_start_pos_dna
                            end_final = orf_end_pos_dna

                        orf_id = f"ORF_{orf_counter}_Strand{strand}_Frame{frame + 1}"
                        orf_pos = f"{start_final + 1}-{end_final}"

                        self.orf_text.insert(tk.INSERT,
                                             f"{orf_id} [{orf_pos}] ({end_final - start_final} bp, {orf_aa_len_protein} aa)\n"
                                             )
                        self.orf_text.insert(tk.INSERT, f"{orf_aa_seq}\n\n", "AA_SEQ")

                        # Zapisanie do eksportu
                        self.last_detected_orfs.append({
                            "id": orf_id,
                            "pos": orf_pos,
                            "seq": orf_aa_seq,
                            "parent_id": self.current_seq_id
                        })
                        orf_counter += 1

        self.orf_text.config(state=tk.DISABLED)

    def export_orfs_to_fasta(self):
        """NOWA FUNKCJA: Zapisuje wykryte ORF do pliku FASTA."""
        if not self.last_detected_orfs:
            messagebox.showwarning("Brak Danych",
                                   "Najpierw uruchom 'Wykryj ORF', aby wygenerować sekwencje do eksportu.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".fasta",
            filetypes=[("FASTA amino acid", "*.fasta *.faa *.fa"), ("All files", "*.*")],
            initialfile=f"orf_export_{self.current_seq_id}.fasta"
        )
        if not filepath:
            return

        try:
            with open(filepath, 'w') as f:
                for orf in self.last_detected_orfs:
                    # Tworzenie nagłówka FASTA
                    header = f">{orf['id']} [position={orf['pos']}] [source={orf['parent_id']}]\n"
                    f.write(header)

                    # Zapis sekwencji z podziałem linii (np. co 70 znaków)
                    seq = orf['seq']
                    for i in range(0, len(seq), 70):
                        f.write(seq[i:i + 70] + '\n')

            messagebox.showinfo("Eksport Zakończony",
                                f"Pomyślnie wyeksportowano {len(self.last_detected_orfs)} sekwencji ORF do:\n{filepath}")

        except Exception as e:
            messagebox.showerror("Błąd Eksportu", f"Nie udało się zapisać pliku: {e}")

    # -----------------------------------------------------------------
    # Metody Eksportu (Analiza Motywów)
    # -----------------------------------------------------------------

    def export_csv(self):
        """Eksportuje wyniki analizy motywów do pliku CSV."""
        if not self.analysis_results:
            messagebox.showwarning("Brak Wyników", "Najpierw przeprowadź analizę sekwencji.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile=f"analysis_{self.current_seq_id}.csv"
        )
        if not filepath:
            return

        data_for_df = []
        seq = self.sequences[self.current_seq_id]
        stats = self.calculate_sequence_stats(str(seq))

        for motif, data in self.analysis_results.items():
            df_row = {
                "Sekwencja ID": self.current_seq_id, "Długość Sekwencji": stats['length'],
                "GC Content": f"{stats['gc']:.2f}%",
                "Motyw": motif, "Suma (Całk.)": data['liczba'],
                "Częstotliwość": f"{data['czestotliwosc']:.6f}",
                "Pozycje (rozd. średnikiem)": ";".join(map(str, data['pozycje'])),
            }
            data_for_df.append(df_row)

        try:
            df = pd.DataFrame(data_for_df)
            df.to_csv(filepath, index=False)
            messagebox.showinfo("Eksport Zakończony", f"Pomyślnie wyeksportowano wyniki CSV do:\n{filepath}")
        except Exception as e:
            messagebox.showerror("Błąd Eksportu CSV", f"Nie udało się zapisać pliku: {e}")

    def export_pdf_report(self):
        """Eksport PDF używając FPDF2 zamiast ReportLab."""
        if not self.analysis_results:
            messagebox.showwarning("Brak Wyników", "Najpierw przeprowadź analizę sekwencji.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf")],
            initialfile=f"report_{self.current_seq_id}.pdf"
        )
        if not filepath: return

        try:
            pdf = FPDF()
            pdf.set_auto_page_break(auto=True, margin=15)
            pdf.add_page()
            pdf.set_font('Helvetica', '', 12)

            pdf.cell(0, 10, "RAPORT Z ANALIZY SEKWENCJI DNA", align="C")
            pdf.ln(10)

            # --- Sekcja Statystyki ---
            seq = self.sequences[self.current_seq_id]
            stats = self.calculate_sequence_stats(str(seq))
            pdf.set_font("Helvetica", '', 12)
            pdf.cell(0, 8, f"Analizowana Sekwencja: {self.current_seq_id}", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.cell(0, 8, f"Dlugosc: {stats['length']:,d} bp", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.cell(0, 8, f"GC Content: {stats['gc']:.2f}%", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.ln(10)

            # --- Tabela Motywów ---
            pdf.cell(60, 8, "Motyw", 1)
            pdf.cell(40, 8, "Suma (Calk.)", 1)
            pdf.cell(50, 8, "Czestotliwosc (%)", 1)
            pdf.ln()
            for motif, data in self.analysis_results.items():
                pdf.cell(60, 8, motif, 1)
                pdf.cell(40, 8, f"{data['liczba']:,d}", 1)
                pdf.cell(50, 8, f"{data['czestotliwosc'] * 100:.4f}", 1)
                pdf.ln()
            pdf.ln(10)

            # --- Wykres Rozmieszczenia Motywów ---
            buf_plot = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
            self.fig.savefig(buf_plot.name, format='PNG', dpi=150)
            pdf.image(buf_plot.name, w=pdf.w * 0.8)
            buf_plot.close()

            # --- Heatmapa (opcjonalnie) ---
            fig_hm = self.generate_comparison_heatmap(return_fig=True)
            if fig_hm:
                buf_hm = tempfile.NamedTemporaryFile(suffix=".png", delete=False)
                fig_hm.savefig(buf_hm.name, format='PNG', dpi=150)
                pdf.add_page()
                pdf.image(buf_hm.name, w=pdf.w * 0.8)
                buf_hm.close()
                plt.close(fig_hm)

            pdf.output(filepath)
            messagebox.showinfo("Eksport Zakończony", f"Pomyślnie wyeksportowano raport PDF do:\n{filepath}")

        except Exception as e:
            messagebox.showerror("Błąd Eksportu PDF", f"Nie udało się utworzyć pliku PDF. Błąd: {e}")

    def _generate_csv_data(self):
        """Pomocnicza funkcja do tworzenia DataFrame dla CSV (używana w 2 miejscach)."""
        if not self.analysis_results:
            return pd.DataFrame()  # Pusty DataFrame

        data_for_df = []
        seq = self.sequences[self.current_seq_id]
        stats = self.calculate_sequence_stats(str(seq))

        for motif, data in self.analysis_results.items():
            df_row = {
                "Sekwencja ID": self.current_seq_id, "Długość Sekwencji": stats['length'],
                "GC Content": f"{stats['gc']:.2f}%",
                "Motyw": motif, "Suma (Całk.)": data['liczba'],
                "Częstotliwość": f"{data['czestotliwosc']:.6f}",
                "Pozycje (rozd. średnikiem)": ";".join(map(str, data['pozycje'])),
            }
            data_for_df.append(df_row)

        return pd.DataFrame(data_for_df)


# ---------------------------------
# Uruchomienie Aplikacji
# ---------------------------------
if __name__ == '__main__':
    root = tk.Tk()
    app = DNAAnalyzerAppPro(root)
    root.mainloop()