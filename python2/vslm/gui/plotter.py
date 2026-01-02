# python2/vslm/gui/plotter.py
import numpy as np
from matplotlib.figure import Figure
from .. import leq

class ResultPlotter:
    """
    Handles drawing analysis results onto a Matplotlib Figure.
    """
    @staticmethod
    def plot(fig: Figure, 
             results: list, 
             mode_id: int, 
             weighting: str, 
             speed: str, 
             leq_interval_txt: str, 
             block_size_ms: float):
        
        fig.clear()
        
        if not results: 
            return

        match mode_id:
            case 1: # LEQ MODE
                ResultPlotter._plot_leq_dashboard(fig, results, weighting, leq_interval_txt, block_size_ms)

            case 0: # LEVEL VS TIME (Lp)
                ResultPlotter._plot_lp_history(fig, results, weighting, speed)

            case 2 | 3: # SPECTRAL (Octave or Third)
                # mode_id 2 is Octave, 3 is Third Octave
                is_third = (mode_id == 3)
                ResultPlotter._plot_spectrum(fig, results, weighting, is_third)

    @staticmethod
    def _plot_leq_dashboard(fig, results, weighting, interval_txt, block_size_ms):
        # 1. Parse Interval
        match interval_txt:
            case "100 ms": interval = 0.1
            case "1 sec": interval = 1.0
            case "10 sec": interval = 10.0
            case "1 min": interval = 60.0
            case "15 min": interval = 900.0
            case "1 hour": interval = 3600.0
            case _: interval = 1.0
        
        # 2. Calculate Stats
        stats = leq.calculate_leq_analysis(results, block_size_ms, interval)
        
        # 3. Top Plot (Time History Step Plot)
        ax1 = fig.add_subplot(2, 1, 1)
        if len(stats.history['time']) > 0:
            t_plot = list(stats.history['time'])
            # Append end point for step plot closure
            t_plot.append(t_plot[-1] + interval)
            l_plot = list(stats.history['leq'])
            l_plot.append(l_plot[-1])
            
            ax1.step(t_plot, l_plot, where='post', color='b', linewidth=1.5)
        
        ax1.set_title(f"LEQ History ({interval_txt} interval, {weighting}-weighted)")
        ax1.set_ylabel("LEQ (dB)")
        ax1.grid(True)
        
        # 4. Bottom Panel (Text Dashboard)
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.axis('off')
        
        # Layout Columns
        col1, col2, col3 = 0.05, 0.35, 0.65
        
        # Header
        ax2.text(0.5, 0.95, f"Overall LEQ: {stats.overall:.1f} dB", 
                 ha='center', fontsize=14, fontweight='bold', color='blue')
        
        # Column 1 (Extremes & Percentiles)
        ax2.text(col1, 0.80, f"Lmax: {stats.max:.1f} dB")
        ax2.text(col1, 0.65, f"Lmin: {stats.min:.1f} dB")
        ax2.text(col1, 0.50, f"L10: {stats.ln[10]:.1f} dB")
        ax2.text(col1, 0.35, f"L50: {stats.ln[50]:.1f} dB")
        ax2.text(col1, 0.20, f"L90: {stats.ln[90]:.1f} dB")
        
        # Column 2 (More Percentiles)
        ax2.text(col2, 0.80, f"L20: {stats.ln[20]:.1f} dB")
        ax2.text(col2, 0.65, f"L30: {stats.ln[30]:.1f} dB")
        ax2.text(col2, 0.50, f"L40: {stats.ln[40]:.1f} dB")
        ax2.text(col2, 0.35, f"L60: {stats.ln[60]:.1f} dB")
        ax2.text(col2, 0.20, f"L80: {stats.ln[80]:.1f} dB")
        
        # Column 3 (Dose)
        ax2.text(col3, 0.80, f"Dose ({stats.dose['standard']})", fontweight='bold')
        ax2.text(col3, 0.65, f"Dose %: {stats.dose['dose']:.1f}%")
        ax2.text(col3, 0.50, f"TWA: {stats.dose['twa']:.1f} dB")
        
        fig.tight_layout()

    @staticmethod
    def _plot_lp_history(fig, results, weighting, speed):
        ax = fig.add_subplot(1, 1, 1)
        t = [r['time'] for r in results]
        l = [r['lp'] for r in results]
        ax.plot(t, l)
        ax.set_title(f"Sound Pressure Level vs Time ({weighting}-weighted, {speed})")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Level (dB)")
        ax.grid(True)

    @staticmethod
    def _plot_spectrum(fig, results, weighting, is_third_octave):
        ax = fig.add_subplot(1, 1, 1)
        freqs = results[0]['band_freqs']
        
        # Energy Average
        energy_sums = np.zeros(len(freqs))
        for r in results:
            pressures = (10**(r['bands']/10.0)) * (20e-6**2)
            energy_sums += pressures
        mean_db = 10 * np.log10((energy_sums / len(results)) / (20e-6**2) + 1e-30)
        
        x = np.arange(len(freqs))
        ax.bar(x, mean_db, color='#2ca02c', alpha=0.8)
        ax.set_xticks(x)
        
        # Labels
        lbls = []
        for f in freqs:
            if f >= 1000: lbls.append(f"{f/1000:.0f}k")
            else: lbls.append(f"{f:.0f}")
        
        # Sparse labels for 1/3 octave to prevent crowding
        if is_third_octave: 
            lbls = [l if i % 3 == 0 else "" for i, l in enumerate(lbls)]
            
        ax.set_xticklabels(lbls, rotation=90)
        ax.set_title(f"Average Spectrum ({weighting}-weighted)")
        ax.set_ylabel("Level (dB)")
        ax.grid(axis='y')