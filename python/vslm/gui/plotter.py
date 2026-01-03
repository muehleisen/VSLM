import matplotlib.pyplot as plt
import numpy as np
from .. import leq
import traceback

class ResultPlotter:
    @staticmethod
    def plot(figure, results, mode_id, weighting, speed, leq_int_txt, block_size_ms, 
             dose_params=None, ref_pressure=20e-6, 
             # NEW parameters
             autoscale=True, ymin=0.0, ymax=120.0):
        
        figure.clear()
        
        if not results:
            ax = figure.add_subplot(111)
            ax.text(0.5, 0.5, "No Data", ha='center', va='center')
            return

        try:
            # Dispatch based on Mode
            if mode_id == 1: # LEQ
                ResultPlotter._plot_leq_dashboard(
                    figure, results, weighting, leq_int_txt, 
                    block_size_ms, dose_params, ref_pressure,
                    autoscale, ymin, ymax
                )
            elif mode_id == 0: # Lp
                ResultPlotter._plot_lp_history(
                    figure, results, weighting, speed,
                    autoscale, ymin, ymax
                )
            elif mode_id in [2, 3]: # Spectrum
                is_third = (mode_id == 3)
                ResultPlotter._plot_spectrum(
                    figure, results, weighting, is_third, ref_pressure,
                    autoscale, ymin, ymax
                )
        except Exception as e:
            # Fallback if specific plotter fails
            ax = figure.add_subplot(111)
            ax.text(0.5, 0.5, f"Plot Error:\n{str(e)}", ha='center', va='center', color='red')
            print(f"Plot Error Traceback:\n{traceback.format_exc()}")

        figure.tight_layout()

    @staticmethod
    def _plot_leq_dashboard(fig, results, weighting, interval_txt, 
                            block_size_ms, dose_params, ref_pressure,
                            autoscale, ymin, ymax):
        
        # Safe Interval Parsing
        intervals = {"100 ms": 0.1, "1 sec": 1.0, "10 sec": 10.0, 
                     "1 min": 60.0, "15 min": 900.0, "1 hour": 3600.0}
        interval = intervals.get(interval_txt, 1.0)
        
        # Calculate Stats
        stats = leq.calculate_leq_analysis(
            results, block_size_ms, interval, dose_params, ref_pressure
        )
        
        # Top Plot (Time History)
        ax1 = fig.add_subplot(2, 1, 1)
        if len(stats.history['time']) > 0:
            t_plot = list(stats.history['time'])
            # Extend last point for step plot
            t_plot.append(t_plot[-1] + interval)
            l_plot = list(stats.history['leq'])
            l_plot.append(l_plot[-1])
            
            ax1.step(t_plot, l_plot, where='post', color='b', linewidth=1.5)
            
            # Apply Scaling
            if not autoscale:
                ax1.set_ylim(ymin, ymax)
            else:
                ax1.autoscale(axis='y')

        ax1.set_title(f"LEQ History ({interval_txt} interval, {weighting}-weighted)")
        ax1.set_ylabel("LEQ (dB)")
        ax1.grid(True)
        
        # Bottom Panel (Text Dashboard)
        ax2 = fig.add_subplot(2, 1, 2)
        ax2.axis('off')
        
        col1, col2, col3 = 0.05, 0.35, 0.65
        
        # Overall
        ax2.text(0.5, 0.95, f"Overall LEQ: {stats.overall:.1f} dB", 
                 ha='center', fontsize=14, fontweight='bold', color='blue')
        
        # Stats Columns
        ax2.text(col1, 0.80, f"Lmax: {stats.max:.1f} dB")
        ax2.text(col1, 0.65, f"Lmin: {stats.min:.1f} dB")
        ax2.text(col1, 0.50, f"L10: {stats.ln[10]:.1f} dB")
        ax2.text(col1, 0.35, f"L50: {stats.ln[50]:.1f} dB")
        ax2.text(col1, 0.20, f"L90: {stats.ln[90]:.1f} dB")
        
        ax2.text(col2, 0.80, f"L20: {stats.ln[20]:.1f} dB")
        ax2.text(col2, 0.65, f"L30: {stats.ln[30]:.1f} dB")
        ax2.text(col2, 0.50, f"L40: {stats.ln[40]:.1f} dB")
        ax2.text(col2, 0.35, f"L60: {stats.ln[60]:.1f} dB")
        ax2.text(col2, 0.20, f"L80: {stats.ln[80]:.1f} dB")
        
        # Dose
        std_name = dose_params.get('name', 'Standard')
        # Handle case where stats.dose might miss a key (robustness)
        dose_val = stats.dose.get('dose', 0.0)
        twa_val = stats.dose.get('twa', 0.0)
        
        ax2.text(col3, 0.80, f"Dose ({std_name})", fontweight='bold')
        ax2.text(col3, 0.65, f"Dose %: {dose_val:.1f}%")
        ax2.text(col3, 0.50, f"TWA: {twa_val:.1f} dB")

    @staticmethod
    def _plot_lp_history(fig, results, weighting, speed, autoscale, ymin, ymax):
        ax = fig.add_subplot(1, 1, 1)
        
        # Extract data safely
        t = [r.get('time', 0) for r in results]
        l = [r.get('lp', -300) for r in results] # Default to silence if key missing
        
        ax.plot(t, l)
        ax.set_title(f"Sound Pressure Level vs Time ({weighting}-weighted, {speed})")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Level (dB)")
        ax.grid(True)
        
        # Apply Scaling
        if not autoscale:
            ax.set_ylim(ymin, ymax)
        else:
            # If data is all -300 (silence), autoscale might look weird. 
            # Force reasonable range if max < 0
            if l and max(l) < 0:
                ax.set_ylim(-100, 100) # Fallback view
            else:
                ax.autoscale(axis='y')

    @staticmethod
    def _plot_spectrum(fig, results, weighting, is_third_octave, ref_pressure, 
                       autoscale, ymin, ymax):
        ax = fig.add_subplot(1, 1, 1)
        
        last_res = results[-1]
        freqs = last_res.get('band_freqs', [])
        
        # Average Spectrum Calculation
        if 'bands' in results[0]:
            energy_sums = np.zeros(len(freqs))
            count = 0
            for r in results:
                if 'bands' in r:
                    # Convert dB to Pressure^2
                    pressures = (10**(r['bands']/10.0)) * (ref_pressure**2)
                    energy_sums += pressures
                    count += 1
            
            if count > 0:
                mean_db = 10 * np.log10((energy_sums / count) / (ref_pressure**2) + 1e-30)
            else:
                mean_db = np.zeros(len(freqs))
        else:
            mean_db = []

        if len(freqs) > 0:
            x = np.arange(len(freqs))
            ax.bar(x, mean_db, color='#2ca02c', alpha=0.8)
            ax.set_xticks(x)
            
            lbls = [f"{f/1000:.0f}k" if f >= 1000 else f"{f:.0f}" for f in freqs]
            if is_third_octave: 
                lbls = [l if i % 3 == 0 else "" for i, l in enumerate(lbls)]
            ax.set_xticklabels(lbls, rotation=90)
            
        ax.set_title(f"Average Spectrum ({weighting}-weighted)")
        ax.set_ylabel("Level (dB)")
        ax.grid(axis='y')
        
        # Apply Scaling
        if not autoscale:
            ax.set_ylim(ymin, ymax)
        else:
            ax.autoscale(axis='y')