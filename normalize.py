import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import minimum_filter, maximum_filter

# ---------------------------------------------------------------------------
# Configuration & Balmer Masking
# ---------------------------------------------------------------------------
BALMER_CATALOG = [
    (6562.81, 95.0), (4861.34, 45.0), (4340.47, 35.0),
    (4101.74, 30.0), (3970.07, 25.0), (3889.05, 22.0),
    (3835.39, 20.0), (3797.90, 18.0), (3770.63, 15.0),
]

def balmer_mask(wave):
    """Generates a boolean mask covering the specified Balmer line regions."""
    mask = np.zeros(len(wave), dtype=bool)
    for center, hw in BALMER_CATALOG:
        mask |= (wave >= center - hw) & (wave <= center + hw)
    return mask
    
# ---------------------------------------------------------------------------
# selected line list for LSD
# ---------------------------------------------------------------------------
lines = np.loadtxt('sirius.mask.10000.40.p00.correct.clean', skiprows=1)
wavelengths = lines[:, 0]

# ---------------------------------------------------------------------------
# Data Loading
# ---------------------------------------------------------------------------
def load_orders(filepath, min_points=200):
    """Parses spectrum file, segmenting by wavelength gaps/order jumps."""
    orders = []
    w_buf, f_buf = [], []
    prev_w = None
    with open(filepath) as fh:
        next(fh); next(fh)
        for line in fh:
            tokens = line.split()
            if len(tokens) < 3 or 'nan' in tokens[1]: continue
            w, f = float(tokens[0]), float(tokens[1])
            if prev_w is not None and (w < prev_w or w - prev_w > 10.0):
                if len(w_buf) >= min_points:
                    orders.append({'wave': np.array(w_buf), 'flux': np.array(f_buf)})
                w_buf, f_buf = [], []
            w_buf.append(w); f_buf.append(f)
            prev_w = w
    if len(w_buf) >= min_points:
        orders.append({'wave': np.array(w_buf), 'flux': np.array(f_buf)})
    return orders

# ---------------------------------------------------------------------------
# Continuum Fitting Pipeline
# ---------------------------------------------------------------------------
def fit_continuum_pipeline(wave, flux, bm, window_width=150, sigma=2.0):
    """
    1. Rolling Ball: Establishes blaze trend.
    2. Windowed Sigma Clipping: Locally rejects absorption/emission lines.
    3. Anchor Pruning: Final check to remove anchors sitting in line profiles.
    """
    telluric_mask = (wave >= 7580.0) & (wave <= 7720.0)
    exclusion_mask = bm | telluric_mask
    mask = ~exclusion_mask
    
    # 1. Rolling Ball blaze trend
    ball_radius = window_width
    blaze_trend = maximum_filter(minimum_filter(flux, size=ball_radius), size=ball_radius)
    norm_flux = flux / blaze_trend
    
    # 2. Windowed Sigma Clipping
    clean_mask = mask.copy()
    for i in range(0, len(wave), window_width // 2):
        win_slice = slice(i, i + window_width)
        local_f = norm_flux[win_slice]
        if len(local_f) == 0: continue
        is_line = np.abs(local_f - np.median(local_f)) > (sigma * np.std(local_f))
        clean_mask[win_slice] &= ~is_line

    # 3. Anchor Pruning
    temp_spline = UnivariateSpline(wave[clean_mask], flux[clean_mask], s=len(wave)*0.05)
    resid = flux[clean_mask] - temp_spline(wave[clean_mask])
    final_mask = clean_mask.copy()
    final_mask[clean_mask] &= (resid > -0.5 * np.std(resid))
    
    final_spline = UnivariateSpline(wave[final_mask], flux[final_mask], s=len(wave)*0.05)
    return final_spline, wave[final_mask], flux[final_mask]

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def plot_results(orders):
    fig = plt.figure(figsize=(14, 12))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1], hspace=0.1)
    ax_raw = fig.add_subplot(gs[0])
    ax_norm = fig.add_subplot(gs[1], sharex=ax_raw)
    ax_resid = fig.add_subplot(gs[2], sharex=ax_raw)

    all_flux = []
    edge_delta = 5.0  # Size of the side zones in Angstroms

    for o in orders:
        all_flux.extend(o['flux'])
        bm = balmer_mask(o['wave'])
        
        cont_func, cw, cf = fit_continuum_pipeline(o['wave'], o['flux'], bm)
        continuum = cont_func(o['wave'])
        
        effective_continuum = continuum.copy()
        
        for center, hw in BALMER_CATALOG:
            w_left_inner = center - hw
            w_right_inner = center + hw
            
            interval_mask = (o['wave'] >= w_left_inner) & (o['wave'] <= w_right_inner)
            if not np.any(interval_mask):
                continue
            
            # Define small side zones just outside the Balmer region, ensuring they don't fall into another Balmer line
            left_zone_mask = (o['wave'] >= w_left_inner - edge_delta) & (o['wave'] < w_left_inner) & (~balmer_mask(o['wave']))
            right_zone_mask = (o['wave'] > w_right_inner) & (o['wave'] <= w_right_inner + edge_delta) & (~balmer_mask(o['wave']))
            
            # Fallbacks if clean unmasked points aren't found in the immediate window
            if not np.any(left_zone_mask):
                left_zone_mask = (o['wave'] >= w_left_inner) & (o['wave'] <= w_left_inner + 2.0)
            if not np.any(right_zone_mask):
                right_zone_mask = (o['wave'] >= w_right_inner - 2.0) & (o['wave'] <= w_right_inner)
                
            c_left = np.median(o['flux'][left_zone_mask]) if np.any(left_zone_mask) else cont_func(w_left_inner)
            c_right = np.median(o['flux'][right_zone_mask]) if np.any(right_zone_mask) else cont_func(w_right_inner)
            
            w_left_pos = np.median(o['wave'][left_zone_mask]) if np.any(left_zone_mask) else w_left_inner
            w_right_pos = np.median(o['wave'][right_zone_mask]) if np.any(right_zone_mask) else w_right_inner
            
            w_interval = o['wave'][interval_mask]
            if w_right_pos != w_left_pos:
                line_continuum = c_left + (w_interval - w_left_pos) * (c_right - c_left) / (w_right_pos - w_left_pos)
            else:
                line_continuum = np.full_like(w_interval, c_left)
                
            effective_continuum[interval_mask] = line_continuum
        
        normalized_flux = o['flux'] / effective_continuum
        
        # Plotting raw flux and identified continuum anchors
        ax_raw.plot(o['wave'], o['flux'], color='blue', lw=0.8, alpha=0.3)
        ax_raw.plot(o['wave'], continuum, color='red', lw=1.0)
        ax_raw.scatter(cw, cf, color='red', marker='x', s=5, alpha=0.5)
        
        # Plot normalized flux with corrected ylim
        ax_norm.plot(o['wave'], normalized_flux, color='black', lw=0.8)
        ax_norm.vlines(wavelengths*10, ymin=0, ymax=1, color='red')
        
        # Residuals
        ax_resid.scatter(cw, cf - cont_func(cw), color='green', s=2, alpha=0.5)

    ax_raw.set_ylim(min(all_flux), max(all_flux))
    ax_norm.set_ylim(0, 1.1)
    
    for ax in [ax_raw, ax_norm, ax_resid]:
        for center, hw in BALMER_CATALOG:
            ax.axvspan(center-hw, center+hw, color='yellow', alpha=0.05)
            
    ax_raw.set_title('Raw Flux with Pruned Windowed Clipping (758-772nm excluded)')
    ax_norm.set_title('Normalized Flux (Safe Side-Zone Linear Pseudo-Continuum)')
    ax_resid.set_title('Residuals (Continuum Anchors Only)')
    ax_resid.set_xlabel('Wavelength (Å)')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    orders = load_orders('new.s_ext')
    plot_results(orders)
