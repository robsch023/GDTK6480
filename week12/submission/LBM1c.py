import numpy as np
import scipy.fft as fft
import matplotlib.pyplot as plt
from git import Repo
import time

# Git revision extraction
timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
repo = Repo('.', search_parent_directories=True)
revsha = repo.head.object.hexsha[:8]
print(revsha)

U0, D, RE = 0.05, 10, 100
REF_FILE = 'week12/submission/pressure_history_centreline.csv'
OWN_FILE = 'week12/submission/pressure_history.csv'
START_REF = 5000
START_OWN = 5000


def determine_frequency(pressure, start, n_pad=2**18):
    pressure = pressure[start:]
    mean_pressure = pressure - np.mean(pressure)
    window = np.hanning(len(mean_pressure))
    windowed_signal = mean_pressure * window
    spectrum = fft.rfft(windowed_signal, n=n_pad)
    freqs = fft.rfftfreq(n_pad)
    spectrum[freqs < 1e-4] = 0  # Remove numerical noise
    freq_max = freqs[np.argmax(np.abs(spectrum))]
    return freq_max


def strouhal(f, U=U0, D=D):
    return f*D/U


def strouhal_empirical(Re):
    # Williamson (1988), valid for 49 < Re < 178
    return -3.3265/Re + 0.1816 + 1.6e-4*Re


def load_trace(fname):
    """Load a (time step, pressure) csv."""
    data = np.genfromtxt(fname, delimiter=',')
    data = data[~np.isnan(data).any(axis=1)]
    return data[:, 0], data[:, 1]


def plot_trace(steps, p, start, title, fname):
    """Full trace with the chosen start marked, plus a zoom on the retained part."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 6), constrained_layout=True)
    a1.plot(steps, p, lw=0.8)
    a1.axvline(start, color='r', ls='--', label=f'start = {start}')
    a1.set_ylabel('Pressure'); a1.set_title(title); a1.legend()
    a2.plot(steps[start:], p[start:], lw=0.8)
    a2.set_xlabel('Iteration'); a2.set_ylabel('Pressure')
    a2.set_title('Retained part')
    fig.text(0.02, 0.995, f"[rev {revsha}]", ha='left', va='top')
    fig.text(0.75, 0.995, timestamp, ha='left', va='top')
    fig.savefig(fname, dpi=200)


def start_sensitivity(p, starts):
    """St vs. choice of start"""
    return [strouhal(determine_frequency(p, s)) for s in starts]


# (i) verification
f_test = 8e-4
test_signal = np.sin(2*np.pi*f_test*np.arange(20000))
f_found = determine_frequency(test_signal, start=0)
print("(i) Verification")
print(f"    expected f = {f_test:.4e}, returned f = {f_found:.4e}, "
      f"error = {abs(f_found - f_test)/f_test*100:.3f} %\n")

# (ii) reference (centreline) trace
steps_ref, p_ref = load_trace(REF_FILE)
plot_trace(steps_ref, p_ref, START_REF, 'Reference trace (centreline probe)',
           'reference_trace.png')
f_ref = determine_frequency(p_ref, START_REF)
St_ref = strouhal(f_ref)
St_will = strouhal_empirical(RE)
print("(ii) Reference trace (centreline)")
print(f"    f = {f_ref:.4e}, St (as measured) = {St_ref:.4f}")
print(f"    Williamson St(Re={RE}) = {St_will:.4f}")
print(f"    Centreline probe sees two dips per cycle, so f_peak = 2f:")
print(f"    St_ref/2 = {St_ref/2:.4f}  "
      f"({(St_ref/2 - St_will)/St_will*100:+.2f} % vs Williamson)\n")

# (iii) own trace
steps_own, p_own = load_trace(OWN_FILE)
plot_trace(steps_own, p_own, START_OWN, 'Own trace (probe at 0.75L, 0.6H)',
           'own_trace.png')
f_own = determine_frequency(p_own, START_OWN)
St_own = strouhal(f_own)
print("(iii) Own trace")
print(f"    f = {f_own:.4e}, St = {St_own:.4f}")
print(f"    Williamson St(Re={RE}) = {St_will:.4f}, "
      f"difference = {(St_own - St_will)/St_will*100:+.2f} %\n")

# sensitivity of St to the start choice (own trace)
starts = np.arange(2000, 10001, 1000)
print("St vs start (own trace):")
for s, St in zip(starts, start_sensitivity(p_own, starts)):
    print(f"    start = {s:6d}  St = {St:.4f}")

# St-Re plot
Re_curve = np.linspace(49, 178, 200)
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(Re_curve, strouhal_empirical(Re_curve), 'k-', label='Williamson (1988)')
ax.plot(RE, St_own, 'ro', ms=8, label=f'LBM, Re = {RE} (St = {St_own:.4f})')
ax.set_xlabel('Re'); ax.set_ylabel('St')
ax.set_xlim(49, 178); ax.grid(alpha=0.3); ax.legend()
fig.text(0.02, 0.995, f"[rev {revsha}]", ha='left', va='top')
fig.text(0.75, 0.995, timestamp, ha='left', va='top')
fig.savefig('St_Re.png', dpi=200)
plt.show()