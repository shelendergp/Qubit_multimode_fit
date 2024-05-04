# MinimizationFunctions.py

import numpy as np
import matplotlib.pyplot as plt
from labmate.acquisition_notebook import AcquisitionAnalysisManager

fig = plt.figure(num = 1) # empty dataset



def remove_lines(figure):
    ax = fig.axes[0]
    try:
        ax.get_legend().remove()
    except:
        pass
    for i,l in enumerate(ax.lines):
        l.remove()
    fig.canvas.draw_idle()


def process_data(data):
    flux_list = data.flux
    frequency = data.frequencies
    I = data.I
    Q = data.Q

    flattened_I = I.flatten()
    flattened_Q = Q.flatten()

    # Convert results into Volts and normalize
    S = (flattened_I + 1j * flattened_Q)
    R = np.abs(S)  # Amplitude
    R_reshaped = R.reshape((len(flux_list), -1))
    mag = R.reshape((len(flux_list), -1))
    phase = np.angle(S)  # Phase
    phase_reshaped = phase.reshape((len(flux_list), -1))

    # Normalize data
    row_sums_R = R_reshaped.sum(axis=0)
    R_reshaped /= row_sums_R[np.newaxis, :]
    mag /= row_sums_R[np.newaxis, :]

    row_sums_phase = phase_reshaped.sum(axis=0)
    phase_reshaped /= row_sums_phase[np.newaxis, :]

    return flux_list, frequency, R_reshaped, mag, phase_reshaped


# Define a function to load and process data
def load_and_process_data(data_directory, file_paths):
    data_list = [AcquisitionAnalysisManager(data_directory).load_file(file_path) for file_path in file_paths]
    return [process_data(data) for data in data_list]



# Define plot_mag function
def plot_mag(mag_rescaled, flux_list, frequency, alpha=1.0, vmin=None, vmax=None):
    plt.pcolormesh(flux_list, frequency, mag_rescaled.T, shading='auto', cmap='viridis', rasterized=True, alpha=alpha, vmin=vmin, vmax=vmax)


# Define plot_phase function
def plot_phase(phase_rescaled, flux_list, frequency, alpha=1.0, vmin=None, vmax=None):
    plt.pcolormesh(flux_list, frequency, phase_rescaled.T, shading='auto', cmap='viridis', rasterized=True, alpha=alpha, vmin=vmin, vmax=vmax)


# Define rescale_data function
def rescale_data(data, rescale=True):
    if rescale:
        return np.array([m - np.mean(m) for m in data])
    else:
        return data


# Define main plotting function
def plot_data(processed_data, alphas, vmin_values, vmax_values, auto_files):
    for i, (flux_list, frequency, R_reshaped, mag, phase_reshaped) in enumerate(processed_data):
        mag_rescaled = rescale_data(mag)
        phase_rescaled = rescale_data(phase_reshaped)
        if i in auto_files:
            # Calculate vmin and vmax automatically based on the data
            vmin_auto = np.min(mag_rescaled)
            vmax_auto = np.max(mag_rescaled)
            plot_mag(phase_rescaled, flux_list, frequency, alpha=alphas[i], vmin=vmin_auto, vmax=vmax_auto)
        else:
            plot_mag(phase_rescaled, flux_list, frequency, alpha=alphas[i], vmin=vmin_values[i], vmax=vmax_values[i])

    cbar = plt.colorbar()
    cbar.set_label("Magnitude (a.u)")
    plt.xlabel('Voltage (V)')
    plt.ylabel('Frequency (Hz)')
    plt.title(r"$R=\sqrt{I^2 + Q^2}$ (normalized_Averaged)")
    plt.tight_layout()
    plt.show()


# Main script
if __name__ == "__main__":
    # Your main script logic here
    pass
