# MinimizationFunctions.py

import numpy as np
from scipy.optimize import minimize
from scqubits import *
import scqubits
scqubits.settings.NUM_CPUS = 4
scqubits.settings.OVERLAP_THRESHOLD  = 1.0 #  In order to look for the indices
n_oscillators = 3
Amp_per_phi0 = 0.1429013 
current_integer_flux = 0.02724955



def current_to_phiext(current_value):
    return (current_value - current_integer_flux) / Amp_per_phi0

def phiext_to_current(phiext):
    return phiext * Amp_per_phi0 + current_integer_flux


def from_parameter_array_to_parameter_dict(parameter_array, n_oscillators):
    keys_list = ['Ej', 'Ec', 'El'] + [f'g{i}' for i in range(1, n_oscillators + 1)] + [f'w{i}' for i in range(1, n_oscillators + 1)]
    parameters = {k: parameter_array[i] for i, k in enumerate(keys_list)}
    return parameters

def rescale_parameter_array(parameters_list):
    non_oscillator_params_length = 3  # Length of non-oscillator parameters
    oscillator_params_length = len(parameters_list) - non_oscillator_params_length
    scaling_factors = [1e-9] * non_oscillator_params_length + [1e-9] * oscillator_params_length
    return np.array(parameters_list * np.array(scaling_factors))

def unrescale_parameter_array(parameters_list):
    non_oscillator_params_length = 3  # Length of non-oscillator parameters
    oscillator_params_length = len(parameters_list) - non_oscillator_params_length
    scaling_factors = [1e9] * non_oscillator_params_length + [1e9] * oscillator_params_length
    return np.array(parameters_list * np.array(scaling_factors))



def plot_theory_lines(fig, parameters, transitions=["03","05","13","23"], **kwargs):
    # Extract the number of oscillators from the parameters
    n_oscillators = (len(parameters) - 3) // 2

    # Generate dynamic title
    title = ', '.join([f'{k} = {v/1e9:.3f} GHz' for k, v in parameters.items()])
    
    
    ### Get parameters from kwargs
    if 'flux_values' in kwargs:
        flux_values = kwargs.get('flux_values')
    else:
        flux_values = np.linspace(-1.0, 1.0, 101)
        
    if 'cutoff' in kwargs:
        cutoff = kwargs.get('cutoff')
    else:
        cutoff = 70
        
    ### calculate the lines
    lines_toplot = {}
    for t in transitions:
        lines_toplot[t] = []
    for f in flux_values:
        fluxonium = Fluxonium(EJ=parameters['Ej']/1e9, EC=parameters['Ec']/1e9, EL=parameters['El']/1e9, flux=f, cutoff=cutoff,)

        # Define oscillators and coupling strengths dynamically
        oscillator_params = {f'w{i}': parameters[f'w{i}'] for i in range(1, n_oscillators + 1)}
        g_strength_params = {f'g{i}': parameters[f'g{i}'] for i in range(1, n_oscillators + 1)}


        # Define oscillators
        oscillators = [Oscillator(oscillator_params[f'w{i}'] / 1e9, l_osc=1.0, truncated_dim=4) for i in range(1, n_oscillators + 1)]

        # Define coupling strengths
        g_strengths = [g_strength_params[f'g{i}'] / 1e9 for i in range(1, n_oscillators + 1)]


        hilbertspace = HilbertSpace([fluxonium, *oscillators])

 

        # Add interactions
        for osc, g_strength in zip(oscillators, g_strengths):
            hilbertspace.add_interaction(
            g_strength=g_strength,
            op1=osc.n_operator,
            # op1=osc.creation_operator,
            op2=fluxonium.n_operator,   
            add_hc=False,
            # add_hc=True
            )



        dressed_hamiltonian = hilbertspace.hamiltonian()
        eigenenergies = dressed_hamiltonian.eigenenergies(0)
        for t in transitions:
            lines_toplot[t].append(1e9*np.array([eigenenergies[int(t[1:])] - eigenenergies[int(t[0])]]))
        

    ### Plotting the lines
    ax = fig.axes[0]
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    current = phiext_to_current(flux_values)
    for t in transitions:
        ax.plot(current, lines_toplot[t], label=t)
        # first_part = str(hilbertspace.bare_index(int(t[0])))
        # second_part = str(hilbertspace.bare_index(int(t[1])))
        # transition_str = first_part + "->" + second_part
        # ax.plot(current, lines_toplot[t], label=transition_str)

    
    
    # ### Plotting the lines
    # ax = fig.axes[0]
    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # labelLines(plt.gca().get_lines(), align=True, backgroundcolor="none",fontsize=6)
    ax.legend(loc='upper right')
    ax.set_title(title)
    fig.canvas.draw_idle()
    
def plot_theory_energies(fig, parameters, levels=["0","1","2","3","4","5",], **kwargs):

    # Extract the number of oscillators from the parameters
    n_oscillators = (len(parameters) - 3) // 2

    # Generate dynamic title
    title = ', '.join([f'{k} = {v/1e9:.3f} GHz' for k, v in parameters.items()])
    
    ### Get parameters from kwargs
    if 'flux_values' in kwargs:
        flux_values = kwargs.get('flux_values')
    else:
        flux_values = np.linspace(-1.0, 1.0, 101)
        
    if 'cutoff' in kwargs:
        cutoff = kwargs.get('cutoff')
    else:
        cutoff = 110
        
    ### calculate the lines
    energies_toplot = {}
    for l in levels:
        energies_toplot[l] = []
    for f in flux_values:
        fluxonium = Fluxonium(EJ=parameters['Ej']/1e9, EC=parameters['Ec']/1e9, EL=parameters['El']/1e9, flux=f, cutoff=cutoff,)
        # Define oscillators and coupling strengths dynamically
        oscillator_params = {f'w{i}': parameters[f'w{i}'] for i in range(1, n_oscillators + 1)}
        g_strength_params = {f'g{i}': parameters[f'g{i}'] for i in range(1, n_oscillators + 1)}


        # Define oscillators
        oscillators = [Oscillator(oscillator_params[f'w{i}'] / 1e9, l_osc=1.0, truncated_dim=2) for i in range(1, n_oscillators + 1)]

        # Define coupling strengths
        g_strengths = [g_strength_params[f'g{i}'] / 1e9 for i in range(1, n_oscillators + 1)]

        hilbertspace = HilbertSpace([fluxonium, *oscillators])

        # Add interactions
        for osc, g_strength in zip(oscillators, g_strengths):
            hilbertspace.add_interaction(
            g_strength=g_strength,
            op1=osc.n_operator,
            # op1=osc.creation_operator,
            op2=fluxonium.n_operator,   
            add_hc=False,
            # add_hc=True
            )



        dressed_hamiltonian = hilbertspace.hamiltonian()
        eigenenergies = dressed_hamiltonian.eigenenergies(0)
        
        for l in levels:
            energies_toplot[l].append(1e9*np.array(eigenenergies[int(l)] ))
            
    hilbertspace.generate_lookup()
    ax = fig.gca() if fig.axes else fig.add_subplot(1, 1, 1)
    for l in levels:
        ax.plot(flux_values, energies_toplot[l], label = l)
        # ax.plot(flux_values, energies_toplot[l], label = str(hilbertspace.generate_lookup(l[0])))
        
    ax.set_xlabel(r'$\phi_{ext}$')
    ax.set_ylabel('Energy (GHz)')
    ax.set_title(title)
    #labelLines(plt.gca().get_lines(), align=True, backgroundcolor="none",fontsize=8)
    ax.legend()
    


 # Implement plot_extracted_points function
def plot_extracted_points(figure, extracted_points):
    ax = figure.axes[0]  # Use the figure argument instead of fig
    flux_values = list(extracted_points.keys())
    for f in flux_values:
        list_transitions = list(extracted_points[f].keys())
        for t in list_transitions:
            ax.plot(f, extracted_points[f][t], marker='x', markersize=8, color="red")
    figure.canvas.draw_idle()

def cost_function(parameters, extracted_points, model, cost_type='absolute', **kwargs):    
    cost = 0
    current_values = list(extracted_points.keys())
    for c in current_values:
        transitions = list(extracted_points[c].keys())
        data = np.array([extracted_points[c][t] for t in transitions])
        fit = model(c, parameters, transitions, **kwargs)*1e9
        if cost_type=='absolute':
            cost += np.sum(np.abs(data - fit))/1e9
        elif cost_type=='relative':
            cost += np.sum(np.abs(data - fit))/data
    return cost

def model(current, parameters, transitions, **kwargs):

    # Extract the number of oscillators from the parameters
    n_oscillators = (len(parameters) - 3) // 2

    
    if 'cutoff' in kwargs:
        cutoff = kwargs.get('cutoff')
    else:
        cutoff = 70
    fluxonium = Fluxonium(EJ=parameters['Ej']/1e9, EC=parameters['Ec']/1e9, EL=parameters['El']/1e9, flux=current_to_phiext(current), cutoff=cutoff,)

    
    # Define oscillators and coupling strengths dynamically
    oscillator_params = {f'w{i}': parameters[f'w{i}'] for i in range(1, n_oscillators + 1)}
    g_strength_params = {f'g{i}': parameters[f'g{i}'] for i in range(1, n_oscillators + 1)}


    # Define oscillators
    oscillators = [Oscillator(oscillator_params[f'w{i}'] / 1e9, l_osc=1.0, truncated_dim=4) for i in range(1, n_oscillators + 1)]

    # Define coupling strengths
    g_strengths = [g_strength_params[f'g{i}'] / 1e9 for i in range(1, n_oscillators + 1)]


    hilbertspace = HilbertSpace([fluxonium, *oscillators])

    # Add interactions
    for osc, g_strength in zip(oscillators, g_strengths):
        hilbertspace.add_interaction(
        g_strength=g_strength,
        op1=osc.n_operator,
        # op1=osc.creation_operator,
        op2=fluxonium.n_operator,
        add_hc=False
        # add_hc=True
        )
    
    
    dressed_hamiltonian = hilbertspace.hamiltonian()
    eigenenergies = dressed_hamiltonian.eigenenergies(0)
    return np.array([eigenenergies[int(t[1:])] - eigenenergies[int(t[0])] for t in transitions])









