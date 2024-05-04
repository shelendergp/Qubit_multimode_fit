# Qubit_multimode_fit
This repository calculates the values of Ej, Ec, El, g, and w for the fluxonium coupled to the multimode cavity from the two-tone data. To use this repository, one has to work with the MainScript.py file and specify the following parameters: Ej, Ec, El, w1, w2, ..., wn, g1, g2, ..., gn, integer_flux_current_per_phi0, and the number of oscillators. With this repository, you can:

(i) Plot the 2D data
(ii) Observe the eigenenergies
(iii) Fit the two-tone spectroscopy data to find the values of Ej, Ec, El, wn, and gn

Installation:

1) Download the source code from GitHub by clicking the Code button on the top right, then unzip the source code folder in a <directory>.
OR
2)  Alternatively, open a terminal, navigate to <directory> where you would like to store the source code of Qubit_multimode_fit, and execute the following commands
git clone https://github.com/shelendergp/Qubit_multimode_fit.git
3)(Optional but recommended) Create a virtual Python environment (version > 3.10). If you are using a Conda terminal, you can create and activate the environment using the following commands:
conda create -n <env_name> python=3.10
conda activate <env_name>
4) Install the required packages:
pip install scqubits
pip install matplotlib-label-lines


