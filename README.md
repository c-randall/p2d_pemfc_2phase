# p2d_pemfc_2phase
Pseudo-2D half-cell model of a PEMFC cathode that includes two-phase transport.

[![DOI](https://zenodo.org/badge/444956566.svg)](https://zenodo.org/badge/latestdoi/444956566)

## Objective(s)
To increase performance and durability metrics in low-cost proton exchange
membrane fuel cells (PEMFCs), an improved understanding of the complex physiochemical
processes within the cathode is needed. Therefore, this model sets out
to identify limiting transport phenomena (i.e. protons, water, or oxygen) within
PEMFCs fabricated with low Pt loadings. The model discretizes both the gas difussion layer
(GDL) and catalyst layer (CL) and solves equations that govern conservation of mass and 
charge. These conservation equations include two-phase (i.e. gas and liquid) transport of
molecular species, ionic transport, and reversible chemical and charge transfer reactions.

Although other groups have modeled similar systems, this model is unique for its inclusion
of structure-property relationships. In other models, transport properties for thin-film Nafion 
are implemented using empirical formulas derived from bulk membranes. These relationships 
fail to capture the more restrictive transport of thin-films caused by their reduced water uptake 
and confinement effects near substrate interfaces. Therefore, we have instead derived
more detailed descriptions for oxygen and proton transport using structural data from neutron 
reflectometry experiments and separate conductivity measurements, both performed using thin-fim
Nafion.

After identifying the largest contributors to performance losses, the model is also used to 
investigate optimal CL microstructures through a parametric study. Some variations in this 
study are simple, such as changing the CL thickness or ionomer loading. In addition, more unique 
designs are included, such as functionally graded CLs. In a graded CL, Pt and/or ionomer loadings
are concentrated to be higher near the GDL and lower near the membrane. 

## Modeling Domains
Figure 1 gives a representation of the included species, physics, and geometry used in the model. 
Most state variables are only tracked in the depth direction of the GDL and CL. Within the thin-film 
Nafion shells found in the CL, state variables are also tracked radially. 

<p align="center"> <img src="https://user-images.githubusercontent.com/39809042/148299226-82796eb7-8c15-4267-89d5-437c974a4b0b.png"> </p>
<p align="center"> Figure 1: Illustration of model species, processes, and geometry. </p>

Note that although the model only discretizes and solves equations for the cathode, the membrane
is modeled as a simple resistor to account for the loss in potential between the anode and cathode. 
The model also outputs a voltage that accounts for the full-cell potential by shifting the solution
according to the equilibrium potential of the anode half-cell reaction.

## Simulation Methods
This model uses a finite volume method in order to conserve mass within the system. 
Although performance at steady state conditions are the output from this model, an 
initial value problem ODE integrator is used along with transient differential 
equations until steady state is reached. Due to this solution method, it is important 
that the user set the "t_sim" variable within "pemfc_runner.py" to a sufficiently 
long time. In order to check that steady state is being reached, the "debug" input 
(also within "pemfc_runner.py") can be set to the value of `1` in order to produce plots 
of the variables within the solution vector against time. If steady state is being 
reached, then each of these plots should reach a constant value by the end of the 
simulation time set by the user.

## Installation Instructions
1. Install [Anaconda](https://www.anaconda.com/distribution/) - make sure to get 
Python 3 syntax.
2. Launch "Anaconda Prompt" once the installation has finished.
3. Type `conda create --name echem --channel cantera/label/dev cantera numpy scipy pandas matplotlib` 
into the terminal of "Anaconda Prompt" to set up an environment named "echem" with the 
needed packages.
4. When prompted, type `y` and press enter to finish setting up the environment. 
Agree to any required pop-up messages.
5. Test the new environment by typing `activate echem` followed by the enter key.
6. Install an editor for Python files. A good option is [Atom](https://atom.io/).
6. Download all of the files from this repository onto your local machine.
7. Follow the operating instructions below to edit and run the model.

## Operating Instructions
1. Open "Anaconda Prompt" and type `activate echem` followed by the enter key.
2. Use `cd` to change into the directory where all of the repository files were 
downloaded to.
3. Once inside the correct directory, run the model by typing in `python pemfc_runner.py` 
and pressing enter.
4. To edit any of the model inputs or options, open the "pemfc_runner.py" file in any 
Python editor (e.g. Atom).
5. After making any desired changes to "pemfc_runner.py", save the file and repeat 
steps 1-3 to rerun the model.

Optional: If you would prefer to use a developer environment (sort of like Matlab) 
instead of the "Anaconda Prompt" terminal, then do the following: open "Anaconda Navigator", 
select "echem" from the dropdown menu labeled "Applications on" near the top of the page, 
and install "spyder" from the tiles on the screen. Once Spyder is installed, the 
"pemfc_runner.py" file can be opened within the program where it can be both edited and 
run without the need for a separate editor and terminal. For more details visit Spyder's 
website [here](https://www.spyder-ide.org/).

## License
This tool is released under the BSD-3 clause license, see LICENSE for details.

## Citing the Model
This model is versioned using Zenodo:

If you use this tool as part of a scholarly work, please cite using:

> C.R. Randall and S.C. DeCaluwe. (2021) P2D PEMFC Model (2 Phase) v1.0 [software]. Zenodo.

A BibTeX entry for LaTeX users is

```TeX
@article{Randall2022,
  title = {Predicted Impacts of Pt and Ionomer Distributions on Low-Pt-Loaded PEMFC Performance},
  volume = {169},
  DOI = {10.1149/1945-7111/ac8cb4},
  number = {9},
  journal = {Journal of The Electrochemical Society},
  publisher = {The Electrochemical Society},
  author = {Corey R. Randall and Steven C. DeCaluwe},
  year = {2022},
  pages = {094512}
}

@misc{P2D-PEMFC-2Phase,
    author = {Corey R. Randall and Steven C. DeCaluwe},
    year = 2021,
    title = {P2D PEMFC Model(2 Phase) v1.0},
    doi = {10.5281/zenodo.5823445},
    url = {https://github.com/c-randall/p2d_pemfc_2phase},
}
```

In both cases, please update the entry with the version used. The DOI for the latest 
version is given in the badge at the top, or alternately <https://doi.org/10.5281/zenodo.5823445> will
take you to the latest version (and generally represents all versions).
