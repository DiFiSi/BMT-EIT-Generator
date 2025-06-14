[CURRENTLY UNDER RENOVATION FOR EVEN FASTER SIMULATION]

Fast generator of cardiac and respiratory mixed electrical impedance tomography (EIT) recordings first introduced in [Fast 4D FEM Model for EIT Source Separation Benchmarking](https://www.degruyterbrill.com/document/doi/10.1515/cdbme-2023-1097/html) - distinguished among best papers in [BMT2023](https://bmt2023.de/).  

This model was designed for computational speed and practicality. It could provide a reliable benchmark for evaluating cardiac and respiratory source separation algorithms as well as allow quick, realistic simulation of various physiological scenarios for developing and training more sophisticated algorithms (e.g., Deep Learning).

## Anatomical Model
The generator employs a simplified yet effective thoracic model represented as an ellipsoidal cylinder, with realistic geometric shapes (ellipsoids) for heart and lungs, enhanced by a 16-electrode system for electrical impedance tomography (EIT) with adjacent-adjacent pattern measurements.

<img src="https://github.com/user-attachments/assets/d7526f03-5b17-48ea-9619-026b2c59737d" width=50% height=50%>

## Modelling of Cardiac and Respiratory Signal Sources
Dynamic cardiac and respiratory signals are created using physiological templates:
* Cardiac: Ventricular and myocardial volume changes via realistic cardiac cycles.
* Respiratory: Templates for spontaneous and mechanical ventilation breathing patterns.
* Coupled Frequencies: Integrates physiological coupling between cardiac and respiratory rhythms for realism.

## Physiological Mixing of Conductivity Sources
Realistic physiological mixing is achieved by:
* Dynamically varying conductivity of lung tissue and blood based on volume and flow.
* Employing spatial distribution of pulmonary vessels with accurate blood conductivity modeling.
* Using Maxwell Garnett’s mixing formula for precise conductivity mixtures.

![sections_final](https://github.com/user-attachments/assets/7005bb90-ed88-4124-956e-1ef93a330f12)

## Folder Structure
Running simulations is possible through the script .\main.m, which calls the generator functions, and logs progress in .\simLog.txt.
* Generator functions are found in .\code (.\code\cond - conductivity mixing, .\code\fem - controlling FEM model, .\code\load - load simulation parameters, .\code\synth -  synthesize time dynamics, .\code\utils - auxiliary functions)
* Third-party auxiliaryy functions are found in .\others
* Simulation results are stored in .\results
* Simulation variables are declared in .\vars
