Fast generator of cardiac and respiratory mixed electrical impedance tomography (EIT) recordings first introduced in [Fast 4D FEM Model for EIT Source Separation Benchmarking](https://www.degruyterbrill.com/document/doi/10.1515/cdbme-2023-1097/html) - distinguished among best papers in [BMT2023](https://bmt2023.de/).  

This model was designed for computational speed and practicality. It could provide a reliable benchmark for evaluating cardiac and respiratory source separation algorithms as well as allow quick, realistic simulation of various physiological scenarios for developing and training more sophisticated algorithms (e.g., Deep Learning).

## Anatomical Model
The generator employs a simplified yet effective thoracic model represented as an ellipsoidal cylinder, with realistic geometric shapes (ellipsoids) for heart and lungs, enhanced by a 16-electrode system for electrical impedance tomography (EIT) with adjacent-adjacent pattern measurements.

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
