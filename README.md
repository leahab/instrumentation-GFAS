#  Codes and Routines for the Ground-based Fog and Aerosol Spectrometer (GFAS)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
## Overview

The **GFAS** (Ground-based Fog and Aerosol Spectrometer) is a newly developed instrument for in-situ single-particle measurements of cloud microphysical properties. Unlike conventional light scattering spectrometers that derive only particle size and number from forward-scattering intensity, the GFAS additionally measures the **change in polarization of backward scattered light**, enabling investigation of particle morphology and composition.

Particle sizing is based on Mie-Lorenz theory. A key methodological contribution of the GFAS is a novel **autocalibration procedure** using a polydisperse spray of pure water droplets coupled with theoretical Mie-Lorenz modeling. This approach yields lower uncertainties than the standard glass bead calibration and is applicable to the polarized backward scattering detectors as well.

The instrument has been validated in both laboratory and field settings. Under laboratory conditions, the GFAS can distinguish water droplets from non-spherical dust particles via polarization measurements. Field measurements further reveal that cloud particle optical properties differ measurably from pure water, suggesting an elevated refractive index in ambient cloud droplets.

This repository is used to collect code developed jointly by the scientific community, starting with the work presented in Haberstock et al. (2026).

> 📄 **Reference:** *Haberstock et al. (2026, to be submitted to AMT)*

---

## Repository Structure

```
GFAS/
│
├── GFAS_reading/               # Functions for reading raw GFAS data files
│   └── GFAS101_read_in             # GFAS (version 101) read in functions for PbP and 1hz data formats
│
├── calibrations/                   # Calibration routines and data processing
│   ├── MieCalculations             # subroutine to calculate theoretical scattering cross sections and plot Figure 2 and 3
│   ├── glass_bead                  # Glass bead calibration, plotting Figure A4, A5, S1,
│   ├── autocalibration_I           # Autocalibration method I
│   └── autocalibration_II_III      # Autocalibration method II and  III, plotting Figure  S4
│   
│
├── PbP_processing/                 # GFAS particle-by-particle (PbP) processing pipeline
│   ├── detector_scaling_factors    # The scaling factors needed to convert the measurement intensity into scattering cross section
│   ├── GFAS101_processing          # subroutines needed for all GFAS raw data to yield reliable results
│   ├── loss_Calc_functions_GFAS101 #subroutines needed for sample loss corrections
│   ├── PbP_level0_processing       # script to run PbP processing routine for PbP raw data to level 1 PbP data
│   └── PbP_to_1hz                  # script to run the PbP processing routine from level 1 to level 2 (binned) data
│
└── plots/                          # Figure generation scripts (Paper I by Haberstock et al., 2026)
│   ├── cloudevents                 # list of start and end times of the cloudevents during the measurement campaign
│   ├── modeled_SCS                 # models the theoretical scattering cross section for different refractive indices
│   ├── ParticleCoincidence         # calculates and plots probability of particles coinciding during measurements, Figure C1
│   ├── Polarization_FB_analysis    # Plotting Figures 6-9 
│   ├── SBO_plots                   # Plotting Figure 5, S5, S6
│   └──...                          # updates might follow
```

### Folder Descriptions

| Folder | Description |
|---|---|
| `GFAS_reading/` | Functions for ingesting and parsing raw GFAS instrument output |
| `calibrations/glass_bead/` | Standard glass bead calibration workflow |
| `calibrations/autocalibration_I/` | First-generation water-spray autocalibration |
| `calibrations/autocalibration_II/` | Second-generation autocalibration with refinements |
| `calibrations/autocalibration_III/` | Third-generation autocalibration, final methodology |
| `P2P_processing/` | End-to-end particle-to-particle data processing pipeline |
| `plots/` | Scripts to reproduce all figures in Paper I by Haberstock et al. (2026) |

---

## Instrument & Manufacturer

The GFAS (model **GFAS-DPOL**) is manufactured by [envea](https://envea.global/). The instrument combines high-sensitivity forward and backscatter particle size measurement with polarization detection and wind speed/direction measurement — designed for fog formation studies, number concentration, liquid water content, extinction coefficient and visibility retrieval in warm and mixed-phase fog.

🔗 **Full technical details:** [DMT GFAS-DPOL product page](https://envea.global/product/ground-based-fog-aerosol-spectrometer-with-polarization-detection/?_gl=1*1oeb8xm*_up*MQ..*_ga*ODg4MTEyNDI5LjE3ODI5ODc0MzQ.*_ga_5XMNKCHKE6*czE3ODI5ODc0MzQkbzEkZzEkdDE3ODI5ODc0MzQkajYwJGwwJGgxNDYwNzU1MTIy*_ga_1KNED64NVJ*czE3ODI5ODc0MzQkbzEkZzEkdDE3ODI5ODc0MzQkajYwJGwwJGgxMzk0MTA2NDQz)


## Getting Started

```bash
git clone https://github.com/SU-air/instrumentation-GFAS.git
cd instrumentation-GFAS
```

*required Python environment: GFAS101.yml*

---

## Citation

If you use the GFAS or this codebase in your work, please cite:

> *Haberstock et al. (2026), [reference will be added soon]*

---

## Contact

*Paul Zieger, Department of Environmental Science, Stockholm University, paul.zieger@aces.su.se*
