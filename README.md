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
│   └── ...
│
├── calibrations/               # Calibration routines and data processing
│   ├── glass_bead/             # Glass bead calibration
│   ├── autocalibration_I/      # Autocalibration method I
│   ├── autocalibration_II/     # Autocalibration method II
│   └── autocalibration_III/    # Autocalibration method III
│
├── P2P_processing/             # GFAS particle-to-particle (P2P) processing pipeline
│   └── ...
│
└── plots/                      # Figure generation scripts (Paper I by Haberstock et al., 2026)
    └── ...
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

The GFAS (model **GFAS-DPOL**) is manufactured by [Droplet Measurement Technologies (DMT)](https://www.dropletmeasurement.com/). The instrument combines high-sensitivity forward and backscatter particle size measurement with polarization detection and wind speed/direction measurement — designed for fog formation studies, number concentration, liquid water content, extinction coefficient and visibility retrieval in warm and mixed-phase fog.

🔗 **Full technical details:** [DMT GFAS-DPOL product page](https://www.dropletmeasurement.com/product/ground-based-fog-aerosol-spectrometer-with-polarization-detection/)


## Getting Started

```bash
git clone https://github.com/SU-air/instrumentation-GFAS.git
cd instrumentation-GFAS
```

*Add installation/dependency instructions here (e.g., required MATLAB version, Python environment, etc.)*

---

## Citation

If you use the GFAS or this codebase in your work, please cite:

> *Haberstock et al. (2026), [reference will be added soon]*

---

## Contact

*Paul Zieger, Department of Environmental Science, Stockholm University, paul.zieger@aces.su.se*
