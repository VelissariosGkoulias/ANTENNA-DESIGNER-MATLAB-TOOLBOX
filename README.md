# RCS Calculation of Custom Geometries using MATLAB Antenna Designer Toolbox

## Overview

This repository contains a MATLAB script (`rcs_analysis.m`) that utilizes the Antenna Designer Toolbox to calculate the Radar Cross Section (RCS) of custom 3D geometries imported from STL files.

## Features

* **STL File Import:** Reads 3D geometry data from standard STL (StereoLithography) files.
* **Custom Antenna Geometry:** Creates custom antenna geometry objects within the MATLAB Antenna Toolbox based on the imported STL data.
* **Monostatic RCS Calculation:** Calculates the monostatic RCS (radar cross section when the transmitter and receiver are at the same location) using the Method of Moments (MoM), also known as Boundary Element Method (BEM), solver.
* **Parameter Configuration:** Allows users to easily define simulation parameters such as operating frequency, incident wave angles (phi and theta), and polarization.
* **Batch Processing:** The script can process multiple STL files in a loop, allowing for comparative analysis of different geometries.
* **Basic Visualization:** Includes an optional section to visualize the imported 3D geometries.
* **RCS Value Output:** Displays the calculated RCS values in dBsm (decibels relative to one square meter) for each processed geometry.
* **Guidance for Further Analysis:** Provides suggestions for extending the analysis, such as varying incident angles and frequency, exploring bistatic RCS, and considering stealth design principles.

## Getting Started

### Prerequisites

* **MATLAB:** You need to have MATLAB installed.
* **Antenna Toolbox:** The MATLAB Antenna Toolbox is required to use the functions for creating custom geometries and calculating RCS.
* **STereoLithography (STL) File Format Support Package:** This package is usually included with MATLAB or can be installed through the Add-Ons Explorer. It provides the `stlread` function used for importing STL files.

### Installation

1.  **Clone the Repository:**
    ```bash
    git clone https://github.com/VelissariosGkoulias/ANTENNA-DESIGNER-MATLAB-TOOLBOX.git
    cd rcs-calculation-matlab
    ```

2.  **Place STL Files:** Ensure that the STL files of the 3D geometries you want to analyze are located in the same directory as the `rcs_analysis.m` script, or update the file paths in the `stlFiles` cell array within the script.

### Usage

1.  **Open MATLAB:** Launch your MATLAB environment.
2.  **Navigate to the Project Directory:** In the MATLAB command window, navigate to the directory where you saved the `rcs_analysis.m` file.
3.  **Edit Simulation Parameters (Optional):** Open the `rcs_analysis.m` script and modify the following parameters in the "**1. Define Simulation Parameters**" section according to your needs:
    * `frequency`: The operating frequency of the radar (in Hz).
    * `incidentAnglePhi`: The incident angle in phi (azimuth, in degrees).
    * `incidentAngleTheta`: The incident angle in theta (elevation, in degrees).
    * `polarization`: The polarization of the incident wave ('HH', 'VV', 'HV', 'VH').
4.  **Specify STL File Paths:** In the "**2. Specify STL File Paths**" section, update the `stlFiles` cell array with the names (or full paths) of your STL files.
5.  **Run the Script:** Execute the script by typing `rcs_analysis` in the MATLAB command window and pressing Enter.

### Output

The script will output the calculated monostatic RCS values (in dBsm) for each of the specified STL files to the MATLAB command window. If enabled, it will also display a figure showing the 3D geometries.

## Further Analysis and Contributions

The script provides a basic framework. For a more in-depth analysis, consider:

* **Varying Incident Angles:** Modify the script to loop through different `incidentAnglePhi` and `incidentAngleTheta` values to study the RCS across various viewing angles.
* **Frequency Sweeps:** Change the `frequency` parameter to observe how RCS changes with different radar frequencies.
* **Bistatic RCS:** Explore the `bistaticrcs` function in the Antenna Toolbox to calculate RCS when the transmitter and receiver are at different locations.
* **Advanced Meshing:** For more accurate results, especially at higher frequencies or for complex geometries, investigate mesh refinement techniques within the Antenna Toolbox (though the solver often handles this automatically).

Contributions to this project are welcome. If you have suggestions for improvements, bug fixes, or new features, please feel free to open an issue or submit a pull request.

## Authors

* **Velissarios Gkoulias** - *Initial work* - (https://github.com/VelissariosGkoulias)

## Extra Material

* https://www.mathworks.com/help/antenna/ug/radar-cross-section-benchmarking.html
* https://github.com/UTAustinCEMGroup/AustinCEMBenchmarks/tree/master/Austin-RCS-Benchmarks/Problem%20III-Almonds/Problem%20Set%20IIIB-Solid%20Resin%20Almonds
* https://github.com/UTAustinCEMGroup/AustinCEMBenchmarks/blob/ba9f243c3210eb3cda78f1bd18971c9c576c53ae/Austin-BioEM-Benchmarks/LICENSE.txt