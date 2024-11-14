Optimal Shield Thickness Calculation for GM Detector

This project calculates the optimal thickness of radiation shielding for a GM (Geiger-Muller) detector. The program determines the best thickness and material based on efficiency, cost, and weight constraints, helping optimize the detector's protection from radioactive sources.

The code simulates photon interactions with materials using the Compton scattering model and attenuation cross-section data, aiming to minimize the radiation penetrating the shielding.
Features

    Materials Supported: Aluminum, Titanium, Iron, Lead
    Cost Optimization: Calculates the cost based on material price per kg
    Weight Constraints: Considers a maximum shield weight (currently 2 kg)
    Photon Energies Simulated: Cs-137, Technetium-99m, GM detector threshold
    Correction Factor for Experimental Data: Reads and applies correction factors from a CSV file

Installation
Prerequisites

Ensure you have Python 3 and the following libraries installed:

    numpy
    matplotlib
    periodictable
    xcom (for cross-section calculations)

Installation

    Clone this repository:

git clone https://github.com/tapravkerlc/Shield-Size-Calc
cd Shield-Size-Calc

Install dependencies:

    pip install numpy matplotlib periodictable xcom

    Ensure you have the correction factors CSV file in the project directory:
        zp1201EkeV.csv

Usage

    Configure Input Parameters:

    Adjust input parameters in the script:
        Material IDs: Set up in materiali list (e.g., [13, 22, 26, 82] for Al, Ti, Fe, Pb).
        Thickness Range: Adjustable via debel.
        Photon Energies: Specified in photon_energies for Cs-137, Tc-99m, etc.

    Run the Simulation:

    python main.py

    View Results:
        Optimal shield thickness and material are displayed.
        The program plots the shield's effectiveness (response) relative to energy levels.
        Two graphs display mass and cost variations by shield thickness.

Functions Overview

    correct_data(): Applies correction factors from a CSV file to measured data.
    compton(): Simulates energy and angle changes due to Compton scattering.
    min_linear_sublist(): Identifies the shield configuration with the highest linearity in response.
    plot(): Plots data showing the material’s attenuation effect.

Example Output

Upon running, the program outputs:

    Optimal shield thickness and cost.
    Plot of shield effectiveness against photon energies.
    Graphs of mass and cost growth with increasing shield thickness.

Future Improvements

    Extend support to additional materials.
    Add more photon sources for simulation.
    Improve performance for larger datasets.

License

MIT License.
