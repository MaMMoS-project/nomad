# Code used to create a upload on Nomad for data from UU

Into the folder called 'dataSets' one should put the simulation dataset to be processed.
 a specific dataset. Each dataset subfolder should have the following structure:

- A `.cif` file containing the crystal structure information.
- A `.README` file with details about the dataset and the expected results.
- A `GS` folder, which contains files related to ground state calculations (e.g., output files).
- A `JiJ` folder, which includes files for exchange interaction calculations (e.g., Jij calculation outputs).
- An `MC` folder, which holds files for the Monte Carlo simulations.

This organization ensures that all relevant data and documentation for each chemical composition are kept together and easy to locate.
