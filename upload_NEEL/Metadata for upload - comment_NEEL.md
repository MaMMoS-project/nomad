# Comment:
The data was created within the project MaMMoS 101135546 (HORIZON-CL4-2023-DIGITAL-EMERGING-01).
Disclaimer: Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or European Health and Digital Executive Agency (HADEA). Neither the European Union nor the granting authority can be held responsible for them.


# References:
https://mammos-project.github.io/
https://neel.cnrs.fr/en/equipes-poles-et-services/micro-et-nanomagnetisme-mnm




################## yaml ###########################
description: "Data from one sample on a compositionally graded NdFeCeB film deposited on a Si wafer at CNRS, within the framework of the EU project MaMMoS.\n The film was sputtered from a 2-part target (1/2 NdFeB + 1/2 CeFeB), at a target-substrate distance of 100 mm and RF power 50 W (film name: NdCeFeB_2-5). The extracted data includes the x- and y- position of the sample on the wafer, the chemical elements measured by EDX and the coercivity measured by MOKE.\n EDX measurements only allow for the quantification of Nd, Ce, and Fe. Therefore the B concentration was assigned to impose an Fe/B ratio of 14, while maintaining the measured ratio between Nd, Ce and Fe, so that Fe+Nd+Ce+B=100%. conc\n . MOKE measurements are performed under pulsed magnetic field. The value of coercivity is extracted from MOKE measurements when the magnetization passes through zero in the reversal curve. MOKE measurements may not have been performed for all samples, so the value of coercivity is missing for some samples.\n\nA juypter-notebook for showing how to handle the data can be found next to the raw data.\n
Additional information like CHADA diagrams can be found on https://mammos-project.github.io/resources.html .
"

## MODA-CHADA
?? Link to CHADA diagram on github?

##  Available fields in yaml
  lab_id: 'CNRS Institut Néel'
  institute: 'CNRS Institut Neel'
  short_name: 'NEEL-Sample--30.0_40.0'
  owner: 'Nora Dempsey, Pierre le Berre, William Rigaut'
  
  # possible fields from the .hdf5 file	
  # fabrication_date: '13/5/2025'
  # operator: 'Pierre'
  # sample_name: NOMAD_test2


