# Monte Carlo simulations related to FRS experiment

Simulations and numerical calculations related to the future hypernuclear experiments with FRS at GSI.
Each directory presents different aspect of the simulations and calculations

## Contents:

### Simulations_frs directory:
  
  Phase space Monte Carlo simulations of the physics processes and decay of hypernuclei. 
  + MomSimu.C :
	ROOT script for the momentum setting split of the decay fragment aimed to be measured in the FRS.
	
  + SimuPhaseSpace.C :
	Study of the produced phase space of the decay particles from the hypernuclei decay. The Monte Carlo
	simulation is based on the generator of the phase space decay ported within ROOT framework. 
	A random walk is also use in the randomization. Momentum, resulting invariant masses and other kinematics 
	are then visualized.
	
	
	
### Mocadi directory:

  Gather the codes for MOCADI simulations (MC simulation program for ion transport through ion optical system).
  Several perl scripts have been developed in order to process a full systematics on the input parameters.
  
  * LoopMocadi.pl : Optimize transmission calculated in the Mocadi simulation for a fragment of interest. 
	- INPUT header\_file.in template\_body.in mocadi\_binary mean\_beam\_energy
      + header\_file.in : load and set the fragment of interest. It corresponds to the start of normal MOCADI input file until the FRAGMENT statement
      + template\_heb.in : rest of the input file of Mocadi 
      + (optional) The path to the mocadi binary. The perl script does not depends on the binary itself.
	  + (optional) Set a fix average beam energy that will set the central B$\rho$ of the magnets. If not set the script will automatically optimize the beam energy between each segment and then run next following segment with the optimized beam energy that maximizes the total transmission and efficiency. 
	
	- It loops over the different save states to optimize the efficiency and the transmission:
      + It creates intermediate input file of Mocadi, then runs it.
      + It read the output file with the results of the simulations (*.out) then it saves the mean energy of the fragments
      + It takes the energy of the fragments of interest and applies it to the template file up to the next save state
      + It runs again Mocadi with updated new header until arriving to the end of the geometry configuration of the magnetic spectrometer. 

  * SystematicMocadi.pl : Produce a systematic study of all fragments through all beam/target set available and several target thickness
    - INPUT header\_file\_template + (parasite option)
    - It loads, from ROOT "RadioNuclides.txt", the set of stable nuclei (c$\tau$ > 10 m) / name of elements
    - A set of target thickness is defined ex: (1, 5, 10, 15 cm)
    - It loads results from Epax calculation which will define the tuple [Beam,Target] I want to run
    - Then for each :
      + Fragment of interest / Beam / Target
      + It creates parasite fragments :
        - (A/Z * (Z-1))-1 $\to$ A
        - Z-1 $\to$ Z
        - check if stable enough otherwise A+1 $\to$ A
		- 3 possible parasites : Z-1 / Z-2 / Z-3
      + It creates the header file to run
      + It runs LoopMocadi.pl
      + Finally it creates in an output the table of results : \[thickness\]\[id\_frag\]\[save_state\] = transmission %
        - \*\_result.dat : Only fragment of interest at the save state : trans % \& cross section
        - \*\_result\_para.dat : Parasite at the last save state : cross section (not the proper need a correction processing)

  * MetaSystematicMocadi.pl : Run the systematics in parallels
    - Multicore processing (threads)
    - it takes the Epax file and split int in N subfile and runs SystematicMocadi.pl over the subfile in parallel

  * ScanningEnergy.pl : Study the efficiency dependency to the beam energy for the fragment of the hypernuclear decay.
	- In those case, Mocadi needs the phase space of the decay fragment from the $^{3}\_{\Lambda}$H, $^{4}\_{\Lambda}H, $^{3}\_{\Lambda}$n:
	  + A function responsible for loading the ROOT file in which GEANT4 simulation output of the precise physics simulations is saved.
	  + The function that is compiled as shared lib is injected inside the MOCADI simulation code and read the option line of the MOCADI input file for loading the ROOT file.
	  + Event by event the initial kinematics of fragment of interest is then set from the output of the GEANT4 simulation results.
    - The efficiency and transmission of the spectrometer is then obtained as a function of the B$\rho$ set through the full phase space of hypernuclear decay fragments. 
  * FuncLoadTree : 	The source of those functions to be injected in MOCADI.
	- FuncLoadTree.cc
	- FuncLoadTreeMC.cc
	- FuncLoadTreeMC_nnL.cc

  * AnaSysEnergy.C : C++ code of the analysis of the efficiencies as a function of the B$\rho$ of the FRS spectrometer.
   
### Epax :
   
  Gathers the codes for the analysis of fragmentation processes described by the Epax model.

  - epax_script.pl : Run the Epax model to obtain the production cross section of all fragment of interest from the beam/target collision and fragmentation mechanism.
    - INPUT A\_target\_max Z\_target\_max A\_projectile\_max Z\_projectile\_max
    - It loads stable elements to remove from the systematic loop the possible unstable nuclei for the target and the stable beam.
	- From proton on proton reaction to the beam A\_projectile\_max Z\_projectile\_max  on a target A\_target\_max Z\_target\_max reaction, the result of the Epax model is written in database of fragment / primary beam / target / production cross section (in mb)
	
  - Extract.C : C++ code (ROOT based) for analysis and data visualization.
  - SimplifyData.C : C++ code (ROOT based) for producing a dataset based on the output of the epax_script code that includes more information: The target density, the differential cross section per cm of target, the number of associated of fragment produced conjointly with the fragment of interest as a function of the cross section differences. 
  
### Ana\_UCLM :

  Analysis code for optimization procedure for finding the experimental optimal design for the future hypernuclear experiment within SuperFRS at FAIR. The production of proton- and neutron-rich hypernuclei can be maximized as a function of the primary beam / production target that will produce the secondary exotic beam that will collide with a carbon target to produce the hypernucleus of interest. Theoretical calculations from Epax and QGSM models, convoluted with the experimental aspects simulated by the MOCADI simulations, the experimental expectations was formulated. The code were used to obtain the results published in C. Rappold, _et al._, "_Examination of experimental conditions for the production of proton-rich and neutron-rich hypernuclei_", Phys. Rev. C __94__ (2016) 044616.
     
### fit TMVA :

  Gather the code based on TMVA framework of ROOT for supervised neutral network and boosted decision tree for modeling the nonlinear function of phase space modulation of the ion-optics of the magnetic spectrometer. Compared with the 3rd order tensorial Taylor expansion of beam trajectory along the magnetic axes, the base idea is to obtain a data driven nonlinear model for the tracking within the full magnetic elements of the spectrometers. Python implementation using _sklearn_ for more advance multivariate and machine learning algorithms is still private for incoming publication.
  
### Other scripts :

  Other C++ codes for data visualization and data analysis of initial parameters and conditions for the FRS experiments.

## Requirements

External: 

  * ROOT v6 ( > 6.08)
  * Epax model (vers. > 3.1)
  * MOCADI (> 3.5)
  * perl5 (> 5.12)
	
The perl packages needed for the perl scripts :

  - IO::CaptureOutput
  - threads
  - Config

Build for C++ code: gcc > 5.1 or clang > 3.8, they are using C++14 standard. 


## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Merge Request
