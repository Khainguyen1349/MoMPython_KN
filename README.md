# MoMPython_KN
Try to translate Matlab code from Makarov book into python and then play with antenna designs and analysis.

As the Jupiter notebooks have some dynamic figures, it is suggested to check the notebook with

https://nbviewer.jupyter.org/

______________________________________________________________________________________

* For example, one can copy the link of the RWG_Bowtie notebook and add to viewer to open

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_Bowtie.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/Bowtie.png" width=40% height=40%>

This is the earliest version of the simulator where one can still find the generation of the antenna mesh, calculation of the RWG basis functions, and the analysis of the antenna performances. Later on, those calculations are placed in libraries [meshlib.py](./meshlib.py) and [rwglib.py](./rwglib.py) for the sake of clarity.

______________________________________________________________________________________

* Simulation of a (Planar) Inverted-F Antenna

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_IFA.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/IFA.png" width=40% height=40%>

A second version with denser mesh is also presented in the project.

______________________________________________________________________________________


* Lumped components analysis on a Dipole antenna @75MHz:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_Dipole_75MHz_wLumpedElement.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/Dipole.png" width=40% height=40%>

The omnidirectional pattern of a typical dipole is found
 
<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/Dipole_75MHz_Pattern.png" width=40% height=40%>

______________________________________________________________________________________

* Directivity improvement with lumped components:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_Directivity_Optim_withLumped.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/DirOpt.png" width=40% height=40%>

By adding the lumped components, the impedance of the antenna is modified a bit but not to much. 

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/DirOpt_S11diff.png">

The radiation pattern of the antenna before adding the lumped components is

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/DirOpt_beforeAddingLumped.png" width=40% height=40%>

and after adding the components

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/DirOpt_afterAddingLumped.png" width=40% height=40%>

______________________________________________________________________________________


* Simulation for a Yagi antenna functioning at 1Ghz:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_Yagi_1GHz.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/Yagi.png" width=40% height=40%>

This antenna has more than 7dBi gain and its radiation pattern is as follow

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/Yagi_1GHz_Pattern.png" width=40% height=40%>

______________________________________________________________________________________


* Analysis of Charateristic Modes of a conductive plate:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_CharacteristicModes.ipynb

Surface currents of the first 4 Characteristic Modes

 <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM1.png" width=40% height=40%>   <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM2.png" width=40% height=40%>  
 <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM3.png" width=40% height=40%>   <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM4.png" width=40% height=40%>  

and its radiation patterns corresponding


 <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM1_pattern.png" width=40% height=40%>   <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM2_pattern.png" width=40% height=40%>  
 <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM3_pattern.png" width=40% height=40%>   <img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/CM4_pattern.png" width=40% height=40%>  

The surface currents and radiation pattern associated are very useful for certain types of design, e.g. Capacitive Coupling Element design for mobile terminals, etc. The generalize eigenvalues problem is solved thank to the work of Benyamin Ghojogh.

______________________________________________________________________________________


* Analysis of Quality-Factor of a conductive plate:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_QualityFactor_868MHz.ipynb

<img src="https://github.com/Khainguyen1349/MoMPython_KN/blob/main/figures/QF.png" width=40% height=40%>

______________________________________________________________________________________


* I have also tried to simulate a Meta-surface:

https://nbviewer.org/github/Khainguyen1349/MoMPython_KN/blob/main/RWG_Yagi_Reflector_1GHz.ipynb

but, well, I failed! It seems that MoM is not a good method to simulate a large sparse structures :D


Anyway have fun trying to simulate your antenna!
