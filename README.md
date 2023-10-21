# wiscobolt
A deterministic, finite-element photon-electron Boltzmann transport solver for problems in radiation therapy.  

Contact the author at myounis@wisc.edu  

Important documents are located in the "Documents" folder. Licensed with GNU General Public License 3.0. The full license can be found at LICENSE/gpl-3.0.txt.  

To understand what wiscobolt aims to do, read the foreword of the **wiscobolt** document.  
To understand how to format your input file, what external files you must generate (i.e., mesh, incident beam spectra), and wiscobolt's output, read the first chapter of the **wiscobolt implementation** document.  

Requirements for wiscobolt (Windows only):  
  &emsp; UnxUtils  
  &emsp; gfortran 12.1.0 or later  
  
In order to run wiscobolt, use the command prompt/shell to navigate to the main wiscobolt directory, and enter the command 'wiscobolt'.  
At the moment, the provided makefile and batch file will not work in Linux.  

Note finally that wiscobolt uses physics data extracted from a number of sources, as well as one other program. The appropriate citations, as well as licensing information, are available in the corresponding folders within the "Physics data" folder.  
