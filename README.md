# wiscobolt
A deterministic, finite-element photon-electron Boltzmann transport solver for problems in radiation therapy.


Contact the author at myounis@wisc.edu


Important documents are located in the "Documents" folder. License information is in the "LICENSE" folder.


To understand what wiscobolt aims to do, read the foreword of the **wiscobolt** document.

To understand how to format your input file, what external files you must generate (i.e., mesh, incident beam spectra), and wiscobolt's output, read the first chapter of the **wiscobolt implementation** document.


Requirements for wiscobolt (Windows only):

  &emsp; UnxUtils
  
  &emsp; gfortran 12.1.0 or later

At the moment, the provided makefile and batch file will not work in Linux.
  
In order to run wiscobolt, use the command prompt/shell to navigate to the main wiscobolt directory, and enter the command 'wiscobolt'.


Note finally that wiscobolt uses physics data extracted from a number of sources, as well as one other program. The appropriate citations, as well as licensing information, are available in the corresponding folders within the "Physics data" folder.
