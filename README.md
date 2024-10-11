# tandem
Python application for automation of tandem aerosol classifier experiments and data inversion

Copyright (C) 2023-24 Cambustion, The University of Alberta, NRC Canada

Original Authors: Jonathan Symonds, Morteza Kiasadegh, Tim Sipkens

See LICENSE for details

## Environments
So far working ok under:  

WSL2 Ubuntu for Windows 11 (note serial ports not supported yet by WSL2)    
Anaconda / Conda for Windows 11  
MacOS Ventura  

## Requires

python 3.8+  
pyqt5 including designer  
pyserial  
numpy  

## Invocation

cd src/  
python tandem.py &  

## Usage

1. From file menu choose "Log raw data to..." and change where raw data file is saved (if required).    
2. For each classifier choose type, then ethernet or serial, then enter IP or serial port. Either use full name / path (e.g. COM1 or /dev/ttyS0) or just the number for Windows or Linux (1 for COM1 or /dev/ttyS1 etc).  
3. For each classifier choose sample flow, sheath flow or resolution, start and end sizes (or masses), number of classes per decade, and delay time between points.  
4. For the CPC, choose type, then ethernet or serial, then enter IP or serial port. Choose the number of seconds to average over.  
5. If you want to do scans with classifier 1 bypassed (1) before the first, and (2) after the last, regular scan then check the box in the second classifier config. Connect the bypass valve box to analogue output **3** of any AAC or CPMA involved.
6. Press scan, wait until finished.  
7. Press invert, TODO. 
8. If you click File->Save Settings (or check "Save Settings on Exit"), then all parameters will be reloaded next time.   

## Directories

*data/*   data files  
*doc/*  documentation  
*res/*  resources (e.g. images)  
*src/*  source files  

## Source files

*tandem.py*   main user interface  
*tandem.ui*   Qt designer file for the UI. To edit: designer tandem.ui &  
*instruments.py*  defines instrument comms and data (see below)  
*mplwidget.py*  needed to embed matplotlib plot in a Qt widget  
*inversion.py*  does the inversion (TODO)  

## Instrument class structure 

Instruments (i.e. their communications and data) are defined in a hierarchical class structure in src/instruments.py.

At the top level a general instrument which can have serial or ethernet comms is defined. This then splits into classifiers and CPCs. CPCs split into Cambustion CPC, TSI CPCs, Magic CPC and then as needed specific model ranges. Classifiers split into Cambustion Classifiers (CPMA or AAC), or TSI DMAs. Cambustion classifiers split into AACs or CPMAs, TSI DMAs into 3080 or 3082. Common features should be defined as high as possible in the hierarchy. When a function is called, the lowest one in the tree for that object will be the one used. 

![Instrument class structure](/doc/instrument%20classes.png)
