# DFE_Wright-Fisher
An individual-based Wright-Fisher model with infinite site genomes and customizable distribution of fitness effects

Parameters

All parameter of the model are defined in the parameter.h file. Parameter values can be changed in parameters.cpp before producing an executable. 
The default parameters are those used in the main text [ref].


How to compile DFE_Wright-Fisher

Simulations were conducted on a high performance computer (Linux). The executable was compiled in Release mode with the compiler gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-36).
To compile the code on a Linux system, the macro LINUX in PolyDisp.h must be set to CLUSTER 1. To compile on a Windows system, CLUSTER in DFE-Wright-Fisher.h must be set to CLUSTER 0.
