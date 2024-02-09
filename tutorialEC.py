'''
The code was created by Angelo Orletti Del Rey for studing Ephaptic Coupling in scientific initiation project at Federal University of ABC (UFABC)

This program creates neurons with one soma with hodking & huxley mechanism and one dendrite with passive mechanism and every cell compartment
have the extracellular mechanism that models extracellular space. This space is represented as a resistor mesh (ohmic medium), the border cells
are the only connected with the ground potential. This connection is built by considering that each resistor have fixed value, the resistors 
that do not connect with any cell are connected to the ground potential. The cells are positioned parallel to each other in a square mesh.
'''
import ephCoupling as ec
'''
The first thing to configure the simulation is to determine the number of cells and the distance between them. 
'''
ec.ncell = 10 #number of cells
ec.cell
