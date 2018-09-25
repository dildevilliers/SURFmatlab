# /fieldFuncs
A variety of functions that manipulates the standard (far)field structure [FFstruct]

Function list:
    example1.m: Example function that does nothing and does not exist

FFstruct contains:
  At least
    .th - [Na x Nf] matrix of theta angles (in rad) between [0,th_max]
    .ph - [Na x Nf] matrix of phi angles (in rad) between [0,ph_max]
    .Eth - [Na x Nf] matrix of Eth values (complex) evaluated at the corresponding (th,ph) row
    .Eph - [Na x Nf] matrix of Eph values (complex) evaluated at the corresponding (th,ph) row 
    .freq - [1 x Nf] vector of frequencies (in Hz)
  Optional (depending on calling function)  
    .pol - 'x' or 'y'  

Na = number of angular samples: length(unique(.th))*length(unique(.ph))
Nf = number of frequency samples    


Example for angles (first columns):
    .th = [0,1,2,3,0,1,2,3,0,1,2,3].'
    .ph = [0,0,0,0,1,1,1,1,2,2,2,2].'


