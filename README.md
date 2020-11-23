# preComputationThis code is an implementation of our precomputation scheme for window $\tau$NAF with widths from 4 to 7; Solinas' pre-computation scheme for window $\tau$NAF with widths from 4 to 5; Hankerson, Menezes, and Vanstone's pre-computation scheme for window $\tau$NAF with widths from 4 to 6; Trost and Xu's pre-computation scheme for window $\tau$NAF with widths from 4 to 6;  and scalar multiplication using these pre-conmputations; and constant-time scalar multiplication using these pre-conmputations. It is compiled by by Microsoft visual studio 2015. 
When one wants to use regular TNAF, one should set the parameter w=2. 
The code is shown in main.cpp,preComputation.cpp and preComputation.h.


Miracl lib  is used to implement  big number arithmetic. Our experiments  
are compiled and executed on Intel R Core TM i7-6567U @3.3 GHZ with Skylake architecture by using C++. The time costs are shown in our paper.  

If one wants to test for more times, please change the value of CYCLES in main.cpp. If one wants to test for more size integers, please change the value of PLambda and using different prepare_basis(...) in main.cpp. 
