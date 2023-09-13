Matlab files for solving model in "Optimal Environmental Policy Under Economic Fluctuations."
Garth Heutel
March 4, 2011
gaheutel@uncg.edu

The folder contains six matlab .m files.

SolveSystem.m solves the first-best social planner problem.  It calls steadystate.m to solve for the steady state values for arbitrary parameter values.

SolveSystem_tax.m solves the decentralized model with an emissions tax.  It calls steadystate_tax.m to solve for the steady state values for arbitrary parameter values.

SolveSystem_tax_AsymInfo.m solves the decentralized tax policy model with asymmetric information between regulators and firms/consumers.  It calls steadystate_tax.m. 
SolveSystem_quota_AsymInfo.m solves the decentralized quota policy model, and it also calls steadystate_tax.m (since the steady-states of those two policies are identical).

All of these programs also call on the Anderson-Moore algorithm programs, available at http://www.bos.frb.org/economic/econbios/fuhrer/matlab.htm.  

Feel free to email or call me if you have any questions: gaheutel@uncg.edu, (336) 334 4872.