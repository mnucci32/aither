# Two reaction, two species dissociating oxygen mechanism
# Rates from Luke, adjusted to be in mol instead of kmol
#
# Luke, E. A. "A Rule-Based Specification System For Computational Fluid 
# Dynamics". Ph. D. Thesis, December 1999.

2 O2 <=> 2 O + O2 : forwardRate=arrhenius(C=1.9e15, eta=-1.5, theta=59500)
O2 + O <=> 3 O : forwardRate=arrhenius(C=6.4e17, eta=-2, theta=59500)