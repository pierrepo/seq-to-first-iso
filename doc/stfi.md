# Theoretical background

## Isotopologues in mass spectrometry

Show a graph of isotopic cluster
![](_static/iso_dist_NC.png)
*Mass spectrum of peptide with sequence "VGEVFINYIQRQNELFQGKLAYLIIDTCLSIVRPNDSKPLDNR"*

## Compute M0 intensity

add formula

.. math::
M_0 =  a(^{12}{\rm C})^{n({\rm C})}
\times a(^{ 1}{\rm H})^{n({\rm H})}
\times a(^{16}{\rm O})^{n({\rm O})}
\times a(^{14}{\rm N})^{n({\rm N})}
\times a(^{32}{\rm S})^{n({\rm S})},

## Compute M1 intensity

add formula

.. math::
\begin{aligned}
M_1 &=
n({\rm C})
\times a(^{12}{\rm C})^{n({\rm C})-1}
\times a(^{13}{\rm C})
\times a(^{ 1}{\rm H})^{n({\rm H})}
\times a(^{16}{\rm O})^{n({\rm O})}
\times a(^{14}{\rm N})^{n({\rm N})}
\times a(^{32}{\rm S})^{n({\rm S})} \\
&+
n({\rm H})
\times a(^{12}{\rm C})^{n({\rm C})}
\times a(^{ 1}{\rm H})^{n({\rm H})-1}
\times a(^{ 2}{\rm H})
\times a(^{16}{\rm O})^{n({\rm O})}
\times a(^{14}{\rm N})^{n({\rm N})}
\times a(^{32}{\rm S})^{n({\rm S})} \\
&+
n({\rm O})
\times a(^{12}{\rm C})^{n({\rm C})}
\times a(^{ 1}{\rm H})^{n({\rm H})}
\times a(^{16}{\rm O})^{n({\rm O})-1}
\times a(^{17}{\rm O})
\times a(^{14}{\rm N})^{n({\rm N})}
\times a(^{32}{\rm S})^{n({\rm S})} \\
&+
n({\rm N})
\times a(^{12}{\rm C})^{n({\rm C})}
\times a(^{ 1}{\rm H})^{n({\rm H})}
\times a(^{16}{\rm O})^{n({\rm O})}
\times a(^{14}{\rm N})^{n({\rm N})-1}
\times a(^{15}{\rm N})
\times a(^{32}{\rm S})^{n({\rm S})} \\
&+
n({\rm S})
\times a(^{12}{\rm C})^{n({\rm C})}
\times a(^{ 1}{\rm H})^{n({\rm H})}
\times a(^{16}{\rm O})^{n({\rm O})}
\times a(^{14}{\rm N})^{n({\rm N})}
\times a(^{32}{\rm S})^{n({\rm S})-1}
\times a(^{33}{\rm S}) \\
\end{aligned}

## Strategy to take into account auxotrophies

(SLIM)
todo Lilian. Add more explanations here.
