#=
=======================================================================
                 GAS CONTINUUM OPACITY/EMISSION MODULE

 In this module various gas continuum opacities and emissivities are
 implemented. Currently we have:
  
  - General free-free extinction and emission for a thermal electron 
    and ion distribution (Gordon & Sorochenko 2002, Kluwer book)

 Planned future continua are:
 
  - Various bound-free continua, which may be linked to a future 
    planned ionization module. 
  - H- opacity (T.L.John 1988, A&A 193, 189-192), consisting of:
    - Free-free part:                       hnu + e- + H --> H + e- 
    - Bound-free part (photo-detachment):   hnu + H-     --> H + e-
  - H2- opacity (W.B.Somerville 1963, ApJ 139, 192), consisting of:
    - Free-free part:                       hnu + e- + H2 --> H2 + e-
    - (bound-free is argued to be negligible)
  - H2+ opacity (Mihalas 1965)
  - He opacity (Peach 1970)
  - He- opacity (Carbon, Gingerich & Latham 1969)
  - H2+H2 and H2+He opacity (Borysov 2002)
 =======================================================================
=#