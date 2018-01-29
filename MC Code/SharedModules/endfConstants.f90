module endfConstants

  use numPrecision

  implicit none
  integer(shortInt), parameter :: histogramInterpolation = 1 ,&
                                  linLinInterpolation    = 2, &
                                  linLogInterpolation    = 3, &
                                  logLinInterpolation    = 4, &
                                  loglogInterpolation    = 5, &
                                  chargedParticleInterpolation = 6
  ! Seperate parameters for interpolationFlags for tabular PDF. Modern ACE files use flags
  ! consistant with ENDF interpolation paramethers, but old MCNP manual(year 2000) specifies
  ! inconsistant flags: histogram = 0 and linLin = 1
  integer(shortInt), parameter :: tabPdfHistogram = 1 ,&
                                  tabPdfLinLin    = 2


  ! List of reaction MT numbers, See Serpent 2 Wiki or ENDF manual for exact details:
  !
  ! Trkov, A., M. Herman, and D. A. Brown. “ENDF-6 Formats Manual.”
  ! Data Formats and Procedures for the Evaluated Nuclear Data Files ENDF/B-VI and ENDF/B-VII,
  ! National Nuclear Data Center Brookhaven National Laboratory, Upton, NY, 2012, 11973–5000.
  !
  ! Dictionary:
  ! N - neutron ; d - deutron; a - alpha particle ; f - fission; p - proton
  !
  ! Syntax:
  ! Example:  N_2N -> (n,2n)
  ! words indicate sum of relevant reactions e.g: (n,fission) is sum (n,f)+(n,nf)+(n,2nf)+(n3nf)

  integer(shortInt), parameter :: N_total       = 1   ,&
                                  N_N_elastic   = 2   ,&
                                  N_N_inelastic = 4   ,&
                                  N_anything    = 5   ,&
                                  N_2Nd         = 11  ,&
                                  N_2N          = 16  ,&
                                  N_3N          = 17  ,&
                                  N_fission     = 18  ,&
                                  N_f           = 19  ,&
                                  N_Nf          = 20  ,&
                                  N_2Nf         = 21  ,&
                                  N_Na          = 22  ,&
                                  N_N3a         = 23  ,&
                                  N_2Na         = 24  ,&
                                  N_3Na         = 25  ,&
                                  N_absorbtion  = 27  ,&
                                  N_Np          = 28  ,&
                                  N_N2a         = 29  ,&
                                  N_2N2a        = 30  ,&
                                  N_Nd          = 32  ,&
                                  N_Nt          = 33  ,&
                                  N_Nhe3        = 34  ,&
                                  N_Nd2a        = 35  ,&
                                  N_Nt2a        = 36  ,&
                                  N_4N          = 37  ,&
                                  N_3Nf         = 38  ,&
                                  N_2Np         = 41  ,&
                                  N_3Np         = 42  ,&
                                  N_N2p         = 44  ,&
                                  N_Npa         = 45  ,&
                                 !N_Nl(:) 51-90
                                 !Inelastic scattering from levels 1-40 is defined at the end
                                  N_Ncont       = 91  ,&
                                  N_disap       = 101 ,&
                                  N_gamma       = 102 ,&
                                  N_p           = 103 ,&
                                  N_d           = 104 ,&
                                  N_t           = 105 ,&
                                  N_he3         = 106 ,&
                                  N_a           = 107 ,&
                                  N_2a          = 108 ,&
                                  N_3a          = 109 ,&
                                  N_2p          = 111 ,&
                                  N_pa          = 112 ,&
                                  N_t2a         = 113 ,&
                                  N_d2a         = 114 ,&
                                  N_pd          = 115 ,&
                                  N_pt          = 116 ,&
                                  N_da          = 117


  integer(shortInt),private    :: i  ! Local, private integer to use array constructor
  integer(shortInt),parameter  :: N_Nl(40)      = [(50+i, i =1,40)]



end module endfConstants