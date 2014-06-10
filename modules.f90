module Particles
  integer, parameter :: DIM=2
  logical, dimension(:,:)  , allocatable :: part

  double precision :: ene_tot, n_tot
end module Particles

module Statistics
  double precision :: ene_ave, ene_aux
  double precision :: n_ave  , n_aux
  
  double precision :: cv, kt

  integer :: total_attempts, accepted
end module Statistics

module Simulation_Control
  integer            :: L=0

  
  double precision :: TempIni, TempFin
  integer          :: NTemp
  double precision :: ChemIni, ChemFin
  integer          :: NChem

  
  integer :: Nsteps, NCorr, NTherm

  character*20 :: SampIn, SampOut

  double precision :: temp
  double precision :: chem

  
  logical :: InitPos
  character*10  :: title
end module Simulation_Control

module Potential
  double precision :: Evall=-0.3, Ebump=0.1
  integer :: Rvall=2, Rbump=7
end module Potential
