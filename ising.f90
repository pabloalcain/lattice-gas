program main
  call Initialize
  call Evolve
  call Terminate
end program main

!!! INITIALIZATION ROUTINES

subroutine Initialize
  implicit none
  call Read_Input
  call Init_Variables
  call Read_Sample
  call Start_Global_Files
  call Start_Screen
end subroutine Initialize

!===> Initialize
subroutine Read_Input
  use Simulation_Control
  use Particles
  implicit none
  integer lis

  open(unit=01, file='input.in', status='old', &
       action='read')
  
9 format(A20,/,&
       I4.4,//,&
       I7.7,/,I7.7,/,I7.7,//,&
       I4.4,/,F7.4,/,F7.4,//,&
       I4.4,/,F7.4,/,F7.4,//,&
       A15) 
  
  read(unit=01, fmt=9, end=200, err=800) SampIn,&
       L,&
       Nsteps, NTherm, NCorr,&
       NTemp, TempIni, TempFin,&
       NChem, ChemIni, ChemFin,&
       SampOut
  close(unit=01)
  
  lis=len_trim(SampIn)
  InitPos = ('None'.ne.SampIn(1:lis))
  if (InitPos) then !! Override L
     open(unit=02, file=SampIn(1:lis), status='old', &
          action='read', err=700)
     read(unit=02,fmt=*, err=600) L
  endif
  
  return   !  successful exit
200 continue
  print*,'Read_Input: FATAL: premature end-of-file in standard input'
  stop
600 continue
  print*,'Read_Sample: FATAL: ',SampIn(1:lis),' is empty?'
  stop
700 continue
  print*,'Read_Sample: FATAL: ',SampIn(1:lis),' not found.'
  stop
800 continue
  print*,'Read_Input: FATAL: read error in standard input'
  stop
end subroutine Read_Input

subroutine Init_Variables
  use Particles
  use Simulation_Control
  implicit none

  allocate ( part(L,L) )
end subroutine Init_Variables

subroutine Read_Sample
  use Particles
  use Simulation_Control
  implicit none
  integer :: i, j

  call randomize

  if (InitPos) then
     do i=1,L
        read(02,*,end=800,err=900) part(i,:)
     enddo
     close(unit=02)
  else
     do i=1,L
        do j=1,L
           part(i,j)=(rand().gt.0.5)
        enddo
     enddo
  endif
  return   !  successful exit

800 continue
  print*,'Read_Sample: FATAL: premature end-of-file at line ',i
  close(unit=1)
  stop   
900 continue
  print*,'Read_Sample: FATAL: read error in ',SampIn(1:len_trim(SampIn))
  close(unit=1)
  stop
end subroutine Read_Sample

subroutine Start_Global_Files
  implicit none
  
  open(unit=30,file='data/stat.dat',status='unknown',&
       action='write')
  write(30,'(3a)') '#Title    Temperature    ',&
       '    Chemical       Energy        CV     ',&
       '       N             Kt'


  open(unit=31,file='data/key.dat',status='unknown',&
       action='write')
  write(31,'(2a)') '   Title          Temperature          L'
end subroutine Start_Global_Files

subroutine Start_Screen
  use Particles
  use Simulation_Control
  implicit none
  character(8)  :: date
  character(10) :: time


  call date_and_time(date,time)
  call system('clear')
  
  write(06, '(a,/,a,/a,////)') '*******************'&
       ,'* Ising MonteCarlo','*******************'
  write(06, '(a,5a,2x,6a)') 'Start:      ', date(1:4),&
       '-', date(5:6), '-', date(7:8), time(1:2), 'h',&
       time(3:4),'m',time(5:6),'s'
  write(06, '(a,i9)')  &
       'Number of steps:',Nsteps
  
  write(06, '(a,i6)') 'Number of particles per side:',L
  write(06, '(a,f7.4,a,f7.4,a,i4,a)') 'Temp sweep from ', Tempini,&
       ' to ',Tempfin,' with a total of ',Ntemp,' steps'
  write(06, '(a,f7.4,a,f7.4,a,i4,a)') 'Chem sweep from ', Chemini,&
       ' to ',Chemfin,' with a total of ',Nchem,' steps'
  if (InitPos) then
     write(06,'(a,a)') "Reading initial config from ", SampIn(1:len_trim(SampIn))  
  else
     write(06,'(a,a)') "Random initial config"
  endif
  write(06,'(a)') 
end subroutine Start_Screen
!<=== Initalize

!!! TIME EVOLUTION ROUTINES

subroutine Evolve
  use Simulation_Control
  implicit none
  integer :: temp_step, chem_step
  double precision :: Interpolate
  
  chemical: do chem_step=1,NChem
     chem=Interpolate(ChemIni,ChemFin,NChem,chem_step)
     temperature: do temp_step=1,Ntemp
        temp=Interpolate(TempIni,TempFin,NTemp,temp_step)
        write(06,'(/,a,i4.4,a,i4.4,a)',advance='no') 'Configuración ', &
             temp_step+(chem_step-1)*Ntemp, ' de ', Ntemp*Nchem,&
             ':   Thermalization...'
        call Start_Setup
        call MonteCarlo_Run
        call End_Setup
     end do temperature
     
  end do chemical
end subroutine Evolve

!===> Evolve
subroutine Start_Setup
  implicit None
  call Reset_Statistics
  call Start_Local_Files
end subroutine Start_Setup

subroutine MonteCarlo_Run
  use Particles
  use Simulation_Control
  use Statistics
  implicit none
  integer :: step, ipart
  
  therm: do step=1,NTherm
     do ipart=1,L**2
        call Trial_Move(int(rand()*L)+1,int(rand()*L)+1)
     enddo
  enddo therm
  
  write(06,'(a,f4.2,a)',advance='no') 'Done!  Acc=',&
       float(accepted)/float(total_attempts),'..............'
  
  evol: do step=1,Nsteps
     do ipart=1,L**2
        call Trial_Move(int(rand()*L)+1,int(rand()*L)+1)
     enddo
     call Update_Information(step)
  enddo evol
end subroutine MonteCarlo_Run

subroutine End_Setup
  call Update_Statistics
  call Update_Global_Files
  call End_Local_Files
end subroutine End_Setup

!------> Start_Setup
subroutine Start_Local_Files
  use Particles
  use Simulation_Control
  implicit none
  integer time_since_epoch
  character(8)  :: date
  character(10) :: time

  character(23) :: evol, fina, traj, info
  call system_clock(time_since_epoch)
  call date_and_time(date,time)
  
  write(title, '(i10.10)') time_since_epoch
  write(evol, '(a9,a10,a4)') 'data/evol', &
       title,'.dat'
  write(fina, '(a9,a10,a4)') 'data/fina', &
       title,'.dat'
  write(traj, '(a9,a10,a4)') 'data/traj', &
       title,'.dat'
  write(info, '(a9,a10,a4)') 'data/info', &
       title,'.dat'
  
  open(unit=40, file=evol, status='new', &
       action='write')
  open(unit=41, file=fina, status='new', &
       action='write')
  open(unit=42, file=traj, status='new', &
       action='write')
  open(unit=43, file=info, status='new', &
       action='write')
  


  write(40, '(a,/,a,/,a)')  '#', &
       '#  Step     Temperature     Chemical      Energy         N',&
       '# ------- --------------- ------------ ------------ -----------'
  
  write(43, '(a,/,a)') '#','# Ising MonteCarlo'
  write(43, '(a,5a,2x,6a)') '# Start:      ', date(1:4),&
       '-', date(5:6), '-', date(7:8), time(1:2), 'h',&
       time(3:4),'m',time(5:6),'s'
  write(43, '(2a)') '# Title:         ',title(1:len_trim(title))
  write(43, '(2a)') '# Input sample:  ',SampIn(1:len_trim(SampIn))
  write(43, '(2a)') '# Output sample: ',SampOut(1:len_trim(SampOut))
  write(43, '(a,i9)')  &
       '# Number of steps:',Nsteps
  write(43, '(a,i6)') '# Number of particles per side:',L
  write(43, '(a,f12.6)') '# Constant T run with T =',temp
  write(43, '(a,f12.6)') '# Constant mu run with mu =',chem
end subroutine Start_Local_Files

subroutine Reset_Statistics
  use Statistics
  use Particles
  use Simulation_Control
  implicit none
  integer :: i,j
  double precision :: energy

  ene_tot=0.d0
  n_tot=0
  do i=1,L
     do j=1,L
        ene_tot=ene_tot+energy(i,j)
        if (part(i,j)) then
           n_tot=n_tot+1
        endif
     enddo
  enddo
  
  ene_tot=ene_tot/2
  
  ene_ave  = 0.d0
  ene_aux  = 0.d0

  n_ave  = 0.d0
  n_aux  = 0.d0

  total_attempts = 0
  accepted = 0
end subroutine Reset_Statistics
!<------Start_Setup

!------> MonteCarlo_Run

subroutine Trial_Move(i, j)
  use Particles
  use Statistics
  use Simulation_Control
  
  implicit none
  integer, intent(in) :: i, j
  double precision :: deltae, deltau
  double precision :: energy, flip_energy
  integer :: deltan
  
  total_attempts=total_attempts+1
  if (part(i,j)) then
     deltae = -energy(i,j)
     deltan = -1
  else
     deltae = flip_energy(i,j)
     deltan = 1
  endif
  deltau=deltae-chem*deltan
  
  if (deltau.lt.0.or.rand().lt.exp(-deltau/temp)) then
     accepted=accepted+1
     part(i,j)=.not.part(i,j)
     ene_tot = ene_tot + deltae
     n_tot   = n_tot   + deltan
  endif

end subroutine Trial_Move

subroutine Update_Information(step)
  use Particles
  use Simulation_Control
  use Statistics
  implicit none
  integer, intent(in) :: step
  integer i
  
  if (mod(step,NCorr).eq.0) then
     ene_ave = ene_ave + ene_tot
     ene_aux = ene_aux + ene_tot**2

     n_ave = n_ave + n_tot
     n_aux = n_aux + n_tot**2
          

     write(40,*) step, temp, chem, ene_tot, n_tot
     
     if (mod(step,NCorr*10).eq.0) then
        write(42, *) L, int(step/(NCorr*10))
        do i=1,L
           write(42,*) part(i,:)
        enddo
     endif
  endif
  
  if (mod(step*1000,Nsteps).eq.0) then
     write(06,'(a,f5.1,a)',advance='no') '\b\b\b\b\b\b',&
          float(step*100)/Nsteps,'%'
  endif
end subroutine Update_Information

!<------ MonteCarlo_Run

!------> End_Setup
subroutine Update_Statistics
  use Statistics
  use Simulation_Control
  implicit none
  ene_ave  = ene_ave   /(Nsteps/NCorr)
  cv       = (ene_aux  /(Nsteps/NCorr) - ene_ave**2)/temp**2

  n_ave    = n_ave     /(Nsteps/NCorr)
  kt       = L**2*(n_aux  /(Nsteps/NCorr) - n_ave**2)/(n_ave**2 * temp)

end subroutine Update_Statistics

subroutine Update_Global_Files
  use Statistics
  use Particles
  use Simulation_Control
  implicit none
  
  write(30,'(1x,a,6f16.6)') title, temp, chem, ene_ave, cv, n_ave, kt
  write(31,'(1x,a,2f16.6,I7.7)') title, temp, chem, L
end subroutine Update_Global_Files

subroutine End_Local_Files
  use Particles
  use Statistics
  use Simulation_Control
  implicit none
  integer i

  character(8)  :: date
  character(10) :: time

  call date_and_time(date,time)
  
  write(40,'(a,/,a,4f14.6)'),'#','# Means', &
       temp, chem, ene_ave, n_ave
  close(unit=40)
  
  write(41, *) L, 1
  do i=1,L
     write(41,*) part(i,:)
  enddo
  close(unit=41)
  
  close(unit=42)
  
  write(43, '(a,5a,2x,6a)') '# End:      ', date(1:4),&
       '-', date(5:6), '-', date(7:8), time(1:2), 'h',&
       time(3:4),'m',time(5:6),'s'
  close(unit=43)

end subroutine End_Local_Files
!<------ End_Setup

!<=== Evolve

!!! TERMINATION ROUTINES

subroutine Terminate
  call End_Screen
  call End_Global_Files
  call Term_Variables
end subroutine Terminate

!===> Terminate
subroutine End_Screen
  write(06,'(/,a)') 'Finalizado'
end subroutine End_Screen

subroutine End_Global_Files
  close(unit=30)
  close(unit=31)
end subroutine End_Global_Files

subroutine Term_Variables
  use Particles
  implicit none
  deallocate(part)
end subroutine Term_Variables
!<=== Terminate

!!! AUXILIAR ROUTINES
double precision function Interpolate(yi,yf,N,i)
  double precision, intent(in) :: yi, yf
  integer         , intent(in) :: N, i

  if (N.gt.1) then
     Interpolate=(yf-yi)*real(i-1)/real(N-1)+yi
  else if (N.eq.1) then
     Interpolate=yi
  else
     Interpolate=0
     write(06,'(a)') "Number of Temperature steps not allowed"
  endif
end function Interpolate

subroutine randomize
  implicit none
  integer*4 time_array(8)
  real random1
  call date_and_time(values=time_array)
  random1=rand(time_array(7)+time_array(8))
end subroutine randomize

function energy(x,y)
  use Potential
  use Particles
  use Simulation_Control
  implicit none
  integer, intent(in) :: x, y
  integer :: i, j, ip, jp
  double precision :: rsq
  
  double precision :: energy
  
  
  energy=0.0
  if (part(x,y)) then
     do ip=x-Rbump,x+Rbump
        do jp=y-Rbump,y+Rbump
           Rsq=(ip-x)**2+(jp-y)**2
           Rsq=(ip-x)**2+(jp-y)**2
           if (Rsq.le.Rbump**2) then
              i=mod(L+ip-1,L)+1
              j=mod(L+jp-1,L)+1
              if (part(i,j).and.i.ne.x.and.j.ne.y) then
                 if (Rsq.le.Rvall**2) then
                    energy=energy+Evall
                 else
                    energy=energy+Ebump
                 endif
              endif
           endif
        enddo
     enddo
  endif
end function energy

function flip_energy(x,y)
  use Potential
  use Particles
  use Simulation_Control
  implicit none
  integer, intent(in) :: x, y
  integer :: i, j, ip, jp
  double precision :: rsq

  double precision :: flip_energy

  flip_energy=0.0
  if (.not.part(x,y)) then
     do ip=x-Rbump,x+Rbump
        do jp=y-Rbump,y+Rbump
           Rsq=(ip-x)**2+(jp-y)**2
           if (Rsq.le.Rbump**2) then
              i=mod(L+ip-1,L)+1
              j=mod(L+jp-1,L)+1
              if (part(i,j).and.i.ne.x.and.j.ne.y) then ! la segunda condición es trivial
                 if (Rsq.le.Rvall**2) then
                    flip_energy=flip_energy+Evall
                 else
                    flip_energy=flip_energy+Ebump
                 endif
              endif
           endif
        enddo
     enddo
  endif
end function flip_energy
