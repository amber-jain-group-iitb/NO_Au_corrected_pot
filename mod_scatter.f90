Module mod_scatter
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Dynamics
integer nsteps_eq,nstep_write,nsteps_ad
integer traj_num
real*8 dtc,tim_eq,tim_ad
integer flag_hop,flag_average,flag_terminate,flag_write,flag_friction,flag_write_vmd
real*8 curr_time,en_avg,ensq_avg
integer n_Etr,n_traj
real*8 E_transl,z_cutoff,adsorbed,E_rot
real*8 NO_com(3)

!! Classical
integer n_Au,n_layers,n_rows,n_columns
real*8 box_length(3),box_pbc_min(2),box_pbc_max(2)
real*8,allocatable,dimension(:,:) :: x_Au,v_Au,acc_Au
real*8,allocatable,dimension(:) :: mass_Au
real*8,dimension(3):: x_O,v_O,acc_O,x_N,v_N,acc_N
real*8 mass_O,mass_N
real*8 pot_en,energy,KE_AU,KE_NO,KE_com,KE_vib,KE_rot
real*8,allocatable :: temp_layer(:)
real*8 temperature,gamma_fric

!! Quantum
integer nquant,state
real*8,allocatable:: d_ij(:,:,:),V_k(:),si_adiab(:,:),si_adiab_old(:,:)

!! Potential parameters
  !! H_00
  real*8 A0,alpha0,B0,beta0,F0,gamma0,r0_NO,e_alpha0_rc,e_beta0_rc
  !! H_11
  real*8 A1,alpha1,B1,beta1,r1_AuN,D,C,z_image,F1,gamma1,r1_NO,phi,Ea,phi_Ea,e_alpha1_rc,e_2beta1_rc_r1,e_beta1_rc_r1,Csq
  !! H_01
  real*8 A2,A3,gamma2,B2,B3,gamma3,V_01_AuO_cutoff,V_01_AuN_cutoff
  !! Misc
  real*8 r_cutoff,a_Au
  !! Au-Au
  real*8 alpha,beta,gamma
  real*8,allocatable::A_mat(:,:,:,:)

!! Potential
integer n_N_Au,n_O_Au
integer,allocatable :: n_Au_Au(:),layer_num(:)
integer,allocatable :: list_N(:),list_O(:),list_Au(:,:)
real*8,allocatable :: dist_N(:),dist_O(:),dist_Au(:,:)
real*8,allocatable :: rij_0(:,:,:)
real*8 r_NO,cos_theta,z_COM,dcostheta_dxN(3)
real*8,allocatable:: H_diab(:,:),dH_diab_Au(:,:,:,:),dH_diab_N(:,:,:),dH_diab_O(:,:,:)
real*8 V_Au_Au
real*8,allocatable:: dV_Au_Au(:,:)

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer,allocatable:: seed(:)
integer nold,cnt_rate
real*8 tim_tot,tim_diag,tim_evolve,tim_wr_out,tim_pot_cl,tim_comp_list,tim_diab
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  integer i,size_seed,seed2(2)
  character st_ch

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  open(10,file="scatter.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) n_layers
  read(10,*) n_rows
  read(10,*) n_columns
  read(10,*) nquant
  read(10,*) tim_eq
  read(10,*) tim_ad
  read(10,*) dtc
  read(10,*) n_traj
  read(10,*) n_Etr
  read(10,*) z_cutoff
  read(10,*) flag_write
  read(10,*) flag_write_vmd
  read(10,*) nstep_write
  read(10,*) temperature
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  if(mod(n_rows,2).ne.0) then
    write(6,*) "n_rows must be even",n_rows
    stop
  endif

  !---------------------------------------------------------- 

  n_Au=n_layers*n_rows*n_columns
  nsteps_eq=nint(tim_eq/dtc)
  nsteps_ad=nint(tim_ad/dtc)

  allocate(x_Au(n_Au,3),v_Au(n_Au,3),acc_Au(n_Au,3),mass_Au(n_Au))
  allocate(list_N(n_Au),list_O(n_Au),dist_N(n_Au),dist_O(n_Au))
  allocate(n_Au_Au(n_Au),layer_num(n_Au),list_Au(n_Au,12),dist_Au(n_Au,12))
  allocate(temp_layer(n_layers))
  allocate(A_mat(3,3,n_Au,n_Au),rij_0(3,n_Au,n_Au))
  allocate(H_diab(nquant,nquant),dH_diab_Au(nquant,nquant,n_Au,3),dH_diab_N(nquant,nquant,3),dH_diab_O(nquant,nquant,3))
  allocate(dV_Au_Au(n_Au,3))
  allocate(V_k(nquant),si_adiab(nquant,nquant),si_adiab_old(nquant,nquant))
  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !---------------------------------------------------------- 

  call set_parameters

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0
  tim_wr_out=0.d0
  tim_diag=0.d0
  tim_pot_cl=0.d0
  tim_comp_list=0.d0
  tim_evolve=0.d0
  tim_diab=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j
  integer c1,c2

  call system_clock(c1)

  if(iflow==1.or.iflow==2) then
    call files(0)
    do i=1,n_Etr
      if(n_Etr>1)E_transl=5.d0+75.d0*(i-1)/dfloat(n_Etr-1)
      if(n_Etr==1)E_transl=10.d0
      E_transl=E_transl*1.66d-21 !! kj_mol to J
      adsorbed=0.d0
      E_rot=0.d0
      do j=1,n_traj
        traj_num=j
        call equilibrate
        call adsorb
      enddo
      write(100,'(4f15.7)')E_transl/1.66d-21,adsorbed/n_traj,E_rot*6.02d20,(dfloat(n_traj)-adsorbed)
    enddo
    
    call system_clock(c2)
    tim_tot=(c2-c1)/real(cnt_rate)
    call files(1)
  endif

end subroutine main
!---------------------------------------------------------- 

subroutine files(nflag)
  implicit none
  integer,intent(in)::nflag

  if(nflag==0) then
    open(10,file="output")
    open(11,file="output_cl_N")
    open(12,file="output_cl_O")
    open(14,file="output_cl_NO")
    open(13,file="output_cl_Au")
    open(100,file="adsorb.out")
  endif

  if(nflag==1) then
    write(10,*) "Total time",tim_tot
    write(10,*) "  Write output",tim_wr_out
    write(10,*) "  Evolve",tim_evolve
    write(10,*) "    compute_list",tim_comp_list
    write(10,*) "    Diag",tim_diag
    write(10,*) "      potential_classical",tim_pot_cl
    write(10,*) "      compute_diab",tim_diab
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(100)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine equilibrate
  !! Equilibrates Au lattice
  implicit none
  integer i

  call set_Au_lattice
  call set_NO_frozen
  call Boltzmann_velocities

  call setup_evolve(1)
  do i=1,10
    call evolve(nsteps_eq/10)
    call Boltzmann_velocities
    call evaluate_variables(0)
    call evaluate_variables(1)
  enddo

end subroutine equilibrate
!-----------------------------------------------------------------  

subroutine adsorb
  implicit none
  integer i

  call set_NO_position
  call setup_evolve(2)
  call evolve(nsteps_ad)
  if(z_com<z_cutoff) adsorbed=adsorbed+1.d0
  if(z_com>z_cutoff) E_rot=E_rot+KE_rot
  if(z_com>z_cutoff) write(6,*) traj_num,KE_NO/wave_to_J,KE_rot/wave_to_J,KE_vib/wave_to_J,KE_com/wave_to_J

end subroutine adsorb
!-----------------------------------------------------------------  

subroutine setup_evolve(nflag)
  implicit none
  integer,intent(in):: nflag

  !! Equilibration setup
  if(nflag==1) then
    state=1;  flag_hop=0;  flag_average=1;  flag_terminate=0;  flag_friction=0
    gamma_fric=10.d0* 80.d0*2*pi*clight
    curr_time=0.d0
  endif

  !! Adsorption setup
  if(nflag==2) then
    state=1;  flag_hop=0;  flag_average=1;  flag_terminate=1;  flag_friction=0
  endif

  call evaluate_variables(0)
  call evaluate_variables(1)

  en_avg=0.d0
  ensq_avg=0.d0

end subroutine setup_evolve
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,iterm
  integer c1,c2

  call system_clock(c1)

  do i=1,nsteps
    call write_output(i,0,nsteps)
    if(flag_average==1) call average(i)
    call evolve_classical(dtc)
    call pbc

    if(flag_hop==1)call hop

    if(flag_terminate==1) call terminate(i,nsteps,iterm)
    if(iterm==1)exit

    curr_time=curr_time+dtc
  enddo
  call write_output(1,1,nsteps)

  call system_clock(c2)
  tim_evolve=tim_evolve+(c2-c1)/real(cnt_rate)

end subroutine evolve
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i

  en_avg=en_avg+energy
  ensq_avg=ensq_avg+energy*energy

end subroutine average
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  real*8,intent(in) :: dt
  real*8 acc_old_N(3),acc_old_O(3),acc_old_Au(n_Au,3)
  integer i
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(n_Au,3),delta_v(n_Au,3)

  if(flag_friction==0) then

    x_N =x_N +v_N*dt +0.5*acc_N*dt*dt
    x_O =x_O +v_O*dt +0.5*acc_O*dt*dt
    do i=1,n_Au
      if(layer_num(i)>1)x_Au(i,:)=x_Au(i,:)+v_Au(i,:)*dt+0.5*acc_Au(i,:)*dt*dt
    enddo

    acc_old_N=acc_N
    acc_old_O=acc_O
    acc_old_Au=acc_Au

    call evaluate_variables(0)
    v_N =v_N +0.5*dt*(acc_old_N+acc_N)
    v_O =v_O +0.5*dt*(acc_old_O+acc_O)
    do i=1,n_Au
      if(layer_num(i)>1)v_Au(i,:)=v_Au(i,:)+0.5*dt*(acc_old_Au(i,:)+acc_Au(i,:))
      if(layer_num(i)==1)v_Au(i,:)=0.d0
    enddo

    call evaluate_variables(1)
  endif

  if(flag_friction==1) then
    gama_dt=gamma_fric*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)

    call stochastic_force(delta_r,delta_v,dt)
    do i=1,n_Au
      if(layer_num(i)>1)x_Au(i,:)=x_Au(i,:)+c1*dt*v_Au(i,:)+c2*dt*dt*acc_Au(i,:)+delta_r(i,:)
    enddo

    acc_old_Au=acc_Au
    call evaluate_variables(0)

    do i=1,n_Au
      if(layer_num(i)>1)v_Au(i,:)=v_Au(i,:)+(c1-c2)*dt*acc_old_Au(i,:)+c2*dt*acc_Au(i,:)+delta_v(i,:)
      if(layer_num(i)==1)v_Au(i,:)=0.d0
    enddo

    call evaluate_variables(1)
  endif

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine hop
  implicit none

end subroutine hop
!-----------------------------------------------------------------  

subroutine terminate(i,nsteps,iterm)
  implicit none
  integer,intent(in)::i,nsteps
  integer,intent(out)::iterm
  integer j

  iterm=0

  if(z_com>11.d-10) then
    iterm=1
    do j=i+1,nsteps
      call write_output(j,0,nsteps)
      call average(j)
      curr_time=curr_time+dtc
    enddo
  endif

end subroutine terminate
!-----------------------------------------------------------------  

subroutine write_output(n,nflag,nsteps)
  implicit none
  integer,intent(in)::nflag,n,nsteps
  integer i
  integer c1,c2

  call system_clock(c1)

  if(nflag==0) then
    if(flag_write==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(3es17.5,i5)')curr_time*1.d15,energy/wave_to_J,KE_AU/(1.5*(n_Au-n_rows*n_columns)*wave_to_J),state
        write(11,'(4f15.5)')curr_time*1.d15,x_N*1.d10
        write(12,'(4f15.5)')curr_time*1.d15,x_O*1.d10
        write(14,'(8f15.5)')curr_time*1.d15,r_NO*1.d10,dacos(cos_theta)*180/pi,z_com*1.d10,&
KE_NO/wave_to_J,KE_COM/wave_to_J,KE_vib/wave_to_J,KE_rot/wave_to_J
        write(15,'(5f15.5)')curr_time*1.d15,temp_layer
        !write(13,'(f15.5)')curr_time*1.d15
        !do i=n_Au,1,-1
        !  write(13,'(3f15.7,i5)')x_Au(i,:)*1.d10,layer_num(i)
        !  if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(13,*)
        !  if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(13,*)
        !enddo
        !write(13,*)
        !write(13,*)
        if(flag_write_vmd==1) then
          write(13,'(i5)') n_Au+2
          write(13,*) curr_time*1.d15
          write(13,'(A,3f15.7)')'N ',x_N(:)*1.d10
          write(13,'(A,3f15.7)')'O ',x_O(:)*1.d10
          do i=1,n_Au
            write(13,'(A,3f15.7)')'Au ',x_Au(i,:)*1.d10!,layer_num(i)
            !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
            !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
          enddo
        endif
      endif
    endif
  endif

  if(nflag==1) then
    write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nsteps))/dfloat(nsteps))/wave_to_J
    if(flag_write==1) then
      write(10,*)"traj num=",traj_num
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
      !write(13,*);write(13,*)
      write(14,*);write(14,*)
    endif
  endif

  call system_clock(c2)
  tim_wr_out=tim_wr_out+(c2-c1)/real(cnt_rate)

end subroutine write_output
!-----------------------------------------------------------------  

subroutine Boltzmann_velocities
  !! Boltzmann velocities for Au atoms
  implicit none
  integer i,j
  real*8 sig_p,rnd,fac,vcom

  do i=1,n_Au
    do j=1,3
      sig_p=dsqrt(kb*temperature*mass_Au(i))
      call gaussian_random_number(rnd)
      if(layer_num(i)>1)v_Au(i,j)=dabs(1.d0/mass_Au(i)*(rnd*sig_p))
      if(layer_num(i)==1)v_Au(i,j)=0.d0
    enddo
  enddo

  !! COM velocity
  do j=1,3
    vcom=sum(v_Au(:,j))/real(n_Au-n_rows*n_columns)
    do i=1,n_Au
      if(layer_num(i)>1)v_Au(i,j)=v_Au(i,j)-vcom
    enddo
  enddo

  !! rescaling
  call evaluate_variables(0)
  call evaluate_variables(1)
  fac=dsqrt(3*(n_Au-n_rows*n_columns)*kb*temperature/(2*KE_Au))
  v_Au=v_Au*fac

  call evaluate_variables(1)

end subroutine Boltzmann_velocities
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j
  real*8 v_com(3),v12(3),r12(3),vr,v_ang_sq,mm,mu
  real*8 x_com(3),rn,ro,mom_inert,fr

  if(flag==0) then
    call NO_values
    call compute_list
    call tise
  endif

  if(flag==1) then
    energy=pot_en
    KE_Au=0.d0
    temp_layer=0.d0
    do i=1,3
      do j=1,n_Au
        KE_Au=KE_Au+0.5*(mass_Au(j)*v_Au(j,i)*v_Au(j,i))
        temp_layer(layer_num(j))=temp_layer(layer_num(j))+0.5*mass_Au(j)*v_Au(j,i)*v_Au(j,i)
      enddo
    enddo
    KE_NO=0.5*mass_N*sum(v_N*v_N)+0.5*mass_O*sum(v_O*v_O)
    temp_layer=temp_layer*2/(3*n_rows*n_columns*kb)
    energy=energy+KE_AU+KE_NO

    mm=mass_N+mass_O
    mu=mass_N*mass_O/mm
    v_com=(mass_N*v_N+mass_O*v_O)/mm
    v12=v_O-v_N
    r12=x_O-x_N
    call pbc_rij(r12)
    vr=sum(v12*r12)/r_NO
    v_ang_sq=(sum(v12*v12)-vr**2)

    KE_COM=0.5*mm*sum(v_com*v_com)
    KE_vib=0.5*mu*vr**2
    KE_rot=0.5*mu*v_ang_sq

!    x_com=(mass_N*x_N+mass_O*x_O)/mm
!    rN=dsqrt(sum((x_com-x_N)*(x_com-x_N)))
!    rO=dsqrt(sum((x_com-x_O)*(x_com-x_O)))
!write(6,*) rN,RO,r_NO
!
!    mom_inert=mass_N*rn**2+mass_O*rO**2
!    fr=sum(v_N*(x_O-x_N))/(r_NO*rN)
!write(6,*) fr
!    fr=sum(v_O*(x_O-x_N))/(r_NO*rO)
!write(6,*) fr
!stop

  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 ens(nquant),vect(nquant,nquant),H_diab_sav(nquant,nquant)
  real*8 delH_dels(nquant,nquant),delH_dels_ad(nquant,nquant)
  real*8 pot_cl,acc_cl_N(3),acc_cl_O(3),acc_cl_Au(n_Au,3)
  integer c1,c2

  call system_clock(c1)

  call compute_diabats
  H_diab_sav=H_diab
  call diag(H_diab_sav,nquant,ens,vect,nquant)

  si_adiab_old=si_adiab
  do i=1,nquant
    si_adiab(:,i)=vect(:,i)
    if(sum(si_adiab(:,i)*si_adiab_old(:,i))<0.d0)si_adiab(:,i)=-si_adiab(:,i)
  enddo

  V_k=ens
  pot_en=V_k(state)
  do i=1,3
    acc_N(i)=-1.d0/mass_N*sum(si_adiab(:,state)*matmul(dH_diab_N(:,:,i),si_adiab(:,state)))
    acc_O(i)=-1.d0/mass_O*sum(si_adiab(:,state)*matmul(dH_diab_O(:,:,i),si_adiab(:,state)))
    do j=1,n_Au
      acc_Au(j,i)=-1.d0/mass_Au(j)*sum(si_adiab(:,state)*matmul(dH_diab_Au(:,:,j,i),si_adiab(:,state)))
    enddo
  enddo

  call potential_classical(pot_cl,acc_cl_N,acc_cl_O,acc_cl_Au)
  pot_en=pot_en+pot_cl
  acc_N=acc_N+acc_cl_N
  acc_O=acc_O+acc_cl_O
  acc_Au=acc_Au+acc_cl_Au

  !d_ij=0.d0
  !d_ij(1,2,1)=delH_dels_ad(1,2)/(V_k(2)-V_k(1))
  !d_ij(2,1,1)=-d_ij(1,2,1)

  call system_clock(c2)
  tim_diag=tim_diag+(c2-c1)/real(cnt_rate)

end subroutine tise
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,acc_cl_N,acc_cl_O,acc_cl_Au)
  implicit none
  real*8,intent(out) :: pot_cl,acc_cl_N(3),acc_cl_O(3),acc_cl_Au(n_Au,3)
  integer i
  integer c1,c2

  call system_clock(c1)

  call gold_gold_coupling

  pot_cl=V_Au_Au
  acc_cl_N=0.d0
  acc_cl_O=0.d0
  do i=1,n_Au
    acc_cl_Au(i,:)=-1.d0/mass_Au(i)*(dV_Au_Au(i,:))
  enddo

  call system_clock(c2)
  tim_pot_cl=tim_pot_cl+(c2-c1)/real(cnt_rate)

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine compute_diabats
  !! computes H_00,H_11,H_01
  !!          acc_00,acc_11,acc_01
  integer i
  integer c1,c2

  call system_clock(c1)

  H_diab=0.d0
  dH_diab_O=0.d0;dH_diab_N=0.d0;dH_diab_Au=0.d0

  call gold_oxygen_coupling
  call gold_nitrogen_coupling
  call oxygen_nitrogen_coupling
  call image_potential

  H_diab(2,2)=H_diab(2,2)+phi_Ea

  H_diab(2,1)=H_diab(1,2)
  dH_diab_N(2,1,:)=dH_diab_N(1,2,:)
  dH_diab_O(2,1,:)=dH_diab_O(1,2,:)
  dH_diab_Au(2,1,:,:)=dH_diab_Au(1,2,:,:)

  call system_clock(c2)
  tim_diab=tim_diab+(c2-c1)/real(cnt_rate)

end subroutine compute_diabats
!-----------------------------------------------------------------  

subroutine gold_oxygen_coupling
  implicit none
  integer i,j
  real*8 tmp1,deriv(3),tmp2(3)

  !! Gold oxygen potential
  do i=1,n_O_au
    j=list_O(i)
    tmp2=(x_O-x_Au(j,:))!/dist_O(i)
    call pbc_rij(tmp2)
    tmp2=tmp2/dist_O(i)

    tmp1=dexp(-alpha0*dist_O(i))
    deriv=A0*tmp1*(-alpha0)*tmp2
    H_diab(1,1)=H_diab(1,1)+A0*(tmp1-e_alpha0_rc)
    dH_diab_O(1,1,:)=dH_diab_O(1,1,:)+deriv
    dH_diab_Au(1,1,j,:)=dH_diab_Au(1,1,j,:)-deriv

    tmp1=dexp(-alpha1*dist_O(i))
    deriv=A1*tmp1*(-alpha1)*tmp2
    H_diab(2,2)=H_diab(2,2)+A1*(tmp1-e_alpha1_rc)
    dH_diab_O(2,2,:)=dH_diab_O(2,2,:)+deriv
    dH_diab_Au(2,2,j,:)=dH_diab_Au(2,2,j,:)-deriv

    tmp1=dexp(gamma2*dist_O(i))
    deriv=A2*A3*gamma2*tmp1*tmp2/(1+A3*tmp1)**2
    H_diab(1,2)=H_diab(1,2)-A2/(1+A3*tmp1)+V_01_AuO_cutoff
    dH_diab_O(1,2,:)=dH_diab_O(1,2,:)+deriv
    dH_diab_Au(1,2,j,:)=dH_diab_Au(1,2,j,:)-deriv
  enddo

end subroutine gold_oxygen_coupling
!-----------------------------------------------------------------  

subroutine gold_nitrogen_coupling
  implicit none
  integer i,j
  real*8 tmp1,deriv(3),tmp2(3)

  do i=1,n_N_Au
    j=list_N(i)
    tmp2=(x_N-x_Au(j,:))!/dist_N(i)
    call pbc_rij(tmp2)
    tmp2=tmp2/dist_N(i)

    tmp1=dexp(-beta0*dist_N(i))
    deriv=B0*tmp1*(-beta0)*tmp2
    H_diab(1,1)=H_diab(1,1)+B0*(tmp1-e_beta0_rc)
    dH_diab_N(1,1,:)=dH_diab_N(1,1,:)+deriv
    dH_diab_Au(1,1,j,:)=dH_diab_Au(1,1,j,:)-deriv

    tmp1=dexp(-2*beta1*(dist_N(i)-r1_AuN))
    deriv=B1*tmp1*(-2*beta1)*tmp2
    H_diab(2,2)=H_diab(2,2)+B1*(tmp1-e_2beta1_rc_r1)
    dH_diab_N(2,2,:)=dH_diab_N(2,2,:)+deriv
    dH_diab_Au(2,2,j,:)=dH_diab_Au(2,2,j,:)-deriv

    tmp1=dexp(-beta1*(dist_N(i)-r1_AuN))
    deriv=B1*tmp1*(-beta1)*tmp2
    H_diab(2,2)=H_diab(2,2) - 2*B1*cos_theta**2*(tmp1-e_beta1_rc_r1)
    dH_diab_N(2,2,:)=dH_diab_N(2,2,:)-2*cos_theta**2*deriv
    dH_diab_N(2,2,:)=dH_diab_N(2,2,:)-4*B1*cos_theta*(tmp1-e_beta1_rc_r1)*dcostheta_dxN
    dH_diab_O(2,2,:)=dH_diab_O(2,2,:)-4*B1*cos_theta*(tmp1-e_beta1_rc_r1)*(-dcostheta_dxN)
    dH_diab_Au(2,2,j,:)=dH_diab_Au(2,2,j,:)+2*cos_theta**2*deriv

    tmp1=dexp(gamma3*dist_N(i))
    deriv=B2*B3*gamma3*tmp1*tmp2/(1+B3*tmp1)**2
    H_diab(1,2)=H_diab(1,2)-B2/(1+B3*tmp1)+V_01_AuN_cutoff
    dH_diab_N(1,2,:)=dH_diab_N(1,2,:)+deriv
    dH_diab_Au(1,2,j,:)=dH_diab_Au(1,2,j,:)-deriv
  enddo

end subroutine gold_nitrogen_coupling
!-----------------------------------------------------------------  

subroutine oxygen_nitrogen_coupling
  implicit none
  real*8 tmp1,tmp2(3),deriv(3)

  tmp2=(x_N-x_O)!/r_NO
  call pbc_rij(tmp2)
  tmp2=tmp2/r_NO

  tmp1=dexp(-gamma0*(r_NO-r0_NO))
  deriv=2*F0*(1.d0-tmp1)*tmp1*gamma0*tmp2
  H_diab(1,1)=H_diab(1,1)+F0*(1.d0-tmp1)*(1.d0-tmp1)
  dH_diab_N(1,1,:)=dH_diab_N(1,1,:)+deriv
  dH_diab_O(1,1,:)=dH_diab_O(1,1,:)-deriv

  tmp1=dexp(-gamma1*(r_NO-r1_NO))
  deriv=2*F1*(1-tmp1)*tmp1*gamma1*tmp2
  H_diab(2,2)=H_diab(2,2)+F1*(1.d0-tmp1)*(1.d0-tmp1)
  dH_diab_N(2,2,:)=dH_diab_N(2,2,:)+deriv
  dH_diab_O(2,2,:)=dH_diab_O(2,2,:)-deriv

end subroutine oxygen_nitrogen_coupling
!-----------------------------------------------------------------  

subroutine gold_gold_coupling
  implicit none
  integer i,j,k
  real*8 tmp(3),rij(3)

  V_Au_Au=0.d0
  dV_Au_Au=0.d0
  do i=1,n_Au
    do j=1,n_Au_Au(i)
      k=list_Au(i,j)
      rij=x_Au(i,:)-x_Au(k,:)
      call pbc_rij(rij)
      rij=rij-rij_0(:,i,k)
      tmp=matmul(A_mat(:,:,i,k),rij)
      V_Au_Au=V_Au_Au+0.5d0*sum(rij*tmp)
      dV_Au_Au(i,:)=dV_Au_Au(i,:)+tmp
      dV_Au_Au(k,:)=dV_Au_Au(k,:)-tmp
    enddo
  enddo

end subroutine gold_gold_coupling
!-----------------------------------------------------------------  

subroutine image_potential
  implicit none
  real*8 tmp1,tmp2

  tmp1=dsqrt(Csq+(z_COM-z_image)**2)

  H_diab(2,2)=H_diab(2,2)-D/tmp1
  dH_diab_N(2,2,3)=dH_diab_N(2,2,3)+(D*(z_COM-z_image)/tmp1**3*(mass_N/(mass_N+mass_O)))
  dH_diab_O(2,2,3)=dH_diab_O(2,2,3)+(D*(z_COM-z_image)/tmp1**3*(mass_O/(mass_N+mass_O)))

end subroutine image_potential
!-----------------------------------------------------------------  

subroutine NO_values
  implicit none
  real*8 r12(3)

  NO_com=(mass_N*x_N+mass_O*x_O)/(mass_N+mass_O)

  r12=x_O-x_N

  r_NO=dsqrt(sum(r12*r12))
  cos_theta=(x_O(3)-x_N(3))/r_NO
  z_COM=NO_COM(3)

  dcostheta_dxN=-(cos_theta)/(r_NO**2)*(-r12)
  dcostheta_dxN(3)=dcostheta_dxN(3)-1.d0/r_NO

end subroutine NO_values
!-----------------------------------------------------------------  

subroutine compute_list
  !! Output - n_N_Au,n_O_Au,list_near_*,dist_near_*
  !! finds the Au atom number,distance within r_cutoff of N,Au
  implicit none
  integer i,j
  real*8 d1,d2
  integer c1,c2

  call system_clock(c1)

  n_N_Au=0;n_O_Au=0
  do i=1,n_Au
    !! N
    d1=distance(i,-1)
    if(d1<r_cutoff) then
      n_N_Au=n_N_Au+1
      list_N(n_N_Au)=i
      dist_N(n_N_Au)=d1
    endif

    !! O
    d1=distance(i,-2)
    if(d1<r_cutoff) then
      n_O_Au=n_O_Au+1
      list_O(n_O_Au)=i
      dist_O(n_O_Au)=d1
    endif
  enddo

  if(n_N_Au>n_Au.or.n_O_Au>n_Au) then
    write(6,*) "problem in compute_list"
    write(6,*) n_Au,n_N_Au,n_O_Au
    stop
  endif

  call system_clock(c2)
  tim_comp_list=tim_comp_list+(c2-c1)/real(cnt_rate)

end subroutine compute_list
!-----------------------------------------------------------------  

function distance(k1,k2)
  implicit none
  integer,intent(in)::k1,k2
  real*8 distance,dx,dy,dz
  real*8 r12(3)

  if(k1<=0.or.k2==0.or.k2<-2.or.k1>n_Au.or.k2>n_Au) then
    write(6,*) "problem in input of distance"
    stop
  endif

  if(k2==-1) then
    r12=x_Au(k1,:)-x_N(:)
    call pbc_rij(r12)
    distance=dsqrt(sum(r12*r12))
  endif

  if(k2==-2) then
    r12=x_Au(k1,:)-x_O(:)
    call pbc_rij(r12)
    distance=dsqrt(sum(r12*r12))
  endif

  if(k2>0) then
    r12=x_Au(k1,:)-x_Au(k2,:)
    call pbc_rij(r12)
    distance=dsqrt(sum(r12*r12))

!    dx=dabs(x_Au(k1,1)-x_Au(k2,1));if(dx>box_length(1)/2.d0)dx=box_length(1)-dx
!    dy=dabs(x_Au(k1,2)-x_Au(k2,2));if(dy>box_length(2)/2.d0)dy=box_length(2)-dy
!    dz=dabs(x_Au(k1,3)-x_Au(k2,3))
!    distance=dsqrt(dx*dx+dy*dy+dz*dz)
  endif

end function distance
!-----------------------------------------------------------------  

subroutine pbc_rij(rij)
  implicit none
  real*8,intent(inout) :: rij(3)
  integer i

  do i=1,2
    if(rij(i)>box_length(i)/2.d0) rij(i)=rij(i)-box_length(i)
    if(rij(i)<-box_length(i)/2.d0) rij(i)=rij(i)+box_length(i)
  enddo

end subroutine pbc_rij
!-----------------------------------------------------------------  

subroutine pbc
  !! Periodic boundary conditions
  !! Places any atom outside the box into the box
  implicit none
  integer i,j

  do i=1,n_Au
    do j=1,2
      if(x_Au(i,j)>box_pbc_max(j))x_Au(i,j)=x_Au(i,j)-box_length(j)
      if(x_Au(i,j)<box_pbc_min(j))x_Au(i,j)=x_Au(i,j)+box_length(j)
    enddo
  enddo

  do j=1,2
    if(NO_com(j)>box_pbc_max(j)) then
      x_N(j)=x_N(j)-box_length(j)
      x_O(j)=x_O(j)-box_length(j)
      call evaluate_variables(0)
      call evaluate_variables(1)
    endif

    if(NO_com(j)<box_pbc_min(j)) then
      x_N(j)=x_N(j)+box_length(j)
      x_O(j)=x_O(j)+box_length(j)
      call evaluate_variables(0)
      call evaluate_variables(1)
    endif
  
    !if(x_N(j)>box_pbc_max(j))x_N(j)=x_N(j)-box_length(j)
    !if(x_N(j)<box_pbc_min(j))x_N(j)=x_N(j)+box_length(j)
    !if(x_O(j)>box_pbc_max(j))x_O(j)=x_O(j)-box_length(j)
    !if(x_O(j)<box_pbc_min(j))x_O(j)=x_O(j)+box_length(j)
  enddo

end subroutine pbc
!-----------------------------------------------------------------  

subroutine pbc_check
  implicit none
  integer error
  integer i,j

  error=0

  do i=1,n_Au
    do j=1,2
      if(x_Au(i,j)>box_pbc_max(j))error=1
      if(x_Au(i,j)<box_pbc_min(j))error=1
    enddo
  enddo

  if(error==1) then
    write(6,*) "pbc check not passed"
    stop
  else
    write(6,*) "pbc check passed"
    stop
  endif

end subroutine pbc_check
!-----------------------------------------------------------------  

subroutine set_Au_lattice
  !!! 111 surface
  implicit none
  integer i
  integer layer,row,column,i_Au
  real*8 r0(3),x0,y0,z0

  x0=a_Au
  y0=a_Au*dsin(pi/3.d0)
  z0=a_Au*dsqrt(2.d0/3.d0)

  box_length(1)=n_columns*x0
  box_length(2)=n_rows*y0
  box_length(3)=(n_layers-1)*z0

  box_pbc_min=-box_length(1:2)/2.d0
  box_pbc_max= box_length(1:2)/2.d0

  write(6,*) "box_lengths="
  write(6,'(3f15.7)') box_length*1.d10
  write(6,*) 

  i_Au=0
  r0(1)=-box_length(1)/2.d0
  r0(2)=-box_length(2)/2.d0
  r0(3)=-box_length(3)

  do layer=1,n_layers
    do row=1,n_rows
      do column=1,n_columns
        i_Au=i_Au+1
        x_Au(i_Au,2)=r0(2);x_Au(i_Au,3)=r0(3)
        x_Au(i_Au,1)=r0(1)+(column-1)*x0
        layer_num(i_Au)=layer
      enddo
      r0(2)=r0(2)+y0
      r0(1)=r0(1)+a_Au*dcos(pi/3.d0)
    enddo
    r0(3)=r0(3)+z0
    r0(1)=-box_length(1)/2.d0!0.d0
    r0(2)=-box_length(2)/2.d0+layer*a_Au*dsqrt(1.d0/3.d0)
  enddo

  if(i_Au.ne.n_Au) then
    write(6,*) "problem in set_Au_lattice,i_Au.ne.n_Au",i_Au,n_Au
  endif

  !x_Au=x_Au-box_length(3)

  call pbc
  call gold_list
  call set_force_matrix

  open(200,file="gold_pos.out")
  write(200,'(i5)') n_Au
  write(200,*) 
  do i=1,n_Au
    write(200,'(A,3f15.7)')'Au ',x_Au(i,:)*1.d10!,layer_num(i)
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
    !if(i>1.and.layer_num(i)-layer_num(i-1)>0) write(200,*) 
  enddo
  write(200,*) 
  close(200)

end subroutine set_Au_lattice
!-----------------------------------------------------------------  

subroutine set_NO_position
  implicit none
  integer i
  real*8 zneq,zN0,vn0
  real*8 rnd
  real*8 orient,en,fr,m_com
  real*8 tmp1

  !! x,y coordinates
  do i=1,2
    call random_number(rnd)
    x_N(i)=2*(rnd-0.5d0)*a_Au
    x_O(i)=2*(rnd-0.5d0)*a_Au
    v_N(i)=0.d0
    v_O(i)=0.d0
  enddo

  !! orientation
  call random_number(rnd)
  orient=rnd*2*pi

  !! equilibrium geometry
  !! v=0
  z_com=10.d-10
  m_com=mass_O*mass_N/(mass_N+mass_O)
  fr=dsqrt(2*F0*gamma0**2/m_com)
  en=0.5*hbar*fr  !! v=0
  zneq=mass_O*r0_NO/(mass_O+mass_N)
  zn0=dsqrt(2*en/m_com/fr**2) * mass_O/(mass_N+mass_O)
  vn0=dsqrt(2*mass_O*en/(mass_N*(mass_N+mass_O)))

  !! vibrational phase (based on SHO)
  call random_number(rnd)
  rnd=rnd*2*pi

  x_N(3)=z_com+(zneq+zn0*dcos(rnd))*dsin(orient)
  x_N(2)=x_N(2)+(zneq+zn0*dcos(rnd))*dcos(orient)
  v_N(3)=vn0*dsin(rnd)*dsin(orient)
  v_N(2)=vn0*dsin(rnd)*dcos(orient)

  x_O(3)=((mass_N+mass_O)*z_com-mass_N*x_N(3))/mass_O
  x_O(2)=((mass_N+mass_O)*x_O(2)-mass_N*x_N(2))/mass_O
  v_O(3)=-mass_N*v_N(3)/mass_O
  v_O(2)=-mass_N*v_N(2)/mass_O

  v_N(3)=v_N(3)-dsqrt(2*E_transl/(mass_O+mass_N))
  v_O(3)=v_O(3)-dsqrt(2*E_transl/(mass_O+mass_N))

end subroutine set_NO_position
!-----------------------------------------------------------------  

subroutine set_NO_frozen
  !! initialized NO for equilibrium simulation
  implicit none

  x_N=0.d0;x_O=0.d0
  x_N(3)=15.d-10;x_O(3)=x_N(3)+r0_NO
  v_N=0.d0;v_O=0.d0

end subroutine set_NO_frozen
!-----------------------------------------------------------------  

subroutine gold_list
  implicit none
  integer i,j
  real*8 d1

  n_Au_Au=0
  do i=1,n_Au-1
    do j=i+1,n_Au
      if(i.ne.j) then
        d1=distance(i,j)
        if(d1<a_Au+0.1d-10) then
          n_Au_Au(i)=n_Au_Au(i)+1
          list_Au(i,n_Au_Au(i))=j
          dist_Au(i,n_Au_Au(i))=d1
        endif
      endif
    enddo
  enddo

end subroutine gold_list
!-----------------------------------------------------------------  

subroutine set_force_matrix
  implicit none
  integer i,j,k
  real*8 A_x(3,3) !! force matrix for atom 1 at origin, and 2 along x-axis
  real*8 A_tmp(3,3)
  real*8 rot_mat_y(3,3),rot_mat_z(3,3)  !! rotate first along y, then along z for general orientation
  real*8 cos_theta_z,sin_theta_z
  real*8 cos_theta_y,sin_theta_y
  real*8 rij(3),ll,rij_z(3)
  real*8 vec(3)

  A_x=0.d0
  A_x(1,1)=beta+gamma
  A_x(2,2)=beta-gamma
  A_x(3,3)=alpha

  do i=1,n_Au
    do j=1,n_Au_Au(i)
      k=list_Au(i,j)
      if(i==k) then
        write(6,*)"problem in set_force_matrix",i,k
        stop
      endif
      rij=x_Au(k,:)-x_Au(i,:)
      call pbc_rij(rij)
      ll=dsqrt(sum(rij*rij))

      vec=0.d0;vec(1)=ll

      !! rotation along z
      sin_theta_z=rij(2)/ll
      cos_theta_z=dsqrt(1-sin_theta_z**2)
      rot_mat_z=0.d0;rot_mat_z(3,3)=1.d0
      rot_mat_z(1,1)=cos_theta_z;rot_mat_z(1,2)=-sin_theta_z
      rot_mat_z(2,1)=sin_theta_z;rot_mat_z(2,2)=cos_theta_z

      A_tmp=matmul(rot_mat_z,matmul(A_x,transpose(rot_mat_z)))
      vec=matmul(rot_mat_z,vec)

      !! rotation along y
      rij_z(1)=dsqrt(rij(1)**2+rij(3)**2);rij_z(2)=rij(2);rij_z(3)=0.d0
      cos_theta_y=rij(1)/dsqrt(rij(1)**2+rij(3)**2)
      sin_theta_y=rij(3)/dsqrt(rij(1)**2+rij(3)**2)
      rot_mat_y=0.d0;rot_mat_y(2,2)=1.d0
      rot_mat_y(1,1)=cos_theta_y;rot_mat_y(1,3)=-sin_theta_y
      rot_mat_y(3,1)=sin_theta_y;rot_mat_y(3,3)=cos_theta_y

      A_tmp=matmul(rot_mat_y,matmul(A_tmp,transpose(rot_mat_y)))
      vec=matmul(rot_mat_y,vec)

      A_mat(:,:,i,k)=A_tmp

      rij_0(:,i,k)=-rij

    enddo
  enddo

end subroutine set_force_matrix
!-----------------------------------------------------------------  

subroutine set_parameters
  !! E in J and distances in m
  implicit none
  real*8 conv

  conv=1.d3/av!1.66d-21

  r_cutoff=10.d-10
  mass_O=16.d0*amu2kg
  mass_N=14.d0*amu2kg
  mass_Au=197.d0*amu2kg

  a_Au=2.95d-10

  !! H_00
  A0=457095.d0*conv
  alpha0=3.7594/1.d-10
  B0=30707*conv
  beta0=3.0082/1.d-10
  F0=638.5*conv
  gamma0=2.743/1.d-10
  r0_NO=1.15077*1.d-10
  
  e_alpha0_rc=dexp(-alpha0*r_cutoff)
  e_beta0_rc =dexp(-beta0*r_cutoff)

  !! H_11
  A1=A0
  alpha1=alpha0
  B1=24.056*conv
  beta1=1.9649/1.d-10
  r1_AuN=2.3491*1.d-10
  D=347.22*conv*1.d-10
  C=1.2423*1.d-10
  z_image=1.1536d-10
  F1=495.98*conv
  gamma1=2.4890/1.d-10
  r1_NO=1.2904d-10
  phi=511.37*conv
  Ea=-0.67540*conv

  phi_Ea=phi-Ea
  e_alpha1_rc=dexp(-alpha1*r_cutoff)
  e_beta1_rc_r1=dexp(-beta1*(r_cutoff-r1_AuN))
  e_2beta1_rc_r1=dexp(-2*beta1*(r_cutoff-r1_AuN))
  Csq=C*C

  !! H_01
  A2=11.842*conv
  A3=0.0061803
  gamma2=1.3693/1.d-10
  B2=50*conv
  B3=0.0047213
  gamma3=2.0194/1.d-10

  V_01_AuO_cutoff=A2/(1+A3*dexp(gamma2*r_cutoff))
  V_01_AuN_cutoff=B2/(1+B3*dexp(gamma3*r_cutoff))

  !! Au-Au
  alpha=-4.94d0!*conv*1.d-3
  beta=17.15d0!*conv*1.d-3
  gamma=19.40d0!*conv*1.d-3

end subroutine set_parameters
!-----------------------------------------------------------------  

subroutine stochastic_force(delr,delv,dt)
  !! Allen Tildesley, page 262
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(n_Au,3),delv(n_Au,3)!f(nclass)
  integer i,j
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  !! Assuming all Au have constant mass

  gdt=gamma_fric*dt
  sig_r=dt*dsqrt(kb*temperature/mass_Au(1) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
  sig_v=dsqrt(kb*temperature/mass_Au(1)*(1-dexp(-2*gdt)))
  sig_rv=(dt*kb*temperature/mass_Au(1)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient


  do i=1,n_Au
    do j=1,3
      call gaussian_random_number(rnd1)
      call gaussian_random_number(rnd2)
      delr(i,j)=sig_r*rnd1
      delv(i,j)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
    enddo
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

subroutine check_acceleration
  implicit none
  integer i,j,nflag
  real*8 delx,en_old,acc_sav_Au(n_Au,3),acc_sav_N(3),acc_sav_O(3)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  call set_Au_lattice
  call set_NO_position

  !x_Au(41,3)=x_Au(41,3)+0.1d-10

  do i=1,n_Au
    do j=1,3
      call random_number(rnd)
      x_Au(i,j)=x_Au(i,j)+(rnd*2-1.d0)*0.1d-10
    enddo
  enddo
!  do i=1,3
!    call random_number(rnd)
!    x_N(i)=(rnd*2-1.d0)*1.d-10
!    call random_number(rnd)
!    x_O(i)=(rnd*2-1.d0)*1.d-10
!  enddo

  call evaluate_variables(0)
  en_old=pot_en
  acc_sav_Au=acc_Au;acc_sav_N=acc_N;acc_sav_O=acc_O

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,3
    x_N(i)=x_N(i)+delx
    call evaluate_variables(0)
    acc_N(i)=-(pot_en-en_old)/delx/mass_N
    write(6,*)"Nitrogen, axes",i
    write(6,*)"Analytical acceleration =",acc_sav_N(i)
    write(6,*)"Numerical acceleration  =",acc_N(i)
    write(6,*)"Error =",(acc_N(i)-acc_sav_N(i))/acc_N(i)*100.d0
    write(6,*)
    x_N(i)=x_N(i)-delx
  enddo

  do i=1,3
    x_O(i)=x_O(i)+delx
    call evaluate_variables(0)
    acc_O(i)=-(pot_en-en_old)/delx/mass_O
    write(6,*)"Oxygen, axes",i
    write(6,*)"Analytical acceleration =",acc_sav_O(i)
    write(6,*)"Numerical acceleration  =",acc_O(i)
    write(6,*)"Error =",(acc_O(i)-acc_sav_O(i))/acc_O(i)*100.d0
    write(6,*)
    x_O(i)=x_O(i)-delx
  enddo

  do j=1,n_Au
    write(6,*)"Au atom ",j
    do i=1,3
      x_Au(j,i)=x_Au(j,i)+delx
      call evaluate_variables(0)
      acc_Au(j,i)=-(pot_en-en_old)/delx/mass_Au(j)
      write(6,*)"Axes ",i
      write(6,*)"Analytical acceleration =",acc_sav_Au(j,i)
      write(6,*)"Numerical acceleration  =",acc_Au(j,i)
      write(6,*)"Error =",(acc_Au(j,i)-acc_sav_Au(j,i))/acc_Au(j,i)*100.d0
      write(6,*)
      x_Au(j,i)=x_Au(j,i)-delx
    enddo
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

subroutine draw_NO_potential
  !! V00_NO,V11_NO
  implicit none
  integer i
  real*8 tmp1,V1,V2

  x_O=0.d0
  x_N=0.d0

  do i=1,100
    x_N(3)=0.8d-10+1.5d-10*(i-1)/99.d0
    call NO_values
    tmp1=dexp(-gamma0*(r_NO-r0_NO))
    V1=F0*(1.d0-tmp1)*(1.d0-tmp1)
    tmp1=dexp(-gamma1*(r_NO-r1_NO))
    V2=F1*(1.d0-tmp1)*(1.d0-tmp1)
    write(20,*) x_N(3),V1/wave_to_J,V2/wave_to_J
  enddo
  stop


end subroutine draw_NO_potential
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!----------------------------------------------------------

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!----------------------------------------------------------


subroutine compute_KE_matrix_dvr(KE_matrix,ngrid,delq,mass)
  !! computes KE matrix in DVR basis
  !! Appendix A of JCP 96, 1982 (1991)
  implicit none
  integer,intent(in) :: ngrid       !! size of DVR grid
  real*8,intent(inout) :: KE_matrix(ngrid,ngrid)
  real*8,intent(in) :: delq         !! step size in position of DVR basis
  real*8,intent(in) :: mass         !! mass
  integer i,j
  real*8 pi,hbar

  pi=dacos(-1.d0);hbar=1.05457266D-34

  KE_matrix=hbar**2/(2*mass*delq**2)
  do i=1,ngrid
    do j=1,ngrid
      KE_matrix(i,j)=KE_matrix(i,j)*(-1.d0)**(i-j)
      if(i==j) then
        KE_matrix(i,j)=KE_matrix(i,j)*pi**2/3.d0
      else
        KE_matrix(i,j)=KE_matrix(i,j)*2.d0/real(i-j)**2
      endif
    enddo
  enddo
end subroutine compute_KE_matrix_dvr
!----------------------------------------------------------

subroutine pes_Fig3
  implicit none
  integer i,j

  x_N(3)=1.6d-10
  x_N(1)=0.d0
  x_N(2)=0.8515916d-10

  x_O=x_N

  state=1

  do i=1,100
    do j=1,100
      x_N(3)=0.8d-10+(i-1)*5.2d-10/99.d0
      x_O(3)=x_N(3)+0.8d-10+(j-1)*1.2d-10/99.d0
      call evaluate_variables(0)
      write(10,'(4f15.7)') x_N(3)*1.d10,(x_O(3)-x_N(3))*1.d10,pot_en*av/1.d3,-si_adiab(2,1)**2
    enddo
    write(10,*)
  enddo
  stop

end subroutine pes_Fig3
!-----------------------------------------------------------------  

End Module mod_scatter
