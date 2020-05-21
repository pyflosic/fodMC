module fodmc_mod
  
contains
subroutine get_guess()

!program fodMC

! Copyright    2019    Kai Trepte
! 
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
! 
! http://www.apache.org/licenses/LICENSE-2.0
! 
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License






! Author: Kai Trepte
! Idea Alex Johnson: Introduce a connectivity matrix to describe how many bonds each atom 
!                    has to each other atom, create FODs accordingly
!                    Input: Connectivity pattern and number of UP and DN lone FODs
! Idea Jakob Krause: Initialize single bonds according to covalent radii of the atoms.
! 
! Version 1
! 12. November 2018
! NEW: No input for H atoms needed anymore
! Version 2
! 21. January 2019
! Change of user input towards a simpler description of the bond pattern. Idea: (center1-center2)-(UP_FODS_in_bond-DN_FODS_in_bond)  for all pairs of atoms which have a bond. Lone FODs are specified afterwards
! 22. January 2019
! Exclude explicit need for X-H bonds. They can be introduced, but automatically it will always be a (1-1) bond
! In addition, lone FODs are now only necessary to place for atoms which actually have lone FODs
! 23. January 2019
! This is a test: Introduce cutoff radius. Set bond order of all pairs of atoms to 1 (1-1). With that, only bonds which are not single bonds need to be specified.
! 05. February 2019
! Rotate double bonds in linear molecules against each other (e.g. Lewis structure of CO2)
! 08. February 2019
! More consistent treatment of initial points for lone-FODs (now scaled with (FOD_UP+FOD_DN)/2)
! Point charge dipole as output (excluding 1s core-FODs)
! Change initial assessment of single bonds via covalent radii of the elements (Thanks to Jakob Krause) - revoke changes from 23. January
! Initial lone FODs are place per atom assuming charge neutrality of the system. Taking the bond information into account.
! 11. February 2019
! Read number of core FODs as input. Necessary to get the initial lone-FOD setting universally applicable.
! 11. March 2019
! Simplify structured. Avoid double counting as much as possible. Few instances where this is possible right now
! Implemented PBCs (periodic boundary conditions). Read in a key word. If it is 'pbc', then read coordinates, cell vectors and bond matrix
! 18. March 2019
! Simplify the specification for different amount of UP and DN bond FODs. Now, the specification (45-31)-(3-2) will automatically assign the correct number of bond FODs between the atoms, regardless of their order
! 19. March 2019
! New determination whether atoms are in a planar or linear environment. Much clearer now
! 13. May 2019
! Introduce 'fix1s' option. Allows to place all 1s FODs at the atomic positions
! 15. May 2019
! Re-introduced the distribution of points on a sphere
! 24. May 2019
! Introduce structural motifs. Use these for atoms and core FODs. Thus, no distribution of 
! points on a sphere necessary for such FODs.
! Further, a Metropolis-like algorithm is now used for the rotation of core vs. valence FODs.

implicit none
! new array structure: 3 real     (coordinates of center)
!		       1 char	  (element symbols)
!                      1 integer  (shell number)
!                      1 integer  (points on the shells)
!                      3 real     (r, theta, phi of the point)
!                      3 real     (x, y, z       of the point)
type global_array
  real(8)              :: center_x_y_z(3)                         ! atomic positions
  character(len=17)    :: elements                                ! element symbols.
  integer              :: n_shells                                ! Number of shells
  integer, allocatable :: n_points(:)                             ! dimensions: n_shells
  real(8), allocatable :: r_theta_phi(:,:,:)                      ! dimensions: n_shells, n_points, 3
  real(8), allocatable :: point_x_y_z(:,:,:)                      ! dimensions: n_shells, n_points, 3
end type global_array
type(global_array), dimension(:), allocatable :: pos1_up          ! global_array(1) - first center, global_array(2) - second center, ...
type(global_array), dimension(:), allocatable :: pos1_dn      
type(global_array), dimension(:), allocatable :: pos2_up          ! 1 - before MC step, 2 - after MC step
type(global_array), dimension(:), allocatable :: pos2_dn      


! Used for assessment of input data and to write output to the screen
integer                        :: number_of_centers ! number of different centers/ different atoms
integer, allocatable           :: element_number(:) ! order number per element. Used for writing CLUSTER
integer, allocatable           :: pseudo_charge(:)  ! number of core-electrons which are not treated explicitely (pseudopotentials)
integer, allocatable           :: core_charge(:)    ! number of core-electrons in total 
real(8)                        :: charge, spin      ! numbers to store charge and spin of the system
real(8), allocatable           :: covalent_R(:)     ! list covalent radii
character(len=2), allocatable  :: element_symbol(:) ! element symbols
character(len=4), allocatable  :: element_count(:)  ! number of specific type of element
logical                        :: periodic          ! keywords determining whether PBCs shall be used or not
logical                        :: fix1s             ! keywords determining whether to fix the 1s FODs at the nuclei
real(8)                        :: cell_a(3), cell_b(3), cell_c(3)   ! unit cell vectors. In case of PBCs
real(8)                        :: len_a, len_b, len_c, V_tot        ! length of unit cell vectors and total volume of the unit cell. In case of PBCs
real(8)                        :: angle_ab, angle_bc,  angle_ac     ! angles between unit cell vectors. In case of PBCs
character(len=1000)            :: bla2
character(len=25)              :: formatBLA
character(len=25)              :: formatBLA2

! Reading database -> need junk
character(len=17)              :: junk, junk2
character(len=8)               :: units             ! either bohr or angstrom . 
real(8)                        :: units_factor      ! conversion factor to bohr. 

! Used for evaluation
real(8)              :: ave_dist1_up_core, ave_dist2_up_core  ! average distances or 1/r. Before and after MC step. UP channel. For core FODs
real(8)              :: ave_dist1_dn_core, ave_dist2_dn_core  ! average distances or 1/r. Before and after MC step. DN channel. For core FODs
real(8)              :: ave_dist1_up, ave_dist2_up            ! average distances or 1/r. Before and after MC step. UP channel. For all or valence FODs
real(8)              :: ave_dist1_dn, ave_dist2_dn            ! average distances or 1/r. Before and after MC step. DN channel. For all or valence FODs
real(8)              :: ave_dist1, ave_dist2                  ! average distances or 1/r. Before and after MC step. Globally
real(8)              :: tmp_dist                              ! temporary distance for evaluation using PBCs
integer              :: cycles                                ! maximum number of steps     
real(8)              :: step_size                             ! maximal step size per point
integer              :: a,b,c,d,e,f,g,h,i,j,t                 ! loop variable
real(8)              :: start_t, finish_t                     ! start and finish time
real(8)              :: d_bond                                ! distance to closest atom
real(8)              :: scale_r                               ! scaling factor for the determination of where to place multiple FODs
integer              :: counter_up, counter_dn, counter       ! counter. To count things
real(8)              :: full_rot(3,3)                         ! rotation matrix  rot_x(3,3), rot_y(3,3), rot_z(3,3)
real(8), parameter   :: pi = 3.14159265358979323846           ! define pi. It is as accurate as needed :)
real(8)              :: dip_x, dip_y, dip_z                   ! point charge dipole in x, y, z
real(8)              :: cent_x, cent_y, cent_z                ! center of molecule/cell. For dipole evaluation
! Used for the connectivity matrix
character(len=50), allocatable :: bond_pattern(:)             ! bond pattern to be read in (NOT including lone pairs)
character(len=15)    :: bond_centers, bond_fods               ! for each individual entry of centers and FODs
integer, allocatable :: con_mat(:,:)                          ! connectivity matrix between all atoms. Order N x N
integer, allocatable :: lone_fods(:,:)                        ! number of UP and DN lone FODs per atom
real(8)              :: bond_center(3)                        ! position of the bond center between two atoms
real(8)              :: bond_vector(3)                        ! vector describing the bond
real(8), allocatable :: tmp_vectors(:,:)                      ! bond vector in case of a planar molecule -> get normal vector by cross product
real(8)              :: tmp_vector(3)                         ! tmeporary vector
real(8)              :: perpendicular_vec(3)                  ! vector perpedicular to the bond vector
integer, allocatable :: bond_count_up(:)                      ! counting the number of bonds already assumed by an atom -> assign valence properly. 
integer, allocatable :: bond_count_dn(:)                      ! counting the number of bonds already assumed by an atom -> assign valence properly.
! Used for Metropolis
real(8)              :: rand_metro                            ! random number for optimization with metropolis algorithm 

! for planarity/linearity
integer, allocatable :: is_planar_linear(:)                   ! determine whether there are planar/linear elements in a molecule. Necessary for bond and lone FOD generation

! for random seed
integer                             :: values(1:8), k
integer, dimension(:), allocatable  :: seed
call date_and_time(values=values)
call random_seed(size=k)
allocate(seed(1:k))
seed(:) = values(8)
call random_seed(put=seed)
! end for random seed


!!!!!!!!!
! START !
!!!!!!!!!
!!!!!!!!!!!!!!!!!!
! Read in values !
!!!!!!!!!!!!!!!!!!
open(unit=17,file='system',status='old',action='read')

read(17,*) number_of_centers                                                    ! read number of centers and allocate the global arrays
allocate(pos1_up(number_of_centers))
allocate(pos2_up(number_of_centers))
allocate(pos1_dn(number_of_centers))
allocate(pos2_dn(number_of_centers))
allocate(element_number(number_of_centers))
allocate(pseudo_charge(number_of_centers))
allocate(core_charge(number_of_centers))

read(17,*) units, junk, junk2                                                   ! read in the unit to be used (angstrom or bohr). and whether PBCs shall be used or not. And whether to fix the 1s FODs at the origin
if ((units == 'angstrom') .or. (units == 'Angstrom')) then
  units_factor = 0.529177D0
else if ((units == 'bohr') .or. (units == 'Bohr')) then
  units_factor = 1.000000D0
else
  write(6,*) 'Specified unit is neither angstrom nor bohr. Please check the second line of the file "system"'
  stop
end if
!
! If PBCs and/or fixing 1s at the nuclear positions
!
if ((junk == 'pbc' .and. junk2 == 'fix1s') .or. (junk == 'fix1s' .and. junk2 == 'pbc')) then
  periodic = .true.
  fix1s = .true.
else if (junk == 'pbc') then
  periodic = .true.
  fix1s = .false.
  backspace(17)                                                                ! If one of the keywords is missing -> need to go one line back up
else if (junk == 'fix1s') then
  periodic = .false.
  fix1s = .true.
  backspace(17)                                                                ! If one of the keywords is missing -> need to go one line back up
else
  periodic = .false.
  fix1s = .false.
  backspace(17)                                                                ! If both keywords are missing -> need to go one line back up
end if

do a = 1, number_of_centers
  read(17,*) pos1_up(a)%elements, pos1_up(a)%center_x_y_z(1:3)                 ! coordinates of atoms and element symbols. Same for all global arrays
  pos2_up(a)%elements = pos1_up(a)%elements
  pos1_dn(a)%elements = pos1_up(a)%elements
  pos2_dn(a)%elements = pos1_up(a)%elements
  pos2_up(a)%center_x_y_z(1:3) = pos1_up(a)%center_x_y_z(1:3)
  pos1_dn(a)%center_x_y_z(1:3) = pos1_up(a)%center_x_y_z(1:3)
  pos2_dn(a)%center_x_y_z(1:3) = pos1_up(a)%center_x_y_z(1:3)
end do
if (periodic) then                                                             ! in case of periodic boundary conditions: Read cell vectors
  read(17,*)                                                                   ! skip line after the atomic coordinates.
  read(17,*) cell_a(1), cell_a(2), cell_a(3)
  read(17,*) cell_b(1), cell_b(2), cell_b(3)
  read(17,*) cell_c(1), cell_c(2), cell_c(3)
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR ATOMS ONLY                         !
! if n = 1 -> read from xx_database_xx.  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (number_of_centers == 1) then
  open(18,file='xx_database_xx',status='old',action='read')
  a = 1                                                                        ! index for the atoms. 1 in this case
  t = 0                                                                        ! counter to abort the next while loop if all required data is read in
  do while (t < 1)
    read(18,*) junk
    if (junk == pos1_up(a)%elements) then                                      ! if the element is found -> read in data
      backspace(18)                                                            ! go up one line
      read(18,*) junk, element_number(a), pseudo_charge(a), core_charge(a)     ! read in the element number, pseudo and core charge
      read(18,*) pos1_up(a)%n_shells, pos1_dn(a)%n_shells                      ! read in the number of shells, for spin up and spin down
     
      pos2_up(a)%n_shells = pos1_up(a)%n_shells                                ! define the number of shells for the second arrays
      pos2_dn(a)%n_shells = pos1_dn(a)%n_shells                                ! 

      allocate(pos1_up(a)%n_points(pos1_up(a)%n_shells))
      allocate(pos1_dn(a)%n_points(pos1_dn(a)%n_shells))
      allocate(pos2_up(a)%n_points(pos2_up(a)%n_shells))
      allocate(pos2_dn(a)%n_points(pos2_dn(a)%n_shells))

      ! store number of points for the different shells
      do b = 1, pos1_up(a)%n_shells                                            ! store the number of points for all shells, spin UP
        read(18,*) junk, pos1_up(a)%n_points(b)
        pos2_up(a)%n_points(b) = pos1_up(a)%n_points(b)
      end do
      do b = 1, pos1_dn(a)%n_shells                                            ! store the number of points for all shells, spin DN
        read(18,*) junk, pos1_dn(a)%n_points(b)
        pos2_dn(a)%n_points(b) = pos1_dn(a)%n_points(b)
      end do

      ! allocate all needed arrays
      allocate(pos1_up(a)%r_theta_phi(pos1_up(a)%n_shells,sum(pos1_up(a)%n_points(:)),3))
      allocate(pos2_up(a)%r_theta_phi(pos2_up(a)%n_shells,sum(pos2_up(a)%n_points(:)),3))
      allocate(pos1_dn(a)%r_theta_phi(pos1_dn(a)%n_shells,sum(pos1_dn(a)%n_points(:)),3))
      allocate(pos2_dn(a)%r_theta_phi(pos2_dn(a)%n_shells,sum(pos2_dn(a)%n_points(:)),3))
      allocate(pos1_up(a)%point_x_y_z(pos1_up(a)%n_shells,sum(pos1_up(a)%n_points(:)),3))
      allocate(pos2_up(a)%point_x_y_z(pos2_up(a)%n_shells,sum(pos2_up(a)%n_points(:)),3))
      allocate(pos1_dn(a)%point_x_y_z(pos1_dn(a)%n_shells,sum(pos1_dn(a)%n_points(:)),3))
      allocate(pos2_dn(a)%point_x_y_z(pos2_dn(a)%n_shells,sum(pos2_dn(a)%n_points(:)),3))

      ! go back to the first entry
      do b = 1, pos1_up(a)%n_shells+pos1_dn(a)%n_shells
        backspace(18)
      end do

      ! store the radii of for the different shells 
      ! INITIALIZE THE POINTS !
      do b = 1, pos1_up(a)%n_shells
        read(18,*) pos1_up(a)%r_theta_phi(b,1,1)
        pos1_up(a)%r_theta_phi(b,1,1) = pos1_up(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
        do c = 1, pos1_up(a)%n_points(b)                                             ! store all other values for this shell accordingly. Add theta and phi.
          pos1_up(a)%r_theta_phi(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
          pos2_up(a)%r_theta_phi(b,c,1:3) = pos1_up(a)%r_theta_phi(b,c,1:3)
          pos1_up(a)%point_x_y_z(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1) + pos1_up(a)%center_x_y_z(1), &
                                             & pos1_up(a)%center_x_y_z(2),     pos1_up(a)%center_x_y_z(3) /)
          pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
        end do
      end do
      do b = 1, pos1_dn(a)%n_shells
        read(18,*) pos1_dn(a)%r_theta_phi(b,1,1)
        pos1_dn(a)%r_theta_phi(b,1,1) = pos1_dn(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
        do c = 1, pos1_dn(a)%n_points(b)                                             ! store all other values for this shell accordingly. Add theta and phi.
          pos1_dn(a)%r_theta_phi(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
          pos2_dn(a)%r_theta_phi(b,c,1:3) = pos1_dn(a)%r_theta_phi(b,c,1:3)
          pos1_dn(a)%point_x_y_z(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1) + pos1_dn(a)%center_x_y_z(1), &
                                             & pos1_dn(a)%center_x_y_z(2),     pos1_dn(a)%center_x_y_z(3) /)
          pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
        end do
      end do
      t = t + 1                                                                ! once everything is done, terminate the while loop
    end if
  end do

!
! Put all 1s core FODs (index (1,1,1:3)) at the origin
!
  if (fix1s) then
    pos1_up(a)%r_theta_phi(1,1,1:3) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /) 
    pos1_up(a)%point_x_y_z(1,1,1:3) = pos1_up(a)%center_x_y_z(1:3)
    pos2_up(a)%r_theta_phi(1,1,1:3) = pos1_up(a)%r_theta_phi(1,1,1:3)
    pos2_up(a)%point_x_y_z(1,1,1:3) = pos1_up(a)%point_x_y_z(1,1,1:3)
    pos1_dn(a)%r_theta_phi(1,1,1:3) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /) 
    pos1_dn(a)%point_x_y_z(1,1,1:3) = pos1_dn(a)%center_x_y_z(1:3)
    pos2_dn(a)%r_theta_phi(1,1,1:3) = pos1_dn(a)%r_theta_phi(1,1,1:3)
    pos2_dn(a)%point_x_y_z(1,1,1:3) = pos1_dn(a)%point_x_y_z(1,1,1:3)
  end if




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FOR MOLECULES -> read core radii from xx_database_xx, rest from system !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else
! create connectivity matrix from bond pattern. Like (1-2)-(2-2) (1-4)-(1-1) (3-4)-(3-2)   ! atom 1 binds to atom 2 with 2-UP and 2DN FODS, 1 to atom 4 with 1-UP and 1-DN FOD and atom 3 to atom 4 with 3-UP and 2-DN FODs
! initialize con_mat and lone_fods with 0 everywhere                                       ! lone FODs will be initialized later
  allocate(con_mat(number_of_centers,number_of_centers))
  allocate(lone_fods(number_of_centers,2))
  allocate(bond_count_up(number_of_centers))
  allocate(bond_count_dn(number_of_centers))
  allocate(covalent_R(number_of_centers))

  con_mat(:,:)     = 0                                                         ! initialize all bonds to be 0
  lone_fods(:,:)   = 0                                                         ! initialize all lone FODs to be 0
  bond_count_up(:) = 1                                                         ! counting UP FODs in bonds for a given atom.
  bond_count_dn(:) = 1                                                         ! counting DN FODs in bonds for a given atom

  ! skip comment line
  read(17,*)
  !
  ! Read number of shells for all atoms
  !
  open(18,file='xx_database_xx',status='old',action='read')
  do a = 1, number_of_centers
    t = 0
    do while (t < 1)
      read(18,*) junk
      if (junk == pos1_up(a)%elements) then                                    ! if the element is found -> read number of shells first 
        backspace(18)                                                          ! go up one line
        read(18,*) junk, element_number(a), pseudo_charge(a), core_charge(a), covalent_R(a)    ! read in the element number, pseudo and core charge and the covalent radius
        read(18,*) pos1_up(a)%n_shells, pos1_dn(a)%n_shells                    ! number of shells, for spin up and spin down
        pos2_up(a)%n_shells = pos1_up(a)%n_shells                              ! define the number of shells for the other arrays
        pos2_dn(a)%n_shells = pos1_dn(a)%n_shells                              !
        allocate(pos1_up(a)%n_points(pos1_up(a)%n_shells))
        allocate(pos1_dn(a)%n_points(pos1_dn(a)%n_shells))
        allocate(pos2_up(a)%n_points(pos2_up(a)%n_shells))
        allocate(pos2_dn(a)%n_points(pos2_dn(a)%n_shells))
        t = t + 1                                                              ! terminate while loop
      end if
    end do
    rewind(18)
  end do
  !
  ! For bonds : simply add one FOD for all bonds. I.e. all bonds are initialized with single bonds
  !
  do a = 1, number_of_centers-1
    do b = a+1, number_of_centers
      if (periodic) then                                                       ! In a peridoic system, take cell vectors into account
        do c = 1, 3
          do d = 1, 3
            do e = 1, 3
              if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                           (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2)) < &
                           (covalent_R(a)+covalent_R(b)+0.20)*units_factor) then
                con_mat(a,b) = 1
                con_mat(b,a) = 1
              end if
            end do
          end do
        end do
      else                                                                     ! For a non-periodic system, just evalaute the distances
        if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:))**2)) < &
                  & ((covalent_R(a)+covalent_R(b)+0.20)*units_factor)) then
          con_mat(a,b) = 1
          con_mat(b,a) = 1
        end if
      end if
    end do
  end do
  !
  ! Read rest of the input pattern regarding the bonding information
  !
  j = 0
  counter = 0 
  loop5: do a = 1, 1000                                                        ! if there is more than one input line -> go through all of them
    read(17,fmt='(A)') bla2
    if (bla2(1:1) == '(') then                                                 ! bond pattern always starts with a '(' 
      j = j + 1                                                                ! count number of entry lines
      do b = 1, len_trim(bla2)+1                                               ! read each string that contains the bonding information. Add one space to read last entry
        if (bla2(b:b) == ' ') then                                             ! if a space if found -> delimiter
          counter = counter + 1                                                ! counter is increased
        end if
      end do
    else
      exit loop5
    end if
  end do loop5
  !
  ! Allocate array according to the counter
  !
  allocate(bond_pattern(counter))
  !
  ! go back to the first input line. j+1 because in the above loop we read the input line and then increase the counter. Thus, it is increase one time too much
  !
  do a = 1, j + 1
    backspace(17)
  end do
   
  counter = 0
  do a = 1, j                                                                  ! go through the number of input lines
    read(17,fmt='(A)') bla2
    c = 0
    do b = 1, len_trim(bla2)+1                                                 ! read each string that contains the bonding information. Add one space to read last entry
      if (bla2(b:b) == ' ') then                                               ! if a space if found -> delimiter
        counter = counter + 1                                                  ! counter is increased
        bond_pattern(counter) = bla2(c+1:b-1)                                  ! store bond pattern
        c = b                                                                  ! increase the index of the string in total. Read all entries individually
      end if
    end do
  end do

  !
  ! CONSTRUCT CON_MAT AND LONE_FODS
  !
  do b = 1, counter                                                            ! for each entry of the bond pattern, iterate through all entries
    do c = 1, len(trim(bond_pattern(b)))                                       ! for each individual entry
      if (bond_pattern(b)(c:c) == '-' .and. bond_pattern(b)(c-1:c-1) == ')' .and. bond_pattern(b)(c+1:c+1) == '(') then   ! delimiter between centers and FOD information
        read(bond_pattern(b)(2:c-2),fmt='(A)') bond_centers
        read(bond_pattern(b)(c+2:len(trim(bond_pattern(b)))-1),fmt='(A)') bond_fods
        !
        ! go through the bond centers and assign variables to them accordingly
        !
        do d = 1, len(trim(bond_centers))
          if (bond_centers(d:d) == '-') then                                   ! found delimiter -> read in the atomic index per center involved in the bond
            read(bond_centers(1:d-1),fmt='(I4)') e                             ! bond center 1
            read(bond_centers(d+1:len(trim(bond_centers))),fmt='(I4)') f       ! bond center 2
          end if
        end do
        !
        ! now, go through the corresponding number of bond FODs (UP, DN) and write con_mat
        !
        do d = 1, len(trim(bond_fods))
          if (bond_fods(d:d) == '-') then                                      ! found delimiter -> read in the number of UP and DN FODs in the bond and write con_mat
            read(bond_fods(1:d-1),fmt='(I4)') g                                ! UP FODs
            read(bond_fods(d+1:len(trim(bond_fods))),fmt='(I4)') h             ! DN FODs
          end if
        end do
        ! 
        ! WRITE CON_MAT. Smaller atom index should always get the number of UP FODs, like (a-b)-(c-d) with a<b and c=UP, d=DN. If pattern is different -> adjust
        !
        if ((e > f) .and. (g < h)) then                                        ! e.g. (46-30)-(1-2)  -> need 1UP and 2DN between these two atoms. Assignment accordingly
          con_mat(e,f) = h
          con_mat(f,e) = g
        else if ((e > f) .and. (g > h)) then                                   ! e.g. (75-31)-(2-1)  -> need 2UP and 1DN between these two atoms. Assignment accordingly
          con_mat(e,f) = h
          con_mat(f,e) = g
        else
          con_mat(e,f) = g
          con_mat(f,e) = h
        end if
      end if
    end do
  end do


! Construct lone FODs from bond information and nuclear charge. 
! Take into account all core charges 
! ATTENTION: This assumes charge neutrality per atom.
!
! Compute core charges
! Get nuclear charge
! Get total number of bond FODs in the bond to the atom
! Assumption: half of these bond FODs are from this atom 
! Thus, if nuclear_charge - core_charge - 1/2*bond_FODs != 0  -> add number of lone FODs
! Number of lone FODs = (nuclear_charge - core_charge - 1/2*bond_FODs)/2. For UP and DN. If odd number, add more UP FODs

  do b = 1, number_of_centers                                                   ! do for all atoms
    if (((sum(con_mat(b,:)) + sum(con_mat(:,b)))/2.0D0 + core_charge(b)) < element_number(b)) then
      lone_fods(b,1) = ceiling((element_number(b)-core_charge(b)-(sum(con_mat(b,:))+sum(con_mat(:,b)))/2.0D0)/2.0D0)            ! UP lone FODs. If uneven number of UP and DN -> give it more UP
      lone_fods(b,2) =   floor((element_number(b)-core_charge(b)-(sum(con_mat(b,:))+sum(con_mat(:,b)))/2.0D0)/2.0D0)            ! DN lone FODs. If uneven number of UP and DN -> give it more UP
    end if
  end do

! GET THE REST OF THE LONE FODs. User input!
! Idea: have a pattern like 1-(1-1), meaning that atom 1 has 1 UP and 1 DN lone FOD
! Read lone FODs per center. 
  loop8: do b = 1, number_of_centers                                            ! do for all atoms
    read(17,fmt='(A)') bla2                                                     ! read lone FOD pattern
    if (len_trim(bla2) == 0) then                                               ! if the string has 0 length -> empty space of input -> END of input
      exit loop8
    else
      loop9: do c = 1, len_trim(bla2)+1                                         ! go through the entire string
        if (bla2(c:c) == '-') then                                              ! delimiter between atomic index and lone FOD information
          read(bla2(1:c-1),fmt='(I4)') i                                        ! atomic index of the corresponding atom
        end if
        if (bla2(c:c) == '(') then                                              ! start of lone FOD information
          do d = c+1, len_trim(bla2)+1                                          ! go through the rest of the string
            if (bla2(d:d) == '-') then                                          ! delimiter between UP-lone and DN-lone FOD
              read(bla2(c+1:d-1),fmt='(I4)') lone_fods(i,1)                     ! UP lone FODs
              read(bla2(d+1:len_trim(bla2)-1),fmt='(I4)') lone_fods(i,2)        ! DN lone FODs
              exit loop9
            end if
          end do
        end if
      end do loop9
    end if
  end do loop8

  !
  ! OUTPUT ON SCREEN
  !
  write(6,*) 'Bonding information and lone FODs'
  do a = 1, number_of_centers
    c = 0
    d = 0
    write(bla2,fmt='(I4)') a
    do b = 1, number_of_centers
      if ((con_mat(a,b) /= 0) .or. (con_mat(b,a) /= 0)) then
        if (a < b) then
          c = c + con_mat(a,b)   ! spin UP in the bond
          d = d + con_mat(b,a)   ! spin DN in the bond
          write(formatBLA,fmt='(I4)')  con_mat(a,b)
          write(formatBLA2,fmt='(I4)') con_mat(b,a)
        else
          c = c + con_mat(b,a)   ! spin UP in the bond
          d = d + con_mat(a,b)   ! spin DN in the bond
          write(formatBLA,fmt='(I4)')  con_mat(b,a)
          write(formatBLA2,fmt='(I4)') con_mat(a,b)
        end if
        write(junk,fmt='(I4)') b
        write(6,*) trim(pos1_up(a)%elements),'(',trim(bla2),') - ',trim(pos1_up(b)%elements),'(',trim(junk),')     ',&
                   '       Bond pattern   ',trim(adjustl(formatBLA)),'-',trim(adjustl(formatBLA2))
      end if
    end do
    write(formatBLA,fmt='(I4)') c
    write(junk,fmt='(I4)') d
    write(6,*) trim(pos1_up(a)%elements),'(',trim(bla2),') : ',trim(adjustl(formatBLA)),' UP and ',trim(adjustl(junk)),' DN - Bond'
    write(formatBLA,fmt='(I4)') lone_fods(a,1)
    write(junk,fmt='(I4)') lone_fods(a,2)
    write(6,*) trim(pos1_up(a)%elements),'(',trim(bla2),') : ',trim(adjustl(formatBLA)),' UP and ',trim(adjustl(junk)),' DN - lone' 
    write(6,*) ' '
  end do
  !
  ! END OUTPUT
  !

  !
  ! Go through all atoms and store number of points
  !
  do a = 1, number_of_centers
    t = 0
    do while (t < 1)
      read(18,*) junk
      if (junk == pos1_up(a)%elements) then                                    ! if the element is found -> read radii
        read(18,*) 
        do b = 1, pos1_up(a)%n_shells                                          ! store the number of points for all shells, spin UP
          read(18,*) junk, pos1_up(a)%n_points(b)
          pos2_up(a)%n_points(b) = pos1_up(a)%n_points(b)
        end do
        do b = 1, pos1_dn(a)%n_shells                                          ! store the number of points for all shells, spin DN
          read(18,*) junk, pos1_dn(a)%n_points(b)
          pos2_dn(a)%n_points(b) = pos1_dn(a)%n_points(b)
        end do
        t = t + 1                                                              ! terminate while loop
      end if
    end do
    !
    ! go back to the first entry in order to read the radii
    !
    do b = 1, pos1_up(a)%n_shells+pos1_dn(a)%n_shells
      backspace(18)
    end do
    !
    ! Get number of valence FODs from con_mat and lone_fods 
    ! UP CHANNEL
    pos1_up(a)%n_points(pos1_up(a)%n_shells) = 0                               ! Initialize new number of valence UP
    do b = 1, number_of_centers
      if (a < b) then
        pos1_up(a)%n_points(pos1_up(a)%n_shells) = pos1_up(a)%n_points(pos1_up(a)%n_shells) + con_mat(a,b)
      else
        continue
      end if
    end do
    pos1_up(a)%n_points(pos1_up(a)%n_shells) = pos1_up(a)%n_points(pos1_up(a)%n_shells) + lone_fods(a,1)
    pos2_up(a)%n_points(pos1_up(a)%n_shells) = pos1_up(a)%n_points(pos1_up(a)%n_shells)

    ! DN CHANNEL
    pos1_dn(a)%n_points(pos1_dn(a)%n_shells) = 0                               ! Initialize new number of valence DN
    do b = 1, number_of_centers
      if (a > b) then
        pos1_dn(a)%n_points(pos1_dn(a)%n_shells) = pos1_dn(a)%n_points(pos1_dn(a)%n_shells) + con_mat(a,b)
      else
        continue
      end if
    end do
    pos1_dn(a)%n_points(pos1_dn(a)%n_shells) = pos1_dn(a)%n_points(pos1_dn(a)%n_shells) + lone_fods(a,2)
    pos2_dn(a)%n_points(pos1_dn(a)%n_shells) = pos1_dn(a)%n_points(pos1_dn(a)%n_shells)
    !
    ! allocate radii, theta, phi and xyz for the FODs
    !
    allocate(pos1_up(a)%r_theta_phi(pos1_up(a)%n_shells,sum(pos1_up(a)%n_points(:)),3))
    allocate(pos2_up(a)%r_theta_phi(pos2_up(a)%n_shells,sum(pos2_up(a)%n_points(:)),3))
    allocate(pos1_dn(a)%r_theta_phi(pos1_dn(a)%n_shells,sum(pos1_dn(a)%n_points(:)),3))
    allocate(pos2_dn(a)%r_theta_phi(pos2_dn(a)%n_shells,sum(pos2_dn(a)%n_points(:)),3))
    allocate(pos1_up(a)%point_x_y_z(pos1_up(a)%n_shells,sum(pos1_up(a)%n_points(:)),3))
    allocate(pos2_up(a)%point_x_y_z(pos2_up(a)%n_shells,sum(pos2_up(a)%n_points(:)),3))
    allocate(pos1_dn(a)%point_x_y_z(pos1_dn(a)%n_shells,sum(pos1_dn(a)%n_points(:)),3))
    allocate(pos2_dn(a)%point_x_y_z(pos2_dn(a)%n_shells,sum(pos2_dn(a)%n_points(:)),3))
    !
    ! store the radii of for the different shells 
    ! INITIALIZE THE POINTS ! Here for core FODs
    do b = 1, pos1_up(a)%n_shells-1
      read(18,*) pos1_up(a)%r_theta_phi(b,1,1)
      pos1_up(a)%r_theta_phi(b,1,1) = pos1_up(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
      do c = 1, pos1_up(a)%n_points(b)                                             ! store all other values for this shell accordingly. Add theta and phi.
        pos1_up(a)%r_theta_phi(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
        pos2_up(a)%r_theta_phi(b,c,1:3) = pos1_up(a)%r_theta_phi(b,c,1:3)
        pos1_up(a)%point_x_y_z(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1) + pos1_up(a)%center_x_y_z(1), &
                                           & pos1_up(a)%center_x_y_z(2),     pos1_up(a)%center_x_y_z(3) /)
        pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
      end do
    end do
    read(18,*)                                                                     ! ignore the valence value here
    do b = 1, pos1_dn(a)%n_shells-1
      read(18,*) pos1_dn(a)%r_theta_phi(b,1,1)
      pos1_dn(a)%r_theta_phi(b,1,1) = pos1_dn(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
      do c = 1, pos1_dn(a)%n_points(b)                                             ! store all other values for this shell accordingly. Add theta and phi.
        pos1_dn(a)%r_theta_phi(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
        pos2_dn(a)%r_theta_phi(b,c,1:3) = pos1_dn(a)%r_theta_phi(b,c,1:3)
        pos1_dn(a)%point_x_y_z(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1) + pos1_dn(a)%center_x_y_z(1), &
                                           & pos1_dn(a)%center_x_y_z(2),     pos1_dn(a)%center_x_y_z(3) /)
        pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
      end do
    end do
  rewind(18)                                                                        ! go back to the top of the database file. Read it out for the next atom
  end do
  ! 
  ! get valence FOD radii. Here for non-bonded atoms. Use xx_database_xx file and use the atomic radii
  !
  do a = 1, number_of_centers
    if (sum(con_mat(a,:)) == 0) then                                           ! No bonds
      t = 0
      do while (t < 1)
        read(18,*) junk
        if (junk == pos1_up(a)%elements) then                                  ! if the element is found -> read radii and number of points
          read(18,*)                                                           ! ignore the number of shells now
          ! UP CHANNEL
          b = pos1_up(a)%n_shells                                              ! store the different radii and number of points in arrays. ONLY VALENCE
          do c = 1, b-1
            read(18,*)                                                         ! skip the radii for the core FODs
          end do
          read(18,*) pos1_up(a)%r_theta_phi(b,1,1)
          pos1_up(a)%r_theta_phi(b,1,1) = pos1_up(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
          do c = 1, pos1_up(a)%n_points(b)                                             ! initialize all lone_fods as specified in the system file.  pos1_up(a)%n_points(b) = lone_fods(a,1)
            pos1_up(a)%r_theta_phi(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
            pos2_up(a)%r_theta_phi(b,c,1:3) = pos1_up(a)%r_theta_phi(b,c,1:3)
            pos1_up(a)%point_x_y_z(b,c,1:3) = (/ pos1_up(a)%r_theta_phi(b,1,1) + pos1_up(a)%center_x_y_z(1), &
                                               & pos1_up(a)%center_x_y_z(2),     pos1_up(a)%center_x_y_z(3) /)
            pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
          end do

          ! DN CHANNEL
          b = pos1_dn(a)%n_shells                                              ! store the different radii and number of points in arrays
          do c = 1, b-1
            read(18,*)                                                         ! skip the radii for the core FODs
          end do
          read(18,*) pos1_dn(a)%r_theta_phi(b,1,1)
          pos1_dn(a)%r_theta_phi(b,1,1) = pos1_dn(a)%r_theta_phi(b,1,1)*units_factor   ! convert the radius to the correct units
          do c = 1, pos1_dn(a)%n_points(b)                                             ! initialize all lone_fods as specified in the system file.  pos1_dn(a)%n_points(b) = lone_fods(a,2)
            pos1_dn(a)%r_theta_phi(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1), real(pi/2.0,8), real(0.0,8) /)
            pos2_dn(a)%r_theta_phi(b,c,1:3) = pos1_dn(a)%r_theta_phi(b,c,1:3)
            pos1_dn(a)%point_x_y_z(b,c,1:3) = (/ pos1_dn(a)%r_theta_phi(b,1,1) + pos1_dn(a)%center_x_y_z(1), &
                                               & pos1_dn(a)%center_x_y_z(2),     pos1_dn(a)%center_x_y_z(3) /)
            pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
          end do
          t = t + 1                                                            ! terminate while loop
        end if
      end do
      rewind(18)                                                               ! go back to the top of the database to read the next element
    !
    ! For bonded atoms: get valence FOD radii from molecular geometry. Evaluate distances to adjacent atoms
    ! use distance between atoms as radii. Only relevant for lone FODs
    else
      d_bond = 10.0
      do b = 1, number_of_centers
        if (a /= b) then
          ave_dist1 = sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:))**2))         ! take shortest distance to neighboring atoms -> use 1/2 of this as radius
          if (periodic) then                                                                          ! In a peridoic system, take cell vectors into account
            do c = 1, 3
              do d = 1, 3
                do e = 1, 3
                  if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                               (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2)) < ave_dist1) then
                    ave_dist1 = sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                         (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2))
                  end if
                end do
              end do
            end do
          end if

          if (ave_dist1 < d_bond) then
            d_bond = ave_dist1
            i = b                                                                                     ! store atom index
          end if
        end if
      end do
      if (pos1_up(a)%elements == 'H' .or. pos1_up(i)%elements == 'H') then                            ! if the distance was taken to a H atoms -> take 1.85*distance/2.0 = 0.875*distance as radius. See further below as well
        d_bond = d_bond*1.85D0
      end if
      
      ! INITIALIZE THE POINTS ! Here for valence UP-FODs
      b = pos1_up(a)%n_shells
      do c = 1, pos1_up(a)%n_points(b)
        pos1_up(a)%r_theta_phi(b,c,1:3) = (/ d_bond/2.0D0, real(pi/2.0,8), real(0.0,8) /)
        pos2_up(a)%r_theta_phi(b,c,1:3) = pos1_up(a)%r_theta_phi(b,c,1:3)
        pos1_up(a)%point_x_y_z(b,c,1:3) = (/ d_bond/2.0D0 + pos1_up(a)%center_x_y_z(1), &
                                           & pos1_up(a)%center_x_y_z(2), pos1_up(a)%center_x_y_z(3) /)
        pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
      end do
      ! INITIALIZE THE POINTS ! Here for valence DN-FODs
      b = pos1_dn(a)%n_shells
      do c = 1, pos1_dn(a)%n_points(b)
        pos1_dn(a)%r_theta_phi(b,c,1:3) = (/ d_bond/2.0D0, real(pi/2.0,8), real(0.0,8) /)
        pos2_dn(a)%r_theta_phi(b,c,1:3) = pos1_dn(a)%r_theta_phi(b,c,1:3)
        pos1_dn(a)%point_x_y_z(b,c,1:3) = (/ d_bond/2.0D0 + pos1_dn(a)%center_x_y_z(1), &
                                           & pos1_dn(a)%center_x_y_z(2), pos1_dn(a)%center_x_y_z(3) /)
        pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
      end do
    end if
  end do

!
! Put all 1s core FODs (index (1,1,1:3)) at the origin
!
  if (fix1s) then
    do a = 1, number_of_centers 
      pos1_up(a)%r_theta_phi(1,1,1:3) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
      pos1_up(a)%point_x_y_z(1,1,1:3) = pos1_up(a)%center_x_y_z(1:3)
      pos2_up(a)%r_theta_phi(1,1,1:3) = pos1_up(a)%r_theta_phi(1,1,1:3)
      pos2_up(a)%point_x_y_z(1,1,1:3) = pos1_up(a)%point_x_y_z(1,1,1:3)
      pos1_dn(a)%r_theta_phi(1,1,1:3) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
      pos1_dn(a)%point_x_y_z(1,1,1:3) = pos1_dn(a)%center_x_y_z(1:3)
      pos2_dn(a)%r_theta_phi(1,1,1:3) = pos1_dn(a)%r_theta_phi(1,1,1:3)
      pos2_dn(a)%point_x_y_z(1,1,1:3) = pos1_dn(a)%point_x_y_z(1,1,1:3)
    end do
  end if

end if

!!!!!!!!!!!!!!!!!!!
! Initialization. !
!!!!!!!!!!!!!!!!!!!
call cpu_time(start_t)
!
! DEFINE NUMBER OF MC CYCLES ACCORDING TO THE NUMBER OF ELECTRONS. MINIMUM IS 5000, MAXIMUM IS 50000
!
cycles = 0
do b = 1, size(pos1_up)
  cycles = cycles + sum(pos1_up(b)%n_points(:))*1000
  cycles = cycles + sum(pos1_dn(b)%n_points(:))*1000
end do
if (cycles < 5000) then
  cycles = 5000
else if (cycles > 50000) then
  cycles = 50000
end if
! DEFINE STEP SIZE
step_size = 0.001D0
! DEFINE SCALING FACTOR FOR MULTIPLE BONDS/LONE FODS
scale_r = 0.60D0

!!!write(6,*) ' '
!!!!write(6,fmt='(A,I7,2X,A,2X,F7.4,2X,A,2X,F7.4)') 'Cycles:',cycles,' Step_size:',step_size,' Scale factor: ',scale_r
!!!write(6,666) 'Cycles:',cycles,' Step_size:',step_size,' Scale factor: ',scale_r
!!!666 FORMAT (A9,I7,2X,A,2X,F7.4,2X,A,2X,F7.4)

close(unit=17) ! close 'system'
close(unit=18) ! close 'xx_database_xx'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If only one center is present -> Atomic guess creation !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (number_of_centers == 1) then
!
! If distribution of points on a sphere is wanted -> change number of MC steps and step size
!
  if (pos1_up(1)%elements == 'POS') then
    write(6,*) ' '
    write(6,*) '!!! POINTS ON A SPHERE WITH ', sum(pos1_up(1)%n_points(:)), ' UP and ',sum(pos1_dn(1)%n_points(:)),' DN'
    write(6,*) ' '
    write(6,*) 'How many Monte-Carlo steps?'
    read(5,*) cycles
    write(6,*) 'What step-size?'
    read(5,*) step_size
!
! ELSE: For an atomic guess generation
!
  else
    a = 1
    charge = 0.0
    spin = 0.0

    charge     = real(element_number(1) - pseudo_charge(1) - &
               & sum(pos1_up(1)%n_points(:)) - sum(pos1_dn(1)%n_points(:)),8)
    spin       = real(sum(pos1_up(1)%n_points(:)) - sum(pos1_dn(1)%n_points(:)),8)

    write(6,*) ' '
    write(6,*) '!!! ATOMIC GUESS WILL BE CREATED FOR  ', pos1_up(1)%elements, '  charge = ',charge,' spin = ',spin
    write(6,*) ' '
  end if

!!!!!!!!!!!!!!!!!!!!!
! CORE OPTIMIZATION !
!!!!!!!!!!!!!!!!!!!!!
! IDEA: Use structural motifs instead of individual optimization. Rotate these via MC (USE METROPOLIS!)
! NEW  struct_motif(n_point_tot,n_point,radius,position_point)
  !
  ! Get UP CHANNEL  
  !
  do b = 1, pos1_up(a)%n_shells                                              ! All shells !
    do c = 1, pos1_up(a)%n_points(b)
      call struct_motif(pos1_up(a)%n_points(b),c,pos1_up(a)%r_theta_phi(b,c,1),&
                        pos1_up(a)%center_x_y_z(1:3),pos1_up(a)%point_x_y_z(b,c,1:3))
      !
      ! invert motif if b is an even number
      !
      if (mod(b,2) == 0) then
        pos1_up(a)%point_x_y_z(b,c,1:3) = -1.0D0*pos1_up(a)%point_x_y_z(b,c,1:3)
      end if

      pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
    end do
    !
    ! Optimize core structure globally (via rotations) with respect to any smaller shell
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ROTATIONS                                              !
    ! Rotate current core against any lower core. UP Channel !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((b > 1) .and. (b <= pos1_up(a)%n_shells)) then                             ! any shell except the 1s. 
      ave_dist1_up = 100000000.0
      do t = 1, cycles
        call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
        do c = 1, pos1_up(a)%n_points(b)
          call rotate_pos(full_rot, pos1_up(a)%point_x_y_z(b,c,1:3), &
                          pos1_up(a)%center_x_y_z(1:3), pos2_up(a)%point_x_y_z(b,c,1:3))
        end do
        !
        ! 1/r calculation
        !
        ave_dist2_up = 0.0                                                         ! 1/r UP
        do c = 1, pos2_up(a)%n_points(b)
          do d = 1, b-1                                                            ! Use any smaller shell structure and optimize the higher shell according to the lower core structure
            do f = 1, pos1_up(a)%n_points(d)                                       ! Unchanged index, thus index 1
              ave_dist2_up  = ave_dist2_up + &
                     1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(b,c,:) - pos1_up(a)%point_x_y_z(d,f,:))**2))
            end do
          end do
        end do
        !
        ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
        ! How: ... Introduce random number rand_metro [0,1]
        !          if (ave_dist2_up < ave_dist1_up) -> take new config
        !          if (ave_dist2_up > ave_dist1_up) -> if (ave_dist2_up-ave_dist1_up) < rand_metro*step_size => still take new config. Else: Take old config.
        !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
        if (ave_dist2_up < ave_dist1_up) then
          pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
          pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
          ave_dist1_up = ave_dist2_up
        else
          call random_number(rand_metro)
          if ((ave_dist2_up-ave_dist1_up) < rand_metro*step_size*0.01) then              ! HERE: NOT 100% good. Moves stuff away from the minimum (e.g. Mg atom) ....
            pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
            pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
            ave_dist1_up = ave_dist2_up
          end if
        end if
      end do  ! end MC cycles
    end if    ! end if - shells
  end do      ! end UP channel
  !
  ! Get DN CHANNEL  
  !
  do b = 1, pos1_dn(a)%n_shells                                              ! All shells !
    do c = 1, pos1_dn(a)%n_points(b)
      call struct_motif(pos1_dn(a)%n_points(b),c,pos1_dn(a)%r_theta_phi(b,c,1),&
                        pos1_dn(a)%center_x_y_z(1:3),pos1_dn(a)%point_x_y_z(b,c,1:3))
      !
      ! invert motif if b is an odd number
      !
      if (mod(b,2) == 1) then
        pos1_dn(a)%point_x_y_z(b,c,1:3) = -1.0D0*pos1_dn(a)%point_x_y_z(b,c,1:3)
      end if
      pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
    end do
    !
    ! Optimize core structure globally (via rotations) with respect to any smaller shell
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ROTATIONS                                              !
    ! Rotate current core against any lower core. DN Channel !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((b > 1) .and. (b <= pos1_dn(a)%n_shells)) then                             ! any shell except the 1s. 
      ave_dist1_dn = 100000000.0
      do t = 1, cycles
        call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
        do c = 1, pos1_dn(a)%n_points(b)
          call rotate_pos(full_rot, pos1_dn(a)%point_x_y_z(b,c,1:3), &
                          pos1_dn(a)%center_x_y_z(1:3), pos2_dn(a)%point_x_y_z(b,c,1:3))
        end do
        !
        ! 1/r calculation
        !
        ave_dist2_dn = 0.0                                                         ! 1/r DN
        do c = 1, pos2_dn(a)%n_points(b)
          do d = 1, b-1                                                            ! Use any smaller shell structure and optimize the higher shell according to the lower core structure
            do f = 1, pos1_dn(a)%n_points(d)                                       ! Unchanged index, thus index 1
              ave_dist2_dn  = ave_dist2_dn + &
                     1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - pos1_dn(a)%point_x_y_z(d,f,:))**2))
            end do
          end do
        end do
        !
        ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
        ! How: ... Introduce random number rand_metro [0,1]
        !          if (ave_dist2_dn < ave_dist1_dn) -> take new config
        !          if (ave_dist2_dn > ave_dist1_dn) -> if (ave_dist2_dn-ave_dist1_dn) < rand_metro*step_size => still take new config. Else: Take old config.
        !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
        if (ave_dist2_dn < ave_dist1_dn) then
          pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
          pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
          ave_dist1_dn = ave_dist2_dn
        else
          call random_number(rand_metro)
          if ((ave_dist2_dn-ave_dist1_dn) < rand_metro*step_size*0.01) then
            pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
            pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
            ave_dist1_dn = ave_dist2_dn
          end if
        end if
      end do  ! end MC cycles
      !
      ! After optimizing the core structure of the DN channel -> optimize it via rotations with respect to the UP channel (only lower core structure as well, and same level in addition)
      ! Rotate all core structure below the current one AND the current one as well
      ! 
      ave_dist1 = 100000000.0
      do t = 1, cycles
        call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
        do d = 1, b                                                                ! all shell up to the current core level
          do c = 1, pos1_dn(a)%n_points(d)                                           ! rotate all dn points up to the current core level (keep symmetry between the lower core structures the same)
            call rotate_pos(full_rot, pos1_dn(a)%point_x_y_z(d,c,1:3), &
                            pos1_dn(a)%center_x_y_z(1:3), pos2_dn(a)%point_x_y_z(d,c,1:3))
          end do
        end do
    ! 1/r calculation
        ave_dist2 = 0.0                                                           ! 1/r between dn and up !!!
        do d = 1, b
          do c = 1, pos2_dn(a)%n_points(d)
            do e = 1, min(b,pos1_up(a)%n_shells)                                   ! This loop for the UP channel goes only to the current core level. If the current core level doesn't exist -> use the maximum available
              do f = 1, pos1_up(a)%n_points(e)
                if (pos2_dn(a)%elements(3:5) == 'ECP' .or. pos2_dn(a)%elements(4:6) == 'ECP') then   ! if ECP are used -> there is no core. Just evaluate the 1/r
                  ave_dist2 = ave_dist2 + &
                        1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(d,c,:) - pos1_up(a)%point_x_y_z(e,f,:))**2))
                else
                  if ((d == 1) .and. (e == 1)) then                                ! avoid evaluating 1s UP and 1s DN (largest contribution) -> makes optimization extremely inefficient !!
                  else
                    ave_dist2  = ave_dist2 + &
                           1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(d,c,:) - pos1_up(a)%point_x_y_z(e,f,:))**2))
                  end if
                end if
              end do
            end do
          end do
        end do
    ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
    ! How: ... Introduce random number rand_metro [0,1]
    !          if (ave_dist2 < ave_dist1) -> take new config
    !          if (ave_dist2 > ave_dist1) -> if (ave_dist2-ave_dist1) < rand_metro*step_size => still take new config. Else: Take old config.
    !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
        if (ave_dist2 < ave_dist1) then
          pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
          pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
          ave_dist1 = ave_dist2
        else
          call random_number(rand_metro)
          if ((ave_dist2-ave_dist1) < rand_metro*step_size*0.01) then
            pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
            pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
            ave_dist1 = ave_dist2
          end if
        end if
      end do  ! end MC cycles
    end if    ! end if - shells
  end do      ! end UP channel
! END NEW

! INCLUDE POS AT SOME POINT. CURRENTLY MISSING


  call cpu_time(finish_t)
  write(6,*) 'up-dn global   ', ave_dist1
  write(6,*) 'Total CPU time ',finish_t - start_t,'s'

!
! If distribution of points on a sphere
!
  if (pos1_up(1)%elements == 'POS') then
    open(unit=19,file='points_on_sphere.xyz',status='unknown',action='write')
    write (junk, '(I8)') size(pos1_up)+sum(pos1_up(1)%n_points(:))+sum(pos1_dn(1)%n_points(:))-1                     ! number of entries in the xyz file. Exclude 'atom'
    write(19,fmt='(A)') adjustl(junk)
    write(19,*) 'points on a sphere'
    d = 1
    do b = 1, pos1_up(d)%n_shells
      do c = 1, pos1_up(d)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'X',pos1_up(d)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
    do b = 1, pos1_dn(d)%n_shells
      do c = 1, pos1_dn(d)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'He',pos1_dn(d)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
    close(unit=19)

!
! For atomic guess generation
!
  else 
!
! GENERATE CLUSTER AND FRMORB FILE for NRLMOL
!
    charge = real(element_number(1) - pseudo_charge(1) - sum(pos1_up(1)%n_points(:)) - sum(pos1_dn(1)%n_points(:)),8)
    spin   = real(sum(pos1_up(1)%n_points(:)) - sum(pos1_dn(1)%n_points(:)),8)
    open(unit=19,file='CLUSTER',status='unknown',action='write')
    write(19,fmt='(A)') 'LDA-PW91*LDA-PW91            (DF TYPE EXCHANGE*CORRELATION)'
    write(19,fmt='(A)') 'NONE                         (TD, OH, IH, X, Y, XY, ... OR GRP)'
    write(19,fmt='(A)') '1                            (NUMBER OF ATOMS)'
    if ((pos1_up(1)%elements(3:5) == 'ECP') .or. (pos1_up(1)%elements(4:6) == 'ECP')) then
      write(19,fmt='(3(F13.8,2X),I3,1X,A)') pos1_up(1)%center_x_y_z(1:3)/units_factor, element_number(1), & 
                   & '  ECP (R, Z, Pseudopotential)'
    else
      write(19,fmt='(3(F13.8,2X),I3,1X,A)') pos1_up(1)%center_x_y_z(1:3)/units_factor, element_number(1), &
                   & '  ALL (R, Z, ALL-ELECTRON)'
    end if
    write(19,fmt='(2(F6.3,2X),19X,A)') charge, spin, '(NET CHARGE AND NET SPIN)'
    close(unit=19)

    open(unit=19,file='FRMORB',status='unknown',action='write')
    write(19,fmt='(I3,5X,I3)') sum(pos1_up(1)%n_points(:)), sum(pos1_dn(1)%n_points(:))
    d = 1
    do b = 1, pos1_up(d)%n_shells
      do c = 1, pos1_up(d)%n_points(b)
        write(19,fmt='(3(F13.8,2X))') pos1_up(d)%point_x_y_z(b,c,1:3)/units_factor
      end do
    end do
    do b = 1, pos1_dn(d)%n_shells
      do c = 1, pos1_dn(d)%n_points(b)
        write(19,fmt='(3(F13.8,2X))') pos1_dn(d)%point_x_y_z(b,c,1:3)/units_factor
      end do
    end do
    close(unit=19)

!
! GENERATE Nuc_FOD.xyz. For PyFLOSIC
!
    open(unit=19,file='Nuc_FOD.xyz',status='unknown',action='write')
    write (junk, '(I8)') size(pos1_up)+sum(pos1_up(1)%n_points(:))+sum(pos1_dn(1)%n_points(:))                       ! number of entries in the xyz file
    write(19,fmt='(A)') adjustl(junk)
    write(19,*) 'angstrom'
    if (pos1_up(1)%elements(2:2) == '_') then
      write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(1)%elements(1:1),pos1_up(1)%center_x_y_z(1:3)*0.529177D0/units_factor ! always in angstrom
    else
      write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(1)%elements(1:2),pos1_up(1)%center_x_y_z(1:3)*0.529177D0/units_factor
    end if
    d = 1
    do b = 1, pos1_up(d)%n_shells
      do c = 1, pos1_up(d)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'X',pos1_up(d)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
    do b = 1, pos1_dn(d)%n_shells
      do c = 1, pos1_dn(d)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'He',pos1_dn(d)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
    close(unit=19)

  end if




































!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! More than one center is present -> Molecular guess creation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Step-wise procedure:  1. get bonds
  !                       2. get lone pairs
  !                       3. optimize the core accordingly
  !
else                                                                           ! more than one atom -> molecular/periodic guess creation
  !
  ! WRITING STUFF TO SCREEN
  !
  allocate(element_symbol(number_of_centers))                                  ! maximum number of element symbols
  do a = 1, number_of_centers
    element_symbol(a) = 'X'                                                    ! dummy argument
  end do
  counter = 0
  do a = 1, size(pos1_up)
    c = 0
    do b = 1, number_of_centers                                                ! if element already in the list -> ignore
      if ((pos1_up(a)%elements(1:1) == element_symbol(b) .and. (pos1_up(a)%elements(2:2) == ' ')) &   ! avoid confusion betweenn e.g. Na and N
     .or. (pos1_up(a)%elements(1:2) == element_symbol(b))) then
        c = c + 1
      end if
    end do
    if (c == 0) then
      counter = counter + 1
      if ((pos1_up(a)%elements(2:2) == '_') .or. (pos1_up(a)%elements(2:2) == ' ')) then
        element_symbol(counter) = pos1_up(a)%elements(1:1)
      else if ((pos1_up(a)%elements(3:3) == '_') .or. (pos1_up(a)%elements(3:3) == ' ')) then
        element_symbol(counter) = pos1_up(a)%elements(1:2)
      end if
    end if
  end do
  allocate(element_count(counter))
  do a = 1, counter
    d = 0
    do b = 1, number_of_centers
      if (((pos1_up(b)%elements(1:1) == element_symbol(a)) .and. &
        & ((pos1_up(b)%elements(2:2) == ' ') .or. (pos1_up(b)%elements(2:2) == '_'))) .or. &       ! Element with a single symbol (like N)
        &  (pos1_up(b)%elements(1:2) == element_symbol(a))) then                                   ! Element with more than one symbol (like Ni)
        d = d + 1
      end if
    end do
    if (d == 1) then
      element_count(a) = ' '
    else
      write(element_count(a),fmt='(I4)') d
    end if
  end do

  charge     = 0.0
  spin       = 0.0
  do a = 1, size(pos1_up)
    charge     = charge + real(element_number(a) - pseudo_charge(a) - &
               & sum(pos1_up(a)%n_points(:)) - sum(pos1_dn(a)%n_points(:)),8)
    spin       = spin   + real(sum(pos1_up(a)%n_points(:)) - sum(pos1_dn(a)%n_points(:)),8)
  end do

  write(6,*) ' '
  if (periodic) then
    write(6,*) '!!! PERIODIC GUESS WILL BE CREATED FOR    ', &
                  & (trim(element_symbol(a)),trim(adjustl(element_count(a))),a=1,counter)
  else
    write(6,*) '!!! MOLECULAR GUESS WILL BE CREATED FOR     ', &
                  & (trim(element_symbol(a)),trim(adjustl(element_count(a))),a=1,counter)
  end if
  write(6,*) ' '
  deallocate(element_symbol)
  deallocate(element_count)
  !
  ! END WRITING STUFF TO SCREEN
  !



!
! Check whether atoms are in a planar or linear environment. Needed later for bonds and lone FODs
!
  allocate(is_planar_linear(size(pos1_up)))
  is_planar_linear(:) = 0                                                       ! initialize all to be 0 (neither planar nor linear)
  if (size(pos1_up) > 2) then                                                   ! ignore diatomics
    !
    ! Go through all atoms
    !
    do a = 1, size(pos1_up)
      counter = 0
      do b = 1, size(pos1_up)
        if (con_mat(a,b) /= 0) then
          !
          ! count the number of bond partners
          !
          counter = counter + 1
        end if
      end do
      !
      ! When an atom has exactly two bonds -> check whether it is in a linear environment or not -> dot product of the bond vectors needs to be -1 or 1 (angle either 0 or 180)
      ! If not -> planar
      !
      if (counter == 2) then
        counter = 0
        allocate(tmp_vectors(2,3))
        do b = 1, size(pos1_up)
          if (con_mat(a,b) /= 0) then
            counter = counter + 1
            tmp_vectors(counter,:) = (pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:))
            if (periodic) then                                                  ! In a peridoic system, take cell vectors into account
              do d = 1, 3
                do e = 1, 3
                  do f = 1, 3
                    if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                 (d-2)*cell_a(:) + (e-2)*cell_b(:) + (f-2)*cell_c(:))**2)) < &
                                 (sqrt(sum(tmp_vectors(counter,:)**2)))) then
                      tmp_vectors(counter,:) = pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                               (d-2)*cell_a(:) + (e-2)*cell_b(:) + (f-2)*cell_c(:)
                    end if
                  end do
                end do
              end do
            end if 
          end if
        end do
        t = 0
        do c = 1, counter-1
          do d = c+1, counter
            if (abs(dot_product(tmp_vectors(c,:),tmp_vectors(d,:)/&
               (sqrt(sum(tmp_vectors(c,:)**2))*sqrt(sum(tmp_vectors(d,:)**2))))) > 0.98) then
              t = t + 1
            end if
          end do
        end do
        if (t == (counter*(counter-1))/2) then
          is_planar_linear(a) = 2
        else
          is_planar_linear(a) = 1
        end if
        deallocate(tmp_vectors)
      end if

      !
      ! For more bonds
      !
      if (counter > 2) then
        allocate(tmp_vectors(counter+(counter*(counter-1)/2),3))                ! store bond vectors (counter) and perpendicular vectors (counter*(counter-1)/2)
        c = 0                                                                   ! counter for the determination of how many bond vectors are needed
        t = 0                                                                   ! checking counter. If it always increased for all pairs of bonds -> do something
        do b = 1, size(pos1_up)
          !
          ! store all bond vectors
          !
          if (con_mat(a,b) /= 0) then
            c = c + 1  
            tmp_vectors(c,:) = (pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:))
            if (periodic) then                                                  ! In a peridoic system, take cell vectors into account
              do d = 1, 3
                do e = 1, 3
                  do f = 1, 3
                    if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                 (d-2)*cell_a(:) + (e-2)*cell_b(:) + (f-2)*cell_c(:))**2)) < &
                                 (sqrt(sum(tmp_vectors(c,:)**2)))) then
                      tmp_vectors(c,:) = pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                         (d-2)*cell_a(:) + (e-2)*cell_b(:) + (f-2)*cell_c(:)
                    end if
                  end do
                end do
              end do
            end if
          end if
        end do
        ! Cross product for each pair of bond vectors -> perpendicular vectors
        ! Need the angles between these perpendicular vectors 
        ! if they all are aligned/counteralinged (angle = 0 or 180 degree) -> planar !
        do b = 1, size(pos1_up)
          if (con_mat(a,b) /= 0) then
            f = 0                                                               ! variable to store the perpendicular vectors
            do c = 1, counter-1
              do d = c+1, counter
                f = f + 1 
                !
                ! cross product
                !
                tmp_vectors(counter+f,1) = tmp_vectors(c,2)*tmp_vectors(d,3) - tmp_vectors(c,3)*tmp_vectors(d,2)
                tmp_vectors(counter+f,2) = tmp_vectors(c,3)*tmp_vectors(d,1) - tmp_vectors(c,1)*tmp_vectors(d,3)
                tmp_vectors(counter+f,3) = tmp_vectors(c,1)*tmp_vectors(d,2) - tmp_vectors(c,2)*tmp_vectors(d,1)
              end do
            end do
          end if 
        end do
        ! 
        ! get dot product. If the perpendicular vectors are aligned -> their angle is either 0 or 180 -> their cosine is either -1 or 1 -> abs(DOT) must be close to 1 
        ! cos alpha = dot/(len_a*len_b)
        !
        do c = 1, counter-1
          do d = c+1, counter
            if (abs(dot_product(tmp_vectors(c+f,:),tmp_vectors(d+f,:)/&
               (sqrt(sum(tmp_vectors(c+f,:)**2))*sqrt(sum(tmp_vectors(d+f,:)**2))))) > 0.98) then
              t = t + 1
            end if
          end do 
        end do
        !
        ! check whether all bonds are in plane
        !
        if (t == (counter*(counter-1))/2) then
          is_planar_linear(a) = 1
        end if
        !
        ! deallocate for next atom
        !
        deallocate(tmp_vectors)
      end if     
    end do
  end if
!
! DONE planarity/linearity check
!



!
! Find centers between atoms. Place as many FODs around these as defined in con_mat. 
! Single bond:  Place 1 UP and 1 DN FOD at center (for X-H, place closer to H (instead of 0.5*bond length, use 0.85*bond length))
! Double bonds: Place 2 UP and 2 DN FODs perpendicular to the bond axis (using optimization routines based on MC)
! Triple bonds: Place 3 UP and 3 DN FODs in a triangle around the bond axis (using optimization routines based on MC)
! Quadruple bonds: 4 UP and 4 DN in a square
! ...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! 1. get bonds !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do a = 1, size(pos1_up)                                                      ! for each atom
    if (sum(con_mat(a,:)) == 0) then                                           ! In case there are no bonds to that atoms -> skip this section
      continue
    else
      do b = 1, size(pos1_up)                                                  ! for each other atom.
        if (con_mat(a,b) /= 0) then                                            ! for atoms which are bonded
          !
          ! get bond center and vector describing the bond
          !
          bond_center(:) = (pos1_up(a)%center_x_y_z(:) + pos1_up(b)%center_x_y_z(:))/2.0D0
          bond_vector(:) =  pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) ! vector from b to a          
          if (periodic) then                                                        ! In a peridoic system, take cell vectors into account
            do c = 1, 3
              do d = 1, 3
                do e = 1, 3
                  if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                               (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2)) < &
                               (sqrt(sum(bond_vector(:)**2)))) then
                    bond_center(:) = (pos1_up(a)%center_x_y_z(:) + pos1_up(b)%center_x_y_z(:) + &
                                     (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))/2.0D0
                    bond_vector(:) = pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                     (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:)
                  end if
                end do
              end do
            end do
          end if
          !
          ! find a vector which is perpendicular to the bond vector
          !
          call per_vector(bond_vector(:),pos1_up(a)%center_x_y_z(:),perpendicular_vec(:))

          !!!!!!!!!!!!!!!!!!!!
          ! SINGLE BOND FODs !
          !!!!!!!!!!!!!!!!!!!!
          if (con_mat(a,b) == 1) then
            if (a < b) then                                                    ! Add UP FOD
              c = pos1_up(a)%n_shells                                          ! valence shell UP
              if (pos1_up(a)%elements(1:1) == 'H') then                        ! for H bonds -> place FOD closer to H (0.85*bond length)
                pos1_up(a)%point_x_y_z(c,bond_count_up(a),:) = pos1_up(b)%center_x_y_z(:) + 0.85D0*bond_vector(:)
                bond_count_up(a) = bond_count_up(a) + 1
              else if (pos1_up(b)%elements(1:1) == 'H') then
                pos1_up(a)%point_x_y_z(c,bond_count_up(a),:) = pos1_up(b)%center_x_y_z(:) + 0.15D0*bond_vector(:)
                bond_count_up(a) = bond_count_up(a) + 1
              else                                                             ! If not H -> place at the center of the bond
                pos1_up(a)%point_x_y_z(c,bond_count_up(a),:) = bond_center(:)
                bond_count_up(a) = bond_count_up(a) + 1
              end if
            else                                                               ! Add DN FOD
              c = pos1_dn(a)%n_shells                                          ! valence shell DN
              if (pos1_up(a)%elements(1:1) == 'H') then                        ! for H bonds -> place FOD closer to H (0.85*bond length)
                pos1_dn(a)%point_x_y_z(c,bond_count_dn(a),:) = pos1_dn(b)%center_x_y_z(:) + 0.85D0*bond_vector(:)
                bond_count_dn(a) = bond_count_dn(a) + 1
              else if (pos1_up(b)%elements(1:1) == 'H') then
                pos1_dn(a)%point_x_y_z(c,bond_count_dn(a),:) = pos1_dn(b)%center_x_y_z(:) + 0.15D0*bond_vector(:)
                bond_count_dn(a) = bond_count_dn(a) + 1
              else                                                             ! If not H -> place at the center of the bond
                pos1_dn(a)%point_x_y_z(c,bond_count_dn(a),:) = bond_center(:)
                bond_count_dn(a) = bond_count_dn(a) + 1
              end if
            end if

            
          !!!!!!!!!!!!!!!!!!!!!!
          ! MULTIPLE BOND FODs !
          !!!!!!!!!!!!!!!!!!!!!!
          ! Take perpendicular vector to the bond vector -> is perpendicular to the bond axis
          ! Use this vector (and a given radius) to place FODs perpendicular to the bond (vector applied to the bond center!)
          ! Rotate these points around the bond vector using MC to minimize their 1/r contribution
          !!!!!!!!!!!!!!
          ! UP CHANNEL !
          !!!!!!!!!!!!!! 
          else if (con_mat(a,b) > 1) then
            if (a < b) then                                                    ! Add UP FOD 
              do f = 1, con_mat(a,b)
                c = pos1_up(a)%n_shells                                        ! valence shell UP
                pos1_up(a)%point_x_y_z(c,bond_count_up(a),:) = &
                & bond_center(:) + perpendicular_vec(:)*pos1_up(a)%r_theta_phi(c,bond_count_up(a),1)*scale_r           ! initial pos of UP channel. To be rotated around the bond axis. scale_r is a factor
                bond_count_up(a) = bond_count_up(a) + 1
              end do
              !
              ! Rotate the point around the bond vector (bond axis).
              !
              ave_dist1 = 100000.0
              do d = 1, cycles
                do f = 1, con_mat(a,b)
                  pos2_up(a)%point_x_y_z(c,bond_count_up(a)-f,:) = pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)      ! store original positions
                  call rotate_around_axis(bond_vector,bond_center,pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:),step_size)
                end do
                !
                ! 1/r calculation.
                !
                ave_dist2 = 0.0
                do e = 1, con_mat(a,b)-1
                  do f = e+1, con_mat(a,b)
                    ave_dist2 = ave_dist2 + 1.0D0/sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - &
                                                        & pos1_up(a)%point_x_y_z(c,bond_count_up(a) - f,:))**2))
                  end do
                end do
                !
                ! In addition, evaluate the 1/r to all atoms the corresponding atom is bonded to (and to the atom itself) -> better symmetry
                !
                do g = 1,size(pos1_up) 
                  if ((con_mat(a,g)) /= 0 .or. (con_mat(g,a) /= 0) .or. (g == a)) then
                    do e = 1, con_mat(a,b)
                      tmp_dist = sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - &
                                                            & pos1_up(g)%center_x_y_z(:))**2))
                      if (periodic) then                                                     ! In a peridoic system, take cell vectors into account
                        do f = 1, 3
                          do h = 1, 3
                            do i = 1, 3
                              if (sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - pos1_up(g)%center_x_y_z(:) + &
                                           (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2)) < tmp_dist) then
                                tmp_dist = sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:)-pos1_up(g)%center_x_y_z(:)+&
                                           (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2))
                              end if
                            end do
                          end do
                        end do
                      end if
                      ave_dist2 = ave_dist2 + 1.0D0/tmp_dist
                    end do
                  end if
                  !
                  ! In addition, use the atoms which the bonded atom is bonded to as well -> more robust. I.e. include all atoms in con_mat(g,:)
                  !
                  do j = 1,size(pos1_up)
                    if ((con_mat(g,j)) /= 0 .or. (con_mat(j,g) /= 0)) then
                      do e = 1, con_mat(a,b)
                        tmp_dist = sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - &
                                                              & pos1_up(j)%center_x_y_z(:))**2))
                        if (periodic) then                                                     ! In a peridoic system, take cell vectors into account
                          do f = 1, 3
                            do h = 1, 3
                              do i = 1, 3
                                if (sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - pos1_up(j)%center_x_y_z(:) + &
                                             (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2)) < tmp_dist) then
                                  tmp_dist = sqrt(sum((pos1_up(a)%point_x_y_z(c,bond_count_up(a)-e,:)-pos1_up(j)%center_x_y_z(:)+&
                                             (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2))
                                end if
                              end do
                            end do
                          end do
                        end if
                        ave_dist2 = ave_dist2 + 1.0D0/tmp_dist
                      end do
                    end if
                  end do
                  ! END NEW
                end do
                ! 
                ! Keep configuration?
                !
                if (ave_dist2 < ave_dist1) then
                  do f = 1, con_mat(a,b)
                    pos2_up(a)%point_x_y_z(c,bond_count_up(a)-f,:) = pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)
                  end do
                  ave_dist1 = ave_dist2
                else
                  do f = 1, con_mat(a,b)
                    pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:) = pos2_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)
                  end do
                end if
              end do
              !
              ! analyze whether atoms is in an linear environemnt.
              ! If so -> minimize 1/r to neighbouring FODs via rotations
              !
              if (is_planar_linear(a) == 2) then
                do d = 1, cycles
                  ave_dist2 = 0.0
                  !
                  ! generate rotation matrix for bond rotation. Rotate positions
                  !
                  call create_rotMat_bond(full_rot,bond_vector,step_size)
                  do e = 1, con_mat(a,b)
                    call rotate_pos(full_rot, pos1_up(a)%point_x_y_z(c,bond_count_up(a)-e,:), &
                                  & bond_center, pos2_up(a)%point_x_y_z(c,bond_count_up(a)-e,:))
                    !
                    ! evaluate 1/r
                    !
                    do f = 1, size(pos1_up)
                      if (con_mat(a,f) /= 0) then
                        i = pos1_up(f)%n_shells
                        do g = 1, pos1_up(f)%n_points(i)
                          ave_dist2 = ave_dist2 + 1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(c,bond_count_up(a) - e,:) - &
                                                       & pos2_up(f)%point_x_y_z(i,g,:))**2))
                        end do
                      end if
                    end do
                  end do
                  !
                  ! Keep configuration?
                  !
                  if (ave_dist2 < ave_dist1) then
                    do f = 1, con_mat(a,b)
                      pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:) = pos2_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)
                    end do
                    ave_dist1 = ave_dist2
                  else
                    do f = 1, con_mat(a,b)
                      pos2_up(a)%point_x_y_z(c,bond_count_up(a)-f,:) = pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)
                    end do
                  end if
                end do
              end if
              !
              ! If same number of UP and DN FODs in the bond -> just pair them right here (regarding the respective other atoms!)
              !
              if (con_mat(a,b) == con_mat(b,a)) then 
                do f = 1, con_mat(a,b)
                  g = pos1_dn(b)%n_shells
                  pos1_dn(b)%point_x_y_z(g,bond_count_dn(b),:) = pos1_up(a)%point_x_y_z(c,bond_count_up(a)-f,:)
                  bond_count_dn(b) = bond_count_dn(b) + 1
                end do
             end if

            !!!!!!!!!!!!!!
            ! DN CHANNEL !
            !!!!!!!!!!!!!!
            else
              if (con_mat(a,b) == con_mat(b,a)) then                            ! same number of UP and DN FODs in the bond -> nothing to do (already paired)
                continue
              else
                do f = 1, con_mat(a,b)
                  c = pos1_dn(a)%n_shells                                       ! valence shell DN
                  pos1_dn(a)%point_x_y_z(c,bond_count_dn(a),:) = &
                  & bond_center(:) + perpendicular_vec(:)*pos1_dn(a)%r_theta_phi(c,bond_count_dn(a),1)*scale_r         ! initial pos of UP channel. To be rotated around the bond axis
                  bond_count_dn(a) = bond_count_dn(a) + 1
                end do
                !
                ! Rotate the point around the bond vector (bond axis). 
                !
                ave_dist1 = 100000.0
                do d = 1, cycles
                  do f = 1, con_mat(a,b)
                    pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:) = pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:)    ! store original pos
                    call rotate_around_axis(bond_vector,bond_center, & 
                                          & pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:),step_size)
                  end do
                  !
                  ! 1/r calculation
                  !
                  ave_dist2 = 0.0
                  do e = 1, con_mat(a,b)-1
                    do f = e+1, con_mat(a,b)
                      ave_dist2 = ave_dist2 + 1.0D0/sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - &
                                              & pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - f,:))**2))
                    end do
                  end do
                  !
                  ! In addition, evaluate the 1/r to all atoms the corresponding atom is bonded to (and to the atom itself) -> better symmetry
                  !
                  do g = 1,size(pos1_dn)
                    if ((con_mat(a,g)) /= 0 .or. (con_mat(g,a) /= 0) .or. (g == a)) then
                      do e = 1, con_mat(a,b)
                        tmp_dist = sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - &
                                                              & pos1_dn(g)%center_x_y_z(:))**2))
                        if (periodic) then                                                     ! In a peridoic system, take cell vectors into account
                          do f = 1, 3
                            do h = 1, 3
                              do i = 1, 3
                                if (sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - pos1_dn(g)%center_x_y_z(:) + &
                                             (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2)) < tmp_dist) then
                                  tmp_dist = sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:)-pos1_dn(g)%center_x_y_z(:)+&
                                             (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2))
                                end if
                              end do
                            end do
                          end do
                        end if
                        ave_dist2 = ave_dist2 + 1.0D0/tmp_dist  
                      end do
                    end if
                    !
                    ! In addition, use the atoms which the bonded atom is bonded to as well -> more robust. I.e. include all atoms in con_mat(g,:)
                    !
                    do j = 1,size(pos1_dn)
                      if ((con_mat(g,j)) /= 0 .or. (con_mat(j,g) /= 0)) then

                        do e = 1, con_mat(a,b)
                          tmp_dist = sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - &
                                                                & pos1_dn(j)%center_x_y_z(:))**2))
                          if (periodic) then                                                     ! In a peridoic system, take cell vectors into account
                            do f = 1, 3
                              do h = 1, 3
                                do i = 1, 3
                                  if (sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - pos1_dn(j)%center_x_y_z(:) + &
                                               (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2)) < tmp_dist) then
                                    tmp_dist = sqrt(sum((pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-e,:)-pos1_dn(j)%center_x_y_z(:)+&
                                               (f-2)*cell_a(:) + (h-2)*cell_b(:) + (i-2)*cell_c(:))**2))
                                  end if
                                end do
                              end do
                            end do
                          end if
                          ave_dist2 = ave_dist2 + 1.0D0/tmp_dist
                        end do
                      end if
                    end do
                    ! END NEW
                  end do
                  !
                  ! Keep configuration? 
                  !
                  if (ave_dist2 < ave_dist1) then
                    do f = 1, con_mat(a,b)
                      pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:) = pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:)
                    end do
                    ave_dist1 = ave_dist2
                  else
                    do f = 1, con_mat(a,b)
                      pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:) = pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:)
                    end do
                  end if
                end do
                !
                ! analyze whether atoms is in an linear environemnt.
                ! If so -> minimize 1/r to neighbouring FODs via rotations
                !
                if (is_planar_linear(a) == 2) then
                  do d = 1, cycles
                    ave_dist2 = 0.0
                    !
                    ! generate rotation matrix for bond rotation. Rotate positions
                    !
                    call create_rotMat_bond(full_rot,bond_vector,step_size)
                    do e = 1, con_mat(a,b)
                      call rotate_pos(full_rot, pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-e,:), &
                                    & bond_center, pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-e,:))
                      !
                      ! evaluate 1/r
                      !
                      do f = 1, size(pos1_dn)
                        if (con_mat(a,f) /= 0) then
                          i = pos1_dn(f)%n_shells
                          do g = 1, pos1_dn(f)%n_points(i)
                            ave_dist2 = ave_dist2 + 1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(c,bond_count_dn(a) - e,:) - &
                                                         & pos2_dn(f)%point_x_y_z(i,g,:))**2))
                          end do
                        end if
                      end do
                    end do
                    !
                    ! Keep configuration?
                    !
                    if (ave_dist2 < ave_dist1) then
                      do f = 1, con_mat(a,b)
                        pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:) = pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:)
                      end do
                      ave_dist1 = ave_dist2
                    else
                      do f = 1, con_mat(a,b)
                        pos2_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:) = pos1_dn(a)%point_x_y_z(c,bond_count_dn(a)-f,:)
                      end do
                    end if
                  end do
                end if
              end if
            end if
          end if
        end if
      end do
    end if
  end do


!
! OVERWRITE ANY EXISTING pos2
!
  do a = 1, size(pos1_up)                         ! for each atom!
    do b = 1, pos1_up(a)%n_shells                 ! UP CHANNEL
      do c = 1, pos1_up(a)%n_points(b)
        pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
      end do
    end do
  end do
  do a = 1, size(pos1_dn)                         ! for each atom!
    do b = 1, pos1_dn(a)%n_shells                 ! DN CHANNEL
      do c = 1, pos1_dn(a)%n_points(b)
        pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
      end do
    end do
  end do





















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! 2. get lone FODs !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Idea: get all bond vectors to a given atom. Sum them up, normalize the resulting vector => In the opposite direction of this vector is were the lone FODs will be (at least in the way how the bond vectors are defined)
! Optimize the lone FODs in the same way as the bond FODs
! If an atom has no bonds -> optimize the number of lone FODs individually
  do a = 1, size(pos1_up)                                                      ! for each atom
    if (sum(lone_fods(a,:)) == 0) then                                         ! if there are no lone FODs -> skip everything else
      continue
    else
      !
      ! If there are no bonds to this atoms -> optimize valence according to that atom
      !
      if (sum(con_mat(a,:)) == 0) then
        !!!!!!!!!!!!!!
        ! UP CHANNEL !
        !!!!!!!!!!!!!!
        !
        ! Rotate the point around the atom. Optimize with MC. 
        !
        ave_dist1_up = 100000.0
        do j = 1, cycles
          b = pos1_up(a)%n_shells
          do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)                                       ! only lone FODs
            call mc_step(pos1_up(a)%r_theta_phi(b,c,1:3), &
                         pos1_up(a)%center_x_y_z(1:3),    &                                                            ! single MC step to change the pos
                         pos2_up(a)%r_theta_phi(b,c,1:3), &
                         pos2_up(a)%point_x_y_z(b,c,1:3), step_size)
          end do
          !
          ! EVAL UP
          ! 1/r calculation. UP channel. Lone FODs against each other only! And 1/r to all atoms to get symmetry right
          ! valence contributions only.
          ave_dist2_up = 0.0
          b = pos2_up(a)%n_shells
          do c = pos2_up(a)%n_points(b)-lone_fods(a,1)+1, pos2_up(a)%n_points(b)
            do f = pos2_up(a)%n_points(b)-lone_fods(a,1)+1, pos2_up(a)%n_points(b)
              if (c == f) then                                                                                         ! exclude evaluation of a point with itself
                continue
              else
                ave_dist2_up = ave_dist2_up + &
                   1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(b,c,:) - pos2_up(a)%point_x_y_z(b,f,:))**2))
              end if
            end do

            do d = 1,size(pos1_up)                                                                                     ! Take the 1/r to all atoms -> better symmetry
              ave_dist2_up = ave_dist2_up + 1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(b,c,:) - &                            ! change '+' to '-' to place the lone FODs towards the other atoms
                                                        & pos1_up(d)%center_x_y_z(:))**2))
            end do
          end do
          !
          ! Keep configuration?
          !
          if (ave_dist2_up < ave_dist1_up) then
            do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
              pos1_up(a)%point_x_y_z(b,c,:) = pos2_up(a)%point_x_y_z(b,c,:)
              pos1_up(a)%r_theta_phi(b,c,:) = pos2_up(a)%r_theta_phi(b,c,:)
            end do
            ave_dist1_up = ave_dist2_up
          else
            do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
              pos2_up(a)%point_x_y_z(b,c,:) = pos1_up(a)%point_x_y_z(b,c,:)                                            ! reset position2
              pos2_up(a)%r_theta_phi(b,c,:) = pos1_up(a)%r_theta_phi(b,c,:)
            end do
          end if
        end do
        ! 
        ! If same number of UP and DN FODs in the lone pair -> just pair them right here
        ! it is ugly right now, but it works
        !
        if (lone_fods(a,1) == lone_fods(a,2)) then
          b = pos1_up(a)%n_shells
          g = pos1_dn(a)%n_shells
          do h = 1, lone_fods(a,1)
            c = pos1_up(a)%n_points(b) - (h-1)
            i = pos1_dn(a)%n_points(b) - (h-1)
            pos1_dn(a)%point_x_y_z(g,i,:) = pos1_up(a)%point_x_y_z(b,c,:)
            pos1_dn(a)%r_theta_phi(g,i,:) = pos1_up(a)%r_theta_phi(b,c,:)
          end do

          !!!!!!!!!!!!!!
          ! DN CHANNEL !
          !!!!!!!!!!!!!!
        else
          !
          ! Rotate the point around the atom. Optimize with MC. 
          !
          ave_dist1_dn = 100000.0
          do j = 1, cycles
            b = pos1_dn(a)%n_shells
            do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)                                     ! only lone FODs
              call mc_step(pos1_dn(a)%r_theta_phi(b,c,1:3), &
                           pos1_dn(a)%center_x_y_z(1:3),    &                                                          ! single MC step to change the pos
                           pos2_dn(a)%r_theta_phi(b,c,1:3), &
                           pos2_dn(a)%point_x_y_z(b,c,1:3), step_size)
            end do
            !
            ! EVAL DN
            ! 1/r calculation. DN channel. Lone FODs against each other only! And 1/r to all atoms to get symmetry right
            ! valence contributions only.
            ave_dist2_dn = 0.0
            b = pos2_dn(a)%n_shells
            do c = pos2_dn(a)%n_points(b)-lone_fods(a,2)+1, pos2_dn(a)%n_points(b)
              do f = pos2_dn(a)%n_points(b)-lone_fods(a,2)+1, pos2_dn(a)%n_points(b)
                if (c == f) then                                                                                       ! exclude evaluation of a point with itself
                  continue
                else
                  ave_dist2_dn = ave_dist2_dn + &
                     1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - pos2_dn(a)%point_x_y_z(b,f,:))**2))
                end if
              end do

              do d = 1,size(pos1_dn)                                                                                    ! Take the 1/r to all atoms -> better symmetry
                ave_dist2_dn = ave_dist2_dn + 1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - &                           ! change '+' to '-' to place the lone FODs towards the other atoms
                                                          & pos1_dn(d)%center_x_y_z(:))**2))
              end do
            end do
            !
            ! Keep configuration?
            !
            if (ave_dist2_dn < ave_dist1_dn) then
              do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                pos1_dn(a)%point_x_y_z(b,c,:) = pos2_dn(a)%point_x_y_z(b,c,:)
                pos1_dn(a)%r_theta_phi(b,c,:) = pos2_dn(a)%r_theta_phi(b,c,:)
              end do
              ave_dist1_dn = ave_dist2_dn
            else
              do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                pos2_dn(a)%point_x_y_z(b,c,:) = pos1_dn(a)%point_x_y_z(b,c,:)                                          ! reset position2
                pos2_dn(a)%r_theta_phi(b,c,:) = pos1_dn(a)%r_theta_phi(b,c,:)
              end do
            end if
          end do
        end if


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! If there are bonds to the specific atom ! 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else
        counter = 0             ! count how many bonds there are
        bond_vector(:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)                                                   ! initialize bond vector
        do b = 1, size(pos1_up)                                                                                        ! for each other atom
          !
          ! get vector describing the bond. For each bond to an atom, as specified in con_mat
          !
          if (con_mat(a,b) /= 0) then                                                                                  ! for different atoms which are bonded
            counter = counter + 1
            !
            ! get vector describing all bonds
            !
            if (periodic) then                                                                                         ! In a periodic system, take cell vectors into account
              tmp_vector(:) = (/ real(100.0,8),real(100.0,8),real(100.0,8) /)                                          ! temporary vector. Make ridiculously big -> minimizue afterwards
              do c = 1, 3
                do d = 1, 3
                  do e = 1, 3
                    if (sqrt(sum((pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                 (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:))**2)) < &
                                 (sqrt(sum(tmp_vector(:)**2)))) then
                      tmp_vector(:) = pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:) + &
                                       (c-2)*cell_a(:) + (d-2)*cell_b(:) + (e-2)*cell_c(:)
                    end if
                  end do
                end do
              end do
              bond_vector(:) = bond_vector(:) + tmp_vector(:) 
            else                                                                                                       ! in case of no periodicity:
              bond_vector(:) = bond_vector(:) + (pos1_up(a)%center_x_y_z(:) - pos1_up(b)%center_x_y_z(:))              ! vector from b to a
            end if
          end if
        end do

        !
        ! Check if the atom is in a planar or linear environment. If so -> lone FODs need to be out-of-plane
        !

        !
        ! Linear
        !
        if (is_planar_linear(a) == 2) then
          !
          ! For any other atoms: Get a vector describing a bond to this atom
          !
          allocate(tmp_vectors(2,size(pos1_up)))
          tmp_vectors(1,:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
          tmp_vectors(2,:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
          do d = 1, size(pos1_up)
            if (con_mat(a,d) /= 0) then
              tmp_vector(:) = pos1_up(a)%center_x_y_z(:) - pos1_up(d)%center_x_y_z(:)
            end if
          end do
          !
          ! Get a vector that is perpendicular to a bond vector (first call to per_vector) AND one that is perpendicular to this vector (second call should do that)
          !
          call per_vector(tmp_vector(:),pos1_up(a)%center_x_y_z(:),tmp_vectors(1,:))
          call per_vector(tmp_vectors(1,:),pos1_up(a)%center_x_y_z(:),tmp_vectors(2,:))
          ! 
          ! assign center for FOD generation
          !
          if (lone_fods(a,1) > 0) then
            f = pos1_up(a)%n_shells
            g = pos1_up(a)%n_points(f)-lone_fods(a,1)+1                                                      ! outermost shell -> get radius of that shell
            bond_center(:) = pos1_up(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*tmp_vectors(2,:)*pos1_up(a)%r_theta_phi(f,g,1)       ! UP channel, scale to number of lone FODs
          else if (lone_fods(a,2) > 0) then
            f = pos1_dn(a)%n_shells
            g = pos1_dn(a)%n_points(f)-lone_fods(a,2)+1
            bond_center(:) = pos1_dn(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*tmp_vectors(2,:)*pos1_dn(a)%r_theta_phi(f,g,1)       ! DN channel
          end if
          deallocate(tmp_vectors)
        ! 
        ! Planar
        !
        else if (is_planar_linear(a) == 1.and.(counter > 2)) then ! more than two bonds!
          !
          ! get cross product of two vectors to adjacent atoms -> cross product will give a perpendicular vector
          ! 
          allocate(tmp_vectors(2,size(pos1_up)))
          tmp_vectors(1,:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
          tmp_vectors(2,:) = (/ real(0.0,8), real(0.0,8), real(0.0,8) /)
          counter = 0
          loop2: do c = 1, size(pos1_up)                                                                     ! for each other atom!
            if (con_mat(a,c) /= 0) then                                                                      ! if there is a bond between atom a and atom b -> get bond vector
              counter = counter + 1
              tmp_vectors(counter,:) = (pos1_up(a)%center_x_y_z(:) - pos1_up(c)%center_x_y_z(:))
            end if
            if (counter == 2) then                                                                           ! two vectors for cross product
              tmp_vector(:) = tmp_vectors(1,:)                                                               ! define bond vector (around which multiple FODs will be rotated)
              !
              ! create normal vector. Cross product -> Perpendicular vector
              !
              perpendicular_vec(1) = tmp_vectors(1,2)*tmp_vectors(2,3) - tmp_vectors(1,3)*tmp_vectors(2,2)
              perpendicular_vec(2) = tmp_vectors(1,3)*tmp_vectors(2,1) - tmp_vectors(1,1)*tmp_vectors(2,3)
              perpendicular_vec(3) = tmp_vectors(1,1)*tmp_vectors(2,2) - tmp_vectors(1,2)*tmp_vectors(2,1)
              exit loop2
            end if
          end do loop2
          !
          ! normalize vector
          !
          perpendicular_vec(:) = perpendicular_vec(:)/sqrt(sum(perpendicular_vec(:)**2)) 
          ! 
          ! assign center for FOD generation
          !
          if (lone_fods(a,1) > 0) then
            f = pos1_up(a)%n_shells
            g = pos1_up(a)%n_points(f)-lone_fods(a,1)+1                                                      ! outermost shell -> get radius of that shell
            bond_center(:) = pos1_up(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*perpendicular_vec(:)*pos1_up(a)%r_theta_phi(f,g,1)   ! UP channel
            if ((lone_fods(a,1) == 1) .and. (lone_fods(a,2) == 0)) then                                      ! in case there is only exactly one lone FODs (none in the other spin channel)
              bond_center(:) = pos1_up(a)%center_x_y_z(:) + perpendicular_vec(:)*pos1_up(a)%r_theta_phi(f,g,1)   
            end if              
          else if (lone_fods(a,2) > 0) then
            f = pos1_dn(a)%n_shells
            g = pos1_dn(a)%n_points(f)-lone_fods(a,2)+1
            bond_center(:) = pos1_dn(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*perpendicular_vec(:)*pos1_dn(a)%r_theta_phi(f,g,1)   ! DN
            if ((lone_fods(a,1) == 0) .and. (lone_fods(a,2) == 1)) then                                      ! in case there is only exactly one lone FODs (none in the other spin channel)
              bond_center(:) = pos1_up(a)%center_x_y_z(:) + perpendicular_vec(:)*pos1_dn(a)%r_theta_phi(f,g,1)   
            end if              
          end if
          deallocate(tmp_vectors)
        !
        ! Neither linear nor planar
        !
        else                                                                                                 ! if less than one bonding partner -> just initialize 'normally'
          ! 
          ! assign center for FOD generation
          !
          bond_vector(:) = bond_vector(:)/sqrt(sum(bond_vector(:)**2))                                       ! normalize bond vector
          if (lone_fods(a,1) > 0) then
            f = pos1_up(a)%n_shells
            g = pos1_up(a)%n_points(f)-lone_fods(a,1)+1                                                      ! outermost shell -> get radius of that shell
            bond_center(:) = pos1_up(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*bond_vector(:)*pos1_up(a)%r_theta_phi(f,g,1)         ! UP Scale to number of lone FODs -> bring multiple lone_FODs closer to the atoms
          else if (lone_fods(a,2) > 0) then
            f = pos1_dn(a)%n_shells
            g = pos1_dn(a)%n_points(f)-lone_fods(a,2)+1
            bond_center(:) = pos1_up(a)%center_x_y_z(:) + &
            & 1.0D0/((lone_fods(a,1)+lone_fods(a,2))/2.0D0)*bond_vector(:)*pos1_dn(a)%r_theta_phi(f,g,1)         ! DN
          end if
        end if

        !!!!!!!!!!!!!!!!!!!!
        ! SINGLE LONE FODs !
        !!!!!!!!!!!!!!!!!!!!
        if (lone_fods(a,1) == 1) then                                                                        ! UP lone FOD
          b = pos1_up(a)%n_shells                                                                            ! valence shell
          c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1                                                        ! Only one point to place
          pos1_up(a)%point_x_y_z(b,c,:) = bond_center(:)
        end if
        if (lone_fods(a,2) == 1) then                                                                        ! DN lone FOD
          b = pos1_dn(a)%n_shells                                                                            ! valence shell
          c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1                                                        ! Only one point to place
          pos1_dn(a)%point_x_y_z(b,c,:) = bond_center(:)
        end if
        !
        ! OVERWRITE ANY EXISTING pos2 -> currently necessary for next optimization
        !
        do b = 1, pos1_up(a)%n_shells                                                                        ! UP CHANNEL
          do c = 1, pos1_up(a)%n_points(b)
            pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
          end do
        end do
        do b = 1, pos1_dn(a)%n_shells                                                                        ! DN CHANNEL
          do c = 1, pos1_dn(a)%n_points(b)
            pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
          end do
        end do


        !!!!!!!!!!!!!!!!!!!!!!
        ! MULTIPLE LONE FODs !
        !!!!!!!!!!!!!!!!!!!!!!
        ! Take perpendicular vector to the bond vector -> is perpendicular to the bond axis
        ! Use this vector (and a given radius) to place FODs perpendicular to the bond (vector applied to the lone-FOD center!)
        ! Rotate these points around the bond vector via MC to minimize their 1/r contribution
        !
        ! find a vector which is perpendicular to the bond vector
        !
        call per_vector(bond_vector(:),pos1_up(a)%center_x_y_z(:),perpendicular_vec(:))

!!!!!!!!!!!!!!
! UP CHANNEL !
!!!!!!!!!!!!!!
        if (lone_fods(a,1) > 1) then                                                                                   ! UP CHANNEL
          b = pos1_up(a)%n_shells                                                                                      ! valence shell UP
          do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)                                       ! UP FODs which are not in a bond, but are LONE -> only consider these FODs
            pos1_up(a)%point_x_y_z(b,c,:) = &
            & bond_center(:) + perpendicular_vec(:)*pos1_up(a)%r_theta_phi(b,c,1)*scale_r                              ! initial pos of UP channel. To be rotated around the bond axis. scale_r is a factor
          end do
          !
          ! Rotate the point around the bond vector (bond axis).
          !
          ave_dist1_up = 100000.0
          do j = 1, cycles
            do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
              pos2_up(a)%point_x_y_z(b,c,:) = pos1_up(a)%point_x_y_z(b,c,:)                                            ! store original pos
              if (is_planar_linear(a) == 0) then                                                                       ! Not planar/linear. Rotate around bond_center
                call rotate_around_axis(bond_vector,bond_center,pos1_up(a)%point_x_y_z(b,c,:),step_size)
              else                                                                                                     ! planar/linear -> rotate around atom
                call rotate_around_axis(bond_vector,pos1_up(a)%center_x_y_z,pos1_up(a)%point_x_y_z(b,c,:),step_size)
              end if
            end do
            !
            ! EVAL UP
            ! 1/r calculation. UP. Lone FODs against the rest
            ! valence contributions only.
            ave_dist2_up = 0.0
            b = pos1_up(a)%n_shells
            do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
              do d = 1, size(pos1_up)                                                                                  ! for all atoms
                if ((con_mat(a,d) /= 0) .or. (con_mat(d,a) /= 0) .or. (a == d)) then                                   ! only evaluate 1/r for neighbouring atoms -> all which have a bond to the atom under observation. 
                                                                                                                       ! and for the same atom too
                  e = pos1_up(d)%n_shells
                  do f = 1, pos1_up(d)%n_points(e)
                    if ((a == d) .and. (b == e) .and. (c == f)) then                                                   ! exclude evaluation of a point with itself
                      continue
                    else
! HERE: periodicity?
                      ave_dist2_up = ave_dist2_up + &
                         1.0D0/sqrt(sum((pos1_up(a)%point_x_y_z(b,c,:) - pos1_up(d)%point_x_y_z(e,f,:))**2))
                    end if
                  end do
                end if
              end do
            end do
            !
            ! Keep configuration?
            !
            if (ave_dist2_up < ave_dist1_up) then
              do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
                pos2_up(a)%point_x_y_z(b,c,:) = pos1_up(a)%point_x_y_z(b,c,:)
              end do
              ave_dist1_up = ave_dist2_up
            else
              do c = pos1_up(a)%n_points(b)-lone_fods(a,1)+1, pos1_up(a)%n_points(b)
                pos1_up(a)%point_x_y_z(b,c,:) = pos2_up(a)%point_x_y_z(b,c,:)
              end do
            end if
          end do
          !
          ! If same number of UP and DN FODs in the lone pair -> just pair them right here (regarding the respective other atoms!)
          ! it is ugly right now, but it works
          !
          if (lone_fods(a,1) == lone_fods(a,2)) then
            b = pos1_up(a)%n_shells
            g = pos1_dn(a)%n_shells
            do h = 1, lone_fods(a,1)
              c = pos1_up(a)%n_points(b) - (h-1)
              i = pos1_dn(a)%n_points(b) - (h-1)
              pos1_dn(a)%point_x_y_z(g,i,:) = pos1_up(a)%point_x_y_z(b,c,:)
            end do
          end if
        end if

!!!!!!!!!!!!!!
! DN CHANNEL !
!!!!!!!!!!!!!!
        if (lone_fods(a,2) > 1) then                                                                                   ! DN CHANNEL
          if (lone_fods(a,1) == lone_fods(a,2)) then                                                                   ! same number of UP and DN lone FODs -> nothing to do (already paired)
            continue
          else
            b = pos1_dn(a)%n_shells                                                                                    ! valence shell DN
            do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)                                     ! DN FODs which are not in a bond, but are LONE -> only consider these FODs
              pos1_dn(a)%point_x_y_z(b,c,:) = &
              & bond_center(:) + perpendicular_vec(:)*pos1_dn(a)%r_theta_phi(b,c,1)*scale_r                            ! initial pos of DN channel. To be rotated around the bond axis. 
            end do
            !
            ! Rotate the point around the bond vector (bond axis). 
            !
            ave_dist1_dn = 100000.0
            do j = 1, cycles
              do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                pos2_dn(a)%point_x_y_z(b,c,:) = pos1_dn(a)%point_x_y_z(b,c,:)                                          ! store original pos
                if (is_planar_linear(a) == 0) then                                                                     ! Not planar/linear. Rotate around bond_center
                  call rotate_around_axis(bond_vector,bond_center,pos1_dn(a)%point_x_y_z(b,c,:),step_size)
                else                                                                                                   ! planar/linear -> rotate around atom
                  call rotate_around_axis(bond_vector,pos1_dn(a)%center_x_y_z,pos1_dn(a)%point_x_y_z(b,c,:),step_size) 
                end if
              end do
              !
              ! EVAL DN
              ! 1/r calculation. DN. Lone FODs against the rest
              ! valence contributions only.
              ave_dist2_dn = 0.0
              b = pos1_dn(a)%n_shells
              do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                do d = 1, size(pos1_dn)                                                                                ! for all atoms
                  if ((con_mat(a,d) /= 0) .or. (con_mat(d,a) /= 0) .or. (a == d)) then                                 ! only evaluate 1/r for meighbouring atoms -> all which have a bond to the atom under observation. 
                                                                                                                       ! and for the same atom too
                    e = pos1_dn(d)%n_shells
                    do f = 1, pos1_dn(d)%n_points(e)
                      if ((a == d) .and. (b == e) .and. (c == f)) then                                                 ! exclude evaluation of a point with itself
                        continue
                      else
! HERE: Periodicity?
                        ave_dist2_dn = ave_dist2_dn + &
                           1.0D0/sqrt(sum((pos1_dn(a)%point_x_y_z(b,c,:) - pos1_dn(d)%point_x_y_z(e,f,:))**2))
                      end if
                    end do
                  end if
                end do
              end do
              !
              ! Keep configuration?
              !
              if (ave_dist2_dn < ave_dist1_dn) then
                do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                  pos2_dn(a)%point_x_y_z(b,c,:) = pos1_dn(a)%point_x_y_z(b,c,:)
                end do
                ave_dist1_dn = ave_dist2_dn
              else
                do c = pos1_dn(a)%n_points(b)-lone_fods(a,2)+1, pos1_dn(a)%n_points(b)
                  pos1_dn(a)%point_x_y_z(b,c,:) = pos2_dn(a)%point_x_y_z(b,c,:)
                end do
              end if
            end do
          end if
        end if 
      end if 
    end if
  end do


!
! OVERWRITE ANY EXISTING pos2 -> currently necessary for next optimization
!
  do a = 1, size(pos1_up)                                                      ! for each atom!
    do b = 1, pos1_up(a)%n_shells                                              ! UP CHANNEL
      do c = 1, pos1_up(a)%n_points(b)
        pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
      end do
    end do
  end do
  do a = 1, size(pos1_dn)                                                      ! for each atom!
    do b = 1, pos1_dn(a)%n_shells                                              ! DN CHANNEL
      do c = 1, pos1_dn(a)%n_points(b)
        pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
      end do
    end do
  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CREATE CORE FODs. Minimize 1/r for core  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do a = 1, size(pos1_up)                                                      ! for each atom!

! NEW Use structural motifs
    !
    ! Get UP CHANNEL  
    !
    do b = 1, pos1_up(a)%n_shells-1                                            ! Only core
      do c = 1, pos1_up(a)%n_points(b)
        call struct_motif(pos1_up(a)%n_points(b),c,pos1_up(a)%r_theta_phi(b,c,1),&
                          pos1_up(a)%center_x_y_z(1:3),pos1_up(a)%point_x_y_z(b,c,1:3))
        pos2_up(a)%point_x_y_z(b,c,1:3) = pos1_up(a)%point_x_y_z(b,c,1:3)
      end do
      !
      ! Optimize core structure globally (via rotations) with respect to any smaller shell
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ROTATIONS                                              !
      ! Rotate current core against any lower core. UP Channel !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((b > 1) .and. (b <= pos1_up(a)%n_shells-1)) then                           ! any shell except the 1s. 
        ave_dist1_up = 100000000.0
        do t = 1, cycles
          call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
          do c = 1, pos1_up(a)%n_points(b)
            call rotate_pos(full_rot, pos1_up(a)%point_x_y_z(b,c,1:3), &
                            pos1_up(a)%center_x_y_z(1:3), pos2_up(a)%point_x_y_z(b,c,1:3))
          end do
          !
          ! 1/r calculation
          !
          ave_dist2_up = 0.0                                                         ! 1/r UP
          do c = 1, pos2_up(a)%n_points(b)
            do d = 1, b-1                                                            ! Use any smaller shell structure and optimize the higher shell according to the lower core structure
              do f = 1, pos1_up(a)%n_points(d)                                       ! Unchanged index, thus index 1
                ave_dist2_up  = ave_dist2_up + &
                       1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(b,c,:) - pos1_up(a)%point_x_y_z(d,f,:))**2))
              end do
            end do
          end do
          !
          ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
          ! How: ... Introduce random number rand_metro [0,1]
          !          if (ave_dist2_up < ave_dist1_up) -> take new config
          !          if (ave_dist2_up > ave_dist1_up) -> if (ave_dist2_up-ave_dist1_up) < rand_metro*step_size => still take new config. Else: Take old config.
          !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
          if (ave_dist2_up < ave_dist1_up) then
            pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
            pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
            ave_dist1_up = ave_dist2_up
          else
            call random_number(rand_metro)
            if ((ave_dist2_up-ave_dist1_up) < rand_metro*step_size*0.01) then
              pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
              pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
              ave_dist1_up = ave_dist2_up
            end if
          end if
        end do  ! end MC cycles
      end if    ! end if - shells
    end do      ! end UP channel
    ! 
    ! Further, rotate all cores with respect to the valence pattern!
    !
    ave_dist1_up = 100000000.0
    do t = 1, cycles
      call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
      do b = 1, pos1_up(a)%n_shells-1 
        do c = 1, pos1_up(a)%n_points(b)
          call rotate_pos(full_rot, pos1_up(a)%point_x_y_z(b,c,1:3), &
                          pos1_up(a)%center_x_y_z(1:3), pos2_up(a)%point_x_y_z(b,c,1:3))
        end do
      end do
      !
      ! 1/r calculation
      !
      ave_dist2_up = 0.0                                                         ! 1/r UP
      do b = 1, pos2_up(a)%n_shells-1
        do c = 1, pos2_up(a)%n_points(b)
          do g = 1, size(pos1_up)                                                ! for all atoms 
            d = pos1_up(g)%n_shells                                                ! Outermost fods
            do f = 1, pos1_up(g)%n_points(d)                                       ! 
              ave_dist2_up  = ave_dist2_up + &
                     1.0D0/sqrt(sum((pos2_up(a)%point_x_y_z(b,c,:) - pos1_up(g)%point_x_y_z(d,f,:))**2))
            end do
          end do
        end do
      end do
      !
      ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
      ! How: ... Introduce random number rand_metro [0,1]
      !          if (ave_dist2_up < ave_dist1_up) -> take new config
      !          if (ave_dist2_up > ave_dist1_up) -> if (ave_dist2_up-ave_dist1_up) < rand_metro*step_size => still take new config. Else: Take old config.
      !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
      if (ave_dist2_up < ave_dist1_up) then
        pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
        pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
        ave_dist1_up = ave_dist2_up
      else
        call random_number(rand_metro)
        if ((ave_dist2_up-ave_dist1_up) < rand_metro*step_size*0.01) then
          pos1_up(a)%point_x_y_z(:,:,:) = pos2_up(a)%point_x_y_z(:,:,:)
          pos1_up(a)%r_theta_phi(:,:,:) = pos2_up(a)%r_theta_phi(:,:,:)
          ave_dist1_up = ave_dist2_up
        end if
      end if
    end do    ! end MC cycle. final end UP channels

    !
    ! Get DN CHANNEL  
    !
    do b = 1, pos1_dn(a)%n_shells-1                                            ! Only core
      do c = 1, pos1_dn(a)%n_points(b)
        call struct_motif(pos1_dn(a)%n_points(b),c,pos1_dn(a)%r_theta_phi(b,c,1),&
                          pos1_dn(a)%center_x_y_z(1:3),pos1_dn(a)%point_x_y_z(b,c,1:3))
        pos2_dn(a)%point_x_y_z(b,c,1:3) = pos1_dn(a)%point_x_y_z(b,c,1:3)
      end do
      !
      ! Optimize core structure globally (via rotations) with respect to any smaller shell
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ROTATIONS                                              !
      ! Rotate current core against any lower core. DN Channel !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((b > 1) .and. (b <= pos1_dn(a)%n_shells-1)) then                           ! any shell except the 1s. 
        ave_dist1_dn = 100000000.0
        do t = 1, cycles
          call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
          do c = 1, pos1_dn(a)%n_points(b)
            call rotate_pos(full_rot, pos1_dn(a)%point_x_y_z(b,c,1:3), &
                            pos1_dn(a)%center_x_y_z(1:3), pos2_dn(a)%point_x_y_z(b,c,1:3))
          end do
          !
          ! 1/r calculation
          !
          ave_dist2_dn = 0.0                                                         ! 1/r DN
          do c = 1, pos2_dn(a)%n_points(b)
            do d = 1, b-1                                                            ! Use any smaller shell structure and optimize the higher shell according to the lower core structure
              do f = 1, pos1_dn(a)%n_points(d)                                       ! Unchanged index, thus index 1
                ave_dist2_dn  = ave_dist2_dn + &
                       1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - pos1_dn(a)%point_x_y_z(d,f,:))**2))
              end do
            end do
          end do
          !
          ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
          ! How: ... Introduce random number rand_metro [0,1]
          !          if (ave_dist2_dn < ave_dist1_dn) -> take new config
          !          if (ave_dist2_dn > ave_dist1_dn) -> if (ave_dist2_dn-ave_dist1_dn) < rand_metro*step_size => still take new config. Else: Take old config.
          !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
          if (ave_dist2_dn < ave_dist1_dn) then
            pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
            pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
            ave_dist1_dn = ave_dist2_dn
          else
            call random_number(rand_metro)
            if ((ave_dist2_dn-ave_dist1_dn) < rand_metro*step_size*0.01) then
              pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
              pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
              ave_dist1_dn = ave_dist2_dn
            end if
          end if
        end do  ! end MC cycles
      end if
    end do
    ! 
    ! Further, rotate all cores with respect to the valence pattern!
    !
    ave_dist1_dn = 100000000.0
    do t = 1, cycles
      call create_rotMat(full_rot, step_size)                                    ! generates a random rotation matrix
      do b = 1, pos1_dn(a)%n_shells-1
        do c = 1, pos1_dn(a)%n_points(b)
          call rotate_pos(full_rot, pos1_dn(a)%point_x_y_z(b,c,1:3), &
                          pos1_dn(a)%center_x_y_z(1:3), pos2_dn(a)%point_x_y_z(b,c,1:3))
        end do
      end do
      !
      ! 1/r calculation
      !
      ave_dist2_dn = 0.0                                                         ! 1/r UP
      do b = 1, pos2_dn(a)%n_shells-1
        do c = 1, pos2_dn(a)%n_points(b)
          do g = 1, size(pos1_dn)                                                ! for all atoms 
            d = pos1_dn(g)%n_shells                                                ! Outermost fods
            do f = 1, pos1_dn(g)%n_points(d)                                       ! 
              ave_dist2_dn  = ave_dist2_dn + &
                      1.0D0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - pos1_dn(g)%point_x_y_z(d,f,:))**2))
            end do
          end do
        end do
      end do
      !
      ! Use Metropolis-like algorithm. Allow f_r to go uphill for the rotations!
      ! How: ... Introduce random number rand_metro [0,1]
      !          if (ave_dist2_up < ave_dist1_up) -> take new config
      !          if (ave_dist2_up > ave_dist1_up) -> if (ave_dist2_up-ave_dist1_up) < rand_metro*step_size => still take new config. Else: Take old config.
      !          ==>> If new 1/r is only slighlty larger AND rand_metro*step_size is small enough
      if (ave_dist2_dn < ave_dist1_dn) then
        pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
        pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
        ave_dist1_dn = ave_dist2_dn
      else
        call random_number(rand_metro)
        if ((ave_dist2_dn-ave_dist1_dn) < rand_metro*step_size*0.01) then
          pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
          pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
          ave_dist1_dn = ave_dist2_dn
        end if
      end if
    end do    ! end MC cycle. final end DN channel
  end do
! END NEW

! For molecules -> not rotation of DN vs UP in the end. Shouldn't be necessary




!!!!
!!!!
!!!!
!!!!
!!!!
!!!!
!!!!      !
!!!!      ave_dist2_dn_core = 0.0
!!!!      do b = 1, pos2_dn(a)%n_shells-1
!!!!        do c = 1, pos2_dn(a)%n_points(b)
!!!!          do d = 1, size(pos1_dn)
!!!!            do e = 1, pos2_dn(d)%n_shells
!!!!              do f = 1, pos2_dn(d)%n_points(e)
!!!!                if ((a == d) .and. (b == e) .and. (c == f)) then               ! do not evaluate 1/r for one and the same point
!!!!                else
!!!!                   ave_dist2_dn_core = ave_dist2_dn_core + &
!!!!                           1.0/sqrt(sum((pos2_dn(a)%point_x_y_z(b,c,:) - pos2_dn(d)%point_x_y_z(e,f,:))**2))
!!!!                end if
!!!!              end do
!!!!            end do
!!!!          end do
!!!!        end do
!!!!      end do
!!!!! HERE: NEW use Metropolis algorithm. Allow f_r to go uphill for the rotations!
!!!!! How: ... Introduce random number rand_metro [0,1]
!!!!!          if (ave_dist2_up < ave_dist1_up) -> take new config
!!!!!          if (ave_dist2_up > ave_dist1_up) -> if (ave_dist2_up-ave_dist1_up) < rand_metro*0.01 => still take new config. Else: Take old config.
!!!!!          ==>> If new 1/r is only slighlty larger AND rand_metro*0.01 is small enough
!!!!      if (ave_dist2_dn_core < ave_dist1_dn_core) then
!!!!        pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
!!!!        pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
!!!!        ave_dist1_dn_core = ave_dist2_dn_core
!!!!      else
!!!!        call random_number(rand_metro)
!!!!        if ((ave_dist2_dn_core-ave_dist1_dn_core) < rand_metro*step_size) then
!!!!          pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
!!!!          pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
!!!!          ave_dist1_dn_core = ave_dist2_dn_core
!!!!        end if
!!!!      end if
!!!!!!!      !
!!!!!!!      ! Keep configuration?
!!!!!!!      !
!!!!!!!      if (ave_dist2_dn_core < ave_dist1_dn_core) then                          ! < if minimizing 1/r, > if maximizing 1/r
!!!!!!!        pos1_dn(a)%point_x_y_z(:,:,:) = pos2_dn(a)%point_x_y_z(:,:,:)
!!!!!!!        pos1_dn(a)%r_theta_phi(:,:,:) = pos2_dn(a)%r_theta_phi(:,:,:)
!!!!!!!        ave_dist1_dn_core = ave_dist2_dn_core
!!!!!!!      end if
!!!!    end do                                                                     ! end MC cycles for core optimization
!!!!  end do                                                                       ! end atoms

!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
!! END OF fodMC !!
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!

! Point charge dipole evaluation !
! Sum r_{FODs} - SUM r_{atoms}
! Calculate with respect to the center of all atoms
! 
  cent_x = 0.0D0
  cent_y = 0.0D0
  cent_z = 0.0D0
  do a = 1, size(pos1_up)
    cent_x = cent_x + pos1_up(a)%center_x_y_z(1)
    cent_y = cent_y + pos1_up(a)%center_x_y_z(2)
    cent_z = cent_z + pos1_up(a)%center_x_y_z(3)
  end do
  cent_x = cent_x/size(pos1_up)
  cent_y = cent_y/size(pos1_up)
  cent_z = cent_z/size(pos1_up)

  dip_x = 0.0D0                                                                   ! initialize with 0
  dip_y = 0.0D0
  dip_z = 0.0D0
  do a = 1, size(pos1_up)                                                       ! for all atoms
    do b = 1, pos1_up(a)%n_shells                                              ! and all UP-FODs
      do c = 1, pos1_up(a)%n_points(b)
        dip_x = dip_x + pos1_up(a)%point_x_y_z(b,c,1) - cent_x
        dip_y = dip_y + pos1_up(a)%point_x_y_z(b,c,2) - cent_y
        dip_z = dip_z + pos1_up(a)%point_x_y_z(b,c,3) - cent_z
      end do
    end do
    do b = 1, pos1_dn(a)%n_shells                                               ! all DN-FODs
      do c = 1, pos1_dn(a)%n_points(b)
        dip_x = dip_x + pos1_dn(a)%point_x_y_z(b,c,1) - cent_x
        dip_y = dip_y + pos1_dn(a)%point_x_y_z(b,c,2) - cent_y
        dip_z = dip_z + pos1_dn(a)%point_x_y_z(b,c,3) - cent_z
      end do
    end do
  end do
  !
  ! Charge and spin evaluated in the beginning of the molecular guess creation
  !
  write(6,fmt='(1X,A,F8.3,2X,A,F8.3,2X,A,3(F12.5,2X))') 'charge = ',charge,' spin = ',spin, &
                         & ' Point charge dipole ', dip_x, dip_y, dip_z
  write(6,*) ' '

  !
  ! Total time
  !
  call cpu_time(finish_t)
  write(6,*) 'Total CPU time was ',finish_t - start_t,'s'

  !
  ! If periodic system: Generate .cif file
  !
  if (periodic) then
    len_a = sqrt(cell_a(1)**2+cell_a(2)**2+cell_a(3)**2)
    len_b = sqrt(cell_b(1)**2+cell_b(2)**2+cell_b(3)**2)
    len_c = sqrt(cell_c(1)**2+cell_c(2)**2+cell_c(3)**2)
    angle_bc = ACOS((cell_b(1)*cell_c(1) + cell_b(2)*cell_c(2) + cell_b(3)*cell_c(3))/(len_b*len_c))*180.0D0/pi
    angle_ac = ACOS((cell_a(1)*cell_c(1) + cell_a(2)*cell_c(2) + cell_a(3)*cell_c(3))/(len_a*len_c))*180.0D0/pi
    angle_ab = ACOS((cell_a(1)*cell_b(1) + cell_a(2)*cell_b(2) + cell_a(3)*cell_b(3))/(len_a*len_b))*180.0D0/pi

    len_a = sqrt(cell_a(1)**2+cell_a(2)**2+cell_a(3)**2)*0.529177D0/units_factor ! always in angstrom
    len_b = sqrt(cell_b(1)**2+cell_b(2)**2+cell_b(3)**2)*0.529177D0/units_factor
    len_c = sqrt(cell_c(1)**2+cell_c(2)**2+cell_c(3)**2)*0.529177D0/units_factor

    open(unit=19,file='Nuc_FOD.cif',status='unknown',action='write')
    write(19,*) '# fodMC guess'
    write(19,*) ' '
    write(19,fmt='(A)') "_chemical_name_common                  'fodMC'"
    write(19,fmt='(A,F13.8)') "_cell_length_a                         ",len_a
    write(19,fmt='(A,F13.8)') "_cell_length_b                         ",len_b
    write(19,fmt='(A,F13.8)') "_cell_length_c                         ",len_c
    write(19,fmt='(A,F9.4)') "_cell_angle_alpha                      ",angle_bc
    write(19,fmt='(A,F9.4)') "_cell_angle_beta                       ",angle_ac
    write(19,fmt='(A,F9.4)') "_cell_angle_gamma                      ",angle_ab
    write(19,fmt='(A)') "_space_group_name_H-M_alt              'P 1'"
    write(19,fmt='(A)') "loop_"
    write(19,fmt='(A)') "   _atom_site_type_symbol"
    write(19,fmt='(A)') "   _atom_site_fract_x"
    write(19,fmt='(A)') "   _atom_site_fract_y"
    write(19,fmt='(A)') "   _atom_site_fract_z"
    do a = 1, size(pos1_up)
      call cart_to_frac(cell_a(1:3),cell_b(1:3),cell_c(1:3),pos1_up(a)%center_x_y_z(1:3),pos2_up(a)%center_x_y_z(1:3))
      if (pos1_up(a)%elements(2:2) == '_') then
        write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(a)%elements(1:1),pos2_up(a)%center_x_y_z(1:3)
      else
        write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(a)%elements(1:2),pos2_up(a)%center_x_y_z(1:3)
      end if
    end do
    do a = 1, size(pos1_up)
      do b = 1, pos1_up(a)%n_shells
        do c = 1, pos1_up(a)%n_points(b)
          call cart_to_frac(cell_a(1:3),cell_b(1:3),cell_c(1:3),pos1_up(a)%point_x_y_z(b,c,1:3),pos2_up(a)%point_x_y_z(b,c,1:3))
          write(19,fmt='(A,3X,3(F13.8,2X))') 'X',pos2_up(a)%point_x_y_z(b,c,1:3)
        end do
      end do
    end do
    do a = 1, size(pos1_dn)
      do b = 1, pos1_dn(a)%n_shells
        do c = 1, pos1_dn(a)%n_points(b)
          call cart_to_frac(cell_a(1:3),cell_b(1:3),cell_c(1:3),pos1_dn(a)%point_x_y_z(b,c,1:3),pos2_dn(a)%point_x_y_z(b,c,1:3))
          write(19,fmt='(A,3X,3(F13.8,2X))') 'He',pos2_dn(a)%point_x_y_z(b,c,1:3)
        end do
      end do
    end do
    close(unit=19)
  end if


  !
  ! GENERATE CLUSTER AND FRMORB FILE for NRLMOL
  !
  charge     = 0.0
  spin       = 0.0
  counter_up = 0
  counter_dn = 0
  do a = 1, size(pos1_up)
    charge     = charge + real(element_number(a) - pseudo_charge(a) - & 
               & sum(pos1_up(a)%n_points(:)) - sum(pos1_dn(a)%n_points(:)),8)
    spin       = spin   + real(sum(pos1_up(a)%n_points(:)) - sum(pos1_dn(a)%n_points(:)),8)
    counter_up = counter_up + sum(pos1_up(a)%n_points(:))
    counter_dn = counter_dn + sum(pos1_dn(a)%n_points(:))
  end do

  open(unit=19,file='CLUSTER',status='unknown',action='write')
  write(19,fmt='(A)') 'LDA-PW91*LDA-PW91            (DF TYPE EXCHANGE*CORRELATION)'
  write(19,fmt='(A)') 'NONE                         (TD, OH, IH, X, Y, XY, ... OR GRP)'
  write(19,fmt='(I3,A)') size(pos1_up),'                            (NUMBER OF ATOMS)'
  do a = 1, size(pos1_up)
    if ((pos1_up(a)%elements(3:5) == 'ECP') .or. (pos1_up(a)%elements(4:6) == 'ECP')) then
      write(19,fmt='(3(F13.8,2X),I3,1X,A)') pos1_up(a)%center_x_y_z(1:3)/units_factor, element_number(a), &
                                            & ' ECP (R, Z, Pseudopotential)'
    else
       write(19,fmt='(3(F13.8,2X),I3,1X,A)') pos1_up(a)%center_x_y_z(1:3)/units_factor, element_number(a), &
                                            & ' ALL (R, Z, ALL-ELECTRON)'
    end if
  end do
  write(19,fmt='(2(F6.3,2X),19X,A)') charge, spin, '(NET CHARGE AND NET SPIN)'
  close(unit=19)

  open(unit=19,file='FRMORB',status='unknown',action='write')
  write(19,fmt='(I3,5X,I3)') counter_up, counter_dn
  do a = 1, size(pos1_up)
    do b = 1, pos1_up(a)%n_shells
      do c = 1, pos1_up(a)%n_points(b)
        write(19,fmt='(3(F13.8,2X))') pos1_up(a)%point_x_y_z(b,c,1:3)/units_factor
      end do
    end do
  end do
  do a = 1, size(pos1_up)
    do b = 1, pos1_dn(a)%n_shells
      do c = 1, pos1_dn(a)%n_points(b)
        write(19,fmt='(3(F13.8,2X))') pos1_dn(a)%point_x_y_z(b,c,1:3)/units_factor
      end do
    end do
  end do
  close(unit=19)

  !
  ! GENERATE Nuc_FOD.xyz. For PyFLOSIC
  !
  open(unit=19,file='Nuc_FOD.xyz',status='unknown',action='write')
  write (junk, '(I4)') size(pos1_up)+counter_up+counter_dn                     ! number of entries in the xyz file
  write(19,fmt='(A)') adjustl(junk)
  write(19,fmt='(A)') 'angstrom'
  do a = 1, size(pos1_up)
    if (pos1_up(a)%elements(2:2) == '_') then
      write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(a)%elements(1:1),pos1_up(a)%center_x_y_z(1:3)*0.529177D0/units_factor ! always in angstrom
    else
      write(19,fmt='(A,3X,3(F13.8,2X))') pos1_up(a)%elements(1:2),pos1_up(a)%center_x_y_z(1:3)*0.529177D0/units_factor
    end if
  end do
  do a = 1, size(pos1_up)
    do b = 1, pos1_up(a)%n_shells
      do c = 1, pos1_up(a)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'X',pos1_up(a)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
  end do
  do a = 1, size(pos1_up)
    do b = 1, pos1_dn(a)%n_shells
      do c = 1, pos1_dn(a)%n_points(b)
        write(19,fmt='(A,3X,3(F13.8,2X))') 'He',pos1_dn(a)%point_x_y_z(b,c,1:3)*0.529177D0/units_factor
      end do
    end do
  end do
  close(unit=19)

  deallocate(bond_pattern)
  deallocate(con_mat)
  deallocate(lone_fods)
  deallocate(bond_count_up)
  deallocate(bond_count_dn)
end if


! Deallocate everything
do a = 1, number_of_centers
  deallocate(pos1_up(a)%n_points)
  deallocate(pos1_up(a)%r_theta_phi)
  deallocate(pos1_up(a)%point_x_y_z)
  deallocate(pos2_up(a)%n_points)
  deallocate(pos2_up(a)%r_theta_phi)
  deallocate(pos2_up(a)%point_x_y_z)
  deallocate(pos1_dn(a)%n_points)                                              
  deallocate(pos1_dn(a)%r_theta_phi)
  deallocate(pos1_dn(a)%point_x_y_z)
  deallocate(pos2_dn(a)%n_points)
  deallocate(pos2_dn(a)%r_theta_phi)
  deallocate(pos2_dn(a)%point_x_y_z)
end do

!end program fodMC
end subroutine get_guess


subroutine mc_step(rthetaphi1, center_xyz, rthetaphi2, xyz2, stepsize)
real(8), intent(in)  :: rthetaphi1(3)
real(8), intent(in)  :: center_xyz(3)
real(8), intent(in)  :: stepsize
real(8), intent(out) :: rthetaphi2(3), xyz2(3)
real(8)              :: rand1, rand2                ! random numbers to change theta and phi
real(8), parameter   :: pi = 3.14159265358979323846 ! define pi

call random_number(rand1)                         ! random number to change theta
call random_number(rand2)                         ! same for phi
rthetaphi2(1) = rthetaphi1(1)
rthetaphi2(2) = rthetaphi1(2) + (2*rand1 - 1)*pi*stepsize 
rthetaphi2(3) = rthetaphi1(3) + (2*rand2 - 1)*2*pi*stepsize
xyz2(1)       = rthetaphi2(1)*sin(rthetaphi2(2))*cos(rthetaphi2(3)) + center_xyz(1)
xyz2(2)       = rthetaphi2(1)*sin(rthetaphi2(2))*sin(rthetaphi2(3)) + center_xyz(2)
xyz2(3)       = rthetaphi2(1)*cos(rthetaphi2(2)) + center_xyz(3)
return
end subroutine mc_step



subroutine create_rotMat(fullRotMat, stepsize)
real(8), intent(in)  :: stepsize
real(8), intent(out) :: fullRotMat(3,3)
real(8)              :: rotmat_x(3,3), rotmat_y(3,3), rotmat_z(3,3)
real(8)              :: rand1, rand2, rand3
real(8), parameter   :: pi = 3.14159265358979323846 ! define pi

call random_number(rand1)                               ! random number for angle 1.    USE THE SAME ANGLES FOR ALL COORDINATES -> ROTATION
call random_number(rand2)                               ! same for angle 2
call random_number(rand3)                               ! same for angle 3
!
! Rotation based on the rotation around the cartesian axes. Three random numbers for three random angles
!
rand1 = (2*rand1 - 1)*stepsize*pi  ! angle in radians. Should be sufficiently small
rand2 = (2*rand2 - 1)*stepsize*pi  ! angle in radians. Should be sufficiently small
rand3 = (2*rand3 - 1)*stepsize*pi  ! angle in radians. Should be sufficiently small

rotmat_x(1,1:3) = (/ real(1.0,8), real(0.0,8), real(0.0,8) /)
rotmat_x(2,1:3) = (/ real(0.0,8), cos(rand1), -sin(rand1) /)
rotmat_x(3,1:3) = (/ real(0.0,8), sin(rand1),  cos(rand1) /)

rotmat_y(1,1:3) = (/ cos(rand2),  real(0.0,8), sin(rand2) /)
rotmat_y(2,1:3) = (/ real(0.0,8), real(1.0,8), real(0.0,8) /)
rotmat_y(3,1:3) = (/ -sin(rand2), real(0.0,8), cos(rand2) /)

rotmat_z(1,1:3) = (/ cos(rand3), -sin(rand3),  real(0.0,8) /)
rotmat_z(2,1:3) = (/ sin(rand3),  cos(rand3),  real(0.0,8) /)
rotmat_z(3,1:3) = (/ real(0.0,8), real(0.0,8), real(1.0,8) /)

fullRotMat = matmul(rotmat_z,matmul(rotmat_y,rotmat_x))
return
end subroutine create_rotMat


subroutine create_rotMat_bond(fullRotMat,axis_vector,stepsize)
real(8), intent(in)  :: stepsize
real(8), intent(in)  :: axis_vector(3)
real(8), intent(out) :: fullRotMat(3,3)
real(8)              :: len_vector
real(8)              :: norm_vector(3)
real(8)              :: rand1
real(8), parameter   :: pi = 3.14159265358979323846 ! define pi

call random_number(rand1)          ! random number for angle
rand1 = (2*rand1 - 1)*stepsize*10.0*pi  ! angle in radians. Should be sufficiently small
! normalize vector
len_vector = sqrt(sum(axis_vector(:)**2))
norm_vector(:) = axis_vector(:)/len_vector
! generate rotation matrix
fullRotMat(1,1) = cos(rand1) + norm_vector(1)**2*(1 - cos(rand1))
fullRotMat(1,2) = norm_vector(1)*norm_vector(2)*(1 - cos(rand1)) - norm_vector(3)*sin(rand1)
fullRotMat(1,3) = norm_vector(1)*norm_vector(3)*(1 - cos(rand1)) + norm_vector(2)*sin(rand1)
fullRotMat(2,1) = norm_vector(2)*norm_vector(1)*(1 - cos(rand1)) + norm_vector(3)*sin(rand1)
fullRotMat(2,2) = cos(rand1) + norm_vector(2)**2*(1 - cos(rand1))
fullRotMat(2,3) = norm_vector(2)*norm_vector(3)*(1 - cos(rand1)) - norm_vector(1)*sin(rand1)
fullRotMat(3,1) = norm_vector(3)*norm_vector(1)*(1 - cos(rand1)) - norm_vector(2)*sin(rand1)
fullRotMat(3,2) = norm_vector(3)*norm_vector(2)*(1 - cos(rand1)) + norm_vector(1)*sin(rand1)
fullRotMat(3,3) = cos(rand1) + norm_vector(3)**2*(1 - cos(rand1))

return
end subroutine create_rotMat_bond


subroutine rotate_pos(fullRotMat, xyz1, center_xyz, xyz2)
real(8), intent(in)  :: fullRotMat(3,3), xyz1(3), center_xyz(3)
real(8), intent(out) :: xyz2(3)
integer              :: a
! rotation at ORIGIN. Move to origin, rotate and move back to center
do a = 1,3
  xyz2(a) = fullRotMat(a,1)*(xyz1(1) - center_xyz(1)) + &
            fullRotMat(a,2)*(xyz1(2) - center_xyz(2)) + &
            fullRotMat(a,3)*(xyz1(3) - center_xyz(3)) + &
            center_xyz(a)
end do
return
end subroutine rotate_pos




subroutine rotate_around_axis(axis_vector, center_of_bond, xyz, stepsize)
real(8), intent(in)    :: axis_vector(3)
real(8), intent(in)    :: center_of_bond(3)         ! shift all pos to the origin
real(8), intent(in)    :: stepsize
real(8), intent(inout) :: xyz(3)
real(8)                :: tmp_xyz(3)
real(8)                :: len_vector
real(8)                :: norm_vector(3)
real(8)                :: full_rot(3,3)
real(8)                :: rand1
real(8), parameter     :: pi = 3.14159265358979323846 ! define pi

call random_number(rand1)          ! random number for angle
rand1 = (2*rand1 - 1)*stepsize*10.0*pi  ! angle in radians. Should be sufficiently small

! normalize vector
len_vector = sqrt(sum(axis_vector(:)**2))
norm_vector(:) = axis_vector(:)/len_vector

full_rot(1,1) = cos(rand1) + norm_vector(1)**2*(1 - cos(rand1))
full_rot(1,2) = norm_vector(1)*norm_vector(2)*(1 - cos(rand1)) - norm_vector(3)*sin(rand1)
full_rot(1,3) = norm_vector(1)*norm_vector(3)*(1 - cos(rand1)) + norm_vector(2)*sin(rand1)
full_rot(2,1) = norm_vector(2)*norm_vector(1)*(1 - cos(rand1)) + norm_vector(3)*sin(rand1)
full_rot(2,2) = cos(rand1) + norm_vector(2)**2*(1 - cos(rand1))
full_rot(2,3) = norm_vector(2)*norm_vector(3)*(1 - cos(rand1)) - norm_vector(1)*sin(rand1)
full_rot(3,1) = norm_vector(3)*norm_vector(1)*(1 - cos(rand1)) - norm_vector(2)*sin(rand1)
full_rot(3,2) = norm_vector(3)*norm_vector(2)*(1 - cos(rand1)) + norm_vector(1)*sin(rand1)
full_rot(3,3) = cos(rand1) + norm_vector(3)**2*(1 - cos(rand1))

tmp_xyz(1) = full_rot(1,1)*(xyz(1)-center_of_bond(1)) + &
           & full_rot(1,2)*(xyz(2)-center_of_bond(2)) + full_rot(1,3)*(xyz(3)-center_of_bond(3))
tmp_xyz(2) = full_rot(2,1)*(xyz(1)-center_of_bond(1)) + &
           & full_rot(2,2)*(xyz(2)-center_of_bond(2)) + full_rot(2,3)*(xyz(3)-center_of_bond(3))
tmp_xyz(3) = full_rot(3,1)*(xyz(1)-center_of_bond(1)) + & 
           & full_rot(3,2)*(xyz(2)-center_of_bond(2)) + full_rot(3,3)*(xyz(3)-center_of_bond(3))

xyz(:) = tmp_xyz(:)+center_of_bond(:)
return
end subroutine rotate_around_axis



subroutine per_vector(bond_vec,center_xyz,per_vec)
real(8), intent(in)    :: bond_vec(3)
real(8), intent(in)    :: center_xyz(3)         
real(8), intent(inout) :: per_vec(3)

! find a vector which is perpendicular to the bond vector

! take atomic position and make cross product between atomic position and bond center ?
per_vec(1) = center_xyz(2)*bond_vec(3) - center_xyz(3)*bond_vec(2)
per_vec(2) = center_xyz(3)*bond_vec(1) - center_xyz(1)*bond_vec(3)
per_vec(3) = center_xyz(1)*bond_vec(2) - center_xyz(2)*bond_vec(1)

!!!!
! NEW TRY FOR PERPENDICULAR VECTOR
!!!!
if (sqrt(sum(per_vec(:)**2)) < 1.0D-6) then                ! if no proper perpendicular vector has been found, use alternative strategy
! find smallest vector component -> build cross product with the axis corresponding to that component     bond_center or bond_vector? might be bond_vector!
  if (abs(bond_vec(1)) <= abs(bond_vec(2)) .and. (abs(bond_vec(1)) <= abs(bond_vec(3)))) then        ! x is smallest
    per_vec(1) = real(0.0,8)
    per_vec(2) = +bond_vec(3)
    per_vec(3) = -bond_vec(2)
  else if (abs(bond_vec(2)) <= abs(bond_vec(1)) .and. (abs(bond_vec(2)) <= abs(bond_vec(3)))) then   ! y is smallest
    per_vec(1) = +bond_vec(3)   ! was - before   ... wouldn't lead to perpendicular vector!!!
    per_vec(2) = real(0.0,8)
    per_vec(3) = -bond_vec(1)
  else                                                                                               ! z is smallest
    per_vec(1) = +bond_vec(2)
    per_vec(2) = -bond_vec(1)
    per_vec(3) = real(0.0,8)
  end if

!!!!
! NEW TRY FOR PERPENDICULAR VECTOR
!!!!
  if (sqrt(sum(per_vec(:)**2)) < 1.0D-6) then           ! if not proper perpendicular vector has been found, use alternative strategy
    if (abs(bond_vec(1)) <= abs(bond_vec(2)) .and. (abs(bond_vec(1)) <= abs(bond_vec(3)))) then        ! x is smallest
      per_vec(1) = real(0.0,8)
      if (abs(bond_vec(2)) >= abs(bond_vec(3))) then                                                   ! y larger than z
        per_vec(2) = real(1.0,8)
        per_vec(3) = -1.0D0*bond_vec(2)/bond_vec(3)
      else
        per_vec(2) = -1.0D0*bond_vec(3)/bond_vec(2)
        per_vec(3) = real(1.0,8)
      end if
    else if (abs(bond_vec(2)) <= abs(bond_vec(1)) .and. (abs(bond_vec(2)) <= abs(bond_vec(3)))) then   ! y is smallest
      per_vec(2) = real(0.0,8)
      if (abs(bond_vec(1)) >= abs(bond_vec(3))) then                                                   ! x larger than z
        per_vec(1) = real(1.0,8)
        per_vec(3) = -1.0D0*bond_vec(1)/bond_vec(3)
      else
        per_vec(1) = -1.0D0*bond_vec(3)/bond_vec(1)
        per_vec(3) = real(1.0,8)
      end if
    else                                                                                               ! z is smallest
      per_vec(3) = real(0.0,8)
      if (abs(bond_vec(1)) >= abs(bond_vec(2))) then                                                   ! x larger than y
        per_vec(1) = real(1.0,8)
        per_vec(2) = -1.0D0*bond_vec(1)/bond_vec(2)
      else
        per_vec(1) = -1.0D0*bond_vec(2)/bond_vec(1)
        per_vec(2) = real(1.0,8)
      end if
    end if
  end if
end if
! normalize perpendicular vector
per_vec(:) = per_vec(:)/sqrt(sum(per_vec(:)**2))

return
end subroutine per_vector


subroutine cart_to_frac(vecA,vecB,vecC,pos_cart,pos_frac)
! Transform cart to fractional. Make sure there are within the unit cell. Transform back to cartesian. Return both cart and frac
real(8), intent(in)    :: vecA(3), vecB(3), vecC(3)
real(8), intent(inout) :: pos_cart(3)
real(8), intent(inout) :: pos_frac(3)
real(8)                :: trans_matrix(3,3)
real(8)                :: determinant
real(8)                :: lenA, lenB, lenC, angleBC, angleAC, angleAB, vol
integer                :: t,f
real(8), parameter     :: pi = 3.14159265358979323846

lenA    = sqrt(vecA(1)**2+vecA(2)**2+vecA(3)**2)
lenB    = sqrt(vecB(1)**2+vecB(2)**2+vecB(3)**2)
lenC    = sqrt(vecC(1)**2+vecC(2)**2+vecC(3)**2)
angleBC = ACOS((vecB(1)*vecC(1) + vecB(2)*vecC(2) + vecB(3)*vecC(3))/(lenB*lenC))*180.0D0/pi
angleAC = ACOS((vecA(1)*vecC(1) + vecA(2)*vecC(2) + vecA(3)*vecC(3))/(lenA*lenC))*180.0D0/pi
angleAB = ACOS((vecA(1)*vecB(1) + vecA(2)*vecB(2) + vecA(3)*vecB(3))/(lenA*lenB))*180.0D0/pi
!!!! Calculate the total unit cell volume using the triple product V = a . (b x c) 
!!!vol     = (vecA(1)*vecB(2)*vecC(3) + vecA(2)*vecB(3)*vecC(1) + vecA(3)*vecB(1)*vecC(2) - &
!!!           vecA(3)*vecB(2)*vecC(1) - vecA(1)*vecB(3)*vecC(2) - vecA(2)*vecB(1)*vecC(3))
! deteminant of the matrix comntaining the cell vectors
determinant = vecA(1)*vecB(2)*vecC(3)+vecB(1)*vecC(2)*vecA(3)+vecC(1)*vecA(2)*vecB(3) - &
              vecA(3)*vecB(2)*vecC(1)-vecB(3)*vecC(2)*vecA(1)-vecC(3)*vecA(2)*vecB(1)
!!!! transformation matrix to get fractional coordinates. It is the inverse of the matrix containing the cell vectors
!!!trans_matrix(1,1) = (vecB(2)*vecC(3)-vecB(3)*vecC(2))/determinant
!!!trans_matrix(1,2) = (vecB(3)*vecC(1)-vecB(1)*vecC(3))/determinant
!!!trans_matrix(1,3) = (vecB(1)*vecC(2)-vecB(2)*vecC(1))/determinant
!!!trans_matrix(2,1) = (vecA(3)*vecC(2)-vecA(2)*vecC(3))/determinant
!!!trans_matrix(2,2) = (vecA(1)*vecC(3)-vecA(3)*vecC(1))/determinant
!!!trans_matrix(2,3) = (vecA(2)*vecC(1)-vecA(1)*vecC(2))/determinant
!!!trans_matrix(3,1) = (vecA(2)*vecB(3)-vecA(3)*vecB(2))/determinant
!!!trans_matrix(3,2) = (vecA(3)*vecB(1)-vecA(1)*vecB(3))/determinant
!!!trans_matrix(3,3) = (vecA(1)*vecB(2)-vecA(2)*vecB(1))/determinant
! transformation matrix to get fractional coordinates. It is the inverse of the matrix containing the cell vectors
trans_matrix(1,1) = (vecB(2)*vecC(3)-vecB(3)*vecC(2))/determinant
trans_matrix(1,2) = (vecA(3)*vecC(2)-vecA(2)*vecC(3))/determinant
trans_matrix(1,3) = (vecA(2)*vecB(3)-vecA(3)*vecB(2))/determinant
trans_matrix(2,1) = (vecB(3)*vecC(1)-vecB(1)*vecC(3))/determinant
trans_matrix(2,2) = (vecA(1)*vecC(3)-vecA(3)*vecC(1))/determinant
trans_matrix(2,3) = (vecA(3)*vecB(1)-vecA(1)*vecB(3))/determinant
trans_matrix(3,1) = (vecB(1)*vecC(2)-vecB(2)*vecC(1))/determinant
trans_matrix(3,2) = (vecA(2)*vecC(1)-vecA(1)*vecC(2))/determinant
trans_matrix(3,3) = (vecA(1)*vecB(2)-vecA(2)*vecB(1))/determinant
! frac = cart*trans_matrix
pos_frac(1) = pos_cart(1)*trans_matrix(1,1) + pos_cart(2)*trans_matrix(2,1) + pos_cart(3)*trans_matrix(3,1)
pos_frac(2) = pos_cart(1)*trans_matrix(1,2) + pos_cart(2)*trans_matrix(2,2) + pos_cart(3)*trans_matrix(3,2)
pos_frac(3) = pos_cart(1)*trans_matrix(1,3) + pos_cart(2)*trans_matrix(2,3) + pos_cart(3)*trans_matrix(3,3)
! make sure that all fractional coordinates are within 0 and 1
do f = 1, 3
  t = 0
  do while (t < 1)
    if (pos_frac(f) > 1) then
      pos_frac(f) = pos_frac(f) - 1
    end if
    if (pos_frac(f) < 0) then
      pos_frac(f) = pos_frac(f) + 1
    end if
    if ((0 <= pos_frac(f)) .and. (pos_frac(f) <= 1)) then
      t = 1
    end if
  end do
end do
! Transfrom back to cartesian. This ensures that all FODs are inside the unit cell
pos_cart(1) = pos_frac(1)*vecA(1) + pos_frac(2)*vecB(1) + pos_frac(3)*vecC(1)
pos_cart(2) = pos_frac(1)*vecA(2) + pos_frac(2)*vecB(2) + pos_frac(3)*vecC(2)
pos_cart(3) = pos_frac(1)*vecA(3) + pos_frac(2)*vecB(3) + pos_frac(3)*vecC(3)

return
end subroutine cart_to_frac


subroutine struct_motif(n_point_tot,n_point,radius,nuc_pos,position_point)     
! make initial structural motifs. For atoms and cores in molecules. Johnson-bodies.    Orient them all along z (to maximize 1/r between them)
! n_point_tot characterizes which structural motif shall be used (e.g. 4 = tetrahedron). 
! n_point specifies the point which shall be assigned
! radius is the correpsonding radiius which shall be used (scaling of original structural motif)
! nuc_pos is the position of the respective atomic position
! position_point is the x,y,z coordinates of the point
integer, intent(in)    :: n_point_tot, n_point
real(8), intent(in)    :: radius
real(8), intent(in)    :: nuc_pos(3)
real(8), intent(inout) :: position_point(3)
real(8), parameter     :: pi = 3.14159265358979323846
real(8)                :: motif_1(3)
real(8)                :: motif_2(2,3)
real(8)                :: motif_3(3,3)
real(8)                :: motif_4(4,3)
real(8)                :: motif_5(5,3)
real(8)                :: motif_6(6,3)
real(8)                :: motif_7(7,3)
real(8)                :: motif_8(8,3)
real(8)                :: motif_9(9,3)
real(8)                :: motif_10(10,3)
real(8)                :: motif_13(13,3)
real(8)                :: motif_15(15,3)
real(8)                :: motif_18(18,3)

! Motif 1 -> single point
motif_1(1:3)    = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
! Motif 2 -> two points in a line
motif_2(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_2(2,1:3) = (/ real(0.0,8), real(0.0,8), real(-1.0,8) /)
! Motif 3 -> equiliteral triangle
motif_3(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_3(2,1:3) = (/ real(sqrt(0.75),8), real(0.0,8), real(-0.5,8) /)
motif_3(3,1:3) = (/ real(-1.0*sqrt(0.75),8), real(0.0,8), real(-0.5,8) /)
! Motif 4 -> tetrahedron
motif_4(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_4(2,1:3) = (/ real(sqrt(8.0/9.0),8), real(0.0,8), real(-1.0/3.0,8) /)
motif_4(3,1:3) = (/ real(-1.0*sqrt(2.0/9.0),8), real(sqrt(2.0/3.0),8), real(-1.0/3.0,8) /)
motif_4(4,1:3) = (/ real(-1.0*sqrt(2.0/9.0),8), real(-1.0*sqrt(2.0/3.0),8), real(-1.0/3.0,8) /)
! Motif 5 -> Triangle + two points in a line (for now)
motif_5(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_5(2,1:3) = (/ real(0.0,8), real(0.0,8), real(-1.0,8) /)
motif_5(3,1:3) = (/ real(1.0,8), real(0.0,8), real(0.0,8) /)
motif_5(4,1:3) = (/ real(-0.5,8), real(sqrt(0.75),8), real(0.0,8) /)
motif_5(5,1:3) = (/ real(-0.5,8), real(-1.0*sqrt(0.75),8), real(0.0,8) /)
! Motif 6 -> Octahedron
motif_6(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_6(2,1:3) = (/ real(0.0,8), real(0.0,8), real(-1.0,8) /)
motif_6(3,1:3) = (/ real(+1.0,8), real(0.0,8), real(0.0,8) /)
motif_6(4,1:3) = (/ real(-1.0,8), real(0.0,8), real(0.0,8) /)
motif_6(5,1:3) = (/ real(0.0,8), real(+1.0,8), real(0.0,8) /)
motif_6(6,1:3) = (/ real(0.0,8), real(-1.0,8), real(0.0,8) /)
! Motif 7 -> Pentagon + 2 points in a line
motif_7(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_7(2,1:3) = (/ real(0.0,8), real(0.0,8), real(-1.0,8) /)
motif_7(3,1:3) = (/ real(+1.0,8), real(0.0,8), real(0.0,8) /)
motif_7(4,1:3) = (/ real(sqrt(1.0 - sin(72.0/360.0*(2.0*pi))),8), real(sin(72.0/360.0*(2.0*pi)),8), real(0.0,8) /)  ! angles in radians
motif_7(5,1:3) = (/ real(sqrt(1.0 - sin(72.0/360.0*(2.0*pi))),8), real(-1.0*sin(72.0/360.0*(2.0*pi)),8), real(0.0,8) /)  ! angles in radians
motif_7(6,1:3) = (/ real(-1.0*sin(54.0/360.0*(2.0*pi)),8), real(sqrt(1.0 - sin(54.0/360.0*(2.0*pi))),8), real(0.0,8) /)  ! angles in radians
motif_7(7,1:3) = (/ real(-1.0*sin(54.0/360.0*(2.0*pi)),8), real(-1.0*sqrt(1.0 - sin(54.0/360.0*(2.0*pi))),8), real(0.0,8) /)  ! angles in radians
! Motif 8 -> twisted cube  (coordinates from optimized points on a sphere) - Symmetrized
motif_8(1,1:3) = (/ real( 0.523842100,8), real( 0.641286608,8), real( 0.5595346325,8) /)
motif_8(2,1:3) = (/ real(-0.523842100,8), real(-0.641286608,8), real( 0.5595346325,8) /)
motif_8(3,1:3) = (/ real(-0.824785975,8), real(-0.082731438,8), real( 0.5595346325,8) /)
motif_8(4,1:3) = (/ real( 0.824785975,8), real( 0.082731438,8), real(-0.5595346325,8) /)
motif_8(5,1:3) = (/ real(-0.641286608,8), real( 0.523842100,8), real(-0.5595346325,8) /)
motif_8(6,1:3) = (/ real( 0.641286608,8), real(-0.523842100,8), real( 0.5595346325,8) /)
motif_8(7,1:3) = (/ real(-0.082731438,8), real( 0.824785975,8), real( 0.5595346325,8) /)
motif_8(8,1:3) = (/ real( 0.082731438,8), real(-0.824785975,8), real(-0.5595346325,8) /)
!!!motif_8(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
!!!motif_8(2,1:3) = (/ real( 0.25873213833899444,8), real( 0.54267827151894577,8), real(-0.79909823199737307,8) /)
!!!motif_8(3,1:3) = (/ real(-0.90657749317642811,8), real( 0.19967266394387995,8), real(-0.37181727150544863,8) /)
!!!motif_8(4,1:3) = (/ real(-0.27516286901850467,8), real( 0.90864087400522586,8), real( 0.31409736914875530,8) /)
!!!motif_8(5,1:3) = (/ real( 0.64785957780839931,8), real(-0.74233833397548044,8), real( 0.17091448403512255,8) /)
!!!motif_8(6,1:3) = (/ real( 0.89976917722989169,8), real( 0.40149697813135576,8), real( 0.17092570010873220,8) /)
!!!motif_8(7,1:3) = (/ real(-0.63141670896554891,8), real(-0.70898843880608209,8), real( 0.31408331144729251,8) /)
!!!motif_8(8,1:3) = (/  real(0.0,8),                 real(-0.60115619469748838,8), real(-0.79910264551689303,8) /)
! Motif 9 -> twisted cube + 1 point at the side (coordinates from optimized points on a sphere) - Symmetrized
motif_9(1,1:3) = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real( 1.00000000000000000,8) /)
motif_9(2,1:3) = (/ real( 0.47353004441203156,8),      real( 0.80594647912770400,8),      real( 0.35527598454680115,8) /)
motif_9(3,1:3) = (/ real(-0.47353004441203156,8),      real(-0.80594647912770400,8),      real( 0.35527598454680115,8) /)
motif_9(4,1:3) = (/ real(-0.27326802860984000,8),      real( 0.82178524771742300,8),      real(-0.50000000000000000,8) /)
motif_9(5,1:3) = (/ real( 0.27326802860984000,8),      real(-0.82178524771742300,8),      real(-0.50000000000000000,8) /)
motif_9(6,1:3) = (/ real(-0.66769965479981240,8),      real(-0.22203384174912840,8),      real(-0.71054777716017660,8) /)
motif_9(7,1:3) = (/ real( 0.66769965479981240,8),      real( 0.22203384174912840,8),      real(-0.71054777716017660,8) /)
motif_9(8,1:3) = (/ real( 0.86187128595732240,8),      real(-0.36187760873760340,8),      real( 0.35527598454680115,8) /)
motif_9(9,1:3) = (/ real(-0.86187128595732240,8),      real( 0.36187760873760340,8),      real( 0.35527598454680115,8) /)
!!!motif_9(1,1:3) = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
!!!motif_9(2,1:3) = (/ real( 0.47353337042217880,8),      real( 0.80594111308729377,8),      real( 0.35528167598555882,8) /)
!!!motif_9(3,1:3) = (/ real(-0.27326896977780724,8),      real( 0.82179506712346861,8),      real(-0.49997692483868250,8) /)
!!!motif_9(4,1:3) = (/ real( 0.27326708744187284,8),      real(-0.82177542831137740,8),      real(-0.50001023758797603,8) /)
!!!motif_9(5,1:3) = (/ real(-0.47352671840188432,8),      real(-0.80595184516811424,8),      real( 0.35526619093734274,8) /)
!!!motif_9(6,1:3) = (/ real(-0.66770232163854615,8),      real(-0.22202296511526887,8),      real(-0.71054866713357567,8) /)
!!!motif_9(7,1:3) = (/ real( 0.86187138391010965,8),      real(-0.36188184700754855,8),      real( 0.35527348730462649,8) /)
!!!motif_9(8,1:3) = (/ real( 0.66769698796107868,8),      real( 0.22204471838298795,8),      real(-0.71054688718677750,8) /)
!!!motif_9(9,1:3) = (/ real(-0.86187118800453522,8),      real( 0.36187337046765822,8),      real( 0.35528258395967677,8) /)
! Motif 10 -> twisted cube + 2 points on a line (coordinates from optimized points on a sphere)  -  Symmetrized
motif_10(1,1:3)  = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real( 1.00000000000000000,8) /)
motif_10(2,1:3)  = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real(-1.00000000000000000,8) /)
motif_10(3,1:3)  = (/ real( 0.57697333050312170,8),      real(-0.69888428215415700,8),      real( 0.42268501450052165,8) /)
motif_10(4,1:3)  = (/ real(-0.57697333050312170,8),      real( 0.69888428215415700,8),      real( 0.42268501450052165,8) /)
motif_10(5,1:3)  = (/ real(-0.90216746645970310,8),      real( 0.08620418087144224,8),      real(-0.42268501450052165,8) /)
motif_10(6,1:3)  = (/ real( 0.90216746645970310,8),      real(-0.08620418087144224,8),      real(-0.42268501450052165,8) /)
motif_10(7,1:3)  = (/ real(-0.08620418087144224,8),      real(-0.90216746645970310,8),      real(-0.42268501450052165,8) /)
motif_10(8,1:3)  = (/ real( 0.08620418087144224,8),      real( 0.90216746645970310,8),      real(-0.42268501450052165,8) /)
motif_10(9,1:3)  = (/ real( 0.69888428215415700,8),      real( 0.57697333050312170,8),      real( 0.42268501450052165,8) /)
motif_10(10,1:3) = (/ real(-0.69888428215415700,8),      real(-0.57697333050312170,8),      real( 0.42268501450052165,8) /)
!!!motif_10(1,1:3)  = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
!!!motif_10(2,1:3)  = (/ real(0.0,8), real(0.0,8), real(-1.0,8) /)
!!!motif_10(3,1:3)  = (/ real( 0.57699367577309602,8),      real(-0.69886911015654829,8),      real( 0.42268223829292056,8) /)
!!!motif_10(4,1:3)  = (/ real(-0.90217449849732423,8),      real( 0.08621114899157805,8),      real(-0.42266868939874275,8) /)
!!!motif_10(5,1:3)  = (/ real(-0.08621837568764129,8),      real(-0.90217781352417903,8),      real(-0.42266013221634058,8) /)
!!!motif_10(6,1:3)  = (/ real(-0.57695603747667668,8),      real( 0.69889954929960751,8),      real( 0.42268326733446515,8) /)
!!!motif_10(7,1:3)  = (/ real( 0.69888368699059067,8),      real( 0.57697858771608423,8),      real( 0.42267873049693550,8) /)
!!!motif_10(8,1:3)  = (/ real(-0.69888478216988115,8),      real(-0.57696502104663006,8),      real( 0.42269542473702360,8) /)
!!!motif_10(9,1:3)  = (/ real( 0.90215859614081839,8),      real(-0.08618446515938097,8),      real(-0.42270806956210127,8) /)
!!!motif_10(10,1:3) = (/ real( 0.08620273364716869,8),      real( 0.90215895767649101,8),      real(-0.42270356396564368,8) /)
! Motif 13 -> coordinates from optimized points on a sphere
motif_13(1,1:3)  = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_13(2,1:3)  = (/ real( 0.87314079655928101,8),      real( 0.22583905110122207,8),      real(-0.43199753899152005,8) /)
motif_13(3,1:3)  = (/ real( 0.47331019631752635,8),      real(-0.68447851527191150,8),      real(-0.55449672701815378,8) /)
motif_13(4,1:3)  = (/ real(-0.57272616035461443,8),      real(-0.64913004087085091,8),      real( 0.50061456748459643,8) /)
motif_13(5,1:3)  = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real(-0.99866669404274189,8) /)
motif_13(6,1:3)  = (/ real(-0.73810512959591690,8),      real( 0.38449286367475988,8),      real(-0.55440603998636773,8) /)
motif_13(7,1:3)  = (/ real( 0.24329192826524626,8),      real(-0.92043327370046246,8),      real( 0.30596018095562438,8) /)
motif_13(8,1:3)  = (/ real( 0.85761906555009948,8),      real(-0.24899254278372324,8),      real( 0.44999139398081345,8) /)
motif_13(9,1:3)  = (/ real(-0.35380101558582833,8),      real( 0.81997125397645609,8),      real( 0.44996887140917685,8) /)
motif_13(10,1:3) = (/ real( 0.11548367234538193,8),      real( 0.89444558512125383,8),      real(-0.43200767001042928,8) /)
motif_13(11,1:3) = (/ real( 0.58973454278262583,8),      real( 0.66831999446098211,8),      real( 0.45338896774300574,8) /)
motif_13(12,1:3) = (/ real(-0.94353615478643682,8),      real( 0.12678712996613012,8),      real( 0.30604666009756598,8) /)
motif_13(13,1:3) = (/ real(-0.57349289610518150,8),      real(-0.64985984969719979,8),      real(-0.49878660136556108,8) /)
! Motif 15 -> coordinates from optimized points on a sphere
motif_15(1,1:3)  = (/ real(0.0,8), real(0.0,8), real(+1.0,8) /)
motif_15(2,1:3)  = (/ real( 0.86122864492398410,8),      real( 0.34073707612429743,8),      real( 0.37707221853225092,8) /)  
motif_15(3,1:3)  = (/ real( 0.19984675019212100,8),      real( 0.82276229074494278,8),      real( 0.53209350303277769,8) /)
motif_15(4,1:3)  = (/ real( 0.65712478746789471,8),      real(-0.50707523322986026,8),      real( 0.55772907539330119,8) /)
motif_15(5,1:3)  = (/ real(-0.58354746983944061,8),      real(-0.63989109624301221,8),      real(-0.50001173696849988,8) /)
motif_15(6,1:3)  = (/ real(-0.19983082189999282,8),      real(-0.82279071464998832,8),      real( 0.53205552362837882,8) /)
motif_15(7,1:3)  = (/ real(-0.65712824203376385,8),      real( 0.50708531804185519,8),      real( 0.55771583349603859,8) /)
motif_15(8,1:3)  = (/ real( 0.24022449427532711,8),      real(-0.34019918461013371,8),      real(-0.90915163594741877,8) /)
motif_15(9,1:3)  = (/ real(-0.26280297747251935,8),      real( 0.93946649073329280,8),      real(-0.21985747300689301,8) /)
motif_15(10,1:3) = (/ real(-0.86121795990793426,8),      real(-0.34073763623566755,8),      real( 0.37709612528159886,8) /)
motif_15(11,1:3) = (/ real(-0.24022459455029721,8),      real( 0.34022177836207107,8),      real(-0.90914316428633912,8) /)
motif_15(12,1:3) = (/ real( 0.58354622357858787,8),      real( 0.63989873349067738,8),      real(-0.50000341256914083,8) /)
motif_15(13,1:3) = (/ real( 0.26282599585490291,8),      real(-0.93945045016497686,8),      real(-0.21989848841533588,8) /)
motif_15(14,1:3) = (/ real( 0.91373833124116843,8),      real(-0.22565420872107333,8),      real(-0.33787934064686809,8) /)
motif_15(15,1:3) = (/ real(-0.91374155692039261,8),      real( 0.22565927536461444,8),      real(-0.33786722934984625,8) /)
! Motif 18 (coordinates from optimized points on a sphere) - Symmetrized
motif_18(1,1:3)  = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real( 1.00000000000000000,8) /)
motif_18(2,1:3)  = (/ real( 0.97909461138088549,8),      real( 0.00000000000000000,8),      real(-0.20340706000146119,8) /)
motif_18(3,1:3)  = (/ real(-0.97909461138088549,8),      real( 0.00000000000000000,8),      real(-0.20340706000146119,8) /)
motif_18(4,1:3)  = (/ real( 0.00000000000000000,8),      real(-0.73767726267225664,8),      real( 0.67515220301689050,8) /)
motif_18(5,1:3)  = (/ real( 0.00000000000000000,8),      real( 0.73767726267225664,8),      real( 0.67515220301689050,8) /)
motif_18(6,1:3)  = (/ real( 0.69232390573735890,8),      real( 0.69232390573735890,8),      real( 0.20340706000146119,8) /)
motif_18(7,1:3)  = (/ real(-0.69232390573735890,8),      real(-0.69232390573735890,8),      real( 0.20340706000146119,8) /)
motif_18(8,1:3)  = (/ real(-0.52161824888978940,8),      real( 0.52161824888978940,8),      real(-0.67515220301689050,8) /)
motif_18(9,1:3)  = (/ real( 0.52161824888978940,8),      real(-0.52161824888978940,8),      real(-0.67515220301689050,8) /)
motif_18(10,1:3) = (/ real( 0.00000000000000000,8),      real( 0.00000000000000000,8),      real(-1.00000000000000000,8) /)
motif_18(11,1:3) = (/ real( 0.00000000000000000,8),      real(-0.97909461138088549,8),      real(-0.20340706000146119,8) /)
motif_18(12,1:3) = (/ real( 0.00000000000000000,8),      real( 0.97909461138088549,8),      real(-0.20340706000146119,8) /)
motif_18(13,1:3) = (/ real(-0.73767726267225664,8),      real( 0.00000000000000000,8),      real( 0.67515220301689050,8) /)
motif_18(14,1:3) = (/ real( 0.73767726267225664,8),      real( 0.00000000000000000,8),      real( 0.67515220301689050,8) /)
motif_18(15,1:3) = (/ real( 0.69232390573735890,8),      real(-0.69232390573735890,8),      real( 0.20340706000146119,8) /)
motif_18(16,1:3) = (/ real(-0.69232390573735890,8),      real( 0.69232390573735890,8),      real( 0.20340706000146119,8) /)
motif_18(17,1:3) = (/ real(-0.52161824888978940,8),      real(-0.52161824888978940,8),      real(-0.67515220301689050,8) /)
motif_18(18,1:3) = (/ real( 0.52161824888978940,8),      real( 0.52161824888978940,8),      real(-0.67515220301689050,8) /)

! Single point
if (n_point_tot == 1) then
  position_point(:) = motif_1(1:3)*radius + nuc_pos(1:3)
end if
! Two points in a line
if (n_point_tot == 2) then
  position_point(:) = motif_2(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Three points in an equiliteral triangle
if (n_point_tot == 3) then
  position_point(:) = motif_3(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Tetrahedron
if (n_point_tot == 4) then
  position_point(:) = motif_4(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Triangle + 2 points in a line
if (n_point_tot == 5) then
  position_point(:) = motif_5(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Octahedron
if (n_point_tot == 6) then
  position_point(:) = motif_6(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Pentagon + 2 points in a line
if (n_point_tot == 7) then
  position_point(:) = motif_7(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Twisted cube
if (n_point_tot == 8) then
  position_point(:) = motif_8(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Twisted cube + 1 point at the side
if (n_point_tot == 9) then
  position_point(:) = motif_9(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Twisted cube + 2 points on a line
if (n_point_tot == 10) then
  position_point(:) = motif_10(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Motif 13
if (n_point_tot == 13) then
  position_point(:) = motif_13(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Motif 15
if (n_point_tot == 15) then
  position_point(:) = motif_15(n_point,1:3)*radius + nuc_pos(1:3)
end if
! Motif 18
if (n_point_tot == 18) then
  position_point(:) = motif_18(n_point,1:3)*radius + nuc_pos(1:3)
end if
return
end subroutine struct_motif


! End module definition
end module fodmc_mod




program fodMC
  USE fodmc_mod
  call get_guess()
end program fodMC 
