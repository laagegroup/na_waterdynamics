MODULE declarations

USE types

implicit none

!general
integer    :: h,i,j,k,ia,ih,io,it,m,o,s,t
integer    :: idum
integer :: ndcd=7,npsf=8,ndef=9,nout=10,ntau=11,npdbc=12,ntauh=13,nin=14,nrep=15,ntrap=16,ngr=25,nsit=26,nstat=27,nhb=28,nrest=29
character(len=6) :: inputformat
!inputformat 1 = namd psf and dcd, 2 = dlpoly HISTORY
real(kind=4),parameter :: pi=3.14159265
real(kind=4) :: rdum
integer,parameter             :: nmax = 800
integer,parameter             :: nmaxdat = 1000

!input
character(len=40) :: dcdfile,psffile,sitedeffile,exclvolfile,fixcrdfile,atomtypefile,trapwatfile  ! input files
character(len=4)  :: vector   ! choice between OH vector or vector orthogonal to water plane
character(len=3)  :: solutesegname,watertype
character(len=4)  :: htype,otype   ! H and O atomtypes in water
integer           :: icentr  ! index of atom approximately at centre of solute
real(kind=4)      :: roughcut,r_oh,invr_oh  ! rouchcut = radius of sphere around icentr including all of solute, roughcut + 5 A used to find first approximation of hydration shell
integer           :: nmaxsol ! maximum number of solvent vectors in system, for memory allocation
integer           :: restart

!dcd
integer      :: nstep ! Number of steps (frames in dcd)
integer      :: natms 
integer      :: dfr ! Number of timesteps between two frames
real(kind=4) :: cell(3)	
real(kind=4),dimension(:),allocatable       :: x,y,z
integer      :: istarttrj,iendtrj ! for only looking at subpart of trajectory
real(kind=4) :: tstarttrj,tendtrj

!dlpoly/lammps
real(kind=4),dimension(:),allocatable :: xtemp,ytemp,ztemp
character(len=4),dimension(:),allocatable :: atomtypetemp,atomtypetemp2,atomtypeschar
integer,dimension(:),allocatable :: atomnumtemp,atomtypesnum
integer :: natmsfix,ntypescorr,natmsfree

!psf
character(len=4),dimension(:),allocatable     :: resname
character(len=4),dimension(:),allocatable     :: atomname
character(len=4),dimension(:),allocatable     :: atomtype
character(len=3),dimension(:),allocatable     :: segname
integer,dimension(:),allocatable              :: resindex            

!site definitions
integer             :: nsitetypes
type(site),dimension(:),allocatable :: sitedata ! data from site definitions file
type(criteria)        :: jump_criterion ! hbond criteria*jump_criterion = hbond criteria for jump (stable state)

!system definition
integer   :: nsites,nhyd,nox
integer,dimension(:),allocatable         :: sitetypeindex ! sitetype for each site atom
integer,dimension(:),allocatable         :: soluteatom ! site atom indices i.e. not including non-site solute atoms
integer,dimension(:),allocatable         :: hwatatom ! water H atom indices
integer,dimension(:),allocatable         :: owatatom ! water O atom indices
integer,dimension(:),allocatable         :: donorlist !index of O for each H in water
integer,dimension(:),allocatable         :: pdonorlist !index of donor for each D-H site in solute
integer,dimension(:),allocatable         :: site_wat ! site index for each water-oxygen or each OH vector
integer,dimension(:),allocatable         :: site_wat_old 
integer,dimension(:),allocatable       :: poptcf,popjcf,popblkjcf,popblktcf ! population of each site
integer,dimension(:),allocatable         :: hbacc
integer,dimension(:),allocatable         :: state
integer,parameter                        :: maxpossh = 30 ! max number of OH vectors per site
logical,dimension(:),allocatable         :: bulk
real(kind=4)                             :: tconvert

!calc
real(kind=4)       :: dx,dy,dz,dr2  ! interatomic distances
real(kind=4)       :: theta         ! angles

! tcf for reorientation
integer                                  :: nstepblock,nblock,nt0max,ntdecorr,ntmax,delt,istep,iblock,t0,&
                                            tdecorr,tmax
real(kind=4)                             :: dt,roh2
character(len=2)                         :: charblock
character(len=4)                         :: charsite
logical                                  :: new
integer,dimension(:),allocatable         :: nt0
integer,dimension(:,:),allocatable       :: time0,site0
real(kind=4),dimension(:,:),allocatable  :: dx0,dy0,dz0,tcf,tcfblk,tcfsq,tcfstdev,norm,normblk
real(kind=4),dimension(:,:,:),allocatable  ::normjblk
real(kind=4),dimension(:),allocatable    :: tau_reor,beta_reor,tau_reor_long,beta_reor_long,tau_reor2
real(kind=4)                             :: tcfdt
integer                                  :: itcfdt
real(kind=4),dimension(:,:),allocatable  :: p2sqblk,chi4blk,chi4
real(kind=4),dimension(:,:),allocatable  :: cjumpsqblk,chijumpblk,chijump
real(kind=4),dimension(:),allocatable    :: amp_chi4,tau_chi4

! tcf (1st order)
real(kind=4),dimension(:,:),allocatable  :: tcf1,tcf1blk,tcf1sq,tcf1stdev

! jump cf
integer           :: at,a0,st,s0
logical,dimension(:,:),allocatable           :: jumped
integer,dimension(:,:),allocatable        :: abs,acc0,state0,gr0
real(kind=4),dimension(:,:,:),allocatable :: jcf,jcfblk,jcfsq,jcfstdev
real(kind=4),dimension(:),allocatable     :: tau_jump,beta_jump,tau_jump_sol,tau_jump_water,tau_jump2
real(kind=4),dimension(:),allocatable     :: beta_jump_sol,beta_jump_wat
real(kind=4),dimension(:),allocatable     :: popinfwater,popinfsol
real(kind=4)                              :: k_jump

! tcf for reorientation between jumps
real(kind=4),dimension(:,:),allocatable   :: tcffr,tcffrblk,tcffrsq,tcffrstdev,normfrblk
real(kind=4),dimension(:),allocatable     :: tau_frame,tau_frame2,beta_frame
real(kind=4),dimension(:,:),allocatable   :: tcffrw,tcffrwblk,tcffrwsq,tcffrwstdev,normfrwblk

! tcf binned by groove width (o4dist)
integer                                   :: no4bins,bo4,njump,naccsites
real(kind=4)                              :: o4distmin,o4distmax,binwidth,o4distave
character(len=4)                          :: charo4bin
integer,dimension(:),allocatable          :: o4bin
integer,dimension(:,:),allocatable        :: o4bin0,jump_trj,bo40,otheroh0
real(kind=4),dimension(:),allocatable     :: o4dist,o4distold
real(kind=4),dimension(:,:),allocatable   :: o4dist_trj
real(kind=4),dimension(:,:,:),allocatable :: tcf_o4,tcfblk_o4,normblk_o4,o4distcollect
real(kind=4),dimension(:,:,:,:),allocatable :: normjblk_o4
real(kind=4),dimension(:,:,:,:,:),allocatable :: jcf_o4,jcfblk_o4,jcfsq_o4,jcfstdev_o4
real(kind=4),dimension(:,:,:),allocatable :: tcfsq_o4,tcfstdev_o4
type(o4indices),dimension(:),allocatable  :: o4
logical,dimension(:,:),allocatable        :: incl_origin
logical,dimension(:,:),allocatable        :: cut
integer,dimension(:),allocatable          :: accatom,accsitetypeindex

! o4dist tcf and distribution
integer                                   :: no4pairs,ncount,firsto4pair
real(kind=4),dimension(:),allocatable     :: var0,o4dist_long,o4dist_s
real(kind=4),dimension(:,:),allocatable   :: o4dist0,var,o4dist_all,o4dist_t

! fit 
integer                                    :: popcut
integer                                    :: nfiteff,ntfit,itfitmin,itfitmax,idata,tfitmin,tfitmax
integer                                    :: ntfit2,itfitmin2,itfitmax2,tfitmin2,tfitmax2,ntfitmax
real*4                                     :: a,b,siga,sigb,chi2,q
real*4,dimension(:),allocatable            :: xfit,yfit,sig

! errors
real(kind=4),dimension(:),allocatable      :: error_tau_jump,error_tau_frame,error_tau_reor
real(kind=4),dimension(:),allocatable      :: error_tau_jump2,error_tau_frame2,error_tau_reor2
real(kind=4),dimension(:),allocatable      :: error_tau_jump_sol,error_tau_jump_water

! distribution
integer             :: nbins,nbinshist
real(kind=4)       :: tdistmin,tdistmax,histmin,histmax ! boundaries of distribution
real(kind=4),dimension(:),allocatable :: histcoord_jump,histcoord_reor,histcoord_reor_long,weight_jcf,weight_tcf
real(kind=4),dimension(:),allocatable :: weight_jcf_sol,weight_jcf_water,weight_tsev
real(kind=4),dimension(:),allocatable :: hist_jump,hist_reor,hist_reor_long
real(kind=4),dimension(:),allocatable :: hist_jump_water,hist_jump_sol,histcoord_jump_water,histcoord_jump_sol
real(kind=4),dimension(:),allocatable :: hist_tsev,histcoord_tsev

! spatial correlation
real(kind=4),dimension(:,:),allocatable   :: vardeltal0
real(kind=4),dimension(:,:,:),allocatable :: vardeltal,stcfblk,stcf,stcfstdev,stcfsq,deltal0

! tsev factor
real(kind=4)          :: r_ts,theta_ts
integer               :: nsolatoms ! total no. of atoms in solute .ne. nsites
real(kind=4),dimension(:),allocatable :: r_ooinit,r_osol
integer,dimension(:),allocatable :: poprho_tsev,solutetypeindex,soluteatomall ! indices of all solute atoms, .ne. soluteatom which only contains those that are defined as sites
real(kind=4)          :: d_ooinit2,d_ooinit,radius2,radius,r_osol2,r_oinitsol2,length,r2,temp
integer               :: is,f,nsolutetypes,n,nh
logical,dimension(:),allocatable          :: neighbour
real(kind=4),dimension(:,:),allocatable   :: bin,pop_ind,ihrho_tsev
real(kind=4),dimension(:),allocatable     :: rho_tsev,irho_tsev,rho_tsevsq,rho_tsevstdev,centrep
real(kind=4),dimension(:,:,:),allocatable :: itrho_tsev
type(exclvol),dimension(:),allocatable    :: exclvoldata

!water neighbours
real(kind=4)                          :: nwatneighboursi
real(kind=4),dimension(:),allocatable :: nwatneighbours,nwatneighbours_blk,nwatneighboursq,nwatneighbourstdev


!misc
logical              :: lbulk
real(kind=4)         :: r0max
!acceptors minor groove
integer ::naccmg
integer,dimension(:),allocatable ::listacc

!groove
integer,parameter  :: ndist=9 !number of O4O4 groove widths
real(kind=4),dimension(ndist) ::groove
integer,dimension(:),allocatable ::groove_ind
real*4 ::gmin,gmax,bing,halfbin
integer ::nbing
integer ::gr,g0,nskip



END MODULE declarations
