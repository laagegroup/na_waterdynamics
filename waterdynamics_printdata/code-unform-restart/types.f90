MODULE types

implicit none

!derived types
type site
  character(len=4) :: siteatomtype,siteresname
  character(len=3) :: sitetype
  integer          :: iacceptor,idonor,ihydrophobe
  real(kind=4)     :: d_ad_acc,d_ah_acc,th_hda_acc,d_ad_don,d_ah_don,th_hda_don,d_hydr
end type site

type criteria
  real(kind=4)     :: factor_ad,factor_ah,factor_hda
end type criteria

type exclvol
  character(len=4) :: soluteatomtype,soluteresname
  real(kind=4)     :: r_excl,r2_excl
end type exclvol

type o4indices
  integer          :: bi,bj,fi,fj
end type o4indices

END MODULE types
