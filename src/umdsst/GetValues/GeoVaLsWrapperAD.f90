module geovals_wrapper_ad

use atlas_module
use iso_c_binding
use kinds
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_locs_mod, only: ufo_locs, ufo_locs_time_mask
use datetime_mod, only: datetime, c_f_datetime

implicit none

contains

subroutine geovals_wrapper_ad_fill( c_key_geovals, c_key_locs, c_t1, c_t2, c_field) &
       bind(c, name="geovals_wrapper_ad_fill_f90")
    integer(c_int),  intent(inout) :: c_key_geovals
    integer(c_int),  intent(inout) :: c_key_locs
    type(c_ptr),        intent(in) :: c_t1
    type(c_ptr),        intent(in) :: c_t2    
    type(c_ptr), value, intent(in) :: c_field

    type(ufo_geovals), pointer :: geovals
    type(ufo_locs),    pointer :: locs
    type(atlas_field)          :: field
    type(datetime)             :: t1, t2
    real(kind_real),    pointer:: field_data(:,:)
    logical, allocatable :: time_mask(:)

    integer :: ivar, nval, i

    ! get fortran version fo the passed in C arguments
    call ufo_geovals_registry%get(c_key_geovals, geovals)
    call ufo_locs_registry%get(c_key_locs, locs)
    field = atlas_field(c_field)
    call field%data(field_data)
    call c_f_datetime(c_t1, t1)
    call c_f_datetime(c_t2, t2)
  
    ! calculate time mask
    call ufo_locs_time_mask(locs, t1, t2, time_mask)

    ! initialize field_data
    ivar=1
    nval=1
    field_data = 0.0   

    ! fill the field_data from geovals, obeying the time masking
    do i=1, size(time_mask)
        if (time_mask(i)) then
            field_data(:,i) = geovals%geovals(ivar)%vals(:,i)
        end if
    end do
end subroutine

end module