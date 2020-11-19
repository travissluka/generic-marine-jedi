module locations_wrapper

use atlas_module
use iso_c_binding
use kinds
use ufo_locs_mod_c, only: ufo_locs_registry
use ufo_locs_mod, only: ufo_locs

implicit none

contains

subroutine locations_wrapper_fill( c_key_loc, c_field) bind(c, name="locations_wrapper_fill_f90")
    integer(c_int),  intent(inout) :: c_key_loc
    type(c_ptr), value, intent(in) :: c_field

    type(ufo_locs),     pointer :: locs
    type(atlas_field)           :: field
    
    real(kind_real),    pointer :: lonlat(:,:)
    integer :: i

    ! get fortran version of the passed in C arguments
    call ufo_locs_registry%get(c_key_loc, locs)
    field = atlas_field(c_field)
    call field%data(lonlat)
    
    ! fill in the lat/lon
    do i = 1, size(lonlat, dim=2)
        lonlat(1,i) = locs%lon(i)
        lonlat(2,i) = locs%lat(i)
    end do

end subroutine

end module