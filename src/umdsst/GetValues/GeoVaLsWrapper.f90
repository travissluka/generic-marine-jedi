module geovals_wrapper

use atlas_module
use iso_c_binding
use kinds
use ufo_geovals_mod_c, only: ufo_geovals_registry
use ufo_geovals_mod, only: ufo_geovals
use ufo_locations_mod
use datetime_mod, only: datetime, c_f_datetime

implicit none

contains

subroutine geovals_wrapper_fill( c_key_geovals, c_locs, c_t1, c_t2, c_field) &
       bind(c, name="geovals_wrapper_fill_f90")
    integer(c_int),  intent(inout) :: c_key_geovals
    type(c_ptr), value, intent(in) :: c_locs
    type(c_ptr),        intent(in) :: c_t1
    type(c_ptr),        intent(in) :: c_t2
    type(c_ptr), value, intent(in) :: c_field

    type(ufo_geovals), pointer :: geovals
    type(ufo_locations)        :: locs
    type(atlas_field)          :: field
    type(datetime)             :: t1, t2
    real(kind_real),    pointer:: field_data(:,:)
    logical(c_bool), allocatable :: time_mask(:)

    integer :: ivar, nval, i

    ! get fortran version fo the passed in C arguments
    call ufo_geovals_registry%get(c_key_geovals, geovals)
    locs = ufo_locations(c_locs)
    field = atlas_field(c_field)
    call field%data(field_data)
    call c_f_datetime(c_t1, t1)
    call c_f_datetime(c_t2, t2)

    ! calculate time mask
    allocate(time_mask(locs%nlocs()))
    call locs%get_timemask(t1, t2, time_mask)

    ! initialize geovals
    ivar=1
    nval=1
    if (.not. geovals%linit) then
        geovals%geovals(ivar)%nval = nval
        allocate(geovals%geovals(ivar)%vals(nval, geovals%geovals(ivar)%nlocs))
        geovals%geovals(ivar)%vals = 0.0
        geovals%linit = .true.
    end if

    ! fill the geovals, obeying the time masking
    do i=1, size(time_mask)
        if (time_mask(i)) then
            geovals%geovals(ivar)%vals(:,i) = field_data(:,i)
        end if
    end do
end subroutine

end module