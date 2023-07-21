module  module_mp_soil_init

      implicit none

      private

      public init_soil_3_mp,initslc

contains

   subroutine init_soil_3_mp (im, tsk , tmn , smois , tslb , &
                                 st_input , sm_input ,    &
                                 zs , dzs , &
                                 st_levels_input , sm_levels_input , &
                                 num_soil_layers , num_st_levels_input , num_sm_levels_input ,  &
                                 num_st_levels_alloc , num_sm_levels_alloc)

      integer , intent(in) :: num_soil_layers , &
                              num_st_levels_input , num_sm_levels_input , &
                              num_st_levels_alloc , num_sm_levels_alloc

      integer , intent(in) :: im

      integer , dimension(1:num_st_levels_input) , intent(inout) :: st_levels_input
      integer , dimension(1:num_sm_levels_input) , intent(inout) :: sm_levels_input

      real , dimension(1:im,1:num_st_levels_alloc,1:1) , intent(inout) :: st_input
      real , dimension(1:im,1:num_sm_levels_alloc,1:1) , intent(inout) :: sm_input

      real , dimension(1:im,1:1) , intent(in) :: tmn
      real , dimension(1:im,1:1) , intent(inout) :: tsk
      real , dimension(num_soil_layers) :: zs , dzs       

      real , dimension(1:im,num_soil_layers,1:1) , intent(out) :: tslb , smois

      real , allocatable , dimension(:) :: zhave

      logical :: debug_print = .false.
      integer :: i , j , l , lout , lin , lwant , lhave, k
      real :: temp

      !  allocate the soil layer array used for interpolating.      

      if ( ( num_st_levels_input .le. 0 ) .or. &
           ( num_sm_levels_input .le. 0 ) ) then
         write (0, fmt='(a)')&
            'no input soil level data (either temperature or moisture, or both are missing).'
      else
           if (debug_print) write(0, fmt='(a)') ' assume non-ruc lsm input'
           allocate ( zhave( max(num_st_levels_input,num_soil_layers)  ) )
      end if
      

         do i = 1 , im
            st_input(i,1,1) = tsk(i,1)
            st_input(i,num_st_levels_input+2,1) = tmn(i,1)
         end do

         do i = 1 , im
            sm_input(i,1,1) = (sm_input(i,2,1)-sm_input(i,3,1))/   &
                              (st_levels_input(2)-st_levels_input(1))*st_levels_input(1)+  &
                              sm_input(i,2,1)
            sm_input(i,num_sm_levels_input+2,1) = sm_input(i,num_sm_levels_input+1,1)
         end do

         zhave(1) = 0.

         DO l = 1 , num_st_levels_input
            zhave(l+1) = st_levels_input(l) / 100.
         END DO

         zhave(num_st_levels_input+2) = 300. / 100.

      !  interpolate between the layers we have (zhave) and those that we want
      !  (zs).

      z_wantt_2 : do lwant = 1 , num_soil_layers
         z_havet_2 : do lhave = 1 , num_st_levels_input +2
            if ( ( zs(lwant) .ge. zhave(lhave  ) ) .and. &
                 ( zs(lwant) .le. zhave(lhave+1) ) ) then
                  do i = 1 , im
                     tslb(i,lwant,1)= ( st_input(i,lhave,1 ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                        st_input(i,lhave+1,1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                ( zhave(lhave+1) - zhave(lhave) )
                  end do
               exit z_havet_2
            end if
         end do z_havet_2
      end do z_wantt_2


      !  here are the levels that we have from the input for moisture.      

           !  interpolate between the layers we have (zhave) and those that we
           !  want (zs).      

         zhave(1) = 0.
         do l = 1 , num_sm_levels_input
            zhave(l+1) = sm_levels_input(l) / 100.
         end do
         zhave(num_sm_levels_input+2) = 300. / 100.

      z_wantm_2 : do lwant = 1 , num_soil_layers
         z_havem_2 : do lhave = 1 , num_sm_levels_input +2
            if ( ( zs(lwant) .ge. zhave(lhave  ) ) .and. &
                 ( zs(lwant) .le. zhave(lhave+1) ) ) then
                  do i = 1 , im
                     smois(i,lwant,1)= ( sm_input(i,lhave,1 ) * ( zhave(lhave+1) - zs   (lwant) ) + &
                                         sm_input(i,lhave+1,1) * ( zs   (lwant  ) - zhave(lhave) ) ) / &
                                                                 ( zhave(lhave+1) - zhave(lhave) )
                  end do
               exit z_havem_2
            end if
         end do z_havem_2
      end do z_wantm_2

      deallocate (zhave)

   end subroutine init_soil_3_mp

  subroutine initslc(im,nzs, isltyp, ivgtyp,    &
                     sh2o, tslb, smois)

   use noahmp_tables

   implicit none

   integer, intent(in)  ::   im
   integer, intent(in)  ::   nzs

   real, dimension(1:im, 1:nzs, 1:1 ), intent(in) :: tslb, smois
   integer, dimension(1:im, 1:1),intent(inout)    :: isltyp,ivgtyp

   real, dimension( 1:im,1:nzs,1:1), intent(out)  :: sh2o

   !-- local

   real, dimension ( 1:nzs ) :: soiliqw

   integer ::  i,j,l,itf,jtf
   real    ::  riw,xlmelt,tln,dqm,ref,psis,qmin,bclh

   integer :: errflag

        riw=900.*1.e-3
        xlmelt=3.35e+5

! for fim

   itf=im  !  min0(ite,ide-1)
   jtf=1   !  min0(jte,jde-1)

   errflag = 0

     do i = 1,im

       if ( isltyp( i,1 ) .lt. 0 ) then
         errflag = 1
         print *, &
         "slc init: out of range isltyp ",i,isltyp( i,1 )
       endif
     enddo ! im loop

     if ( errflag .eq. 1 ) then
       print *, "slc out of range value of isltype"
     endif

     do i=1,im

       ! in zobler classification isltyp=0 for water. statsgo classification
       ! has isltyp=14 for water

       if (isltyp(i,1) == 0) isltyp(i,1)=14
     
         dqm    = smcmax_table   (isltyp(i,1)) - &
                  smcdry_table   (isltyp(i,1))
         ref    = smcref_table   (isltyp(i,1))
         psis   = - psisat_table (isltyp(i,1))
         qmin   = smcdry_table   (isltyp(i,1))
         bclh   = bexp_table     (isltyp(i,1))

!        mavail(i,j) = max(0.00001,min(1.,(smois(i,1,j)-qmin)/(ref-qmin)))

         do l=1,nzs
           !-- for land points initialize soil ice ?
           tln=log(tslb(i,l,1)/273.15)
          
           if(tln.lt.0.) then
             soiliqw(l)=(dqm+qmin)*(xlmelt*                        &
            (tslb(i,l,1)-273.15)/tslb(i,l,1)/9.81/psis)            &
            **(-1./bclh)
            !**(-1./bclh)-qmin
            soiliqw(l)=max(0.,soiliqw(l))
            soiliqw(l)=min(soiliqw(l),smois(i,l,1))
            sh2o(i,l,1)=soiliqw(l)
!           smfr3d(i,l,j)=(smois(i,l,j)-soiliqw(l))/riw  
           else
             sh2o(i,l,1)=smois(i,l,1)
           endif
         enddo

       enddo

       end subroutine initslc
end module module_mp_soil_init
