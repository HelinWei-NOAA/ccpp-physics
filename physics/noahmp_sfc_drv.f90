!>  \file sfc_noahmp_drv.f
!!  This file contains the NoahMP land surface scheme driver.

!>\defgroup NoahMP_LSM NoahMP LSM Model
!! \brief This is the NoahMP LSM driver module, with the functionality of 
!! preparing variables to run the NoahMP LSM subroutine noahmp_sflx(), calling NoahMP LSM and post-processing
!! variables for return to the parent model suite including unit conversion, as well
!! as diagnotics calculation.

!> This module contains the CCPP-compliant NoahMP land surface model driver.
      module noahmpdrv

      implicit none

      private

      public :: noahmpdrv_init, noahmpdrv_run, noahmpdrv_finalize

      contains

!> \ingroup NoahMP_LSM
!! \brief This subroutine is called during the CCPP initialization phase and calls set_soilveg() to 
!! initialize soil and vegetation parameters for the chosen soil and vegetation data sources.
!! \section arg_table_noahmpdrv_init Argument Table
!! \htmlinclude noahmpdrv_init.html
!!
      subroutine noahmpdrv_init(me, isot, ivegsrc, nlunit, errmsg,errflg)
        
        use set_soilveg_mod,  only: set_soilveg
        
        implicit none
      
        integer,              intent(in)  :: me, isot, ivegsrc, nlunit
        character(len=*),     intent(out) :: errmsg
        integer,              intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ivegsrc /= 1) then
          errmsg = 'The NOAHMP LSM expects that the ivegsrc physics '// &
                  'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if
        if (isot /= 1) then
          errmsg = 'The NOAHMP LSM expects that the isot physics '// &
                  'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if

        !--- initialize soil vegetation
        call set_soilveg(me, isot, ivegsrc, nlunit)
      
      end subroutine noahmpdrv_init

      subroutine noahmpdrv_finalize
      end subroutine noahmpdrv_finalize

!> \ingroup NoahMP_LSM
!! \brief This subroutine is the main CCPP entry point for the NoahMP LSM.
!! \section arg_table_noahmpdrv_run Argument Table
!! \htmlinclude noahmpdrv_run.html
!!
!! \section general_noahmpdrv NoahMP Driver General Algorithm
!!  @{
!!    - Initialize CCPP error handling variables.
!!    - Set a flag to only continue with each grid cell if the fraction of land is non-zero.
!!    - This driver may be called as part of an iterative loop. If called as the first "guess" run, 
!!        save land-related prognostic fields to restore.
!!    - Initialize output variables to zero and prepare variables for input into the NoahMP LSM.
!!    - Call transfer_mp_parameters() to fill a derived datatype for input into the NoahMP LSM.
!!    - Call noahmp_options() to set module-level scheme options for the NoahMP LSM.
!!    - If the vegetation type is ice for the grid cell, call noahmp_options_glacier() to set 
!!        module-level scheme options for NoahMP Glacier and call noahmp_glacier().
!!    - For other vegetation types, call noahmp_sflx(), the entry point of the NoahMP LSM.
!!    - Set output variables from the output of noahmp_glacier() and/or noahmp_sflx().
!!    - Call penman() to calculate potential evaporation.
!!    - Calculate the surface specific humidity and convert surface sensible and latent heat fluxes in W m-2 from their kinematic values.
!!    - If a "guess" run, restore the land-related prognostic fields.
!                                                                       !
!-----------------------------------
  subroutine noahmpdrv_run                                          &
!...................................
!  ---  inputs:
    ( im, km, itime, ps, u1, v1, t1, q1, soiltyp, vegtype,       &
      sigmaf, dlwflx, dswsfc, snet, delt, tg3, cm, ch,           &
      prsl1, prslki, zf, dry, wind, slopetyp,                    &
      shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
      idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc, iopt_frz,   &
      iopt_inf, iopt_rad, iopt_alb, iopt_snf, iopt_tbot,         &
      iopt_stc, xlatin, xcoszin, iyrlen, julian,                 &
      rainn_mp, rainc_mp, snow_mp, graupel_mp, ice_mp,           &
      con_hvap, con_cp, con_jcal, rhoh2o, con_eps, con_epsm1,    &
      con_fvirt, con_rd, con_hfus,                               &

!  ---  in/outs:
      weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
      canopy, trans, tsurf, zorl,                                &

! --- Noah MP specific

      snowxy, tvxy, tgxy, canicexy, canliqxy, eahxy, tahxy, cmxy,&
      chxy, fwetxy, sneqvoxy, alboldxy, qsnowxy, wslakexy, zwtxy,&
      waxy, wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy,   &
      rtmassxy, stmassxy, woodxy, stblcpxy, fastcpxy, xlaixy,    &
      xsaixy, taussxy, smoiseq, smcwtdxy, deeprechxy, rechxy,    &
      albdvis, albdnir,  albivis,  albinir,emiss,                &

!  ---  outputs:
      sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
      cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
      smcwlt2, smcref2, wet1, t2mmp, q2mp, errmsg, errflg)     

  use machine ,   only : kind_phys
  use funcphys,   only : fpvs

  use module_sf_noahmplsm
  use module_sf_noahmp_glacier
  use noahmp_tables, only : isice_table, co2_table, o2_table,       &
                            isurban_table,smcref_table,smcdry_table,   &
                            smcmax_table,co2_table,o2_table,           &
                            saim_table,laim_table

  implicit none
      
  real(kind=kind_phys), parameter :: a2      = 17.2693882
  real(kind=kind_phys), parameter :: a3      = 273.16
  real(kind=kind_phys), parameter :: a4      = 35.86
  real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)
      
  real, parameter                 :: undefined       =  -1.e36   ! TODO change to smaller value

  integer, parameter              :: nsoil   = 4   ! hardwired to Noah
  integer, parameter              :: nsnow   = 3   ! max. snow layers

  real(kind=kind_phys), save  :: zsoil(nsoil)
  data zsoil  / -0.1, -0.4, -1.0, -2.0 /

!
!  ---  CCPP interface fields (in call order)
!

  integer                                 , intent(in)    :: im         ! horiz dimension and num of used pts
  integer                                 , intent(in)    :: km         ! vertical soil layer dimension
  integer                                 , intent(in)    :: itime      ! NOT USED
  real(kind=kind_phys), dimension(im)     , intent(in)    :: ps         ! surface pressure [Pa]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: u1         ! u-component of wind [m/s]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: v1         ! u-component of wind [m/s]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: t1         ! layer 1 temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: q1         ! layer 1 specific humidity [kg/kg]
  integer             , dimension(im)     , intent(in)    :: soiltyp    ! soil type (integer index)
  integer             , dimension(im)     , intent(in)    :: vegtype    ! vegetation type (integer index)
  real(kind=kind_phys), dimension(im)     , intent(in)    :: sigmaf     ! areal fractional cover of green vegetation
  real(kind=kind_phys), dimension(im)     , intent(in)    :: dlwflx     ! downward longwave radiation [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: dswsfc     ! downward shortwave radiation [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: snet       ! total sky sfc netsw flx into ground[W/m2]
  real(kind=kind_phys)                    , intent(in)    :: delt       ! time interval [s]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: tg3        ! deep soil temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: cm         ! surface exchange coeff for momentum [-]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: ch         ! surface exchange coeff heat & moisture[-]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: prsl1      ! sfc layer 1 mean pressure [Pa]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: prslki     ! to calculate potential temperature
  real(kind=kind_phys), dimension(im)     , intent(in)    :: zf         ! height of bottom layer [m]
  logical             , dimension(im)     , intent(in)    :: dry        ! = T if a point with any land
  real(kind=kind_phys), dimension(im)     , intent(in)    :: wind       ! wind speed [m/s]
  integer             , dimension(im)     , intent(in)    :: slopetyp   ! surface slope classification
  real(kind=kind_phys), dimension(im)     , intent(in)    :: shdmin     ! min green vegetation coverage [fraction]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: shdmax     ! max green vegetation coverage [fraction]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: snoalb     ! upper bound on max albedo over deep snow
  real(kind=kind_phys), dimension(im)     , intent(inout) :: sfalb      ! mean surface albedo [fraction]
  logical             , dimension(im)     , intent(in)    :: flag_iter  !
  logical             , dimension(im)     , intent(in)    :: flag_guess !
  integer                                 , intent(in)    :: idveg      ! option for dynamic vegetation
  integer                                 , intent(in)    :: iopt_crs   ! option for canopy stomatal resistance
  integer                                 , intent(in)    :: iopt_btr   ! option for soil moisture factor for stomatal resistance
  integer                                 , intent(in)    :: iopt_run   ! option for runoff and groundwater
  integer                                 , intent(in)    :: iopt_sfc   ! option for surface layer drag coeff (ch & cm)
  integer                                 , intent(in)    :: iopt_frz   ! option for supercooled liquid water (or ice fraction)
  integer                                 , intent(in)    :: iopt_inf   ! option for frozen soil permeability
  integer                                 , intent(in)    :: iopt_rad   ! option for radiation transfer
  integer                                 , intent(in)    :: iopt_alb   ! option for ground snow surface albedo
  integer                                 , intent(in)    :: iopt_snf   ! option for partitioning  precipitation into rainfall & snowfall
  integer                                 , intent(in)    :: iopt_tbot  ! option for lower boundary condition of soil temperature
  integer                                 , intent(in)    :: iopt_stc   ! option for snow/soil temperature time scheme (only layer 1)
  real(kind=kind_phys), dimension(im)     , intent(in)    :: xlatin     ! latitude
  real(kind=kind_phys), dimension(im)     , intent(in)    :: xcoszin    ! cosine of zenith angle
  integer                                 , intent(in)    :: iyrlen     ! year length [days]
  real(kind=kind_phys)                    , intent(in)    :: julian     ! julian day of year
  real(kind=kind_phys), dimension(im)     , intent(in)    :: rainn_mp   ! microphysics non-convective precipitation [mm]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: rainc_mp   ! microphysics convective precipitation [mm]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: snow_mp    ! microphysics snow [mm]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: graupel_mp ! microphysics graupel [mm]
  real(kind=kind_phys), dimension(im)     , intent(in)    :: ice_mp     ! microphysics ice/hail [mm]
  real(kind=kind_phys)                    , intent(in)    :: con_hvap   ! latent heat condensation [J/kg]
  real(kind=kind_phys)                    , intent(in)    :: con_cp     ! specific heat air [J/kg/K] 
  real(kind=kind_phys)                    , intent(in)    :: con_jcal   ! joules per calorie (not used)
  real(kind=kind_phys)                    , intent(in)    :: rhoh2o     ! density of water [kg/m^3]
  real(kind=kind_phys)                    , intent(in)    :: con_eps    ! Rd/Rv 
  real(kind=kind_phys)                    , intent(in)    :: con_epsm1  ! Rd/Rv - 1
  real(kind=kind_phys)                    , intent(in)    :: con_fvirt  ! Rv/Rd - 1
  real(kind=kind_phys)                    , intent(in)    :: con_rd     ! gas constant air [J/kg/K]
  real(kind=kind_phys)                    , intent(in)    :: con_hfus   ! lat heat H2O fusion  [J/kg]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: weasd      ! water equivalent accumulated snow depth [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: snwdph     ! snow depth [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tskin      ! ground surface skin temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tprcp      ! total precipitation [m]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: srflag     ! snow/rain flag for precipitation
  real(kind=kind_phys), dimension(im,km)  , intent(inout) :: smc        ! total soil moisture content [m3/m3]
  real(kind=kind_phys), dimension(im,km)  , intent(inout) :: stc        ! soil temp [K]
  real(kind=kind_phys), dimension(im,km)  , intent(inout) :: slc        ! liquid soil moisture [m3/m3]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: canopy     ! canopy moisture content [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: trans      ! total plant transpiration [m/s]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tsurf      ! surface skin temperature [after iteration]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: zorl       ! surface roughness [cm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: snowxy     ! actual no. of snow layers
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tvxy       ! vegetation leaf temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tgxy       ! bulk ground surface temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: canicexy   ! canopy-intercepted ice [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: canliqxy   ! canopy-intercepted liquid water [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: eahxy      ! canopy air vapor pressure [Pa]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: tahxy      ! canopy air temperature [K]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: cmxy       ! bulk momentum drag coefficient [m/s]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: chxy       ! bulk sensible heat exchange coefficient [m/s]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: fwetxy     ! wetted or snowed fraction of the canopy [-]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: sneqvoxy   ! snow mass at last time step[mm h2o]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: alboldxy   ! snow albedo at last time step [-]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: qsnowxy    ! snowfall on the ground [mm/s]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: wslakexy   ! lake water storage [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: zwtxy      ! water table depth [m]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: waxy       ! water in the "aquifer" [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: wtxy       ! groundwater storage [mm]
  real(kind=kind_phys), dimension(im,-2:0), intent(inout) :: tsnoxy     ! snow temperature [K]
  real(kind=kind_phys), dimension(im,-2:4), intent(inout) :: zsnsoxy    ! snow/soil layer depth [m]
  real(kind=kind_phys), dimension(im,-2:0), intent(inout) :: snicexy    ! snow layer ice [mm]
  real(kind=kind_phys), dimension(im,-2:0), intent(inout) :: snliqxy    ! snow layer liquid water [mm]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: lfmassxy   ! leaf mass [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: rtmassxy   ! mass of fine roots [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: stmassxy   ! stem mass [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: woodxy     ! mass of wood (incl. woody roots) [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: stblcpxy   ! stable carbon in deep soil [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: fastcpxy   ! short-lived carbon, shallow soil [g/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: xlaixy     ! leaf area index [m2/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: xsaixy     ! stem area index [m2/m2]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: taussxy    ! snow age factor [-]
  real(kind=kind_phys), dimension(im,1:4) , intent(inout) :: smoiseq    ! eq volumetric soil moisture [m3/m3]
  real(kind=kind_phys), dimension(im)     , intent(inout) :: smcwtdxy   ! soil moisture content in the layer to the water table when deep
  real(kind=kind_phys), dimension(im)     , intent(inout) :: deeprechxy ! recharge to the water table when deep
  real(kind=kind_phys), dimension(im)     , intent(inout) :: rechxy     ! recharge to the water table
  real(kind=kind_phys), dimension(im)     , intent(out)   :: albdvis    ! albedo - direct  visible [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: albdnir    ! albedo - direct  NIR     [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: albivis    ! albedo - diffuse visible [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: albinir    ! albedo - diffuse NIR     [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: emiss      ! sfc lw emissivity [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: sncovr1    ! snow cover over land [fraction]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: qsurf      ! specific humidity at sfc [kg/kg]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: gflux      ! soil heat flux [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: drain      ! subsurface runoff [mm/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: evap       ! total latent heat flux [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: hflx       ! sensible heat flux [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: ep         ! potential evaporation [mm/s?]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: runoff     ! surface runoff [mm/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: cmm        ! cm*U [m/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: chh        ! ch*U*rho [kg/m2/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: evbs       ! direct soil evaporation [m/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: evcw       ! canopy water evaporation [m/s]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: sbsno      ! sublimation/deposit from snopack [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: snowc      ! fractional snow cover [-]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: stm        ! total soil column moisture content [mm]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: snohf      ! snow/freezing-rain latent heat flux [W/m2]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: smcwlt2    ! dry soil moisture threshold [m3/m3]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: smcref2    ! soil moisture threshold [m3/m3]
  real(kind=kind_phys), dimension(im)     , intent(out)   :: wet1       ! normalized surface soil saturated fraction
  real(kind=kind_phys), dimension(im)     , intent(out)   :: t2mmp      ! combined T2m from tiles
  real(kind=kind_phys), dimension(im)     , intent(out)   :: q2mp       ! combined q2m from tiles
  character(len=*)    ,                     intent(out)   :: errmsg
  integer             ,                     intent(out)   :: errflg

!
!  ---  some new options, hard code for now
!

  integer    :: iopt_rsf  = 4 ! option for surface resistance
  integer    :: iopt_soil = 1 ! option for soil parameter treatment
  integer    :: iopt_pedo = 1 ! option for pedotransfer function
  integer    :: iopt_crop = 0 ! option for crop model
  integer    :: iopt_gla  = 2 ! option for glacier treatment

!
!  ---  guess iteration fields - target for removal
!

  real(kind=kind_phys), dimension(im)       :: weasd_old
  real(kind=kind_phys), dimension(im)       :: snwdph_old
  real(kind=kind_phys), dimension(im)       :: tskin_old
  real(kind=kind_phys), dimension(im)       :: canopy_old
  real(kind=kind_phys), dimension(im)       :: tprcp_old
  real(kind=kind_phys), dimension(im)       :: srflag_old
  real(kind=kind_phys), dimension(im)       :: snow_old
  real(kind=kind_phys), dimension(im)       :: tv_old
  real(kind=kind_phys), dimension(im)       :: tg_old
  real(kind=kind_phys), dimension(im)       :: canice_old
  real(kind=kind_phys), dimension(im)       :: canliq_old
  real(kind=kind_phys), dimension(im)       :: eah_old
  real(kind=kind_phys), dimension(im)       :: tah_old
  real(kind=kind_phys), dimension(im)       :: fwet_old
  real(kind=kind_phys), dimension(im)       :: sneqvo_old
  real(kind=kind_phys), dimension(im)       :: albold_old
  real(kind=kind_phys), dimension(im)       :: qsnow_old
  real(kind=kind_phys), dimension(im)       :: wslake_old
  real(kind=kind_phys), dimension(im)       :: zwt_old
  real(kind=kind_phys), dimension(im)       :: wa_old
  real(kind=kind_phys), dimension(im)       :: wt_old
  real(kind=kind_phys), dimension(im)       :: lfmass_old
  real(kind=kind_phys), dimension(im)       :: rtmass_old
  real(kind=kind_phys), dimension(im)       :: stmass_old
  real(kind=kind_phys), dimension(im)       :: wood_old
  real(kind=kind_phys), dimension(im)       :: stblcp_old
  real(kind=kind_phys), dimension(im)       :: fastcp_old
  real(kind=kind_phys), dimension(im)       :: xlai_old
  real(kind=kind_phys), dimension(im)       :: xsai_old
  real(kind=kind_phys), dimension(im)       :: tauss_old
  real(kind=kind_phys), dimension(im)       :: smcwtd_old
  real(kind=kind_phys), dimension(im)       :: rech_old
  real(kind=kind_phys), dimension(im)       :: deeprech_old
  real(kind=kind_phys), dimension(im,   km) :: smc_old
  real(kind=kind_phys), dimension(im,   km) :: stc_old
  real(kind=kind_phys), dimension(im,   km) :: slc_old
  real(kind=kind_phys), dimension(im,   km) :: smoiseq_old
  real(kind=kind_phys), dimension(im,-2: 0) :: tsno_old  
  real(kind=kind_phys), dimension(im,-2: 0) :: snice_old
  real(kind=kind_phys), dimension(im,-2: 0) :: snliq_old 
  real(kind=kind_phys), dimension(im,-2:km) :: zsnso_old
  real(kind=kind_phys), dimension(im,-2:km) :: tsnso_old

!
!  ---  local inputs to noah-mp and glacier subroutines; listed in order in noah-mp call
!
                                                                            ! intent
  integer                                          :: i_location            ! in    | grid index
  integer                                          :: j_location            ! in    | grid index (not used in ccpp)
  real (kind=kind_phys)                            :: latitude              ! in    | latitude [radians]
  integer                                          :: year_length           ! in    | number of days in the current year
  real (kind=kind_phys)                            :: julian_day            ! in    | julian day of year [floating point]
  real (kind=kind_phys)                            :: cosine_zenith         ! in    | cosine solar zenith angle [-1,1]
  real (kind=kind_phys)                            :: timestep              ! in    | time step [sec]
  real (kind=kind_phys)                            :: spatial_scale         ! in    | spatial scale [m] (not used in noah-mp)
  real (kind=kind_phys)                            :: atmosphere_thickness  ! in    | thickness of lowest atmo layer [m] (not used in noah-mp)
  integer                                          :: soil_levels           ! in    | soil levels
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_interface_depth  ! in    | soil layer-bottom depth from surface [m]
  integer                                          :: max_snow_levels       ! in    | maximum number of snow levels
  real (kind=kind_phys)                            :: vegetation_frac       ! in    | vegetation fraction [0.0-1.0]
  real (kind=kind_phys)                            :: max_vegetation_frac   ! in    | annual maximum vegetation fraction [0.0-1.0]
  integer                                          :: vegetation_category   ! in    | vegetation category
  integer                                          :: ice_flag              ! in    | ice flag (1->ice)
  integer                                          :: surface_type          ! in    | surface type flag 1->soil; 2->lake
  integer                                          :: crop_type             ! in    | crop type category
  real (kind=kind_phys), dimension(       1:nsoil) :: eq_soil_water_vol     ! in    | (opt_run=5) equilibrium soil water content [m3/m3]
  real (kind=kind_phys)                            :: temperature_forcing   ! in    | forcing air temperature [K]
  real (kind=kind_phys)                            :: air_pressure_surface  ! in    | surface air pressure [Pa]
  real (kind=kind_phys)                            :: air_pressure_forcing  ! in    | forcing air pressure [Pa]
  real (kind=kind_phys)                            :: uwind_forcing         ! in    | forcing u-wind [m/s]
  real (kind=kind_phys)                            :: vwind_forcing         ! in    | forcing v-wind [m/s]
  real (kind=kind_phys)                            :: spec_humidity_forcing ! in    | forcing mixing ratio [kg/kg]
  real (kind=kind_phys)                            :: cloud_water_forcing   ! in    | cloud water mixing ratio [kg/kg] (not used in noah-mp)
  real (kind=kind_phys)                            :: sw_radiation_forcing  ! in    | forcing downward shortwave radiation [W/m2]
  real (kind=kind_phys)                            :: radiation_lw_forcing  ! in    | forcing downward longwave radiation [W/m2]
  real (kind=kind_phys)                            :: precipitation_forcing ! in    | total precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_convective     ! in    | convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_non_convective ! in    | non-convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_sh_convective  ! in    | shallow convective precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_snow           ! in    | snow precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_graupel        ! in    | graupel precipitation [mm/s]
  real (kind=kind_phys)                            :: precip_hail           ! in    | hail precipitation [mm/s]
  real (kind=kind_phys)                            :: temperature_soil_bot  ! in    | soil bottom boundary condition temperature [K]
  real (kind=kind_phys)                            :: co2_air               ! in    | atmospheric co2 concentration [Pa]
  real (kind=kind_phys)                            :: o2_air                ! in    | atmospheric o2 concentration [Pa]
  real (kind=kind_phys)                            :: foliage_nitrogen      ! in    | foliage nitrogen [%] [1-saturated]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_ice_frac_old     ! in    | snow ice fraction at last timestep [-]
  real (kind=kind_phys)                            :: forcing_height        ! inout | forcing height [m]
  real (kind=kind_phys)                            :: snow_albedo_old       ! inout | snow albedo at last time step (class option) [-]
  real (kind=kind_phys)                            :: snow_water_equiv_old  ! inout | snow water equivalent at last time step [mm]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil) :: temperature_snow_soil ! inout | snow/soil temperature [K]
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_liquid_vol       ! inout | volumetric liquid soil moisture [m3/m3]
  real (kind=kind_phys), dimension(       1:nsoil) :: soil_moisture_vol     ! inout | volumetric soil moisture (ice + liq.) [m3/m3]
  real (kind=kind_phys)                            :: temperature_canopy_air! inout | canopy air tmeperature [K]
  real (kind=kind_phys)                            :: vapor_pres_canopy_air ! inout | canopy air vapor pressure [Pa]
  real (kind=kind_phys)                            :: canopy_wet_fraction   ! inout | wetted or snowed fraction of canopy (-)
  real (kind=kind_phys)                            :: canopy_liquid         ! inout | canopy intercepted liquid [mm]
  real (kind=kind_phys)                            :: canopy_ice            ! inout | canopy intercepted ice [mm]
  real (kind=kind_phys)                            :: temperature_leaf      ! inout | leaf temperature [K]
  real (kind=kind_phys)                            :: temperature_ground    ! inout | grid ground surface temperature [K]
  real (kind=kind_phys)                            :: spec_humidity_surface ! inout | surface specific humidty [kg/kg]
  real (kind=kind_phys)                            :: snowfall              ! inout | land model partitioned snowfall [mm/s]
  real (kind=kind_phys)                            :: rainfall              ! inout | land model partitioned rainfall [mm/s]
  integer                                          :: snow_levels           ! inout | active snow levels [-]
  real (kind=kind_phys), dimension(-nsnow+1:nsoil) :: interface_depth       ! inout | layer-bottom depth from snow surf [m]
  real (kind=kind_phys)                            :: snow_depth            ! inout | snow depth [m]
  real (kind=kind_phys)                            :: snow_water_equiv      ! inout | snow water equivalent [mm]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_level_ice        ! inout | snow level ice [mm]
  real (kind=kind_phys), dimension(-nsnow+1:    0) :: snow_level_liquid     ! inout | snow level liquid [mm]
  real (kind=kind_phys)                            :: depth_water_table     ! inout | depth to water table [m]
  real (kind=kind_phys)                            :: aquifer_water         ! inout | water storage in aquifer [mm]
  real (kind=kind_phys)                            :: saturated_water       ! inout | water in aquifer+saturated soil [mm]
  real (kind=kind_phys)                            :: lake_water            ! inout | lake water storage (can be neg.) [mm]
  real (kind=kind_phys)                            :: leaf_carbon           ! inout | leaf mass [g/m2]
  real (kind=kind_phys)                            :: root_carbon           ! inout | mass of fine roots [g/m2]
  real (kind=kind_phys)                            :: stem_carbon           ! inout | stem mass [g/m2]
  real (kind=kind_phys)                            :: wood_carbon           ! inout | mass of wood (incl. woody roots) [g/m2]
  real (kind=kind_phys)                            :: soil_carbon_stable    ! inout | stable soil carbon [g/m2]
  real (kind=kind_phys)                            :: soil_carbon_fast      ! inout | short-lived soil carbon [g/m2]
  real (kind=kind_phys)                            :: leaf_area_index       ! inout | leaf area index [-]
  real (kind=kind_phys)                            :: stem_area_index       ! inout | stem area index [-]
  real (kind=kind_phys)                            :: cm_noahmp             ! inout | grid momentum drag coefficient [m/s]
  real (kind=kind_phys)                            :: ch_noahmp             ! inout | grid heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: snow_age              ! inout | non-dimensional snow age [-]
  real (kind=kind_phys)                            :: grain_carbon          ! inout | grain mass [g/m2]
  real (kind=kind_phys)                            :: growing_deg_days      ! inout | growing degree days [-]
  integer                                          :: plant_growth_stage    ! inout | plant growing stage [-]
  real (kind=kind_phys)                            :: soil_moisture_wtd     ! inout | (opt_run=5) soil water content between bottom of the soil and water table [m3/m3]
  real (kind=kind_phys)                            :: deep_recharge         ! inout | (opt_run=5) recharge to or from the water table when deep [m]
  real (kind=kind_phys)                            :: recharge              ! inout | (opt_run=5) recharge to or from the water table when shallow [m] (diagnostic)
  real (kind=kind_phys)                            :: z0_total              !   out | weighted z0 sent to coupled model [m]
  real (kind=kind_phys)                            :: sw_absorbed_total     !   out | total absorbed solar radiation [W/m2]
  real (kind=kind_phys)                            :: sw_reflected_total    !   out | total reflected solar radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_total     !   out | total net lw rad [W/m2]  [+ to atm]
  real (kind=kind_phys)                            :: sensible_heat_total   !   out | total sensible heat [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: ground_heat_total     !   out | ground heat flux [W/m2]   [+ to soil]
  real (kind=kind_phys)                            :: latent_heat_canopy    !   out | canopy evaporation heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: latent_heat_ground    !   out | ground evaporation heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: transpiration_heat    !   out | transpiration heat flux [W/m2] [+ to atm]
  real (kind=kind_phys)                            :: evaporation_canopy    !   out | canopy evaporation [mm/s]
  real (kind=kind_phys)                            :: transpiration         !   out | transpiration [mm/s]
  real (kind=kind_phys)                            :: evaporation_soil      !   out | soil surface evaporation [mm/s]
  real (kind=kind_phys)                            :: temperature_radiative !   out | surface radiative temperature [K]
  real (kind=kind_phys)                            :: temperature_bare_grd  !   out | bare ground surface temperature [K]
  real (kind=kind_phys)                            :: temperature_veg_grd   !   out | below_canopy ground surface temperature [K]
  real (kind=kind_phys)                            :: temperature_veg_2m    !   out | vegetated 2-m air temperature [K]
  real (kind=kind_phys)                            :: temperature_bare_2m   !   out | bare ground 2-m air temperature [K]
  real (kind=kind_phys)                            :: spec_humidity_veg_2m  !   out | vegetated 2-m air specific humidity [K]
  real (kind=kind_phys)                            :: spec_humidity_bare_2m !   out | bare ground 2-m air specfic humidity [K]
  real (kind=kind_phys)                            :: runoff_surface        !   out | surface runoff [mm/s] 
  real (kind=kind_phys)                            :: runoff_baseflow       !   out | baseflow runoff [mm/s]
  real (kind=kind_phys)                            :: par_absorbed          !   out | absorbed photosynthesis active radiation [W/m2]
  real (kind=kind_phys)                            :: photosynthesis        !   out | total photosynthesis [umol CO2/m2/s] [+ out]
  real (kind=kind_phys)                            :: sw_absorbed_veg       !   out | solar radiation absorbed by vegetation [W/m2]
  real (kind=kind_phys)                            :: sw_absorbed_ground    !   out | solar radiation absorbed by ground [W/m2]
  real (kind=kind_phys)                            :: snow_cover_fraction   !   out | snow cover fraction on the ground [-]
  real (kind=kind_phys)                            :: net_eco_exchange      !   out | net ecosystem exchange [g/m2/s CO2]
  real (kind=kind_phys)                            :: global_prim_prod      !   out | global primary production [g/m2/s C]
  real (kind=kind_phys)                            :: net_prim_prod         !   out | net primary productivity [g/m2/s C]
  real (kind=kind_phys)                            :: vegetation_fraction   !   out | vegetation fraction [0.0-1.0]
  real (kind=kind_phys)                            :: albedo_total          !   out | total surface albedo [-]
  real (kind=kind_phys)                            :: snowmelt_out          !   out | snowmelt out bottom of pack [mm/s]
  real (kind=kind_phys)                            :: snowmelt_shallow      !   out | shallow snow melt [mm]
  real (kind=kind_phys)                            :: snowmelt_shallow_1    !   out | additional shallow snow melt [mm]
  real (kind=kind_phys)                            :: snowmelt_shallow_2    !   out | additional shallow snow melt [mm]
  real (kind=kind_phys)                            :: rs_sunlit             !   out | sunlit leaf stomatal resistance [s/m]
  real (kind=kind_phys)                            :: rs_shaded             !   out | shaded leaf stomatal resistance [s/m]
  real (kind=kind_phys), dimension(1:2)            :: albedo_direct         !   out | direct vis/nir albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_diffuse        !   out | diffuse vis/nir albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_direct_snow    !   out | direct vis/nir snow albedo [-]
  real (kind=kind_phys), dimension(1:2)            :: albedo_diffuse_snow   !   out | diffuse vis/nir snow albedo [-]
  real (kind=kind_phys)                            :: canopy_gap_fraction   !   out | between canopy gap fraction [-]
  real (kind=kind_phys)                            :: incanopy_gap_fraction !   out | within canopy gap fraction for beam [-]
  real (kind=kind_phys)                            :: ch_vegetated          !   out | vegetated heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_bare_ground        !   out | bare-ground heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: emissivity_total      !   out | grid emissivity [-]
  real (kind=kind_phys)                            :: sensible_heat_grd_veg !   out | below-canopy ground sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: sensible_heat_leaf    !   out | leaf-to-canopy sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: sensible_heat_grd_bar !   out | bare ground sensible heat flux [W/m2]
  real (kind=kind_phys)                            :: latent_heat_grd_veg   !   out | below-canopy ground evaporation heat flux [W/m2]
  real (kind=kind_phys)                            :: latent_heat_grd_bare  !   out | bare ground evaporation heat flux [W/m2]
  real (kind=kind_phys)                            :: ground_heat_veg       !   out | below-canopy ground heat flux [W/m2]
  real (kind=kind_phys)                            :: ground_heat_bare      !   out | bare ground heat flux [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_grd_veg   !   out | below-canopy ground absorbed longwave radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_leaf      !   out | leaf absorbed longwave radiation [W/m2]
  real (kind=kind_phys)                            :: lw_absorbed_grd_bare  !   out | bare ground net longwave radiation [W/m2]
  real (kind=kind_phys)                            :: latent_heat_trans     !   out | transpiration [W/m2]
  real (kind=kind_phys)                            :: latent_heat_leaf      !   out | leaf evaporation [W/m2]
  real (kind=kind_phys)                            :: ch_leaf               !   out | leaf exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_below_canopy       !   out | below-canopy exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_vegetated_2m       !   out | 2-m vegetated  heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: ch_bare_ground_2m     !   out | 2-m bare-ground heat exchange coefficient [m/s]
  real (kind=kind_phys)                            :: precip_frozen_frac    !   out | precipitation snow fraction [-]
  real (kind=kind_phys)                            :: precip_adv_heat_veg   !   out | precipitation advected heat - vegetation net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_grd_v !   out | precipitation advected heat - below-canopy net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_grd_b !   out | precipitation advected heat - bare ground net [W/m2]
  real (kind=kind_phys)                            :: precip_adv_heat_total !   out | precipitation advected heat - total [W/m2)
  real (kind=kind_phys)                            :: snow_sublimation      !   out | snow sublimation [W/m2]
  real (kind=kind_phys)                            :: lai_sunlit            !   out | sunlit leaf area index [m2/m2]
  real (kind=kind_phys)                            :: lai_shaded            !   out | shaded leaf area index [m2/m2]
  real (kind=kind_phys)                            :: leaf_air_resistance   !   out | leaf boundary layer resistance [s/m]

!
!  ---  local variable
!

  integer :: soil_category(nsoil)
  integer :: slope_category
  integer :: soil_color_category

  real (kind=kind_phys) :: spec_humidity_sat      ! saturation specific humidity
  real (kind=kind_phys) :: vapor_pressure_sat     ! saturation vapor pressure
  real (kind=kind_phys) :: latent_heat_total        ! total latent heat flux [W/m2]
  real (kind=kind_phys) :: density                ! air density
  real (kind=kind_phys) :: virtual_temperature    ! used for penman calculation and density
  real (kind=kind_phys) :: potential_evaporation  ! used for penman calculation
  real (kind=kind_phys) :: potential_temperature  ! used for penman calculation
  real (kind=kind_phys) :: penman_radiation       ! used for penman calculation
  real (kind=kind_phys) :: dqsdt                  ! used for penman calculation
  real (kind=kind_phys) :: precip_freeze_frac_in  ! used for penman calculation
 
  logical               :: is_snowing             ! used for penman calculation
  logical               :: is_freeze_rain         ! used for penman calculation
  integer :: i, k
      
!
!  --- local derived constants:
!
      
  type(noahmp_parameters) :: parameters

!
!  --- end declaration
!     

!
!  --- Initialize CCPP error handling variables
!
  errmsg = ''
  errflg = 0

!
!  --- save land-related prognostic fields for guess run  TARGET FOR REMOVAL
!
  do i = 1, im
    if (dry(i) .and. flag_guess(i)) then
      weasd_old(i)   = weasd(i)
      snwdph_old(i)  = snwdph(i)
      tskin_old(i)   = tskin(i)
      canopy_old(i)  = canopy(i)
      tprcp_old(i)   = tprcp(i)
      srflag_old(i)  = srflag(i)
      snow_old(i)    = snowxy(i)
      tv_old(i)      = tvxy(i)
      tg_old(i)      = tgxy(i)
      canice_old(i)  = canicexy(i)
      canliq_old(i)  = canliqxy(i)
      eah_old(i)     = eahxy(i)
      tah_old(i)     = tahxy(i)
      fwet_old(i)    = fwetxy(i)
      sneqvo_old(i)  = sneqvoxy(i) 
      albold_old(i)  = alboldxy(i)
      qsnow_old(i)   = qsnowxy(i)
      wslake_old(i)  = wslakexy(i)
      zwt_old(i)     = zwtxy(i)
      wa_old(i)      = waxy(i)
      wt_old(i)      = wtxy(i)
      lfmass_old(i)  = lfmassxy(i)
      rtmass_old(i)  = rtmassxy(i)
      stmass_old(i)  = stmassxy(i)
      wood_old(i)    = woodxy(i)
      stblcp_old(i)  = stblcpxy(i)
      fastcp_old(i)  = fastcpxy(i)
      xlai_old(i)    = xlaixy(i)
      xsai_old(i)    = xsaixy(i)
      tauss_old(i)   = taussxy(i)
      smcwtd_old(i)  = smcwtdxy(i)
      rech_old(i)    = rechxy(i)
      deeprech_old(i) = deeprechxy(i)

      do k = 1, km
        smc_old(i,k)     = smc(i,k)
        stc_old(i,k)     = stc(i,k)
        slc_old(i,k)     = slc(i,k)
        smoiseq_old(i,k) = smoiseq(i,k)
      end do

      do k = -2, 0
        tsno_old(i,k)  = tsnoxy(i,k)
        snice_old(i,k) = snicexy(i,k)
        snliq_old(i,k) = snliqxy(i,k)
      end do

      do k = -2, km
        zsnso_old (i,k) = zsnsoxy(i,k)
      end do

    end if  ! dry(i) .and. flag_guess(i)

  end do  ! im _old loop

  do i = 1, im

    if (flag_iter(i) .and. dry(i)) then

!
!  --- variable checks and derived fields
!

      if(vegtype(i) == isice_table ) then
        if(weasd(i) < 0.1) then
          weasd(i)  = 0.1
        end if
      end if                                       

!
!  --- noah-mp input variables (except snow_ice_frac_old done later)
!

      i_location            = i
      j_location            = -9999
      latitude              = xlatin(i)
      year_length           = iyrlen
      julian_day            = julian
      cosine_zenith         = xcoszin(i)
      timestep              = delt
      spatial_scale         = -9999.0
      atmosphere_thickness  = -9999.0
      soil_levels           = km
      soil_interface_depth  = zsoil
      max_snow_levels       = nsnow
      vegetation_frac       = sigmaf(i)
      max_vegetation_frac   = shdmax(i)
      vegetation_category   = vegtype(i)
      surface_type          = 1
      crop_type             = 0
      eq_soil_water_vol     = smoiseq(i,:) ! only need for run=5
      temperature_forcing   = t1(i) 
      air_pressure_surface  = ps(i)
      air_pressure_forcing  = prsl1(i)
      uwind_forcing         = u1(i)
      vwind_forcing         = v1(i)

      spec_humidity_forcing = max(q1(i), 1.e-8)                            ! specific humidity at level 1 (kg/kg)
      virtual_temperature   = temperature_forcing * &
                               (1.0 + con_fvirt * spec_humidity_forcing)   ! virtual temperature
      vapor_pressure_sat    = fpvs( temperature_forcing )                  ! sat. vapor pressure at level 1 (Pa)
      spec_humidity_sat     = con_eps*vapor_pressure_sat / &
                              (prsl1(i) + con_epsm1*vapor_pressure_sat)    ! sat. specific humidity at level 1 (kg/kg)
      spec_humidity_sat     = max(spec_humidity_sat, 1.e-8)                ! lower limit sat. specific humidity (kg/kg)
      spec_humidity_forcing = min(spec_humidity_sat,spec_humidity_forcing) ! limit specific humidity at level 1 (kg/kg)

      cloud_water_forcing   = -9999.0
      sw_radiation_forcing  = dswsfc(i)
      radiation_lw_forcing  = dlwflx(i)
      precipitation_forcing = rhoh2o * tprcp(i) / delt !1000.0 * tprcp(i) / delt
      precip_convective     = rainc_mp(i)
      precip_non_convective = rainn_mp(i)
      precip_sh_convective  = 0.
      precip_snow           = snow_mp(i)
      precip_graupel        = graupel_mp(i)
      precip_hail           = ice_mp(i)
      temperature_soil_bot  = tg3(i)
      co2_air               = co2_table * air_pressure_forcing
      o2_air                = o2_table  * air_pressure_forcing
      foliage_nitrogen      = 1.0
      
!
!  --- noah-mp inout variables
!

      forcing_height               = zf(i)
      snow_albedo_old              = alboldxy(i)
      snow_water_equiv_old         = sneqvoxy(i)
      temperature_snow_soil(-2: 0) = tsnoxy(i,:)
      temperature_snow_soil( 1:km) = stc(i,:)
      soil_liquid_vol              = slc(i,:)
      soil_moisture_vol            = smc(i,:)
      temperature_canopy_air       = tahxy(i)
      vapor_pres_canopy_air        = (air_pressure_surface*spec_humidity_forcing)/(0.622+spec_humidity_forcing) 
      ! TODO recalculated? set to eahxy
      canopy_wet_fraction          = fwetxy(i)
      canopy_liquid                = canliqxy(i)
      canopy_ice                   = canicexy(i)
      temperature_leaf             = tvxy(i)
      temperature_ground           = tgxy(i)
      spec_humidity_surface        = undefined  ! doesn't need inout; should be out
      snowfall                     = qsnowxy(i) ! doesn't need inout; should be out
!      rainfall              ! doesn't need inout; should be out TODO
      snow_levels                  = nint(snowxy(i))
      interface_depth              = zsnsoxy(i,:)
      snow_depth                   = snwdph(i) * 0.001         ! convert from mm to m
      snow_water_equiv             = weasd(i)
          if (snow_water_equiv /= 0.0 .and. snow_depth == 0.0) then  !TODO this should be done elsewhere
            snow_depth = 10.0 * snow_water_equiv /1000.0
          endif
      snow_level_ice               = snicexy(i,:)
      snow_level_liquid            = snliqxy(i,:)
      depth_water_table            = zwtxy(i)
      aquifer_water                = waxy(i)
      saturated_water              = waxy(i)     ! why not wta !!! TODO
      lake_water                   = wslakexy(i)
      leaf_carbon                  = lfmassxy(i)
      root_carbon                  = rtmassxy(i)
      stem_carbon                  = stmassxy(i)
      wood_carbon                  = woodxy(i)
      soil_carbon_stable           = stblcpxy(i)
      soil_carbon_fast             = fastcpxy(i)
      leaf_area_index              = xlaixy(i)
      stem_area_index              = xsaixy(i)
      cm_noahmp                    = cmxy(i)
      ch_noahmp                    = chxy(i)
      snow_age                     = taussxy(i)
!      grain_carbon          ! new variable
!      growing_deg_days      ! new variable
!      plant_growth_stage    ! new variable
      soil_moisture_wtd            = smcwtdxy(i)
      deep_recharge                = deeprechxy(i)
      recharge                     = rechxy(i)
      
      snow_ice_frac_old = 0.0
      do k = snow_levels+1, 0
        if(snow_level_ice(k) > 0.0 ) &
          snow_ice_frac_old(k) = snow_level_ice(k) /(snow_level_ice(k)+snow_level_liquid(k)) 
      end do

!
!  --- some outputs for atm model?
!
      density = air_pressure_forcing / (con_rd * virtual_temperature)
      chh(i) = ch(i)  * wind(i) * density
      cmm(i) = cm(i)  * wind(i)
!
!  --- noah-mp additional variables
!

      soil_category       = soiltyp(i)
      slope_category      = slopetyp(i)
      soil_color_category = 4

      call transfer_mp_parameters(vegetation_category,soil_category, &
                        slope_category,soil_color_category,crop_type,parameters)

      call noahmp_options(idveg ,iopt_crs, iopt_btr , iopt_run, iopt_sfc,  &
                                 iopt_frz, iopt_inf , iopt_rad, iopt_alb,  &
                                 iopt_snf, iopt_tbot, iopt_stc,            &
			         iopt_rsf, iopt_soil, iopt_pedo, iopt_crop )

      if ( vegetation_category == isice_table )  then

        if (precipitation_forcing > 0.0) then
          if (srflag(i) > 0.0) then                              ! TODO rain/snow flag, one condition is enough?
            snowfall = srflag(i) * precipitation_forcing/10.0    ! still use rho water?
          endif
        endif

        ice_flag = -1
        temperature_soil_bot = min(temperature_soil_bot,263.15)

        call noahmp_options_glacier(iopt_alb, iopt_snf, iopt_tbot, iopt_stc, iopt_gla )

        call noahmp_glacier (                                                                      &
          i_location           ,1                    ,cosine_zenith        ,nsnow                , &
          nsoil                ,timestep             ,                                             &
          temperature_forcing  ,air_pressure_forcing ,uwind_forcing        ,vwind_forcing        , &
          spec_humidity_forcing,sw_radiation_forcing ,precipitation_forcing,radiation_lw_forcing , &
          temperature_soil_bot ,forcing_height       ,snow_ice_frac_old    ,zsoil                , &
          snowfall             ,snow_water_equiv_old ,snow_albedo_old      ,                       &
          cm_noahmp            ,ch_noahmp            ,snow_levels          ,snow_water_equiv     , &
          soil_moisture_vol    ,interface_depth      ,snow_depth           ,snow_level_ice       , &
          snow_level_liquid    ,temperature_ground   ,temperature_snow_soil,soil_liquid_vol      , &
          snow_age             ,spec_humidity_surface,sw_absorbed_total    ,sw_reflected_total   , &
          lw_absorbed_total    ,sensible_heat_total  ,latent_heat_ground   ,ground_heat_total    , &
          temperature_radiative,evaporation_soil     ,runoff_surface       ,runoff_baseflow      , &
          sw_absorbed_ground   ,albedo_total         ,snowmelt_out         ,snowmelt_shallow     , &
          snowmelt_shallow_1   ,snowmelt_shallow_2   ,temperature_bare_2m  ,spec_humidity_bare_2m, &
          emissivity_total     ,precip_frozen_frac   ,ch_bare_ground_2m    ,snow_sublimation     , &
#ifdef CCPP
          albedo_direct        ,albedo_diffuse       ,errmsg               ,errflg )
#else
          albedo_direct        ,albedo_diffuse  )
#endif

#ifdef CCPP
        if (errflg /= 0) return
#endif

!
! set some non-glacier fields over the glacier
!

        snow_cover_fraction    = 1.0
        temperature_leaf       = undefined  
        canopy_ice             = undefined
        canopy_liquid          = undefined
        vapor_pres_canopy_air  = undefined
        temperature_canopy_air = undefined
        canopy_wet_fraction    = undefined
        lake_water             = undefined
        depth_water_table      = undefined
        aquifer_water          = undefined
        saturated_water        = undefined
        leaf_carbon            = undefined
        root_carbon            = undefined
        stem_carbon            = undefined
        wood_carbon            = undefined
        soil_carbon_stable     = undefined
        soil_carbon_fast       = undefined
        leaf_area_index        = undefined
        stem_area_index        = undefined
        soil_moisture_wtd      = 0.0
        recharge               = 0.0
        deep_recharge          = 0.0
        eq_soil_water_vol      = soil_moisture_vol
        transpiration_heat     = undefined
        latent_heat_canopy     = undefined
        z0_total               = 0.002
        latent_heat_total      = latent_heat_ground
        t2mmp(i)               = temperature_bare_2m
        q2mp(i)                = spec_humidity_bare_2m

      else  ! not glacier

        ice_flag = 0 

        call noahmp_sflx (parameters                                          , &
          i_location            ,j_location            ,latitude              , &
          year_length           ,julian_day            ,cosine_zenith         , &
          timestep              ,spatial_scale         ,atmosphere_thickness  , &
          soil_levels           ,soil_interface_depth  ,max_snow_levels       , &
          vegetation_frac       ,max_vegetation_frac   ,vegetation_category   , &
          ice_flag              ,surface_type          ,crop_type             , &
          eq_soil_water_vol     ,temperature_forcing   ,air_pressure_forcing  , &
          air_pressure_surface  ,uwind_forcing         ,vwind_forcing         , &
          spec_humidity_forcing ,cloud_water_forcing   ,sw_radiation_forcing  , &
          radiation_lw_forcing  ,precip_convective     , &
          precip_non_convective ,precip_sh_convective  ,precip_snow           , &
          precip_graupel        ,precip_hail           ,temperature_soil_bot  , &
          co2_air               ,o2_air                ,foliage_nitrogen      , &
          snow_ice_frac_old     ,                                               &
          forcing_height        ,snow_albedo_old       ,snow_water_equiv_old  , &
          temperature_snow_soil ,soil_liquid_vol       ,soil_moisture_vol     , &
          temperature_canopy_air,vapor_pres_canopy_air ,canopy_wet_fraction   , &
          canopy_liquid         ,canopy_ice            ,temperature_leaf      , &
          temperature_ground    ,spec_humidity_surface ,snowfall              , &
          rainfall              ,snow_levels           ,interface_depth       , &
          snow_depth            ,snow_water_equiv      ,snow_level_ice        , &
          snow_level_liquid     ,depth_water_table     ,aquifer_water         , &
          saturated_water       ,                                               &
          lake_water            ,leaf_carbon           ,root_carbon           , &
          stem_carbon           ,wood_carbon           ,soil_carbon_stable    , &
          soil_carbon_fast      ,leaf_area_index       ,stem_area_index       , &
          cm_noahmp             ,ch_noahmp             ,snow_age              , &
          grain_carbon          ,growing_deg_days      ,plant_growth_stage    , &
          soil_moisture_wtd     ,deep_recharge         ,recharge              , &
          z0_total              ,sw_absorbed_total     ,sw_reflected_total    , &
          lw_absorbed_total     ,sensible_heat_total   ,ground_heat_total     , &
          latent_heat_canopy    ,latent_heat_ground    ,transpiration_heat    , &
          evaporation_canopy    ,transpiration         ,evaporation_soil      , &
          temperature_radiative ,temperature_bare_grd  ,temperature_veg_grd   , &
          temperature_veg_2m    ,temperature_bare_2m   ,spec_humidity_veg_2m  , &
          spec_humidity_bare_2m ,runoff_surface        ,runoff_baseflow       , &
          par_absorbed          ,photosynthesis        ,sw_absorbed_veg       , &
          sw_absorbed_ground    ,snow_cover_fraction   ,net_eco_exchange      , &
          global_prim_prod      ,net_prim_prod         ,vegetation_fraction   , &
          albedo_total          ,snowmelt_out          ,snowmelt_shallow      , &
          snowmelt_shallow_1    ,snowmelt_shallow_2    ,rs_sunlit             , &
          rs_shaded             ,albedo_direct         ,albedo_diffuse        , &
          albedo_direct_snow    ,albedo_diffuse_snow   ,canopy_gap_fraction   , &
          incanopy_gap_fraction ,ch_vegetated          ,ch_bare_ground        , &
          emissivity_total      ,sensible_heat_grd_veg ,sensible_heat_leaf    , &
          sensible_heat_grd_bar ,latent_heat_grd_veg   ,latent_heat_grd_bare  , &
          ground_heat_veg       ,ground_heat_bare      ,lw_absorbed_grd_veg   , &
          lw_absorbed_leaf      ,lw_absorbed_grd_bare  ,latent_heat_trans     , &
          latent_heat_leaf      ,ch_leaf               ,ch_below_canopy       , &
          ch_vegetated_2m       ,ch_bare_ground_2m     ,precip_frozen_frac    , &
          precip_adv_heat_veg   ,precip_adv_heat_grd_v ,precip_adv_heat_grd_b , &
          precip_adv_heat_total ,snow_sublimation      ,lai_sunlit            , &
#ifdef CCPP
          lai_shaded            ,leaf_air_resistance   ,                        &
          errmsg                ,errflg                )
#else
          lai_shaded            ,leaf_air_resistance   )
#endif
        
#ifdef CCPP
        if (errflg /= 0) return
#endif

        latent_heat_total  = latent_heat_canopy + latent_heat_ground + transpiration_heat

        t2mmp(i) = temperature_veg_2m   * vegetation_fraction + &
                  temperature_bare_2m   * (1-vegetation_fraction) 
         q2mp(i) = spec_humidity_veg_2m * vegetation_fraction + &
                  spec_humidity_bare_2m * (1-vegetation_fraction)

      endif          ! glacial split ends

!
!  --- noah-mp inout and out variables
!

      tsnoxy    (i,:) = temperature_snow_soil(-2: 0)
      stc       (i,:) = temperature_snow_soil( 1:km)
      hflx      (i)   = sensible_heat_total !note unit change below
      evap      (i)   = latent_heat_total   !note unit change below
      evbs      (i)   = latent_heat_ground
      evcw      (i)   = latent_heat_canopy
      trans     (i)   = transpiration_heat
      gflux     (i)   = -1.0*ground_heat_total                          ! opposite sign to be consistent with noah
      snohf     (i)   = snowmelt_out * con_hfus         ! only snow that exits pack
      sbsno     (i)   = snow_sublimation

      cmxy      (i)   = cm_noahmp
      chxy      (i)   = ch_noahmp
      zorl      (i)   = z0_total * 100.0  ! convert to cm

      smc       (i,:) = soil_moisture_vol
      slc       (i,:) = soil_liquid_vol
      snowxy    (i)   = float(snow_levels)
      weasd     (i)   = snow_water_equiv
      snicexy   (i,:) = snow_level_ice
      snliqxy   (i,:) = snow_level_liquid
      snwdph    (i)   = snow_depth * 1000.0       ! convert from mm to m
      canopy    (i)   = canopy_ice + canopy_liquid  ! TODO check units
      canliqxy  (i)   = canopy_liquid
      canicexy  (i)   = canopy_ice
      zwtxy     (i)   = depth_water_table
      waxy      (i)   = aquifer_water
      wtxy      (i)   = saturated_water
      qsnowxy   (i)   = snowfall
      drain     (i)   = runoff_baseflow
      runoff    (i)   = runoff_surface

      lfmassxy  (i)   = leaf_carbon
      rtmassxy  (i)   = root_carbon
      stmassxy  (i)   = stem_carbon
      woodxy    (i)   = wood_carbon
      stblcpxy  (i)   = soil_carbon_stable
      fastcpxy  (i)   = soil_carbon_fast
      xlaixy    (i)   = leaf_area_index
      xsaixy    (i)   = stem_area_index

      snowc     (i)   = snow_cover_fraction
      sncovr1   (i)   = snow_cover_fraction
      ! TODO check this eqn shoud lthis be q2 not q1 why two con_cp
      qsurf     (i)   = q1(i)  + evap(i) / (con_hvap / con_cp * density * con_cp * ch(i) * wind(i))     
      tskin     (i)   = temperature_radiative
      tsurf     (i)   = temperature_radiative
      tvxy      (i)   = temperature_leaf
      tgxy      (i)   = temperature_ground
      tahxy     (i)   = temperature_canopy_air
      eahxy     (i)   = vapor_pres_canopy_air
      emiss     (i)   = emissivity_total

      if(albedo_total > 0.0) then
        sfalb(i)      = albedo_total
        albdvis(i)    = albedo_direct(1)
        albdnir(i)    = albedo_direct(2)
        albivis(i)    = albedo_diffuse(1)
        albinir(i)    = albedo_diffuse(2)
      end if

      zsnsoxy   (i,:) = interface_depth

      wslakexy  (i)   = lake_water         ! not active
      fwetxy    (i)   = canopy_wet_fraction
      taussxy   (i)   = snow_age
      alboldxy  (i)   = snow_albedo_old
      sneqvoxy  (i)   = snow_water_equiv_old

      smcwtdxy  (i)   = soil_moisture_wtd  ! only need for run=5
      deeprechxy(i)   = deep_recharge      ! only need for run=5
      rechxy    (i)   = recharge           ! only need for run=5
      smoiseq   (i,:) = eq_soil_water_vol  ! only need for run=5; listed as in

      stm       (i)   = (0.1*soil_moisture_vol(1) + &
                         0.3*soil_moisture_vol(2) + &
                         0.6*soil_moisture_vol(3) + &        ! clean up and use depths above
                         1.0*soil_moisture_vol(4))*1000.0    ! unit conversion from m to kg m-2

      wet1   (i) = soil_moisture_vol(1) /  smcmax_table(soil_category(1))
      smcwlt2(i) = smcdry_table(soil_category(1))   !!!change to wilt?
      smcref2(i) = smcref_table(soil_category(1))

!      
!  --- change units for output
!
      hflx(i) = hflx(i) / density / con_cp
      evap(i) = evap(i) / density / con_hvap

!
!  --- calculate potential evaporation using noah code
!
      potential_temperature = temperature_forcing * prslki(i)
      virtual_temperature   = temperature_forcing * (1.0 + 0.61*spec_humidity_forcing)
      penman_radiation      = sw_absorbed_total + radiation_lw_forcing
      dqsdt                 = spec_humidity_sat * a23m4/(temperature_forcing-a4)**2

      precip_freeze_frac_in = srflag(i)
      is_snowing            = .false.          
      is_freeze_rain        = .false.
      if (precipitation_forcing > 0.0) then
        if (precip_freeze_frac_in > 0.0) then                   ! rain/snow flag, one condition is enough?
          is_snowing = .true.
        else
          if (temperature_forcing <= 275.15) is_freeze_rain = .true.  
        end if
      end if
      
      call penman (temperature_forcing, air_pressure_forcing , ch_noahmp            , &
                   virtual_temperature, potential_temperature, precipitation_forcing, &
                   penman_radiation   , ground_heat_total    , spec_humidity_forcing, &
                   spec_humidity_sat  , potential_evaporation, is_snowing           , &
                   is_freeze_rain     , precip_freeze_frac_in, dqsdt                , &
                   emissivity_total   , snow_cover_fraction  )

      ep(i) = potential_evaporation

    end if ! flag_iter(i) .and. dry(i)

  end do ! im loop

!
!  --- restore land-related prognostic fields for guess run  TARGET FOR REMOVAL
!

      do i = 1, im
        if (dry(i) .and. flag_guess(i)) then
          weasd(i)      = weasd_old(i)
          snwdph(i)     = snwdph_old(i)
          tskin(i)      = tskin_old(i)
          canopy(i)     = canopy_old(i)
          tprcp(i)      = tprcp_old(i)
          srflag(i)     = srflag_old(i)
          snowxy(i)     = snow_old(i)
          tvxy(i)       = tv_old(i)
          tgxy(i)       = tg_old(i)
          canicexy(i)   = canice_old(i)
          canliqxy(i)   = canliq_old(i)
          eahxy(i)      = eah_old(i)
          tahxy(i)      = tah_old(i)
          fwetxy(i)     = fwet_old(i)
          sneqvoxy(i)   = sneqvo_old(i)
          alboldxy(i)   = albold_old(i)
          qsnowxy(i)    = qsnow_old(i)
          wslakexy(i)   = wslake_old(i)
          zwtxy(i)      = zwt_old(i)
          waxy(i)       = wa_old(i)
          wtxy(i)       = wt_old(i)
          lfmassxy(i)   = lfmass_old(i)
          rtmassxy(i)   = rtmass_old(i)
          stmassxy(i)   = stmass_old(i)
          woodxy(i)     = wood_old(i)
          stblcpxy(i)   = stblcp_old(i)
          fastcpxy(i)   = fastcp_old(i)
          xlaixy(i)     = xlai_old(i)
          xsaixy(i)     = xsai_old(i)
          taussxy(i)    = tauss_old(i)
          smcwtdxy(i)   = smcwtd_old(i)
          rechxy(i)     = rech_old(i)
          deeprechxy(i) = deeprech_old(i)

          do k = 1, km
            smc(i,k)     = smc_old(i,k)
            stc(i,k)     = stc_old(i,k)
            slc(i,k)     = slc_old(i,k)
            smoiseq(i,k) = smoiseq_old(i,k)
          end do

          do k = -2,0
            tsnoxy(i,k)  = tsno_old(i,k)
            snicexy(i,k) = snice_old(i,k)
            snliqxy(i,k) = snliq_old(i,k)
          end do

          do k = -2, km
            zsnsoxy(i,k) =  zsnso_old(i,k)
          end do
       
        else
            tskin(i) = tsurf(i)    

        end if
      end do

      return

      end subroutine noahmpdrv_run
!> @}
!-----------------------------------

!> \ingroup NoahMP_LSM
!! \brief This subroutine fills in a derived data type of type noahmp_parameters with data
!! from the module \ref noahmp_tables.
      subroutine transfer_mp_parameters (vegtype,soiltype,slopetype,    &
                                       soilcolor,croptype,parameters)
     
        use noahmp_tables
        use module_sf_noahmplsm
      
        implicit none
      
        integer, intent(in)    :: vegtype
        integer, intent(in)    :: soiltype(4)
        integer, intent(in)    :: slopetype
        integer, intent(in)    :: soilcolor
        integer, intent(in)    :: croptype
          
        type (noahmp_parameters), intent(out) :: parameters
          
        real    :: refdk
        real    :: refkdt
        real    :: frzk
        real    :: frzfact
        integer :: isoil
      
        parameters%iswater   =  iswater_table
        parameters%isbarren  =  isbarren_table
        parameters%isice     =  isice_table
        parameters%iscrop    =  iscrop_table
        parameters%eblforest =  eblforest_table
      
!-----------------------------------------------------------------------&
        parameters%urban_flag = .false.
        if( vegtype == isurban_table .or. vegtype == 31                 &
     &         .or.vegtype  == 32 .or. vegtype == 33) then
           parameters%urban_flag = .true.
        endif
      
!------------------------------------------------------------------------------------------!
! transfer veg parameters
!------------------------------------------------------------------------------------------!
      
        parameters%ch2op  =  ch2op_table(vegtype)       !maximum intercepted h2o per unit lai+sai (mm)
        parameters%dleaf  =  dleaf_table(vegtype)       !characteristic leaf dimension (m)
        parameters%z0mvt  =  z0mvt_table(vegtype)       !momentum roughness length (m)
        parameters%hvt    =    hvt_table(vegtype)       !top of canopy (m)
        parameters%hvb    =    hvb_table(vegtype)       !bottom of canopy (m)
        parameters%den    =    den_table(vegtype)       !tree density (no. of trunks per m2)
        parameters%rc     =     rc_table(vegtype)       !tree crown radius (m)
        parameters%mfsno  =  mfsno_table(vegtype)       !snowmelt m parameter ()
	parameters%scffac = scffac_table(vegtype)       !snow cover factor
        parameters%saim   =   saim_table(vegtype,:)     !monthly stem area index, one-sided
        parameters%laim   =   laim_table(vegtype,:)     !monthly leaf area index, one-sided
        parameters%sla    =    sla_table(vegtype)       !single-side leaf area per kg [m2/kg]
        parameters%dilefc = dilefc_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%dilefw = dilefw_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%fragr  =  fragr_table(vegtype)       !fraction of growth respiration  !original was 0.3 
        parameters%ltovrc = ltovrc_table(vegtype)       !leaf turnover [1/s]
      
        parameters%c3psn  =  c3psn_table(vegtype)       !photosynthetic pathway: 0. = c4, 1. = c3
        parameters%kc25   =   kc25_table(vegtype)       !co2 michaelis-menten constant at 25c (pa)
        parameters%akc    =    akc_table(vegtype)       !q10 for kc25
        parameters%ko25   =   ko25_table(vegtype)       !o2 michaelis-menten constant at 25c (pa)
        parameters%ako    =    ako_table(vegtype)       !q10 for ko25
        parameters%vcmx25 = vcmx25_table(vegtype)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
        parameters%avcmx  =  avcmx_table(vegtype)       !q10 for vcmx25
        parameters%bp     =     bp_table(vegtype)       !minimum leaf conductance (umol/m**2/s)
        parameters%mp     =     mp_table(vegtype)       !slope of conductance-to-photosynthesis relationship
        parameters%qe25   =   qe25_table(vegtype)       !quantum efficiency at 25c (umol co2 / umol photon)
        parameters%aqe    =    aqe_table(vegtype)       !q10 for qe25
        parameters%rmf25  =  rmf25_table(vegtype)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
        parameters%rms25  =  rms25_table(vegtype)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%rmr25  =  rmr25_table(vegtype)       !root maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%arm    =    arm_table(vegtype)       !q10 for maintenance respiration
        parameters%folnmx = folnmx_table(vegtype)       !foliage nitrogen concentration when f(n)=1 (%)
        parameters%tmin   =   tmin_table(vegtype)       !minimum temperature for photosynthesis (k)
      
        parameters%xl     =     xl_table(vegtype)       !leaf/stem orientation index
        parameters%rhol   =   rhol_table(vegtype,:)     !leaf reflectance: 1=vis, 2=nir
        parameters%rhos   =   rhos_table(vegtype,:)     !stem reflectance: 1=vis, 2=nir
        parameters%taul   =   taul_table(vegtype,:)     !leaf transmittance: 1=vis, 2=nir
        parameters%taus   =   taus_table(vegtype,:)     !stem transmittance: 1=vis, 2=nir
      
        parameters%mrp    =    mrp_table(vegtype)       !microbial respiration parameter (umol co2 /kg c/ s)
        parameters%cwpvt  =  cwpvt_table(vegtype)       !empirical canopy wind parameter
      
        parameters%wrrat  =  wrrat_table(vegtype)       !wood to non-wood ratio
        parameters%wdpool = wdpool_table(vegtype)       !wood pool (switch 1 or 0) depending on woody or not [-]
        parameters%tdlef  =  tdlef_table(vegtype)       !characteristic t for leaf freezing [k]
      
        parameters%nroot  =  nroot_table(vegtype)       !number of soil layers with root present
        parameters%rgl    =    rgl_table(vegtype)       !parameter used in radiation stress function
        parameters%rsmin  =     rs_table(vegtype)       !minimum stomatal resistance [s m-1]
        parameters%hs     =     hs_table(vegtype)       !parameter used in vapor pressure deficit function
        parameters%topt   =   topt_table(vegtype)       !optimum transpiration air temperature [k]
        parameters%rsmax  =  rsmax_table(vegtype)       !maximal stomatal resistance [s m-1]
      
!------------------------------------------------------------------------------------------!
! transfer rad parameters
!------------------------------------------------------------------------------------------!
      
         parameters%albsat    = albsat_table(soilcolor,:)
         parameters%albdry    = albdry_table(soilcolor,:)
         parameters%albice    = albice_table
         parameters%alblak    = alblak_table               
         parameters%omegas    = omegas_table
         parameters%betads    = betads_table
         parameters%betais    = betais_table
         parameters%eg        = eg_table
      
!------------------------------------------------------------------------------------------!
! Transfer crop parameters
!------------------------------------------------------------------------------------------!

      if(croptype > 0) then
        parameters%pltday    =    pltday_table(croptype)    ! planting date
        parameters%hsday     =     hsday_table(croptype)    ! harvest date
        parameters%plantpop  =  plantpop_table(croptype)    ! plant density [per ha] - used?
        parameters%irri      =      irri_table(croptype)    ! irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
        parameters%gddtbase  =  gddtbase_table(croptype)    ! base temperature for gdd accumulation [c]
        parameters%gddtcut   =   gddtcut_table(croptype)    ! upper temperature for gdd accumulation [c]
        parameters%gdds1     =     gdds1_table(croptype)    ! gdd from seeding to emergence
        parameters%gdds2     =     gdds2_table(croptype)    ! gdd from seeding to initial vegetative 
        parameters%gdds3     =     gdds3_table(croptype)    ! gdd from seeding to post vegetative 
        parameters%gdds4     =     gdds4_table(croptype)    ! gdd from seeding to intial reproductive
        parameters%gdds5     =     gdds5_table(croptype)    ! gdd from seeding to pysical maturity 
        parameters%c3c4      =      c3c4_table(croptype)    ! photosynthetic pathway:  1. = c3 2. = c4
        parameters%aref      =      aref_table(croptype)    ! reference maximum co2 assimulation rate 
        parameters%psnrf     =     psnrf_table(croptype)    ! co2 assimulation reduction factor(0-1) (e.g.pests, weeds)
        parameters%i2par     =     i2par_table(croptype)    ! fraction of incoming solar radiation to photosynthetically active radiation
        parameters%tassim0   =   tassim0_table(croptype)    ! minimum temperature for co2 assimulation [c]
        parameters%tassim1   =   tassim1_table(croptype)    ! co2 assimulation linearly increasing until temperature reaches t1 [c]
        parameters%tassim2   =   tassim2_table(croptype)    ! co2 assmilation rate remain at aref until temperature reaches t2 [c]
        parameters%k         =         k_table(croptype)    ! light extinction coefficient
        parameters%epsi      =      epsi_table(croptype)    ! initial light use efficiency
        parameters%q10mr     =     q10mr_table(croptype)    ! q10 for maintainance respiration
        parameters%foln_mx   =   foln_mx_table(croptype)    ! foliage nitrogen concentration when f(n)=1 (%)
        parameters%lefreez   =   lefreez_table(croptype)    ! characteristic t for leaf freezing [k]
        parameters%dile_fc   =   dile_fc_table(croptype,:)  ! coeficient for temperature leaf stress death [1/s]
        parameters%dile_fw   =   dile_fw_table(croptype,:)  ! coeficient for water leaf stress death [1/s]
        parameters%fra_gr    =    fra_gr_table(croptype)    ! fraction of growth respiration
        parameters%lf_ovrc   =   lf_ovrc_table(croptype,:)  ! fraction of leaf turnover  [1/s]
        parameters%st_ovrc   =   st_ovrc_table(croptype,:)  ! fraction of stem turnover  [1/s]
        parameters%rt_ovrc   =   rt_ovrc_table(croptype,:)  ! fraction of root tunrover  [1/s]
        parameters%lfmr25    =    lfmr25_table(croptype)    ! leaf maintenance respiration at 25c [umol co2/m**2  /s]
        parameters%stmr25    =    stmr25_table(croptype)    ! stem maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%rtmr25    =    rtmr25_table(croptype)    ! root maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%grainmr25 = grainmr25_table(croptype)    ! grain maintenance respiration at 25c [umol co2/kg bio/s]
        parameters%lfpt      =      lfpt_table(croptype,:)  ! fraction of carbohydrate flux to leaf
        parameters%stpt      =      stpt_table(croptype,:)  ! fraction of carbohydrate flux to stem
        parameters%rtpt      =      rtpt_table(croptype,:)  ! fraction of carbohydrate flux to root
        parameters%grainpt   =   grainpt_table(croptype,:)  ! fraction of carbohydrate flux to grain
        parameters%bio2lai   =   bio2lai_table(croptype)    ! leaf are per living leaf biomass [m^2/kg]
      end if

!------------------------------------------------------------------------------------------!
! transfer global parameters
!------------------------------------------------------------------------------------------!
      
         parameters%co2          =          co2_table
         parameters%o2           =           o2_table
         parameters%timean       =       timean_table
         parameters%fsatmx       =       fsatmx_table
         parameters%z0sno        =        z0sno_table
         parameters%ssi          =          ssi_table
         parameters%swemx        =        swemx_table
         parameters%tau0         =         tau0_table
         parameters%grain_growth = grain_growth_table
         parameters%extra_growth = extra_growth_table
         parameters%dirt_soot    =    dirt_soot_table
         parameters%bats_cosz    =    bats_cosz_table
         parameters%bats_vis_new = bats_vis_new_table
         parameters%bats_nir_new = bats_nir_new_table
         parameters%bats_vis_age = bats_vis_age_table
         parameters%bats_nir_age = bats_nir_age_table
         parameters%bats_vis_dir = bats_vis_dir_table
         parameters%bats_nir_dir = bats_nir_dir_table
         parameters%rsurf_snow   =   rsurf_snow_table
         parameters%rsurf_exp    =    rsurf_exp_table
         parameters%snow_emis    =    snow_emis_table
      
! ----------------------------------------------------------------------
!  transfer soil parameters
! ----------------------------------------------------------------------
      
      do isoil = 1, size(soiltype)
        parameters%bexp(isoil)   = bexp_table   (soiltype(isoil))
        parameters%dksat(isoil)  = dksat_table  (soiltype(isoil))
        parameters%dwsat(isoil)  = dwsat_table  (soiltype(isoil))
        parameters%psisat(isoil) = psisat_table (soiltype(isoil))
        parameters%quartz(isoil) = quartz_table (soiltype(isoil))
        parameters%smcdry(isoil) = smcdry_table (soiltype(isoil))
        parameters%smcmax(isoil) = smcmax_table (soiltype(isoil))
        parameters%smcref(isoil) = smcref_table (soiltype(isoil))
        parameters%smcwlt(isoil) = smcwlt_table (soiltype(isoil))
      end do
          
      parameters%f1     = f1_table(soiltype(1))
      parameters%refdk  = refdk_table
      parameters%refkdt = refkdt_table

! ----------------------------------------------------------------------
! transfer genparm parameters
! ----------------------------------------------------------------------
          parameters%csoil  = csoil_table
          parameters%zbot   = zbot_table
          parameters%czil   = czil_table
      
          frzk   = frzk_table
          parameters%kdt    = parameters%refkdt * parameters%dksat(1) / parameters%refdk
          parameters%slope  = slope_table(slopetype)
      
          if(parameters%urban_flag)then  ! hardcoding some urban parameters for soil
             parameters%smcmax = 0.45 
             parameters%smcref = 0.42 
             parameters%smcwlt = 0.40 
             parameters%smcdry = 0.40 
             parameters%csoil  = 3.e6
          endif
      
      ! adjust frzk parameter to actual soil type: frzk * frzfact
      
!-----------------------------------------------------------------------&
          if(soiltype(1) /= 14) then
            frzfact = (parameters%smcmax(1) / parameters%smcref(1)) * (0.412 / 0.468)
            parameters%frzx = frzk * frzfact
          end if
      
       end subroutine transfer_mp_parameters

!> \ingroup NoahMP_LSM
!! \brief This subroutine uses a pedotransfer method to calculate soil properties.
SUBROUTINE PEDOTRANSFER_SR2006(nsoil,sand,clay,orgm,parameters)

  use module_sf_noahmplsm
  use noahmp_tables
        
  implicit none
        
  integer,                    intent(in   ) :: nsoil     ! number of soil layers
  real, dimension( 1:nsoil ), intent(inout) :: sand
  real, dimension( 1:nsoil ), intent(inout) :: clay
  real, dimension( 1:nsoil ), intent(inout) :: orgm
    
  real, dimension( 1:nsoil ) :: theta_1500t
  real, dimension( 1:nsoil ) :: theta_1500
  real, dimension( 1:nsoil ) :: theta_33t
  real, dimension( 1:nsoil ) :: theta_33
  real, dimension( 1:nsoil ) :: theta_s33t
  real, dimension( 1:nsoil ) :: theta_s33
  real, dimension( 1:nsoil ) :: psi_et
  real, dimension( 1:nsoil ) :: psi_e
    
  type(noahmp_parameters), intent(inout) :: parameters
  integer :: k

  do k = 1,4
    if(sand(k) <= 0 .or. clay(k) <= 0) then
      sand(k) = 0.41
      clay(k) = 0.18
    end if
    if(orgm(k) <= 0 ) orgm(k) = 0.0
  end do
        
  theta_1500t =   sr2006_theta_1500t_a*sand       &
                + sr2006_theta_1500t_b*clay       &
                + sr2006_theta_1500t_c*orgm       &
                + sr2006_theta_1500t_d*sand*orgm  &
                + sr2006_theta_1500t_e*clay*orgm  &
                + sr2006_theta_1500t_f*sand*clay  &
                + sr2006_theta_1500t_g

  theta_1500  =   theta_1500t                      &
                + sr2006_theta_1500_a*theta_1500t  &
                + sr2006_theta_1500_b

  theta_33t   =   sr2006_theta_33t_a*sand       &
                + sr2006_theta_33t_b*clay       &
                + sr2006_theta_33t_c*orgm       &
                + sr2006_theta_33t_d*sand*orgm  &
                + sr2006_theta_33t_e*clay*orgm  &
                + sr2006_theta_33t_f*sand*clay  &
                + sr2006_theta_33t_g

  theta_33    =   theta_33t                              &
                + sr2006_theta_33_a*theta_33t*theta_33t  &
                + sr2006_theta_33_b*theta_33t            &
                + sr2006_theta_33_c

  theta_s33t  =   sr2006_theta_s33t_a*sand      &
                + sr2006_theta_s33t_b*clay      &
                + sr2006_theta_s33t_c*orgm      &
                + sr2006_theta_s33t_d*sand*orgm &
                + sr2006_theta_s33t_e*clay*orgm &
                + sr2006_theta_s33t_f*sand*clay &
                + sr2006_theta_s33t_g

  theta_s33   = theta_s33t                       &
                + sr2006_theta_s33_a*theta_s33t  &
                + sr2006_theta_s33_b

  psi_et      =   sr2006_psi_et_a*sand           &
                + sr2006_psi_et_b*clay           &
                + sr2006_psi_et_c*theta_s33      &
                + sr2006_psi_et_d*sand*theta_s33 &
                + sr2006_psi_et_e*clay*theta_s33 &
                + sr2006_psi_et_f*sand*clay      &
                + sr2006_psi_et_g
 
  psi_e       =   psi_et                        &
                + sr2006_psi_e_a*psi_et*psi_et  &
                + sr2006_psi_e_b*psi_et         &
                + sr2006_psi_e_c
    
  parameters%smcwlt = theta_1500
  parameters%smcref = theta_33
  parameters%smcmax =   theta_33    &
                      + theta_s33            &
                      + sr2006_smcmax_a*sand &
                      + sr2006_smcmax_b

  parameters%bexp   = 3.816712826 / (log(theta_33) - log(theta_1500) )
  parameters%psisat = psi_e
  parameters%dksat  = 1930.0 * (parameters%smcmax - theta_33) ** (3.0 - 1.0/parameters%bexp)
  parameters%quartz = sand
    
! Units conversion
    
  parameters%psisat = max(0.1,parameters%psisat)     ! arbitrarily impose a limit of 0.1kpa
  parameters%psisat = 0.101997 * parameters%psisat   ! convert kpa to m
  parameters%dksat  = parameters%dksat / 3600000.0   ! convert mm/h to m/s
  parameters%dwsat  = parameters%dksat * parameters%psisat *parameters%bexp / parameters%smcmax  ! units should be m*m/s
  parameters%smcdry = parameters%smcwlt
  
! Introducing somewhat arbitrary limits (based on SOILPARM) to prevent bad things
  
  parameters%smcmax = max(0.32 ,min(parameters%smcmax,             0.50 ))
  parameters%smcref = max(0.17 ,min(parameters%smcref,parameters%smcmax ))
  parameters%smcwlt = max(0.01 ,min(parameters%smcwlt,parameters%smcref ))
  parameters%smcdry = max(0.01 ,min(parameters%smcdry,parameters%smcref ))
  parameters%bexp   = max(2.50 ,min(parameters%bexp,               12.0 ))
  parameters%psisat = max(0.03 ,min(parameters%psisat,             1.00 ))
  parameters%dksat  = max(5.e-7,min(parameters%dksat,              1.e-5))
  parameters%dwsat  = max(1.e-6,min(parameters%dwsat,              3.e-5))
  parameters%quartz = max(0.05 ,min(parameters%quartz,             0.95 ))
    
 END SUBROUTINE PEDOTRANSFER_SR2006

!-----------------------------------------------------------------------&

!> \ingroup NoahMP_LSM
!! brief Calculate potential evaporation for the current point. Various
!! partial sums/products are also calculated and passed back to the
!! calling routine for later use.
      subroutine penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,ssoil,     &
     &                   q2,q2sat,etp,snowng,frzgra,ffrozp,             &
     &                   dqsdt2,emissi_in,sncovr)
 
! etp is calcuated right after ssoil

! ----------------------------------------------------------------------
! subroutine penman
! ----------------------------------------------------------------------
      implicit none
      logical, intent(in)     :: snowng, frzgra
      real, intent(in)        :: ch, dqsdt2,fdown,prcp,ffrozp,          &
     &                           q2, q2sat,ssoil, sfcprs, sfctmp,       &
     &                           t2v, th2,emissi_in,sncovr
      real, intent(out)       :: etp
      real                    :: epsca,flx2,rch,rr,t24
      real                    :: a, delta, fnet,rad,rho,emissi,elcp1,lvs

      real, parameter :: elcp = 2.4888e+3, lsubc = 2.501000e+6,cp = 1004.6
      real, parameter :: lsubs = 2.83e+6, rd = 287.05, cph2o = 4.1855e+3
      real, parameter :: cpice = 2.106e+3, lsubf   = 3.335e5  
      real, parameter :: sigma = 5.6704e-8

! ----------------------------------------------------------------------
! executable code begins here:
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! prepare partial quantities for penman equation.
! ----------------------------------------------------------------------
        emissi=emissi_in
!       elcp1  = (1.0-sncovr)*elcp  + sncovr*elcp*lsubs/lsubc
        lvs    = (1.0-sncovr)*lsubc + sncovr*lsubs

      flx2 = 0.0
      delta = elcp * dqsdt2
!     delta = elcp1 * dqsdt2
      t24 = sfctmp * sfctmp * sfctmp * sfctmp
       rr = t24 * 6.48e-8 / (sfcprs * ch) + 1.0
!     rr = emissi*t24 * 6.48e-8 / (sfcprs * ch) + 1.0
      rho = sfcprs / (rd * t2v)

! ----------------------------------------------------------------------
! adjust the partial sums / products with the latent heat
! effects caused by falling precipitation.
! ----------------------------------------------------------------------
      rch = rho * cp * ch
      if (.not. snowng) then
         if (prcp >  0.0) rr = rr + cph2o * prcp / rch
      else
! ---- ...  fractional snowfall/rainfall
        rr = rr + (cpice*ffrozp+cph2o*(1.-ffrozp))                      &
     &       *prcp/rch
      end if

! ----------------------------------------------------------------------
! include the latent heat effects of frzng rain converting to ice on
! impact in the calculation of flx2 and fnet.
! ----------------------------------------------------------------------
!      fnet = fdown - sigma * t24- ssoil
      fnet = fdown -  emissi*sigma * t24- ssoil
      if (frzgra) then
         flx2 = - lsubf * prcp
         fnet = fnet - flx2
! ----------------------------------------------------------------------
! finish penman equation calculations.
! ----------------------------------------------------------------------
      end if
      rad = fnet / rch + th2- sfctmp
       a = elcp * (q2sat - q2)
!     a = elcp1 * (q2sat - q2)
      epsca = (a * rr + rad * delta) / (delta + rr)
       etp = epsca * rch / lsubc
!     etp = epsca * rch / lvs

! ----------------------------------------------------------------------
      end subroutine penman

      end module noahmpdrv
