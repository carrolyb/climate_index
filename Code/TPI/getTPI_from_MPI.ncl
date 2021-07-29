; Analysis done to MPI-ESM1.1 (1850-2020) to obtain TPI index
; Authors: Bin Yao(1), Yangyang Xu(2), Andrew E. Dessler(2), Chao Liu(3)
;          1. Nanjing Joint Institute for Atmospheric Sciences, Nanjing, China
;          2. Texas A&M University, College Station, TX, USA
;          3. Nanjing University of Information Science & Technology, Nanjing, China
; Last updated in Jul. 2021.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;just 1 run for lat and lon
newf = addfile("../Data/MPIprecip_hist/precip_hist_0001.nc","r")
lat  = newf->lat
lon  = newf->lon

files_hist  = systemfunc("\ls ../Data/MPIprecip_hist/precip_hist*.nc")
files_rcp85 = systemfunc("\ls ../Data/MPIprecip_RCP85/precip_rcp85*.nc")

yrStrt_hist =1850
yrLast_hist =2005
yrStrt_rcp85=2006
yrLast_rcp85=2020

yrStrt_all  =1850
yrLast_all  =2020

;;;;;;;;;
fils_hist   =addfiles(files_hist,"r")
fils_rcp85  =addfiles(files_rcp85,"r")

TIME_hist   = fils_hist[:]->time
TEMP_hist   = tointeger(TIME_hist)
yy_hist     = TEMP_hist/10000
mmdd_hist   = TEMP_hist-yy_hist*10000
mm_hist     = mmdd_hist/100
dd_hist     = mmdd_hist-mm_hist*100
hh_hist     = dd_hist
mn_hist     = dd_hist
sc_hist     = dd_hist
hh_hist     = 0
mn_hist     = 0
sc_hist     = 0
units       = "hours since 1850-01-01"
TIMEE_hist  = cd_inv_calendar(yy_hist,mm_hist,dd_hist,hh_hist,mn_hist,sc_hist,units,0)
TIMEE_hist!0="TIMEE_hist"
YYYY_hist   = cd_calendar(TIMEE_hist,-1)/100
iYYYY_hist  = ind(YYYY_hist .ge. yrStrt_hist .and. YYYY_hist .le. yrLast_hist) ;get the index of the time (month * time * ensemble)

TIME_rcp85  = fils_rcp85[:]->time
TEMP_rcp85  = tointeger(TIME_rcp85)
yy_rcp85    = TEMP_rcp85/10000
mmdd_rcp85  = TEMP_rcp85-yy_rcp85*10000
mm_rcp85    = mmdd_rcp85/100
dd_rcp85    = mmdd_rcp85-mm_rcp85*100
hh_rcp85    = dd_rcp85
mn_rcp85    = dd_rcp85
sc_rcp85    = dd_rcp85
hh_rcp85    = 0
mn_rcp85    = 0
sc_rcp85    = 0
units       = "hours since 2006-01-01"
TIMEE_rcp85 = cd_inv_calendar(yy_rcp85,mm_rcp85,dd_rcp85,hh_rcp85,mn_rcp85,sc_rcp85,units,0)
TIMEE_rcp85!0="TIMEE_rcp85"
YYYY_rcp85  = cd_calendar(TIMEE_rcp85,-1)/100
iYYYY_rcp85 = ind(YYYY_rcp85 .ge. yrStrt_rcp85 .and. YYYY_rcp85 .le. yrLast_rcp85) ;get the index of the time (month * time * ensemble)


TS_long_hist     = fils_hist[:]->temp2
TS_long_rcp85    = fils_rcp85[:]->temp2

TS_monthly_hist  = TS_long_hist(iYYYY_hist,:,:) ; read only 1850-2005
TS_monthly_rcp85 = TS_long_rcp85(iYYYY_rcp85,:,:) ; read only 2006-2020

; annual mean data ; remove seasonal cycle
TS_hist    = month_to_annual(TS_monthly_hist,1)
TS_rcp85   = month_to_annual(TS_monthly_rcp85,1)

time_hist  = ispan(yrStrt_hist,yrLast_hist,1)
time_rcp85 = ispan(yrStrt_rcp85,yrLast_rcp85,1)
time       = ispan(yrStrt_all,yrLast_all,1)


; reshpae hist to 4D (ensemble, time, lat, lon)
TSdim_hist          = dimsizes(TS_hist)
timelength_hist     = dimsizes(time_hist)
TSreshape_hist      = reshape(TS_hist, (/TSdim_hist(0)/timelength_hist,timelength_hist,TSdim_hist(1),TSdim_hist(2)/) )
TSreshape_hist!0    = "ens"
TSreshape_hist!1    = "time_hist"
TSreshape_hist&time_hist = time_hist
TSreshape_hist!2    = "lat"
TSreshape_hist&lat  = lat
TSreshape_hist!3    = "lon"
TSreshape_hist&lon  = lon

; reshpae rcp85 to 4D (ensemble, time, lat, lon)

TSdim_rcp85         = dimsizes(TS_rcp85)
timelength_rcp85    = dimsizes(time_rcp85)
TSreshape_rcp85     = reshape(TS_rcp85, (/TSdim_rcp85(0)/timelength_rcp85,timelength_rcp85,TSdim_rcp85(1),TSdim_rcp85(2)/) )
TSreshape_rcp85!0   = "ens"
TSreshape_rcp85!1   = "time_rcp85"
TSreshape_rcp85&time_rcp85 = time_rcp85
TSreshape_rcp85!2   = "lat"
TSreshape_rcp85&lat = lat
TSreshape_rcp85!3   = "lon"
TSreshape_rcp85&lon = lon

; reshpae hist_rcp85 to 4D (ensemble, time, lat, lon)
ensnumber           = 100
timelength          = dimsizes(time)
TSreshape           = new((/ensnumber, timelength,dimsizes(lat),dimsizes(lon)/),"float")
TSreshape!0         = "ens"
TSreshape!1         = "time"
TSreshape&time      = time
TSreshape!2         = "lat"
TSreshape&lat       = lat
TSreshape!3         = "lon"
TSreshape&lon       = lon


do i=0,ensnumber-1
TSreshape(i,:,:,:) = array_append_record(TSreshape_hist(i,:,:,:),TSreshape_rcp85(i,:,:,:),0)
end do

TSreshape!0    = "ens"
TSreshape!1    = "time"
TSreshape&time = time

TSreshape_all  = TSreshape

;******************************* mask TS -> SST ******************************************
;  only ocean (0) or land (1)
a          = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata     = a->LSMASK
lsm        = landsea_mask(lsdata,TSreshape&lat,TSreshape&lon)
typ        = 0

TSreshape  = mask(TSreshape,lsm.eq.typ,True) ; mask to ocean temp

; ******************************* smooth using low pass filter ****************************
TSreshape_smoothed=TSreshape
TSreshape_smoothed_all=TSreshape_all

wgts_decade    = (/1,6,19,42,71,96,106,96,71,42,19,6,1/)*1.0  ;to make it float
wgts_decade    = wgts_decade/sum(wgts_decade)
kopt           = 1  ;  kopt=1 reflective end pts ;or; kopt=0  missing end pts
TSdummy3       = wgt_runave_n_Wrap(TSreshape, wgts_decade, kopt,1)
TSdummy3_all   = wgt_runave_n_Wrap(TSreshape_all, wgts_decade, kopt,1)

TSreshape_smoothed     = TSdummy3
TSreshape_smoothed_all = TSdummy3_all

;;;****************************** detrending *****************************************
TSreshape_smoothed_detrended     = TSreshape_smoothed
TSreshape_smoothed_detrended_all = TSreshape_smoothed_all

;calculating trend of ens avg (3d)
TSreshape_smoothed_detrended_ensavg     = dim_avg_n_Wrap(TSreshape_smoothed_detrended,0)
TSreshape_smoothed_detrended_ensavg_all = dim_avg_n_Wrap(TSreshape_smoothed_detrended_all,0)

TSreshape_smoothed_detrended            = TSreshape_smoothed_detrended-conform(TSreshape_smoothed_detrended,TSreshape_smoothed_detrended_ensavg,(/1,2,3/))
TSreshape_smoothed_detrended_all        = TSreshape_smoothed_detrended_all-conform(TSreshape_smoothed_detrended_all,TSreshape_smoothed_detrended_ensavg_all,(/1,2,3/))

; ****************************************  before  analysis **************************
; re-order coordinate just for  analysis
TSreshape_smoothed_detrended_reorder     = TSreshape_smoothed_detrended(ens|:, lat|:,lon|:,time|:)
TSreshape_smoothed_detrended_all_reorder = TSreshape_smoothed_detrended_all(ens|:, lat|:,lon|:,time|:)

; weighting
rad   = 4.*atan(1.)/180.
clat  = lat
wgt   = doubletofloat(sqrt( cos(rad*clat) ) )

TSreshape_smoothed_detrended_reorder      = TSreshape_smoothed_detrended_reorder*conform(TSreshape_smoothed_detrended_reorder, wgt, 1)/sum(wgt)
TSreshape_smoothed_detrended_all_reorder = TSreshape_smoothed_detrended_all_reorder*conform(TSreshape_smoothed_detrended_all_reorder, wgt, 1)/sum(wgt)


; creating a few variables to store output
; TPI time series
tpi_ts   = new((/ensnumber, timelength/),"float")
tpi_ts!0 = "ens"
tpi_ts!1 = "time"

; TPI time series (not standardized)
tpi_ts_Notstandardized   = new((/ensnumber, timelength/),"float")
tpi_ts_Notstandardized!0 = "ens"
tpi_ts_Notstandardized!1 = "time"

; TPI pattern over the global domain
TPI3d      = new((/ensnumber, dimsizes(TSreshape_smoothed_detrended_reorder&lat), dimsizes(TSreshape_smoothed_detrended_reorder&lon)/),"float")
TPI3d!0    = "ens"
TPI3d!1    = "lat"
TPI3d!2    = "lon"
TPI3d&lat  = TSreshape_smoothed_detrended_reorder&lat
TPI3d&lon  = TSreshape_smoothed_detrended_reorder&lon

; TPI pattern over the global coverage
TPI3d_all      = new((/ensnumber, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
TPI3d_all!0    = "ens"
TPI3d_all!1    = "lat"
TPI3d_all!2    = "lon"
TPI3d_all&lat  = TSreshape_smoothed_detrended_all_reorder&lat
TPI3d_all&lon  = TSreshape_smoothed_detrended_all_reorder&lon

; TPI pattern over the global coverage: regCoef-Notstand
TPI3d_all_reg_notstand      = new((/ensnumber, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
TPI3d_all_reg_notstand!0    = "ens"
TPI3d_all_reg_notstand!1    = "lat"
TPI3d_all_reg_notstand!2    = "lon"
TPI3d_all_reg_notstand&lat  = TSreshape_smoothed_detrended_all_reorder&lat
TPI3d_all_reg_notstand&lon  = TSreshape_smoothed_detrended_all_reorder&lon

; TPI pattern over the global coverage: regCoef-stand
TPI3d_all_reg_stand     = new((/ensnumber, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
TPI3d_all_reg_stand!0   = "ens"
TPI3d_all_reg_stand!1   = "lat"
TPI3d_all_reg_stand!2   = "lon"
TPI3d_all_reg_stand&lat = TSreshape_smoothed_detrended_all_reorder&lat
TPI3d_all_reg_stand&lon = TSreshape_smoothed_detrended_all_reorder&lon

; *****************************  start  ****************
do i=0, ensnumber-1

TST      = TSreshape_smoothed_detrended_reorder(i,:,:,:)

TS1      = TST({lat|25:45},{lon|140:220}, time|:)
TS2      = TST({lat|-10:10},{lon|170:270}, time|:)
TS3      = TST({lat|-50:-15},{lon|150:200}, time|:)

SSTA1    = dim_avg_n_Wrap(TS1,(/0,1/))
SSTA2    = dim_avg_n_Wrap(TS2,(/0,1/))
SSTA3    = dim_avg_n_Wrap(TS3,(/0,1/))
;=====================================================================
;Calculate TPI
;====================================================================
TPI      = SSTA2-(SSTA1+SSTA3)/2

tpi_ts(i,:)                    = TPI(:)
tpi_ts_Notstandardized(i,:)    = tpi_ts(i,:)
tpi_ts(i,:)                    = dim_standardize(tpi_ts(i,:),0)

TPI3d(i,:,:)                   = (/ escorc(tpi_ts(i,:),TSreshape_smoothed_detrended_reorder(i,:,:,:) ) /)
TPI3d_all(i,:,:)               = (/ escorc(tpi_ts(i,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /)
TPI3d_all_reg_notstand(i,:,:)  = (/ regCoef(tpi_ts_Notstandardized(i,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /)
TPI3d_all_reg_stand(i,:,:)     = (/ regCoef(tpi_ts(i,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /)

end do
;========================  output ========================
filename = "TPI_annualSST_filtering_localdetrending_mpi_merged_all"

system("/bin/rm -rf ./Data/"+filename+".nc" ) ;
anew=addfile("./Data/"+filename+".nc","c")

anew->tpi_ts                 = tpi_ts
anew->tpi_ts_Notstandardized = tpi_ts_Notstandardized
anew->TPI3d                  = TPI3d
anew->TPI3d_all              = TPI3d_all
anew->TPI3d_all_reg_notstand = TPI3d_all_reg_notstand
anew->TPI3d_all_reg_stand    = TPI3d_all_reg_stand


end
