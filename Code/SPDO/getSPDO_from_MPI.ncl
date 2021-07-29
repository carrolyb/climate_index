; Analysis done to MPI-ESM1.1 (1850-2020) to obtain SPDO index
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
lat = newf->lat
lon = newf->lon

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

TS_long_hist    = fils_hist[:]->temp2
TS_long_rcp85   = fils_rcp85[:]->temp2

TS_monthly_hist  = TS_long_hist(iYYYY_hist,:,:) ; read only 1850-2005
TS_monthly_rcp85 = TS_long_rcp85(iYYYY_rcp85,:,:) ; read only 2006-2020


; annual mean data ; remove seasonal cycle
TS_hist    = month_to_annual(TS_monthly_hist,1)
TS_rcp85   = month_to_annual(TS_monthly_rcp85,1)

time_hist  = ispan(yrStrt_hist,yrLast_hist,1)
time_rcp85 = ispan(yrStrt_rcp85,yrLast_rcp85,1)
time       = ispan(yrStrt_all,yrLast_all,1)


; reshpae hist to 4D (ensemble, time, lat, lon)

TSdim_hist           = dimsizes(TS_hist)
timelength_hist      = dimsizes(time_hist)
TSreshape_hist       = reshape(TS_hist, (/TSdim_hist(0)/timelength_hist,timelength_hist,TSdim_hist(1),TSdim_hist(2)/) )
TSreshape_hist!0     = "ens"
TSreshape_hist!1     = "time_hist"
TSreshape_hist&time_hist = time_hist
TSreshape_hist!2     = "lat"
TSreshape_hist&lat   = lat
TSreshape_hist!3     = "lon"
TSreshape_hist&lon   = lon

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
TSreshape(i,:,:,:)=array_append_record(TSreshape_hist(i,:,:,:),TSreshape_rcp85(i,:,:,:),0)
end do

TSreshape!0    = "ens"
TSreshape!1    = "time"
TSreshape&time = time

TSreshape_all  = TSreshape
;******************************* mask TS -> SST ******************************************
;  only ocean (0) or land (1)
a      = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
lsdata = a->LSMASK
lsm    = landsea_mask(lsdata,TSreshape&lat,TSreshape&lon)
typ    = 0

TSreshape = mask(TSreshape,lsm.eq.typ,True) ; mask to ocean temp

; ******************************* smooth using low pass filter ****************************
TSreshape_smoothed     = TSreshape
TSreshape_smoothed_all = TSreshape_all

wgts_decade  = (/1,6,19,42,71,96,106,96,71,42,19,6,1/)*1.0  ;to make it float
wgts_decade  = wgts_decade/sum(wgts_decade)
kopt         = 1  ;  kopt=1 reflective end pts ;or; kopt=0  missing end pts
TSdummy3     = wgt_runave_n_Wrap(TSreshape, wgts_decade, kopt,1)
TSdummy3_all = wgt_runave_n_Wrap(TSreshape_all, wgts_decade, kopt,1)

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

; ****************************************  before EOF analysis **************************
; re-order coordinate just for EOF analysis
;  SST for EOF analysis
TSreshape_smoothed_detrended_reorder     = TSreshape_smoothed_detrended(ens|:, lat|:,lon|:,time|:)
TSreshape_smoothed_detrended_all_reorder = TSreshape_smoothed_detrended_all(ens|:, lat|:,lon|:,time|:)

; weighting
rad   = 4.*atan(1.)/180.
clat  = lat
wgt   = doubletofloat(sqrt( cos(rad*clat) ) )

TSreshape_smoothed_detrended_reorder      = TSreshape_smoothed_detrended_reorder*conform(TSreshape_smoothed_detrended_reorder, wgt, 1)/sum(wgt)
;TSreshape_smoothed_detrended_all_reorder = TSreshape_smoothed_detrended_all_reorder*conform(TSreshape_smoothed_detrended_all_reorder, wgt, 1)/sum(wgt)

; masking
; Q: S Pacific (SPDO)
latMin  = -60
latMax  = 0
lonMin  = 110
lonMax  = 290


TSreshape_smoothed_detrended_reorder_mask = TSreshape_smoothed_detrended_reorder(ens|:, {lat|latMin:latMax},{lon|lonMin:lonMax},time|:)

; creating a few variables to store EOF output

; number of EOFs
neof = 3

; PC time series
eof_ts = new((/ensnumber,neof,timelength/),"float")
eof_ts!0 = "ens"
eof_ts!1 = "evn"
eof_ts!2 = "time"

; PC time series (not standardized)
eof_ts_Notstandardized = new((/ensnumber,neof, timelength/),"float")
eof_ts_Notstandardized!0 = "ens"
eof_ts_Notstandardized!1 = "evn"
eof_ts_Notstandardized!2 = "time"

; Variance explained
eof_pcvar = new((/ensnumber,neof/),"float")
eof_pcvar!0 = "ens"
eof_pcvar!1 = "evn"

; EOF pattern over the small domain
 eof4d = new((/ensnumber, neof, dimsizes(TSreshape_smoothed_detrended_reorder_mask&lat),dimsizes(TSreshape_smoothed_detrended_reorder_mask&lon)/),"float")
eof4d!0 = "ens"
eof4d!1 = "evn"
eof4d!2 = "lat"
eof4d!3 = "lon"
eof4d&lat=TSreshape_smoothed_detrended_reorder_mask&lat
eof4d&lon=TSreshape_smoothed_detrended_reorder_mask&lon

; pattern over the global domain
SPDO4d = new((/ensnumber,neof, dimsizes(TSreshape_smoothed_detrended_reorder&lat), dimsizes(TSreshape_smoothed_detrended_reorder&lon)/),"float")
SPDO4d!0 = "ens"
SPDO4d!1 = "evn"
SPDO4d!2 = "lat"
SPDO4d!3 = "lon"
SPDO4d&lat=TSreshape_smoothed_detrended_reorder&lat
SPDO4d&lon=TSreshape_smoothed_detrended_reorder&lon

; pattern over the global coverage
SPDO4d_all = new((/ensnumber,neof, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
SPDO4d_all!0 = "ens"
SPDO4d_all!1 = "evn"
SPDO4d_all!2 = "lat"
SPDO4d_all!3 = "lon"
SPDO4d_all&lat=TSreshape_smoothed_detrended_all_reorder&lat
SPDO4d_all&lon=TSreshape_smoothed_detrended_all_reorder&lon

; pattern over the global coverage: regCoef-Notstand
SPDO4d_all_reg_notstand = new((/ensnumber,neof, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
SPDO4d_all_reg_notstand!0 = "ens"
SPDO4d_all_reg_notstand!1 = "evn"
SPDO4d_all_reg_notstand!2 = "lat"
SPDO4d_all_reg_notstand!3 = "lon"
SPDO4d_all_reg_notstand&lat=TSreshape_smoothed_detrended_all_reorder&lat
SPDO4d_all_reg_notstand&lon=TSreshape_smoothed_detrended_all_reorder&lon

; pattern over the global coverage: regCoef-stand
SPDO4d_all_reg_stand = new((/ensnumber,neof, dimsizes(TSreshape_smoothed_detrended_all_reorder&lat), dimsizes(TSreshape_smoothed_detrended_all_reorder&lon)/),"float")
SPDO4d_all_reg_stand!0 = "ens"
SPDO4d_all_reg_stand!1 = "evn"
SPDO4d_all_reg_stand!2 = "lat"
SPDO4d_all_reg_stand!3 = "lon"
SPDO4d_all_reg_stand&lat=TSreshape_smoothed_detrended_all_reorder&lat
SPDO4d_all_reg_stand&lon=TSreshape_smoothed_detrended_all_reorder&lon

; ***************************  start of EOF ********************************

optEOF = True

;*** two ways to do EOF ****
; 1 = correlation EOF ; 
; 0 = covariance EOF
optEOF@jopt = 1

do i =0, ensnumber-1

y              = eofunc_Wrap(TSreshape_smoothed_detrended_reorder_mask(i,:,:,:),neof, optEOF)
eof4d(i,:,:,:) = y
eof_pcvar(i,:) = y@pcvar

eof_ts(i,:,:)                 = eofunc_ts_Wrap(TSreshape_smoothed_detrended_reorder_mask(i,:,:,:), eof4d(i,:,:,:), False)
eof_ts_Notstandardized(i,:,:) = eof_ts(i,:,:)
eof_ts(i,:,:)                 = dim_standardize(eof_ts(i,:,:),0)

; regression (regCoef) or correlation (escorc)
; regress/correlated SST (just as it was done before EOF)
SPDO4d(i,:,:,:)                  = (/ escorc(eof_ts(i,:,:),TSreshape_smoothed_detrended_reorder(i,:,:,:) ) /) 
SPDO4d_all(i,:,:,:)              = (/ escorc(eof_ts(i,:,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /) 
SPDO4d_all_reg_notstand(i,:,:,:) = (/ regCoef(eof_ts_Notstandardized(i,:,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /) 
SPDO4d_all_reg_stand(i,:,:,:)    = (/ regCoef(eof_ts(i,:,51:167),TSreshape_smoothed_detrended_all_reorder(i,:,:,51:167) ) /) 

end do  

;========================  output ========================
filename="SPDO_annualSST_filtering_localdetrending_mpi_merged_all"

system("/bin/rm -f ./Data/"+filename+".nc" ) ;
anew=addfile("./Data/"+filename+".nc","c")

eof4d!0     = "ens"
eof4d!1     = "evn"
eof4d!2     = "latS" ; renaming the coordinates
eof4d!3     = "lonS"
eof4d&latS  = TSreshape_smoothed_detrended_reorder_mask&lat
eof4d&lonS  = TSreshape_smoothed_detrended_reorder_mask&lon
anew->eof4d = eof4d 

;time series
anew->eof_ts                 = eof_ts
anew->eof_ts_Notstandardized = eof_ts_Notstandardized
anew->eof_pcvar              = eof_pcvar 

; larger domain
anew->SPDO4d                  = SPDO4d 
anew->SPDO4d_all              = SPDO4d_all 
anew->SPDO4d_all_reg_notstand = SPDO4d_all_reg_notstand 
anew->SPDO4d_all_reg_stand    = SPDO4d_all_reg_stand 

end
