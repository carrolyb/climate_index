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

; switch sign for all data
switchsign = 0

; swap or not
swap = 1 ; or 0

;; plot which EOF 
; EOF1 => nn=0
; EOF2 => nn=1
nn       = 0
whichEOF = tostring(nn+1)

; read data
filename = "SPDO_annualSST_filtering_localdetrending_mpi_merged_all"

anew                   = addfile("./Data/"+filename+".nc","rw")
eof_ts                 = anew->eof_ts
eof_ts_Notstandardized = anew->eof_ts_Notstandardized

eof_pcvar               = anew->eof_pcvar
SPDO4d                  = anew->SPDO4d
SPDO4d_all              = anew->SPDO4d_all
SPDO4d_all_reg_notstand = anew->SPDO4d_all_reg_notstand
SPDO4d_all_reg_stand    = anew->SPDO4d_all_reg_stand
eof4d                   = anew->eof4d

numofens = dimsizes(SPDO4d(:,0,0,0))
numofeof = 3

if (switchsign .eq. 1 )

eof_ts                 = -eof_ts
eof_ts_Notstandardized = -eof_ts_Notstandardized

eof4d                   = -eof4d
SPDO4D                  = -SPDO4d
SPDO4D_all              = -SPDO4d_all
SPDO4D_all_reg_notstand = -SPDO4d_all_reg_notstand
SPDO4D_all_reg_stand    = -SPDO4d_all_reg_stand

end if


; ens avg
eof                     = dim_avg_n_Wrap(eof4d,0)
SPDO                    = dim_avg_n_Wrap(SPDO4d,0)
SPDO_all                = dim_avg_n_Wrap(SPDO4d_all,0)
SPDO_all_reg_notstand   = dim_avg_n_Wrap(SPDO4d_all_reg_notstand,0)
SPDO_all_reg_stand      = dim_avg_n_Wrap(SPDO4d_all_reg_stand,0)

if (swap .eq. 1 )

;; 1) checking the pattern to be of the same sign

re = pattern_cor(conform(eof4d, eof,(/1,2,3/)), eof4d, 1.0, 0)    


do i=0, numofens-1
do j=0, 1 ; numofEOF

if (re(i,j) .lt. 0)

eof_ts(i,j,:)                    = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)    = -eof_ts_Notstandardized(i,j,:)

eof4d(i,j,:,:)                   = -eof4d(i,j,:,:)
SPDO4d(i,j,:,:)                  = -SPDO4d(i,j,:,:)
SPDO4d_all(i,j,:,:)              = -SPDO4d_all(i,j,:,:)
SPDO4d_all_reg_notstand(i,j,:,:) = -SPDO4d_all_reg_notstand(i,j,:,:)
SPDO4d_all_reg_stand(i,j,:,:)    = -SPDO4d_all_reg_stand(i,j,:,:)
end if

end do
end do

; re do ens avg
eof                    = dim_avg_n_Wrap(eof4d,0)
SPDO                   = dim_avg_n_Wrap(SPDO4d,0)
SPDO_all               = dim_avg_n_Wrap(SPDO4d_all,0)
SPDO_all_reg_notstand  = dim_avg_n_Wrap(SPDO4d_all_reg_notstand,0)
SPDO_all_reg_stand     = dim_avg_n_Wrap(SPDO4d_all_reg_stand,0)

;; 2) swap EOF1 and 2 if needed 

re0=new( (/numofens/),typeof(eof4d) )
re1=new( (/numofens/),typeof(eof4d) )

do i=0, numofens-1
re0(i)=pattern_cor(SPDO(0,:,:), SPDO4d(i,0,:,:),1.0,0) 
re1(i)=pattern_cor(SPDO(0,:,:), SPDO4d(i,1,:,:),1.0,0)  
end do

do i=0, numofens-1

if (abs(re0(i)) .lt. abs(re1(i)))

dummy=eof_ts(i,0,:)
eof_ts(i,0,:)=eof_ts(i,1,:)
eof_ts(i,1,:)=dummy
delete(dummy)

dummy=eof_ts_Notstandardized(i,0,:)
eof_ts_Notstandardized(i,0,:)=eof_ts_Notstandardized(i,1,:)
eof_ts_Notstandardized(i,1,:)=dummy
delete(dummy)

dummy=eof4d(i,0,:,:)
eof4d(i,0,:,:)=eof4d(i,1,:,:)
eof4d(i,1,:,:)=dummy
delete(dummy)

dummy=SPDO4d(i,0,:,:)
SPDO4d(i,0,:,:)=SPDO4d(i,1,:,:)
SPDO4d(i,1,:,:)=dummy
delete(dummy)

dummy=SPDO4d_all(i,0,:,:)
SPDO4d_all(i,0,:,:)=SPDO4d_all(i,1,:,:)
SPDO4d_all(i,1,:,:)=dummy
delete(dummy)

dummy=SPDO4d_all_reg_notstand(i,0,:,:)
SPDO4d_all_reg_notstand(i,0,:,:)=SPDO4d_all_reg_notstand(i,1,:,:)
SPDO4d_all_reg_notstand(i,1,:,:)=dummy
delete(dummy)

dummy=SPDO4d_all_reg_stand(i,0,:,:)
SPDO4d_all_reg_stand(i,0,:,:)=SPDO4d_all_reg_stand(i,1,:,:)
SPDO4d_all_reg_stand(i,1,:,:)=dummy
delete(dummy)

dummy=eof_pcvar(i,0)
eof_pcvar(i,0)=eof_pcvar(i,1)
eof_pcvar(i,1)=dummy
delete(dummy)

end if
end do 

eof                   = dim_avg_n_Wrap(eof4d,0)
PDO                   = dim_avg_n_Wrap(SPDO4d,0)
PDO_all               = dim_avg_n_Wrap(SPDO4d_all,0)
PDO_all_reg_notstand  = dim_avg_n_Wrap(SPDO4d_all_reg_notstand,0)
PDO_all_reg_stand     = dim_avg_n_Wrap(SPDO4d_all_reg_stand,0)

;**********************************************************************
do i=0, numofens-1
do j=0, 1 ; numofEOF

eof_ts(i,j,:)                    = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)    = -eof_ts_Notstandardized(i,j,:)

eof4d(i,j,:,:)                   = -eof4d(i,j,:,:)
SPDO4d(i,j,:,:)                  = -SPDO4d(i,j,:,:)
SPDO4d_all(i,j,:,:)              = -SPDO4d_all(i,j,:,:)
SPDO4d_all_reg_notstand(i,j,:,:) = -SPDO4d_all_reg_notstand(i,j,:,:)
SPDO4d_all_reg_stand(i,j,:,:)    = -SPDO4d_all_reg_stand(i,j,:,:)

end do
end do

eof                              = dim_avg_n_Wrap(eof4d,0)
SPDO                             = dim_avg_n_Wrap(SPDO4d,0)
SPDO_all                         = dim_avg_n_Wrap(SPDO4d_all,0)
SPDO_all_reg_notstand            = dim_avg_n_Wrap(SPDO4d_all_reg_notstand,0)
SPDO_all_reg_stand               = dim_avg_n_Wrap(SPDO4d_all_reg_stand,0)
;***************************************************************************

;; re do 1 just to make sure

;;  1) checking the pattern to be of the same sign

re = pattern_cor(conform(eof4d, eof, (/1,2,3/)), eof4d, 1.0, 0)  

do i=0, numofens-1
do j=0,1 ; numofEOF

if (re(i,j) .lt. 0)

eof_ts(i,j,:)                    = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)    = -eof_ts_Notstandardized(i,j,:)
eof4d(i,j,:,:)                   = -eof4d(i,j,:,:)
SPDO4d(i,j,:,:)                  = -SPDO4d(i,j,:,:)
SPDO4d_all(i,j,:,:)              = -SPDO4d_all(i,j,:,:)
SPDO4d_all_reg_notstand(i,j,:,:) = -SPDO4d_all_reg_notstand(i,j,:,:)
SPDO4d_all_reg_stand(i,j,:,:)    = -SPDO4d_all_reg_stand(i,j,:,:)

end if

end do
end do


; ens avg
eof                   = dim_avg_n_Wrap(eof4d,0)
SPDO                  = dim_avg_n_Wrap(SPDO4d,0)
SPDO_all              = dim_avg_n_Wrap(SPDO4d_all,0)
SPDO_all_reg_notstand = dim_avg_n_Wrap(SPDO4d_all_reg_notstand,0)
SPDO_all_reg_stand    = dim_avg_n_Wrap(SPDO4d_all_reg_stand,0)

end if

;**********************************************************************
system("/bin/rm -f "+"./Data_Swap/"+filename+"_Swap.nc" )

aanew=addfile("./Data_Swap/"+filename+"_Swap.nc","c")
aanew->eof_ts                  = eof_ts
aanew->eof_ts_Notstandardized  = eof_ts_Notstandardized
aanew->eof_pcvar               = eof_pcvar  
aanew->eof4d                   = eof4d
aanew->SPDO4d                  = SPDO4d
aanew->SPDO4d_all              = SPDO4d_all
aanew->SPDO4d_all_reg_notstand = SPDO4d_all_reg_notstand
aanew->SPDO4d_all_reg_stand    = SPDO4d_all_reg_stand
aanew->eof                     = eof
aanew->SPDO                    = SPDO
aanew->SPDO_all                = SPDO_all
aanew->SPDO_all_reg_notstand   = SPDO_all_reg_notstand
aanew->SPDO_all_reg_stand      = SPDO_all_reg_stand

end
