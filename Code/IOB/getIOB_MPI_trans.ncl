; Analysis done to MPI-ESM1.1 (1850-2020) to obtain IOB index
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
nn=0
whichEOF = tostring(nn+1)

; read data
filename="IOB_annualSST_filtering_localdetrending_mpi_merged_all"

anew  = addfile("./Data/"+filename+".nc","rw")
eof_ts                 = anew->eof_ts
eof_ts_Notstandardized = anew->eof_ts_Notstandardized

eof_pcvar              = anew->eof_pcvar
IOB4d                  = anew->IOB4d
IOB4d_all              = anew->IOB4d_all
IOB4d_all_reg_notstand = anew->IOB4d_all_reg_notstand
IOB4d_all_reg_stand    = anew->IOB4d_all_reg_stand
eof4d                  = anew->eof4d

numofens               = dimsizes(IOB4d(:,0,0,0))
numofeof               = 3

if (switchsign .eq. 1 )

eof_ts                 = -eof_ts
eof_ts_Notstandardized = -eof_ts_Notstandardized

eof4d                  = -eof4d
IOB4D                  = -IOB4d
IOB4D_all              = -IOB4d_all
IOB4D_all_reg_notstand = -IOB4d_all_reg_notstand
IOB4D_all_reg_stand    = -IOB4d_all_reg_stand

end if

; ens avg
eof                    = dim_avg_n_Wrap(eof4d,0)
IOB                    = dim_avg_n_Wrap(IOB4d,0)
IOB_all                = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand   = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand      = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)

if (swap .eq. 1 )

;; 1) checking the pattern to be of the same sign
re = pattern_cor(conform(eof4d, eof,(/1,2,3/)), eof4d, 1.0, 0)    

do i=0, numofens-1
do j=0, 1 ; numofEOF

if (re(i,j) .lt. 0)

eof_ts(i,j,:)                   = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)   = -eof_ts_Notstandardized(i,j,:)

eof4d(i,j,:,:)                  = -eof4d(i,j,:,:)
IOB4d(i,j,:,:)                  = -IOB4d(i,j,:,:)
IOB4d_all(i,j,:,:)              = -IOB4d_all(i,j,:,:)
IOB4d_all_reg_notstand(i,j,:,:) = -IOB4d_all_reg_notstand(i,j,:,:)
IOB4d_all_reg_stand(i,j,:,:)    = -IOB4d_all_reg_stand(i,j,:,:)
end if

end do
end do

; re do ens avg
eof                  = dim_avg_n_Wrap(eof4d,0)
IOB                  = dim_avg_n_Wrap(IOB4d,0)
IOB_all              = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand    = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)

;; 2) swap EOF1 and 2 if needed 

re0=new( (/numofens/),typeof(eof4d) )
re1=new( (/numofens/),typeof(eof4d) )

do i=0, numofens-1
re0(i)=pattern_cor(IOB(0,:,:), IOB4d(i,0,:,:),1.0,0)  
re1(i)=pattern_cor(IOB(0,:,:), IOB4d(i,1,:,:),1.0,0)  
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

dummy=IOB4d(i,0,:,:)
IOB4d(i,0,:,:)=IOB4d(i,1,:,:)
IOB4d(i,1,:,:)=dummy
delete(dummy)

dummy=IOB4d_all(i,0,:,:)
IOB4d_all(i,0,:,:)=IOB4d_all(i,1,:,:)
IOB4d_all(i,1,:,:)=dummy
delete(dummy)

dummy=IOB4d_all_reg_notstand(i,0,:,:)
IOB4d_all_reg_notstand(i,0,:,:)=IOB4d_all_reg_nostand(i,1,:,:)
IOB4d_all_reg_notstad(i,1,:,:)=dummy
delete(dummy)

dummy=IOB4d_all_reg_stand(i,0,:,:)
IOB4d_all_reg_stand(i,0,:,:)=IOB4d_all_reg_stand(i,1,:,:)
IOB4d_all_reg_stand(i,1,:,:)=dummy
delete(dummy)

dummy=eof_pcvar(i,0)
eof_pcvar(i,0)=eof_pcvar(i,1)
eof_pcvar(i,1)=dummy
delete(dummy)

end if
end do 

eof                  = dim_avg_n_Wrap(eof4d,0)
IOB                  = dim_avg_n_Wrap(IOB4d,0)
IOB_all              = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand    = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)

;**********************************************************************
do i=0, numofens-1
do j=0, 1 ; numofEOF

eof_ts(i,j,:)                   = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)   = -eof_ts_Notstandardized(i,j,:)

eof4d(i,j,:,:)                  = -eof4d(i,j,:,:)
IOB4d(i,j,:,:)                  = -IOB4d(i,j,:,:)
IOB4d_all(i,j,:,:)              = -IOB4d_all(i,j,:,:)
IOB4d_all_reg_notstand(i,j,:,:) = -IOB4d_all_reg_notstand(i,j,:,:)
IOB4d_all_reg_stand(i,j,:,:)    = -IOB4d_all_reg_stand(i,j,:,:)

end do
end do

eof                             = dim_avg_n_Wrap(eof4d,0)
IOB                             = dim_avg_n_Wrap(IOB4d,0)
IOB_all                         = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand            = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand               = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)
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
IOB4d(i,j,:,:)                   = -IOB4d(i,j,:,:)
IOB4d_all(i,j,:,:)               = -IOB4d_all(i,j,:,:)
IOB4d_all_reg_notstand(i,j,:,:)  = -IOB4d_all_reg_notstand(i,j,:,:)
IOB4d_all_reg_stand(i,j,:,:)     = -IOB4d_all_reg_stand(i,j,:,:)

end if

end do
end do


; ens avg
eof                              = dim_avg_n_Wrap(eof4d,0)
IOB                              = dim_avg_n_Wrap(IOB4d,0)
IOB_all                          = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand             = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand                = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)

end if

;**********************************************************************
do i=0, numofens-1
do j=0, 1 ; numofEOF

eof_ts(i,j,:)                    = -eof_ts(i,j,:)
eof_ts_Notstandardized(i,j,:)    = -eof_ts_Notstandardized(i,j,:)

eof4d(i,j,:,:)                   = -eof4d(i,j,:,:)
IOB4d(i,j,:,:)                   = -IOB4d(i,j,:,:)
IOB4d_all(i,j,:,:)               = -IOB4d_all(i,j,:,:)
IOB4d_all_reg_notstand(i,j,:,:)  = -IOB4d_all_reg_notstand(i,j,:,:)
IOB4d_all_reg_stand(i,j,:,:)     = -IOB4d_all_reg_stand(i,j,:,:)

end do
end do

eof                              = dim_avg_n_Wrap(eof4d,0)
IOB                              = dim_avg_n_Wrap(IOB4d,0)
IOB_all                          = dim_avg_n_Wrap(IOB4d_all,0)
IOB_all_reg_notstand             = dim_avg_n_Wrap(IOB4d_all_reg_notstand,0)
IOB_all_reg_stand                = dim_avg_n_Wrap(IOB4d_all_reg_stand,0)
;***************************************************************************

; write back to another nc file
system("/bin/rm -f "+"./Data_Swap/"+filename+"_Swap.nc" )

aanew        = addfile("./Data_Swap/"+filename+"_Swap.nc","c")
aanew->eof_ts                 = eof_ts
aanew->eof_ts_Notstandardized = eof_ts_Notstandardized
aanew->eof_pcvar              = eof_pcvar  
aanew->eof4d                  = eof4d
aanew->IOB4d                  = IOB4d
aanew->IOB4d_all              = IOB4d_all
aanew->IOB4d_all_reg_notstand = IOB4d_all_reg_notstand
aanew->IOB4d_all_reg_stand    = IOB4d_all_reg_stand
aanew->eof                    = eof
aanew->IOB                    = IOB
aanew->IOB_all                = IOB_all
aanew->IOB_all_reg_notstand   = IOB_all_reg_notstand
aanew->IOB_all_reg_stand      = IOB_all_reg_stand

end
