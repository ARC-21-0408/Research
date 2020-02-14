pro spec,wv,spec,ssky,Z=z,WVLO=wvlo,WVHI=wvhi,YMIN=ymin, $
         YMAX=ymax,SMOOTHIE=smoothie,SPECTRUM=spectrum,ITER=iter,LABEL=label,NOZ=noz    
   if n_params() eq 0 then begin
     print,"spec,wv,spec,ssky[,Z,WVLO,WVHI,YMIN,YMAX,SMOOTHIE,'SPECTRUM',/ITER,'LABEL',/NOZ]"
     return
   endif

;;
;; writing to ps file:
;;
; ps_open,'junk',/encapsulated,/portrait,/ps_fonts,/color
; device,xs=9.5,ys=3.5,/inch,/times
; device,xs=7.5,ys=6.5,/inch,/times
;;

!p.charsize=2.0
!p.charthick=1.4
;loadct,2
loadct,0
;invct

if not keyword_set(spectrum) then spectrum='coadd.dat'


if (Not Keyword_Set(iter)) then readcol,spectrum,pix,wv,spec,ssky
;if (Not Keyword_Set(iter)) then readcol,spectrum,wv,spec,ssky
ny=n_elements(wv)
print,'Number of columns is ',ny
iarr=findgen(ny)
if n_elements(wvlo) le 0 then wvlo=4000.
if n_elements(wvhi) le 0 then wvhi=9000.
if not keyword_set(label) then label=' '
yy  = fltarr(n_elements(spec)-1)
if n_elements(smoothie) le 0 then yy=spec

;; is the array spec even or odd?  Check before calling rebin.
yyt=fltarr(n_elements(spec)-1)
wvt=fltarr(n_elements(wv)-1)
evenodd = n_elements(spec) mod 2

if evenodd eq 1 then begin
   if n_elements(smoothie) gt 0 then begin
      yyt(*) = spec(0:n_elements(spec)-2)
      wvt(*) = wv(0:n_elements(wv)-2)
      yy=rebin(yyt,n_elements(yyt)/2)
      wv=rebin(wvt,n_elements(wvt)/2)
      help,spec,yyt,yy,wvt,wv
   endif
endif


;if n_elements(smoothie) gt 0 then yy=smooth(spec,smoothie)
;print,size(spec)



pixb=(where(abs(wv-wvlo) eq min(abs(wv-wvlo))))[0]
pixe=(where(abs(wv-wvhi) eq min(abs(wv-wvhi))))[0]
lopix=pixb
hipix=pixe
if (pixe lt pixb) then lopix=pixe
if (pixe lt pixb) then hipix=pixb
;if n_elements(ymin) le 0 then ymin=min(yy(lopix:hipix))
if n_elements(ymin) le 0 then ymin=0
if n_elements(ymax) le 0 then ymax=max(yy(lopix:hipix))
if n_elements(z) le 0 then z=0.22
z1=z+1.
erase

plot,wv,yy,xr=[wvlo,wvhi],yr=[ymin,ymax],xstyle=1,ystyle=1,xtickname=$
;   (make_array(30,value=' ')),ytit=textoidl('F_\lambda'),$
   (make_array(30,value=' ')),ytit='Counts',$
    pos=[.1,0.35,0.95,0.9],psym=10,xtick_get = xv,$
    /noerase,/nodata
;*!d.n_colors

;oplot,wv,yy,95*!d.n_colors
oplot,wv,yy,psym=10

;;
;; draw an axis matching observed wavelengths to rest wavelengths
;; and pixel position.
;;

nlab=n_elements(xv)
rlab=fltarr(nlab)
plab=fltarr(nlab)
for i=0,nlab-1 do begin
   plab(i)=(where(abs(wv-xv(i)) eq min(abs(wv-xv(i)))))[0]
endfor 

plab=fix(plab)
rlab=fix(xv/z1)
rlab=strtrim(rlab,2)
plab=strtrim(plab,2)

axis, xaxis=1, xtickv=xv, xtickname=rlab,xticks=nlab+1,xminor=10,chars=1.2
;axis,0,ymin,xaxis=0,xtickv=xv,xtickname=plab,xticks=nlab+1,$
;   xminor=10,chars=1.2

;;
;; Plot and label spectral features.
;;

oplot,[wvlo,wvhi],[0,0]

low=0.
hi=1.
barlo=0.75
barhi=0.80
zeq='Z =            '
strput,zeq,z,3
zeq=strcompress(zeq)

plot,[1215.67*z1,1215.67*z1],[low,hi],th=1,linestyle=1,/noerase,$
xstyle=5,ystyle=4,position=[0.1,0.1,0.95,.950],xr=[wvlo,wvhi],$
yr=[0,1],/nodata,title=label
yout=!y.crange(0)+0.84*(!y.crange(1)-!y.crange(0))

if (Keyword_set(noz)) then begin
      awv=7600.
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[7600.,7600.],[low,.93],th=1,linestyle=1
      xyouts,7600.,yout,'Ab',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
endif
if (Keyword_set(noz)) then begin
      awv=6868.
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[6868.,6868],[low,.93],th=1,linestyle=1 
      xyouts,6868.,yout,'Bb',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
endif

if (NOt Keyword_Set(noz)) then begin

plot,[1215.67*z1,1215.67*z1],[low,hi],th=1,linestyle=1,/noerase,$
xstyle=5,ystyle=4,position=[0.1,0.1,0.95,0.90],xr=[wvlo,wvhi],$
yr=[0,1],/nodata

yout10=!y.crange(0)+.95*(!y.crange(1)-!y.crange(0))
xyouts,wvhi-400,yout10,zeq,size=0.85,alignment=0.5

wvlab = [912.,1026.,1215.67,$
   1260.,1302.,1335.,1394.,$
;   1403.,1428.,
   1527,$
   1608.,1671.,1855.,1863.,1909.,$
   2326.,2374.,2424,2799.,3727.,3835.,$
   3933.,3968.,4102.,4304.,$
   4861.,5167.,5173.,5184.,5896,$
   5889,5986,6563,7040,7680,8190,8520]

linlab = textoidl(['Lylim','Ly\beta','Ly\alpha',$
   'SiII','SiII','CII','SiIV',$
;   'SiIV','CIII',$
   'SiII',$
   'FeII','AlII','AlIII','.','CIII',$
   'CII','FeII','NeIV','MgII','OII','H\eta',$
   'K','H','H\delta','Gb',$
   'H\beta',' ','MgI',' ',' ',$
   'NaD','.','H\alpha','TiO','KI','Na','Cs'])

offset=50
bigoff=100
yout2=!y.crange(0)+0.90*(!y.crange(1)-!y.crange(0))
yout3=!y.crange(0)+0.78*(!y.crange(1)-!y.crange(0))

FOR i=0,n_elements(wvlab)-1 DO BEGIN
   awv = wvlab(i)*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[wvlab(i)*z1,wvlab(i)*z1],[low,hi],th=1,linestyle=1 
      xyouts,wvlab(i)*z1,yout2,linlab(i),size=1.0,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   ENDIF
ENDFOR

;; special cases here.   

      awv=1240.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[1240.*z1,1240.*z1],[low,hi],th=1,linestyle=1
      xyouts,1240.*z1,yout,'NV',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=1294.543*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'SiIII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=1296.33*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'CIII/SiIII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=1304.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[1304.*z1,1304.*z1],[low,hi],th=1,linestyle=1 
      xyouts,1304.*z1,yout,'OI',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
  endif

      awv=1323.929*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'CII/NIII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=1343.354*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'OIV',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=1403.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[1403.*z1,1403.*z1],[low,hi],th=1,linestyle=1
   endif

   
      awv=1417.237*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'SiIII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif

      awv=1427.85*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'CIII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif

      awv=1501.76*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1 
      xyouts,awv,yout,'SV',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif

      awv=1550.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[1548.*z1,1548.*z1],[low,hi],th=1,linestyle=1 
      oplot,[1552.*z1,1552.*z1],[low,hi],th=1,linestyle=1  
      xyouts,1550.*z1,yout,'CIV',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
      awv=1640.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[1640.*z1,1640.*z1],[low,hi],th=1,linestyle=1
      xyouts,1640.*z1,yout,'HeII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=2344.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[2344.*z1,2344.*z1],[low,hi],th=1,linestyle=1  
      xyouts,2344.*z1,yout,'FeII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
  endif

awv=2600.*z1
IF (awv GT wvlo AND awv LT wvhi) THEN begin
oplot,[awv,awv],[low,hi],th=1,linestyle=1  
xyouts,awv,yout,textoidl('FeII'),size=1.,$
alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
endif
awv=2587.*z1
IF (awv GT wvlo AND awv LT wvhi) THEN begin
oplot,[awv,awv],[low,hi],th=1,linestyle=1
endif

      awv=2853.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[2853.*z1,2853.*z1],[low,hi],th=1,linestyle=1  
      xyouts,2853.*z1,yout,'MgI',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
      awv=3798.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[3798.*z1,3798.*z1],[low,hi],th=1,linestyle=1  
      xyouts,3798.*z1,yout,textoidl('H\theta'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
  endif
      awv=3346.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1  
      xyouts,awv,yout,textoidl('NeV'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
  endif

      awv=3426.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1  
      xyouts,awv,yout,textoidl('NeV'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
  endif

      awv=3889.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[3889.*z1,3889.*z1],[low,hi],th=1,linestyle=1  
      xyouts,3889.*z1,yout,textoidl('H\zeta'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
      awv=3970.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[3970.*z1,3970.*z1],[low,hi],th=1,linestyle=1
      xyouts,3970.*z1,yout,textoidl('H\epsilon'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
      awv=4340.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[4340.*z1,4340.*z1],[low,hi],th=1,linestyle=1 
      xyouts,4340.*z1,yout,textoidl('H\gamma'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=4959.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[4959.*z1,4959.*z1],[low,hi],th=1,linestyle=1
   endif
      awv=4983.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      xyouts,4983.*z1,yout2,textoidl('OIII'),size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=5007.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[5007.*z1,5007.*z1],[low,hi],th=1,linestyle=1
   endif
    awv=5876.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1
      xyouts,awv,yout,'HeI',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
      awv=6875.
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[6875.,6875.],[low,hi],th=1,linestyle=1
      xyouts,6875.,yout,'Bb',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax] 
   endif
      awv=7600.
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[7600.,7600.],[low,hi],th=1,linestyle=1
      xyouts,7600.,yout,'Ab',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
     awv=6548.*z1
     awv2=6583.4*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1
      oplot,[awv2,awv2],[low,hi],th=1,linestyle=1
      xyouts,6560*z1,yout,'NII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[6716.*z1,6716.*z1],[low,hi],th=1,linestyle=1
      oplot,[6730.*z1,6730.*z1],[low,hi],th=1,linestyle=1
      xyouts,6722.*z1,yout,'SII',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
    awv=6300.*z1
   IF (awv GT wvlo AND awv LT wvhi) THEN begin
      oplot,[awv,awv],[low,hi],th=1,linestyle=1
      xyouts,awv,yout,'OI',size=1.,$
       alignment=0.5,clip=[wvlo,ymin,wvhi,ymax]
   endif
endif
;;
;; Decide upon the sky spectrum upper and lower limits.
;; Plot out the sky spectrum.

pixsky=(where(abs(wv-6300) eq min(abs(wv-6300))))[0]
pixoffsky=(where(abs(wv-6100) eq min(abs(wv-6100))))[0]
sskymax=ssky(pixsky)+5000.
;sskymax=ssky(pixsky)-10000.

;if (ymax lt 5000) then sskymin=ssky(pixoffsky)-2400.
;if (ymax lt 5000) then sskymin=ssky(pixoffsky)-2500.
;if (ymax ge 5000) then sskymin=ssky(pixoffsky)-4000.

if ((abs(6300.-wv(pixsky)) gt 20.) or (pixsky lt 20)) then sskymax=7000.
;if ((abs(6100.-wv(pixoffsky)) gt 20.) or (pixoffsky lt 20)) then sskymin=0.

;print,'Sky spectrum minimum',pixoffsky,wv(pixoffsky),ssky(pixoffsky)
;print,'Sky spectrum maximum',pixsky,wv(pixsky),ssky(pixsky)
;print,'Sky spectrum limits',sskymin,sskymax
sskymin=0
;sskymax=5.e4

plot,wv,ssky,xr=[wvlo,wvhi],yr=[sskymin,sskymax],xstyle=9,ystyle=1,$
pos=[0.1,0.1,0.95,0.35],xtit='Wavelength ('+string(197b)+')',$
/noerase,yticks=2,yminor=4,linestyle=0,/nodata, xchars=2
;*!d.n_colors

oplot,wv,ssky
;oplot,wv,ssky,0 
;oplot,wv,ssky,95*!d.n_colors

;device,/close
;ps_close
end
