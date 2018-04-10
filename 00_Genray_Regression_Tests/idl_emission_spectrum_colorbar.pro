
set_plot,'ps'
device,/helvetica,/bold,/color

emis=ncdf_restore('genray.nc')

print, 'nfreq', emis.nfreq
print, 'nrayelt_emis_nc[0]', emis.nrayelt_emis_nc[0]
print, 'jx_kin', emis.jx_kin

j1=0
;j1=1
j2=emis.nfreq-1

print, 'j1' ,j1
print, 'j2' ,j2

 plot, emis.kinetic_energy_kev[0:emis.jx_kin-1], $ 
 title='kinetic energy grid [KeV] ',$
 xtitle='number of grid points ', $
 xstyle=1

;-------------------------------------------------------------------
; plot contours of the emission flux at the plasma edge
; (kinetic_energy, frequency)
;---------------------------------------------------------------------

;BH   Clip max plotted intensity by factor clip_max [=1. for no clipping]
clip_max=1.0

;---black contours

  max_intensity=0.0 
  for ifreq=0,emis.nfreq-1 do begin
    for j=0,j2 do begin
     if ( emis.wi_0_x_nc[j,ifreq] gt max_intensity) then $
        max_intensity = emis.wi_0_x_nc[j,ifreq]
    endfor
  endfor

  max_intensity=clip_max*max_intensity

  title_frq=string('emission flux spectrum at the plasma edge', format='(a)')
    
;---------color contours
 
 levels=fltarr(10)
  for l=0,9 do begin
    levels[l]=(l+1)*(max_intensity/10)
  endfor

  print, 'levels',levels

  c_labels = [0,1,0,1,0,1,0,1,0,1]

  nlevels =10

  ncolors= nlevels+1
  bottom=1

;  c_levels=[0.,levels,max_intensity]
  c_levels=[0.,levels ]
  print, '  c_levels',  c_levels

  c_labels=[0, replicate(1,nlevels),0]
  print, 'c_labels', c_labels
  c_colors= indgen(ncolors)+bottom
  print, 'c_colors',c_colors
  loadct, 33, ncolors = ncolors, bottom =bottom
 
  contour, emis.wi_0_x_nc[0:emis.jx_kin-1,0:emis.nfreq-1], $
  emis.kinetic_energy_kev[*], emis.wfreq_nc[*],$
  xrange=[0.,400.],yrange=[100.,500.],$
  levels=c_levels, c_colors=c_colors, /fill,$
  title= title_frq,$
  ytitle='frequency [GHz]', xstyle=1, $
  xtitle='kinetic energy [KeV]', ystyle=1

  contour, emis.wi_0_x_nc[0:emis.jx_kin-1,0:emis.nfreq-1], $
  emis.kinetic_energy_kev[*], emis.wfreq_nc[*],$
  levels=c_levels, c_labels=c_labels, /overplot

  c_labels = [0,1,0,1,0,1,0,1,0,1]
  c_labels=[0, replicate(1,nlevels),0]
  print, 'c_labels', c_labels
  c_colors= indgen(ncolors)+bottom
  print, 'c_colors',c_colors
  loadct, 33, ncolors = ncolors, bottom =bottom
 
!P.Multi = [0,2,0,0,1]

  contour, emis.wi_0_x_nc[0:emis.jx_kin-1,0:emis.nfreq-1], $
  emis.kinetic_energy_kev[*], emis.wfreq_nc[*],$
  xrange=[0.,400.],yrange=[100.,500.],$
  levels=c_levels, c_colors=c_colors, /fill,$
  title= title_frq,$
  ytitle='frequency [GHz]', xstyle=1, $
  xtitle='kinetic energy [KeV]', ystyle=1,$
  pos=[0.1,0.1,.7,.9]

  print, 'max_intensity', max_intensity
 
  LoadCT, 33  
  colorbar_fanning, Color=FSC_Color('firebrick'),$
  range=[0.,max_intensity],$
  divisions=11,$
  posision =[0.88, 0,10, 0.95, 0.90],/Vertical, $
  Format ='(e10.2)', Title='flux values', Right=right

print, 'normal end of program' 
ps_close

end
