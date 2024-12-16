reset
set term x11
set term pos eps enhanced color
set output 'poteam-all.eps'
set encoding iso_8859_1
#unset key
set key top right
#set multiplot

#set origin 0,0
set size 0.4,0.5
set xr [2:5] 
set yr [-0.008:0.003]
set xtics 1
set mxtics 2
set ytics 0.002 

set style line 1 lt 1 lw 1
set style line 2 lt 2 dt 2 lw 1
set style line 3 lt 3 dt 3 lw 1
set style line 4 lt 4 dt 4 lw 1
set style line 5 lt 5 dt 5 lw 1
set style line 6 lt 6 dt 6 lw 1
set style line 7 lt 7 dt 7 lw 1
set style line 8 lt 8 dt 8 lw 1
set style line 9 lt 9 dt 9 lw 1


set yl '{/Times-Italic v}({/Times-Italic r})'
set xl '{/Times-Italic r} (\305)'

pl 'BRET_Co-EAM.DAT' u 1:2 t 'Co' w l ls 1,'BRET_Cr-EAM.DAT' u 1:2 t 'Cr' w l ls 2,'BRET_Cu-EAM.DAT' u 1:2 t 'Cu' w l ls 3, 'BRET_Fe-EAM.DAT' u 1:2 t 'Fe' w l ls 4,'BRET_Mn-EAM.DAT' u 1:2 t 'Mn' w l ls 5,'BRET_Ni-EAM.DAT' u 1:2 t 'Ni' w l ls 6,'BRET_Ti-EAM.DAT' u 1:2 t 'Ti' w l ls 7,'BRET_Sc-EAM.DAT' u 1:2 t 'Sc' w l ls 8,'BRET_V-EAM.DAT' u 1:2 t 'V' w l ls 9 
