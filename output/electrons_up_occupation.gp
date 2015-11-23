set encoding iso_8859_1
set mxtics; set mytics
unset key
set style data lines
set terminal pdf enhanced color
set output "plot/electrons_up_occupation.pdf"
set title "occupation (up)" font "Sans,24" enhanced
plot "output/electrons_up_occupation.txt"
