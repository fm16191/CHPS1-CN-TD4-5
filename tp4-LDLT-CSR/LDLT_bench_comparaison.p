set term svg
set size ratio 0.8
set key right bottom Right samplen 2 spacing .8 height 3 font ",10"
set logscale y
set xtics scale default 0,25
set output "LDLT_bench.svg"
set xlabel "Matrix size"
set ylabel offset 2 "Average execution time"
set title font "Helvetica Bold, 10" #offset -5
set style line 1 lt 1 linecolor rgb "#1800FF" lw 2 pt 1
set style line 2 lt 1 linecolor rgb "#00BAE6" lw 2 pt 1
set style line 3 lt 1 linecolor rgb "red" lw 0.5 pt 1
set style line 4 lt 1 linecolor rgb "#004700" lw 2 pt 1
set style line 5 lt 1 linecolor rgb "#00FF01" lw 2 pt 1


set ytics nomirror
set y2tics border scale default 0,10
set logscale y2
set y2label offset -2 "Average speedup" 

set title "Difference of execution time in logarithmic scale according to the size of the matrix"
plot "LDLT3b.dat" using 1:2 ls 1 title "myLDLT3b" with lines, "LDLT3b.dat" using 1:2:4:3 ls 3 notitle with errorbars,\
"LDLT1b.dat" using 1:2 ls 2 title "myLDLT1b" with lines, "LDLT1b.dat" using 1:2:4:3 ls 3 title "dispersion" with errorbars,\
"< paste LDLT3b.dat LDLT1b.dat" using 1:($2/$9) smooth bezier ls 4 title "average speedup" with lines axis x1y2

#set key right center Right
set ytics mirror
unset y2tics
unset y2label
set title "Difference of precision in logarithmic scale according to the size of the matrix"
set output "LDLT_bench_precision.svg"
set ylabel "Average precision"
plot "LDLT3b.dat" using 1:5 ls 4 title "myLDLT3b" with lines, "LDLT3b.dat" using 1:5:7:6 ls 3 notitle with errorbars,\
"LDLT1b.dat" using 1:5 ls 5 title "myLDLT1b" with lines, "LDLT1b.dat" using 1:5:7:6 ls 3 title "dispersion" with errorbars