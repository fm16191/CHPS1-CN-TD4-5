set term svg
set size ratio 0.8
set key left top Right samplen 2 spacing .8 height 3 font ",10"
#set logscale y
set xtics scale default 0,100
set output "data/CSR_bench.svg"
set xlabel "Matrix size"
set ylabel offset 2 "Average execution time"
set title font "Helvetica Bold, 10" #offset -5
set style line 1 lt 1 linecolor rgb "#0d42ba" lw 2 pt 1
set style line 2 lt 1 linecolor rgb "#0d98ba" lw 2 pt 1
set style line 3 lt 1 linecolor rgb "#0dba86" lw 2 pt 1

set style line 4 lt 1 linecolor rgb "red" lw 2 pt 1
set style line 5 lt 1 linecolor rgb "#00FF01" lw 2 pt 1

set style line 9 lt 1 linecolor rgb "red" lw 0.5 pt 1 # red dispersion

set ytics nomirror
set y2tics border scale default 0,50
#set logscale y2
set y2label offset -2 "Average speed difference"

set title "Difference of execution time in logarithmic scale according to the size of the matrix"
plot "data/CSR.dat" using 1:2 ls 1 title "reference mult" with lines, "data/CSR.dat" using 1:2:4:3 ls 3 notitle with errorbars,\
"data/CSR.dat" using 1:8 ls 2 title "mCSRv" with lines, "data/CSR.dat" using 1:8:10:9 ls 9 notitle with errorbars,\
"data/CSR.dat" using 1:5 ls 3 title "csmtCSR" with lines, "data/CSR.dat" using 1:5:7:6 ls 9 title "dispersion" with errorbars,\
"data/CSR.dat" using 1:($8/$2) smooth bezier ls 4 title "speed mCSRv / ref" with lines axis x1y2

unset logscale y
set ytics mirror
unset y2tics
unset y2label
set title "Difference of precision according to the size of the matrix"
set output "data/CSR_bench_precision.svg"
set ylabel "Average precision difference"

plot "data/CSR.dat" using 1:11 ls 5 title "precision difference" with lines, "data/CSR.dat" using 1:11:13:12 ls 9 title "dispersion" with errorbars

#i, mean_t1, min(t(1,:)), max(t(1,:)), mean_t2, min(t(2,:)), max(t(2,:)), mean_norm, min(n(1,:)), max(n(1,:))