# Second-order perturbation
Author: [Hitoshi Mori](mailto:h.mori.pro.xyz.will.inge@gmail.com)


## How to use sparse-ir-fortran
Please refer to [README of sparse-ir-fortran](https://github.com/SpM-lab/sparse-ir-fortran) for how to link the library to your program and a list of functions. If you have not downloaded `sparse-ir-fortran` modules, please download the project from the repository [sparse-ir-fortran](https://github.com/SpM-lab/sparse-ir-fortran), and install some required Python modules and the FFTW3 library for Fortran compiler before starting the tutorial.


## Building the tutorial code
1. Generate a data file and a fortran source file for $\Lambda = 10^5$ and $\epsilon = 10^{-7}$ on the sparse-ir-fortran directory.
```bash
cd /somewhere/sparse-ir-fortran/
python3 dump.py 1e+5 1e-7 ir_nlambda5_ndigit7.dat
python3 mk_preset.py --nlambda 5 --ndigit 7 > sparse_ir_preset.f90
```

2. Download the source code {download}`2nd_order_pert_theory.f90 <2nd_order_pert_theory.f90>` and then copy all the f90 files and the data file from the `sparse-ir-fortran/` directory to a working directory.
```bash
> cd /your/preferred/dir/
> cp /somewhere/sparse-ir-fortran/*f90 ./
> cp /somewhere/sparse-ir-fortran/*dat ./
```

3. Create Makefile and edit it. The example is shown here.
```
FC = gfortran

LDLIBS = -lblas -llapack
LDFLAGS = -L/usr/lib/x86_64-linux-gnu
F90FLAGS = -std=f2003 -Wall -fbounds-check -g -fcheck=array-temps,bounds,do,mem,pointer,recursion
IFLAGS = -I/usr/include
FFT_LIBS = -lfftw3

OBJS = sparse_ir.o sparse_ir_io.o sparse_ir_preset.o

.SUFFIXES: .f90

%.o: %.f90
	$(COMPILE.f) $(OUTPUT_OPTION) $< $(F90FLAGS) $(IFLAGS)

%.mod: %.f90 %.o
	@:

second_order_perturbation_fort.o: sparse_ir.mod sparse_ir_io.mod sparse_ir_preset.mod

.PHONY: second_order_perturbation_fort
second_order_perturbation_fort: 2nd_order_pert_theory.o $(OBJS)
	$(LINK.f) $^ $(FFT_LIBS) $(LDLIBS) $(LDFLAGS) -o $@.x

.PHONY: tutorial
tutorial: second_order_perturbation_fort

.PHONY: clean
clean:
	rm -f *.o *.mod
```

4. You can build `second_order_perturbation_fort.x` as follows:
```bash
> make tutorial
```

## Running the code "`second_order_perturbation_fort.x`"
Before running "`second_order_perturbation_fort.x`", please create the input file `input.in` as follows:
```fortran
&input
    lambda = 1.0D5,
    beta = 1.0D3,
    eps = 1.0D-7,
    nk_lin = 256,
    hubbardu = 2.0D0,
    use_preset = .true. ! .false. if you want to pull data from ir_nlambda*_ndigit**.dat
/
```

Run the program. You can obtain the following output:
```bash
> ./second_order_perturbation_fort.x < input.in
 cond (tau):    334.877581832775
 cond (matsu):    248.851825144482
 Shape of k1: (         256         256 )
 Shape of k2: (         256         256 )
 Shape of ek_: (         256         256 )
 Shape of kp: (           2       65536 )
 Shape of ek: (       65536 )
```
The first two rows show the condition numbers of SVD of the IR-basis sets for imaginary time and Matsubara frequency. When obtaining the IR-basis functions for imaginary time, sampling points of imaginary time is given as $\tau \in [0,\beta]$. In the tutorial of Sec. 8, the Python script does not use both endpoints as sampling points, but `sparse-ir-fortran` always takes them into the sampling points. That is the reason why the condition number for imaginary time, `cond (tau)`, is different from the one obtained by the Python script.
<!-- When the smallest singular value is quite small compared with the largest one, it follows that solving this fitting problem using SVD is not a good way. This is because it corresponds to the fact that the pseudo-inverse matrix $\Sigma^+$ obtained by SVD can be sensitive to numerical errors.
To judge whether a given fitting problem (that is, a given matrix) is suitable to be solved using SVD, one can check the ratio of the largest eigenvalue to the smallest one, which is a quantity known as the "condition number". When employing double precision, the problem can be considered to be stably solved using SVD as long as it does not reach a rather large value such as $10^{14}$.-->

We can plot the results by running `gnuplot` as 
```
gnuplot script.plt
```
where `script.plt` is:
```gnuplot
set terminal pdfcairo color enhanced font "Times,20" size 5,3.5
set colorsequence classic
set lmargin 10
set rmargin 5
set bmargin 3
set samples 200

set xrange[-580:580]
set xtics 200
set format x "% g"
set xlabel "ν" offset 0, 0.5
set yrange[-0.14:0.14]
set ytics 0.05
set format y "% g"
set ylabel "Im G(iν, k)" offset 1, 0
set nologscale
set key right top
set out "second_order_g_k_freq_fort.pdf"
plot "./second_order_gkf_gamma.txt" u 2:4 w lp lw 2 lt 1 pt 2 lc 3 title "k = Γ"

set xrange[-4:74]
set xtics 20
set format x "% g"
set xlabel "l" offset 0, 0.5
set yrange[1E-9:2E0]
set ytics 1E2
set format y "%.0E"
set ylabel "|G(l, k)|" offset 1, 0
set logscale y
set key right top
set out "second_order_abs-g_k_l_fort.pdf"
plot "./second_order_gkl_gamma.txt" u 1:5 w l lw 2 lt 1 lc 3 title "k = Γ", \
     "./second_order_gkl_gamma.txt" u 1:2 w l lw 2 lt 1 dt (10,10) lc rgb "orange" title "Singular values"

set xrange[-50:1050]
set xtics 200
set format x "% g"
set xlabel "τ" offset 0, 0.5
set yrange[-1.05:0.05]
set ytics 0.2
set format y "% g"
set ylabel "Re G(τ, k)" offset 1, 0
set nologscale
set key right bottom
set out "second_order_g_k_tau_fort.pdf"
plot "./second_order_gkt_gamma.txt" u 1:2 w lp lw 2 lt 1 pt 2 lc 3 title "k = Γ", \
     "./second_order_gkt_m.txt" u 1:2 w lp lw 2 lt 1 pt 1 lc rgb "orange" title "k = M"

set xrange[-50:1050]
set xtics 200
set format x "% g"
set xlabel "τ" offset 0, 0.5
set yrange[-0.51:0.01]
set ytics 0.1
set format y "% g"
set ylabel "Re G(τ, r)" offset 1, 0
set nologscale
set key center bottom
set out "second_order_g_r_tau_fort.pdf"
plot "./second_order_grt_r_0.txt" u 1:2 w l lw 2 lt 1 lc 3 title "r = (0, 0)"

set xrange[-50:1050]
set xtics 200
set format x "% g"
set xlabel "τ" offset 0, 0.5
set yrange[-0.51:0.01]
set ytics 0.1
set format y "% g"
set ylabel "Re Σ(τ, r)" offset 1, 0
set nologscale
set key center bottom
set out "second_order_sigma_r_tau_fort.pdf"
plot "./second_order_srt_r_0.txt" u 1:2 w lp lw 2 lt 1 pt 2 lc 3 title "r = (0, 0)"

set xrange[-4:74]
set xtics 20
set format x "% g"
set xlabel "l" offset 0, 0.5
set yrange[1E-8:1E0]
set ytics 1E2
set format y "%.0E"
set ylabel "|Σ(l, r)|" offset 1, 0
set logscale y
set key right top
set out "second_order_abs-sigma_r_l_fort.pdf"
plot "./second_order_srl_r_0.txt" u 1:4 w l lw 2 lt 1 lc 3 title "r = (0, 0)"

set xrange[-4:74]
set xtics 20
set format x "% g"
set xlabel "l" offset 0, 0.5
set yrange[1E-10:5E0]
set ytics 1E-9, 1E+2, 1E-1
set format y "%.0E"
set ylabel "|Σ(l, r)|" offset 1, 0
set logscale y
set key right top
set out "second_order_sigma_k_l_fort.pdf"
plot "./second_order_skl_gamma.txt" u 1:5 w l lw 2 lt 1 lc 3 title "k = Γ", \
     "./second_order_skl_gamma.txt" u 1:2 w l lw 2 lt 1 dt (10,10) lc rgb "orange" title "Singular values"

set xrange[-580:580]
set xtics 200
set format x "% g"
set xlabel "ν" offset 0, 0.5
set yrange[-0.14:0.14]
set ytics 0.05
set format y "% g"
set ylabel "Im Σ(iν, k)" offset 1, 0
set nologscale
set key right top
set out "second_order_sigma_k_freq_fort.pdf"
plot "./second_order_skf_gamma.txt" u 2:4 w lp lw 2 lt 1 pt 2 lc 3 title "k = Γ"
```

Note that the `sparse-ir-fortran` interface can evaluate functions only on sampling frequency points from a set of expansion coefficients. This tutorial does not evaluate the self-energy on a finer mesh of frequency.