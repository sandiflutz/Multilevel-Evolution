set term post eps enha color 20
set output'averInvXt_N10000_Ty100_Kh500_net0_Bv0.001000_cost0.0100_sigma0.05_mu0.010000_mig0.000010_eu.eps'

set title'{/=15 Euler: L=100, K_h=500, # of types=100}'

set xl'{/=25 time}'
set yl'{/=25 Average Investment}'

#set xr[1:200000]

#set xtics 100000

plot '../averInvXt_N10000_Ty100_Kh500_net0_Bv0.001000_cost0.0100_sigma0.05_mu0.010000_mig0.000010_eu.dat' u 1:2 w l lw 3 lc 1 t''
