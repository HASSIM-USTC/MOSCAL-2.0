gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=2 -DBSM -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o bsm_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=2 -DBOSE_QUAD -DNORMAL -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o bose_quad_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=4 -DSPARSE -DFERMI_LINEAR -Wall -Iinclude -funroll-loops main_omp_thermo.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o thermo_fermi_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=2 -DBOSE_QUAD -DNORMAL -Wall -Iinclude -funroll-loops main_omp_thermo.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o thermo_bose_quad_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=4 -DFERMI_QUAD -DSPARSE -DNORMAL -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o fermi_quad_4.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=2 -DFERMI_QUAD -DSPARSE -DNORMAL -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o fermi_quad_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=2 -DBOSE_LINEAR -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o bose_2.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=4 -DSPARSE -DFERMI_LINEAR -DNORMAL -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o fermi_4.out
# # #
gcc -march=native -O3 -std=c++17 -DSTD -DNSYS=4 -DSPARSE -DTEMPLATE -DNORMAL -Wall -Iinclude -funroll-loops main_threads_omp_equ_corr.cpp -Wl,--start-group -lstdc++ -lm -ldl -pthread -fopenmp -Wl,--end-group -o template.out
# # #
mv *.out ../
