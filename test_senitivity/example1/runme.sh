gfortran sensitivity.f90 -o sensitivity
sensitivity<<EOF
3
in_spetral.txt
out_kernel_amplitude
out_kernel_phase
5500
EOF
#parameters inputed into sensitivity are:
#velocity
#spetral file
#amplitude kernel output file
#phase kernel output file
#epicenter distance 

echo "=calculation end="
python3 draw2D out_kernel_amplitude  -o out_kernel_amplitude -c dc -x "lotitude(degree)" -y "latitude(degree)" -t "kernel(/km^-2)"
echo "=draw end="
