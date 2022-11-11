gfortran sensitivity.f90 -o sensitivity
sensitivity<<EOF
3
in_spetral.txt
out_kernel_amplitude
out_kernel_phase
554.963
EOF
#parameters inputed into sensitivity are:
#velocity
#spetral file
#amplitude kernel output file
#phase kernel output file
#epicenter distance 
echo "--kernel calculation end--"
python3 draw2D out_kernel_amplitude  -o out_kernel_amplitude -x "lotitude(degree)" -y "latitude(degree)" -t "kernel(/km^-2)"
python3 draw2D out_kernel_phase  -o out_kernel_phase -x "lotitude(degree)" -y "latitude(degree)" -t "kernel(/km^-2)"
echo "--draw end--"


echo "--integration begen--"
gfortran plus.f90 -o plus
plus<<EOF
in_phvel_velocity.txt
out_kernel_amplitude 
104.435   28.1069   99.1422   26.4552   554.963
out_rotation_kernel_amplitude
out_dlnA
EOF
#phase velocity inputed
#kernel inputed
#source longitude, source latitude, station longitude, station latitude, distance
#kernel outputed
#dlnA outputed
plus<<EOF
in_phvel_velocity.txt
out_kernel_phase
104.435   28.1069   99.1422   26.4552   554.963
out_rotation_kernel_phase
out_dtheta
EOF

echo "--calculation end--"
python3 draw2D in_phvel_velocity.txt -o out_phvel  -x "longitude(degree)" -y "latitude(degree)" -t "phase velocity(km/s)" -c big
python3 draw2D out_rotation_kernel_amplitude -o out_rotation_kernel_amplitude -x "longitude(degree)" -y "latitude(degree)" -t "kernel(/km^-2)" 
python3 draw2D out_rotation_kernel_phase -o out_rotation_kernel_phase -x "longitude(degree)" -y "latitude(degree)" -t "kernel(/km^-2)" 
