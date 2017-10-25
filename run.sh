./output
avconv -i "./plots/u velocity%04d.png" -r 30 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" "homogeneous.mp4"
make
./output
avconv -i "./plots/u velocity%04d.png" -r 30 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" "heterogenous.mp4"
