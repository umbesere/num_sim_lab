for i in 0 1 2 3 4 5 6 7 8 9
do
./MolDyn_NVE.exe
source restart.sh
echo $i
done
