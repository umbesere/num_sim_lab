steps=($(seq 5 1 5))
mus=($(seq 0.75 0.005 0.85))
sigmas=($(seq 0.55 0.005 0.64))
for step in ${steps[*]}
do
    for mu in ${mus[*]}
    do
        for sigma in  ${sigmas[*]}
        do
            ./es08_2.exe $mu $sigma $step
        done
    done
done