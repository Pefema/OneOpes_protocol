
#folders with the replicas
array=(`seq 0 0`)
#echo "${array[2]}"
rm */delta*

compute_fes() {
  cd $i
  rm fes*

  cp ~/Scripts/funnel_FES_from_Reweighting.py .
  head -n 200005 COLVAR.* > COLVAR
  ./funnel_FES_from_Reweighting.py --skiprows 50000 --sigma 0.01 --bias opes.bias --colvar COLVAR --cv cyl.z --bin 200 --temp 300 --min 0.3 --max 1.8 --rfunnel 0.2 --uat 1.5 --bat 0.8 --blocks 3 --outfile fes_blocks.dat
}

for i in ${array[@]}
do
  compute_fes &
done
wait

#append all the deltaF in one file
paste */delta* > deleteme;

#calculate average and stdev of all replicas in time
awk '{sum = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; print NR*5000+10000, sum}' deleteme > temp1;
awk '{sum = 0; sum2 = 0; for (i = 1; i <= NF; i++) sum += $i; sum /= NF; for (i = 1; i <= NF; i++) sum2 += ($i-sum)^2; sum2 /= NF; print sum2}' deleteme > temp2;

paste temp1 temp2 > deltaFall.dat;
rm temp* deleteme 
