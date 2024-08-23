
#folders with the replicas
array=(`seq 0 0`)
#echo "${array[2]}"
rm */delta*

compute_fes() {
  cd $i
  rm fes*

  awk '{
      if ($1 ~ /^#/) {
          print $0;  # Print comment lines as is
      } else {
          $3 = -$3;  # Multiply the third column by -1
          if ($3 >= 0) {
              $3 = " " $3;  # Add a space if the result is non-negative
          }
          print $0;  # Print the modified line
      }
  }' COLVAR.* > COLVAR_reversed

  cp ~/Scripts/funnel_FES_from_Reweighting.py .
  head -n 250005 COLVAR.* > COLVAR
  head -n 250005 COLVAR_reversed > COLVAR_reversed_capped
  ./funnel_FES_from_Reweighting.py --skiprows 100000 --sigma 0.01 --bias opes.bias --colvar COLVAR --cv cyl.z --bin 200 --temp 300 --min -1.5 --max 1.5 --rfunnel 0.2 --uat 1.4 --bat 0.9 --blocks 3 --outfile fes_blocks.dat
  ./funnel_FES_from_Reweighting.py --skiprows 100000 --sigma 0.01 --bias opes.bias --colvar COLVAR_reversed_capped --cv cyl.z --bin 200 --temp 300 --min -1.5 --max 1.5 --rfunnel 0.2 --uat 1.4 --bat 0.9 --blocks 3 --outfile fes_blocks_reversed.dat
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
