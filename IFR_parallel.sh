name=mm10_liver_embryo_d11.5_PE_bio1_ENCLB441LCB_1
thread=12

min_lfp=5
max_lfp=15
step_lfp=2
min_lsh=50
max_lsh=200
step_lsh=20
method=BH
p_cutoff=0.05

awk '{print $1"\t"$2-50"\t"$3+50"\t"$4}' "step3.3_peakcall_"$name"_peaks.narrowPeak" |\
	awk '{if(($3-$2)>1000) {printf $1"\t"$2"\t"; printf "%.0f",($2+$3)/2; print "\t"$4".1"; printf $1"\t"; printf "%.0f",($2+$3)/2; print "\t"$3"\t"$4".2"} else print $0}' |\
	awk '{if(($3-$2)>1000) {printf $1"\t"$2"\t"; printf "%.0f",($2+$3)/2; print "\t"$4".1"; printf $1"\t"; printf "%.0f",($2+$3)/2; print "\t"$3"\t"$4".2"} else print $0}' > temp.peak

split -n l/$thread temp.peak
rm temp.peak

[ -f temp.txt ] && rm list
for file in `ls xa*`
do
	intersectBed -a "step4.2_insertion_site_"$name".bedGraph" -b $file -wa -wb | uniq > $file".bed"
	echo "Rscript /home/chengl/R/ATAC-seq_wellington.R "$file".bed IFR_"$file".txt" $min_lfp $max_lfp $step_lfp $min_lsh $max_lsh $step_lsh $method $p_cutoff >> list
	rm $file
done

cat list | parallel

rm xa*bed list
cat IFR*txt | sed "s/\"//g" | sort -k1,1V -k2,2n | awk '{print $1"\t"$2"\t"$3"\t""found_IFR_"NR"\t"$4"\t"".""\t"$5}' > "found_IFR_"$name".bed"
rm IFR*txt
