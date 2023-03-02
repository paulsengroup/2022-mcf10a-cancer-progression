set -e
set -u
set -x


name=$(echo "$0" | cut -f 1 -d '.')
awk -v filename="$0" -v bin_size_arg="$1" '{ bin_size = bin_size_arg; for (counter_bin1 = $2; counter_bin1 < $3; counter_bin1 += bin_size) {
    if (counter_bin1 + bin_size > $3){second_boundry = $3;} else {second_boundry = counter_bin1 + bin_size;}
    for (counter_bin2 = $5; counter_bin2 < $6; counter_bin2 += bin_size){
        if (counter_bin2 + bin_size > $6){fourth_boundry = $6;} else {fourth_boundry = counter_bin2 + bin_size;}
        print $1  "\t" counter_bin1 "\t" second_boundry "\t" $4 "\t" counter_bin2 "\t" fourth_boundry;
    }
}
}' "$0" > tmp.bedpe

cut -f 1,2,3,4,5,6 tmp.bedpe | awk '{if(!(($1 == $4) && ($2 == $5) && ($3 == $6))) print $0}' > tmp1.bedpe
cut -f 1,2,3 tmp1.bedpe > tmp1_left.bedpe
cut -f 4,5,6 tmp1.bedpe  > tmp1_right.bedpe
cat tmp1_left.bedpe tmp1_right.bedpe | sort -u | bedtools sort > tmp1_beads.bedpe
bedtools complement -L -i tmp1_beads.bedpe -g "$(sort -k 1 "$2")" | cat - tmp1_beads.bedpe | bedtools sort -g "$(sort -k 1 "$2")" > tmp1_beads_complemented1.bedpe
bedtools makewindows -b  tmp1_beads_complemented1.bedpe -w "$1" > tmp1_beads_complemented.bedpe
gunzip -c "$3" | grep acen | bedtools pairtobed -a tmp1.bedpe -b stdin -type neither > tmp1_noncen.bedpe
gunzip -c "$4" | grep acen | bedtools pairtobed -a tmp1_noncen.bedpe -b stdin -type neither > tmp1_noncen_nongap.bedpe
makeGtrack.py tmp1_noncen_nongap.bedpe tmp1_beads_complemented.bedpe > "$name""_""$1"_t.gtrack

bedtools intersect -c -a "$name""_""$1"_t.gtrack -b "$5" | \
awk '{if($8 >= 1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t1\t" $7; else  print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t.\t" $7}' \
| bedtools sort > "$name""_""$1"_LADs_t.gtrack


make_diploid_gtrack.py "$name""_""$1"_LADs_t.gtrack > "$name""_""$1"_LADs_diploid_t.gtrack
awk 'BEGIN{print "###seqid\tstart\tend\tid\tradius\tperiphery\tedges"}1' "$name""_""$1"_LADs_diploid_t.gtrack > "$name""_""$1"_LADs_diploid.gtrack


rm tmp.bedpe tmp1.bedpe tmp1_left.bedpe tmp1_right.bedpe tmp1_beads.bedpe tmp1_beads_complemented1.bedpe tmp1_beads_complemented.bedpe tmp1_noncen.bedpe \
tmp1_noncen_nongap.bedpe "$name""_""$1"_t.gtrack "$name""_""$1"_LADs_diploid_t.gtrack "$name""_""$1"_LADs_t.gtrack