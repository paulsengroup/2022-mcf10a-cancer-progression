#!/usr/bin/env bash

set -e
set -u
set -x
set -o pipefail

name=$(echo "$1" | cut -f 1 -d '.')

# Save path to current working dir
wd="$(readlink -f "$PWD")"

# Create a temporary folder
tmpdir="$(mktemp -d)"

# This will remove the temporary folder and cd back to the original folder upon exit (even in case of errors)
# shellcheck disable=SC2064
trap "rm -rf '$tmpdir'; cd '$PWD'" EXIT

# Copy input files to the temp folder
cp "$1" "${@:3:5}" "$tmpdir/"

cd "$tmpdir"
awk -v filename="$1" -v bin_size_arg="$2" '{ bin_size = bin_size_arg; for (counter_bin1 = $2; counter_bin1 < $3; counter_bin1 += bin_size) {
    if (counter_bin1 + bin_size > $3){second_boundry = $3;} else {second_boundry = counter_bin1 + bin_size;}
    for (counter_bin2 = $5; counter_bin2 < $6; counter_bin2 += bin_size){
        if (counter_bin2 + bin_size > $6){fourth_boundry = $6;} else {fourth_boundry = counter_bin2 + bin_size;}
        print $1 "\t" counter_bin1 "\t" second_boundry "\t" $4 "\t" counter_bin2 "\t" fourth_boundry;
    }
}
}' "$1" > tmp.bedpe

cut -f 1,2,3,4,5,6 tmp.bedpe | awk '{if(!(($1 == $4) && ($2 == $5) && ($3 == $6))) print $0}' > tmp1.bedpe
cut -f 1,2,3 tmp1.bedpe > tmp1_left.bedpe
cut -f 4,5,6 tmp1.bedpe > tmp1_right.bedpe
cat tmp1_left.bedpe tmp1_right.bedpe | sort -u | bedtools sort -g "$3" > tmp1_beads.bedpe
bedtools complement -L -i tmp1_beads.bedpe -g "$3" | cat - tmp1_beads.bedpe | bedtools sort -g "$3" > tmp1_beads_complemented1.bedpe
bedtools makewindows -b tmp1_beads_complemented1.bedpe -w "$2" > tmp1_beads_complemented.bedpe

# I don't think this is working as intended, grepping for acen will discard all transolcations and assembly gaps
#zcat "$4" "$5" | grep acen | bedtools pairtobed -a tmp1.bedpe -b stdin -type neither > tmp1_noncen_nongap.bedpe
#cat "$7" | grep acen | bedtools pairtobed -a tmp1_noncen_nongap.bedpe -b stdin -type neither > tmp1_noncen_nongap_.bedpe

makeGtrack.py tmp1.bedpe tmp1_beads_complemented.bedpe > "$name""_""$2"_t.gtrack

# TODO filter using the unifed blacklist ("$5")

zcat "$4" |
bedtools intersect -c -a "$name""_""$2"_t.gtrack -b stdin | \
awk '{if($8 >= 1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t1\t" $7; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t.\t" $7}' \
| bedtools sort > "$name""_""$2"_LADs_t.gtrack

make_diploid_gtrack.py "$name""_""$2"_LADs_t.gtrack > "$name""_""$2"_LADs_diploid_t.gtrack
awk 'BEGIN{print "###seqid\tstart\tend\tid\tradius\tperiphery\tedges"}1' "$name""_""$2"_LADs_diploid_t.gtrack > "$name""_""$2"_LADs_diploid.gtrack

# TODO instead of generating unique file names based on the inputs, consider printing the final output to stdout
# and let the caller decide what the output should be called
mv ./*diploid.gtrack "$wd/"
