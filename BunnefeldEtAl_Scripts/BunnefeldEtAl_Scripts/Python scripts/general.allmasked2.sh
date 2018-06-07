#!/bin/bash
flist="$2 $3 $4 $5 $6 $7"
for fvalue in $flist
do
echo $1$fvalue
./../../SpeciesVariants/Remove_Masked_From_Individuals.sh $1$fvalue $8
done
