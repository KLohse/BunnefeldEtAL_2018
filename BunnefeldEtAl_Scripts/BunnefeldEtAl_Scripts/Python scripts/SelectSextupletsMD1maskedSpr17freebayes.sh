#Annotation amended 170418
#Selection script 
#Hardcoded for 6 individuals
#Produces bed2 file and vcf files (amended here for only snps)

#Takes as input bed files per individual and filtered vcfs (variants, snps, indels)
#Can be amended for different syntax and/or numbers of individuals 

##Creates subsampled variant files across SIX individuals.
#Here you create a file with contig name and the x:y positions that overlap across each individual's bed files (e.g. sites 1-10 and 21-25 are callable across all individuals)
~/../software/bedtools2-master/bin/multiIntersectBed -i <(grep -E 'CALLABLE|LOW_COVERAGE' $1$2.CallableLoci.mbq10mmq20md2.masked.bed) <(grep -E 'CALLABLE|LOW_COVERAGE' $1$3.CallableLoci.mbq10mmq20md2.masked.bed) <(grep -E 'CALLABLE|LOW_COVERAGE' $1$4.CallableLoci.mbq10mmq20md2.masked.bed) <(grep -E 'CALLABLE|LOW_COVERAGE' $1$5.CallableLoci.mbq10mmq20md2.masked.bed) <(grep -E 'CALLABLE|LOW_COVERAGE' $1$6.CallableLoci.mbq10mmq20md2.masked.bed) <(grep -E 'CALLABLE|LOW_COVERAGE' $1$7.CallableLoci.mbq10mmq20md2.masked.bed)| grep "1,2,3,4,5,6" | awk '{print $0,$3 - $2}'| sed 's/ /\t/g' | cut -f1-3,12 > $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.CallableLoci.mbq10mmq20md1.masked.CALLABLE.bed2

bcftools view -R $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.CallableLoci.mbq10mmq20md1.masked.CALLABLE.bed2 $1.freebayes.C1.b10.m20.decompose.snps.Q20.vcf.gz > $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf
~/../software/vcftools_0.1.12b/bin/vcf-sort -c $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf > $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf.tmp
mv $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf.tmp $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf
#rm $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.vcf.gz

bgzip $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf
tabix -p vcf $1$2.$1$3.$1$4.$1$5.$1$6.$1$7.calibrated_haploid.filtered_snps.CallableLoci.md1.masked.freebayes.vcf.gz

