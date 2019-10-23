#!/bin/bash

# # phenotype unrelated sample selection based PIHAT cutoff # #

PhenotypeList=$1 # List of phenotype positive IIDs 

PlinkGenome=$2 # plink IDB .genome file with Phi-Hat values on 10th column

PhenoPrefix=$3 # phenotype prefix to add

PIHATcutoff=0.125 # PI_HAT cutoff 

echo "$PhenotypeList $PlinkGenome $PhenoPrefix $PIHATcutoff" ;

# # get pheno specific IDB pairs and sort them highest to lowest PI_HAT values (column #10) 
fgrep -wf $PhenotypeList $PlinkGenome | awk -v OFS="\t" '{$1=$1 ; print}' | sort -k10,10rn  >${PhenoPrefix}.genome ; 

# # touch the empty intermediate file
 >${PhenoPrefix}.genome.FlipCase.tsv; 

# #
cp ${PhenoPrefix}.genome ${PhenoPrefix}.genome.tmp ; 

cat ${PhenotypeList} | while read L ; do 

             fgrep -v "$L" ${PhenoPrefix}.genome.tmp >${PhenoPrefix}.genome.tmp1 ; # iterative reduction of the file

	     # report run progress
             vnum="$(wc -l  ${PhenoPrefix}.genome.tmp1 | cut -d" " -f1 )"; 
             echo -ne "$vnum\r" ;
     
             # left align the pheno PTs (flip FID1 and FID2 if right )    
             awk -v pheno="$L" 'BEGIN{FS=OFS="\t"} { 
                                 FID="" ; IID="" ; 
                                 if ($1==pheno){print pheno, $0} else if ($3==pheno){
                                 FID=$1; IID=$2; $1=$3 ; $2=$4; $3=FID; $4=IID; print pheno, $0} }' ${PhenoPrefix}.genome.tmp >> ${PhenoPrefix}.genome.FlipCase.tsv ; 
             # move back  
             mv ${PhenoPrefix}.genome.tmp1 ${PhenoPrefix}.genome.tmp ; 
     done  

echo -ne "\n";
echo "For two sides pheno positives, select right ones and remove" ;

# For two sides pheno positives, select right ones and remove # 
    awk -F"\t" 'NR==FNR{a[$1]=$4; next} ($4 in a){print "^"$4}' ${PhenotypeList} ${PhenoPrefix}.genome.FlipCase.tsv >${PhenoPrefix}.right.pheno.txt ; 

grep -vf ${PhenoPrefix}.right.pheno.txt ${PhenoPrefix}.genome.FlipCase.tsv | tee >(cut -f1 | sort -u > ${PhenoPrefix}.Cases_selected.txt) >( awk -v pihat="$PIHATcutoff" 'BEGIN{FS=OFS="\t"} ($11>=pihat)' | cut -f4 | sort -u >${PhenoPrefix}.Case.relatives2Remove.txt) > /dev/null 


# remove all intermediate files 
# rm ${PhenoPrefix}.genome.FlipCase.tsv ${PhenoPrefix}.genome.tmp ${PhenoPrefix}.right.pheno.txt # #
