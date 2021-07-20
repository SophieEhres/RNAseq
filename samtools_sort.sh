#!/bin/bash

cleandir="${SCRATCH}/RNAseq/cleanbam"
logdir="${SCRATCH}/scripts/RNAseq/logs/sort"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/sort"

mkdir -p ${logdir}
mkdir -p ${jobdir}

for file in $(ls ${cleandir} | grep -v "bai" | grep -v "sorted" ) ; do

    name=$(echo ${file} | cut -d "." -f1)
    quickcheck=$(samtools quickcheck ${cleandir}/${name}_sorted.bam 2>&1)

    if [ -f ${cleandir}/${name}_sorted.bam ]; then
            
        if [ -z ${quickcheck} ] ; then
        echo -e "${name} correctly sorted \n"
        todo="FALSE"
        
        else
        echo -e "${name} removing failed bam \n"
        rm ${cleandir}/${name}_sorted.bam
        todo="TRUE"
        
        fi

    else
    todo="TRUE"

    fi

    if [ "${todo}" == "TRUE" ]; then

    echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=sort_${name}
#SBATCH --output=${logdir}/${name}_sort.log
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

samtools sort -@ 10 ${cleandir}/${file} > ${cleandir}/${name}_sorted.bam
" > ${jobdir}/${name}_sort.sh

echo -e "submitting for ${name} \n"
sbatch ${jobdir}/${name}_sort.sh
    fi

done