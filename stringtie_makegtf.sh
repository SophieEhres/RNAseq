#!/bin/bash


cleandir="${SCRATCH}/RNAseq/cleanbam"
logdir="${SCRATCH}/scripts/RNAseq/logs/denovo"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/denovo"
denovodir="${SCRATCH}/RNAseq/denovo"

annotation="${SCRATCH}/genomes/human/hg38/Homo_sapiens.GRCh38.104.chr_mod.gtf"

mkdir -p ${logdir}
mkdir -p ${jobdir}
mkdir -p ${denovodir}

for file in $(ls ${cleandir} | grep -e "sorted" ); do
    name=$(echo ${file} | rev | cut -d "_" -f3- | rev )
    
    if [ -f ${denovodir}/${name}.gtf ]; then
        echo "already made de novo transcriptome for ${name}"
    else

    echo "#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name=stringtie_${name}
#SBATCH --output=${logdir}/${name}_stringtieGtf.log
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

start=`date +%s`


stringtie -G ${annotation} -p 10 -o ${denovodir}/${name}.gtf ${cleandir}/${file}

end=`date +%s`

runtime=$((end-start))

echo -e '${runtime}' > ${logdir}/${name}_stringtieGtf.log

" > ${jobdir}/${name}_denovo.sh

    echo -e "submitting for ${name} \n"
    sbatch ${jobdir}/${name}_denovo.sh

    fi

done
