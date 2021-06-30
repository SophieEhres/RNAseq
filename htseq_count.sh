#!/bin/bash

cleandir="${SCRATCH}/RNAseq/cleanbam"
logdir="${SCRATCH}/scripts/RNAseq/logs/count"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/count"
countdir="${SCRATCH}/RNAseq/counts"
gtf="${SCRATCH}/genomes/human/hg38/Homo_sapiens.GRCh38.104.chr_mod.gtf"

mkdir -p ${jobdir}
mkdir -p ${logdir}
mkdir -p ${countdir}


for file in $(ls ${cleandir}/*_clean.bam | head -n 1 | xargs -n 1 basename) ; do
    name=$(echo ${file} | rev | cut -d "_" -f2- | rev)

    if [[ -f ${countdir}/${name}_countMatrix.txt ]]; then
        echo "${name} already counted"
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=countmatrix_${name}
#SBATCH --output=${logdir}/${name}_count.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=60G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

module load scipy-stack

samtools sort -@ 20 ${cleandir}/${file} > ${cleandir}/sorted_${file}
samtools index -@ 20 ${cleandir}/sorted_${file} 

htseq-count -f bam \
-s yes -m union -r pos \
${cleandir}/sorted_${file} ${gtf} > ${countdir}/${name}_countMatrix.txt


" > ${jobdir}/${name}_count.sh
    
echo "submitting dup job for ${name}"

sbatch --account=def-sauvagm ${jobdir}/${name}_count.sh

    fi

done