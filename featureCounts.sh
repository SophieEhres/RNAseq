#!/bin/bash

cleandir="${SCRATCH}/RNAseq/cleanbam"
logdir="${SCRATCH}/scripts/RNAseq/logs/count"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/count"
countdir="${SCRATCH}/RNAseq/counts"
gtf="${SCRATCH}/genomes/human/hg38/Homo_sapiens.GRCh38.104.chr_mod.gtf"
#mod has added "chr" to chromosomes 


mkdir -p ${jobdir}
mkdir -p ${logdir}
mkdir -p ${countdir}


for file in $(ls ${cleandir}/*_clean.bam | head -n 1 | xargs -n 1 basename) ; do
    name=$(echo ${file} | rev | cut -d "_" -f2- | rev)

    if [[ -f ${countdir}/${name}_countMatrix.txt ]]; then
        echo "${name} already counted"
    else

echo "#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=countmatrix_${name}
#SBATCH --output=${logdir}/${name}_count.log
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

module load scipy-stack


featureCounts \
    -T 10 -t "exon" \
    -p -F GTF \
    -a ${gtf} \
    -o ${countdir}/${name}_counts.txt \
    ${cleandir}/${file}


" > ${jobdir}/${name}_count.sh
    
echo "submitting count job for ${name}"

sbatch --account=def-sauvagm ${jobdir}/${name}_count.sh

    fi

done