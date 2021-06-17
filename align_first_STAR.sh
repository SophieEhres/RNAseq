#!/bin/bash

trimdir="${SCRATCH}/RNAseq/trim"
index="${SCRATCH}/genomes/human/hg38/STAR_index"
alignfirstdir="${SCRATCH}/RNAseq/align_first"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/align_first"
logdir="${SCRATCH}/scripts/RNAseq/logs/align_first"


mkdir -p ${alignfirstdir}
mkdir -p ${jobdir}
mkdir -p ${logdir}

names=$(ls ${trimdir} | grep -v "unpaired" | rev | cut -d "_" -f3- | rev | uniq)


for name in ${names}; do

files=$(ls ${trimdir}/*.gz | grep -v "unpaired" | grep -e "${name}" | tr '\n' ' ')
logfile="star_alignfirst_${name}.log"

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=staralign_firstpass_${name}
#SBATCH --ntasks=12
#SBATCH --output=staralign_firstpass.out
#SBATCH --mem-per-cpu=20G
#SBATCH --output=${logdir}/${logfile}

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >> ${logdir}/${logfile}

STAR --runMode alignReads \
--genomeDir ${index} \
--runThreadN 12 \
--readFilesIn ${files} \
--outSAMtype BAM Unsorted \
--outFileNamePrefix ${alignfirstdir}/alignmentfirstpass_${name}_ \
--readFilesCommand zcat \
--outSAMattrRGline ID:"${name}"	PL:"ILLUMINA" CN:"IRCM" 2>> ${logdir}/${logfile}
" > ${jobdir}/${name}_alignfirst.sh

done

for file in $(ls ${jobdir}); do
echo "submitting align job for ${name}"
    sbatch ${jobdir}/${file}
done

