#!/bin/bash

trimdir="${SCRATCH}/RNAseq/trim"
index="${SCRATCH}/genomes/human/hg38/STAR_index"
alignfirstdir="${SCRATCH}/RNAseq/align_first"
alignseconddir="${SCRATCH}/RNAseq/align_second"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/align_second"
logdir="${SCRATCH}/scripts/RNAseq/logs/align_second"


mkdir -p ${alignseconddir}
mkdir -p ${jobdir}
mkdir -p ${logdir}

names=$(ls ${trimdir} | grep -v "unpaired" | rev | cut -d "_" -f3- | rev | uniq)

gawk '$6==1 || ($6==0 && ($7>2))' ${alignfirstdir}/*SJ.out.tab > ${alignfirstdir}/ALL_SJ.filt.tab


for name in ${names}; do

files=$(ls ${trimdir}/*.gz | grep -v "unpaired" | grep -e "${name}" | tr '\n' ' ')
logfile="star_alignfirst_${name}.log"

if [ -f ${alignseconddir}/alignmentsecondpass_${name}_*.bam ] ; then
    
    echo "${name} already aligned"
    
elif [ -d ${aligseconddir}/alignmentsecondpass_${name}*tmp ] ; then
    echo "restarting ${name} job"
    ls -d ${alignseconddir}/* | grep -e "${name}" | xargs rm -rf

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=staralign_secondpass_${name}
#SBATCH --ntasks=12
#SBATCH --output=staralign_secondpass.out
#SBATCH --mem-per-cpu=40G
#SBATCH --output=${logdir}/${logfile}

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >> ${logdir}/${logfile}

STAR --runMode alignReads \
--genomeDir ${index} \
--runThreadN 12 \
--readFilesIn ${files} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${alignseconddir}/alignmentsecondpass_${name}_ \
--readFilesCommand zcat \
--sjdbFileChrStartEnd ${alignfirstdir}/ALL_SJ.filt.tab \
--outSAMattrRGline ID:"${name}"	PL:"ILLUMINA" CN:"IRCM" 2>> ${logdir}/${logfile}
" > ${jobdir}/${name}_alignsecond.sh

    echo "submitting align job for ${name}"
    sbatch ${jobdir}/${file}

else

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=staralign_secondpass_${name}
#SBATCH --ntasks=12
#SBATCH --output=staralign_secondpass.out
#SBATCH --mem-per-cpu=40G
#SBATCH --output=${logdir}/${logfile}

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >> ${logdir}/${logfile}

STAR --runMode alignReads \
--genomeDir ${index} \
--runThreadN 12 \
--readFilesIn ${files} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${alignseconddir}/alignmentsecondpass_${name}_ \
--readFilesCommand zcat \
--sjdbFileChrStartEnd ${alignfirstdir}/ALL_SJ.filt.tab \
--outSAMattrRGline ID:"${name}"	PL:"ILLUMINA" CN:"IRCM" 2>> ${logdir}/${logfile}
" > ${jobdir}/${name}_alignsecond.sh


echo "submitting align job for ${name}"
    sbatch ${jobdir}/${file}


fi

done

