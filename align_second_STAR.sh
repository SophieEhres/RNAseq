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

if [ -f ${alignfirstdir}/ALL_SJ.filt.tab ] ; then
    echo "splice junctions already merged, moving to second pass alignment"
else
    gawk '$6==1 || ($6==0 && ($7>2))' ${alignfirstdir}/*SJ.out.tab > ${alignfirstdir}/ALL_SJ.filt.tab
fi

for name in ${names}; do

files=$(ls ${trimdir}/*.gz | grep -v "unpaired" | grep -e "${name}" | tr '\n' ' ')
logfile="star_alignfirst_${name}.log"

echo -e "${name}"
if [ -d ${alignseconddir}/alignmentsecondpass_${name}*_STARtmp ] ; then
    echo "restarting ${name} job"
    ls -d ${alignseconddir}/* | grep -e "${name}" | xargs rm -rf
    
echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=staralign_secondpass_${name}
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=60G
#SBATCH --output=${logdir}/${logfile}

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >> ${logdir}/${logfile}

STAR --runMode alignReads \
--genomeDir ${index} \
--runThreadN 12 \
--readFilesIn ${files} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${alignseconddir}/alignmentsecondpass_${name} \
--readFilesCommand zcat \
--sjdbFileChrStartEnd ${alignfirstdir}/ALL_SJ.filt.tab \
--outSAMattrRGline ID:"${name}"	PL:"ILLUMINA" CN:"IRCM" 2>> ${logdir}/${logfile}
" > ${jobdir}/${name}_alignsecond.sh

    echo "submitting align job for ${name}"
    sbatch ${jobdir}/${name}_alignsecond.sh

elif [ -f ${alignseconddir}/alignmentsecondpass_${name}_*.bam ] ; then
    
    echo "${name} already aligned"

else

echo "#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=staralign_secondpass_${name}
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=40G
#SBATCH --output=${logdir}/${logfile}

now=$(date +"%Y-%m-%d-%H-%M")
echo "${now}" >> ${logdir}/${logfile}

STAR --runMode alignReads \
--genomeDir ${index} \
--runThreadN 12 \
--readFilesIn ${files} \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${alignseconddir}/alignmentsecondpass_${name} \
--readFilesCommand zcat \
--sjdbFileChrStartEnd ${alignfirstdir}/ALL_SJ.filt.tab \
--outSAMattrRGline ID:"${name}"	PL:"ILLUMINA" CN:"IRCM" 2>> ${logdir}/${logfile}
" > ${jobdir}/${name}_alignsecond.sh


echo "submitting align job for ${name}"
    sbatch ${jobdir}/${name}_alignsecond.sh


fi

done

