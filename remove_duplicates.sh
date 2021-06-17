#!/bin/bash

markdupdir="${SCRATCH}/RNAseq/markdup"
alignseconddir="${SCRATCH}/RNAseq/align_second"
logdir="${SCRATCH}/scripts/RNAseq/out_files/dup"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/dup"
metricdir="${SCRATCH}/RNAseq/dupmetrics"
cleandir="${SCRATCH}/RNAseq/cleanbam"

mkdir -p ${logdir}
mkdir -p ${markdupdir}
mkdir -p ${jobdir}
mkdir -p ${metricdir}
mkdir -p ${cleandir}


for name in $(ls ${alignseconddir}/*.bam | xargs -n 1 basename | rev | cut -d "_" -f2- | rev | cut -d "_" -f2-); do

    file=${alignseconddir}/alignmentsecondpass_${name}_*.bam

    if [ -f ${cleandir}/${name}_clean.bam ]; then
        echo "${name} already removed duplicates"
    else

echo "#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=removeduplicates_${name}
#SBATCH --output=${logdir}/${name}_dup.log
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_dup.log
module load picard java/11.0.2


java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates --version >> ${logdir}/${name}_dup.log
java -Xmx20G -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates -I ${file} -M ${metricdir}/${name}_dupmetrics.txt \
-O ${markdupdir}/${name}_markdup.bam \
--REMOVE_DUPLICATES false

samtools view -m 8G -@ 10 -F 1804 -b -h ${markdupdir}/${name}_markdup.bam > ${cleandir}/${name}_clean.bam


" > ${jobdir}/${name}_dup.sh
    
echo "submitting dup job for ${name}"

#sbatch --account=def-sauvagm ${jobdir}/${name}_dup.sh

    fi

done