#!/bin/bash

markdupdir="${SCRATCH}/RNAseq/markdup"
alignseconddir="${SCRATCH}/RNAseq/align_second"
logdir="${SCRATCH}/scripts/RNAseq/logs/dup"
jobdir="${SCRATCH}/scripts/RNAseq/job_files/dup"
metricdir="${SCRATCH}/RNAseq/dupmetrics"
cleandir="${SCRATCH}/RNAseq/cleanbam"

mkdir -p ${logdir}
mkdir -p ${markdupdir}
mkdir -p ${jobdir}
mkdir -p ${metricdir}
mkdir -p ${cleandir}


for name in $(ls ${alignseconddir}/*.bam | xargs -n 1 basename | rev | cut -d "_" -f2- | rev | cut -d "_" -f2-); do

    file=$(ls ${alignseconddir}/*.bam | grep -e "${name}")

    quickcheck=$(samtools quickcheck ${cleandir}/${name}_clean.bam 2>&1)

    if [ -f ${cleandir}/${name}_clean.bam ] ; then
        echo "${name} already removed duplicates"
        if [ -z ${quickcheck} ] ; then
        echo -e "correctly removed duplicates \n"
        else
        echo -e "removing failed bam \n"
        rm ${cleandir}/${name}_clean.bam
        fi
    else

echo "#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --job-name=removeduplicates_${name}
#SBATCH --output=${logdir}/${name}_dup.log
#SBATCH --ntasks=20
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=sophie.ehresmann@gmail.com
#SBATCH --account=def-sauvagm

echo "$(date +"%Y-%m-%d-%H-%M")" >> ${logdir}/${name}_dup.log
module load picard java/11.0.2
if [ -f ${markdupdir}/${name}_markdup.bam ]; then
    
    echo "already marked duplicates"

else

java -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates --version >> ${logdir}/${name}_dup.log
java -Xmx20G -jar /cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/picard/2.23.3/picard.jar MarkDuplicates -I ${file} -M ${metricdir}/${name}_dupmetrics.txt \
-O ${markdupdir}/${name}_markdup.bam \
--REMOVE_DUPLICATES false

fi

samtools view -m 10G -@ 20 -F 1804 -b -h ${markdupdir}/${name}_markdup.bam > ${cleandir}/${name}_clean.bam


" > ${jobdir}/${name}_dup.sh
    
echo -e "submitting dup job for ${name} \n"

sbatch --account=def-sauvagm ${jobdir}/${name}_dup.sh

    fi

done