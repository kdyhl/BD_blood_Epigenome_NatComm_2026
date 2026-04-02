#!/bin/bash
#$ -l os=RedHat7

source ~/.bashrc
source ~/.my.bashrc
export LD_LIBRARY_PATH=$LIBRARY_PATH:/broad/uge/8.4.3/lib/lx-amd64 



start=`date +%s`


#argument
TASK=${SGE_TASK_ID}


# file names
FILE_MOTIF="~/lhou.compbio/data/Epimap/motifs/motif_forHomer/collated_motifs_homerFormat.pval_1e-4.10_$(($TASK-1)).txt"
dir_ou=${dir_ou}/motif_$TASK


if [ $MOD == "rerun" ] && [ -f ${dir_ou}/knownResults.html ]; then
	echo "[STATUS] File already exists at: ${dir_ou}/knownResults.html"
else
	echo "[STATUS] Running matches for $FILE_MOTIF"
	if [ $bg == "randomBG" ]; then
		cmd="findMotifsGenome.pl $fg $genome $dir_ou -nomotif -size given -mknown $FILE_MOTIF"
	else
		cmd="findMotifsGenome.pl $fg $genome $dir_ou -nomotif -size given -bg $bg -mknown $FILE_MOTIF"
	fi
	echo "$cmd"
	bash -c "time $cmd"

	end=`date +%s`
	runtime=$((end-start))
	echo "Finished motif matches ($TASK) sucessfully in $runtime seconds."
fi
    


