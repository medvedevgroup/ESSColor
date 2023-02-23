a=$(cat stat_nkmer_ess)
b=$(cat stat_nkmer_jc)

sub=$(($b - $a))
if [[ $a = $b ]]
        then
                echo "validation successful";		echo "Success" > validate1
        else
		if [[ $sub -eq 1 ]]
		then
			cat jc_matrix.tsv | head -n -1 > jc_fixed ; cp jc_fixed jc_matrix.tsv
			echo "validation successful";           echo "Success" > validate1

		else 
	                echo "Fail"
		fi
        fi
