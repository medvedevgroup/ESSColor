kmer_order_file=$1
kmtricks_file=$2

#kmer_order_file="kmer_id.txt"
#kmtricks_file="jc_matrix.tsv"

# assuming kmtricks file is sorterd

sort -n -k 1 -o serial_sorted_by_kmer $kmer_order_file #At this point, both are sorted by canonical kmer, assumption
echo "1st sort done"
paste -d' ' serial_sorted_by_kmer $kmtricks_file > merged_kmtricks.txt
echo "join  done"
sort -n -k 2 -o merged_kmtricks_sorted.txt merged_kmtricks.txt
cat merged_kmtricks_sorted.txt | cut -f4 -d" " > col_bitmatrix
echo "final sort  done"


