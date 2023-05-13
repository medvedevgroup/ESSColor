
# rule all:
#     input:
#         "esscolor/mega.essd",
#         "esscolor/meta.txt",
#         "matrix.txt"
        
# rule decompress_color:
#     input:
#         "esscolor/mega.essd",
#         "esscolor/meta.txt",
#     output:
#         "matrix.txt"
#     shell:
#         "essColorAuxMatrixDecompress -i esscolor/; mv dec_ess_color ../matrix.txt; mv mega.essd ../simplitigs.fa; cd ../ "

# rule ess_kmers_decompress:
#     input:
#         "esscolor/mega.essc"
#     output:
#         "esscolor/mega.essd"
#     shell:
#         "cd esscolor; essAuxMFCompressD mega.essc; essColorAuxKmersDecompress -i mega.essc.d 1; mv kmers.ess mega.essd; cd ../"

# rule decompress_tar:
#     input:
#         "esscolor.tar.gz"
#     output:
#         "meta.txt",
#         "frequency_sorted",
#         "mega.essc",
#         "rrr_main",
#         "rrr_local_table",
#         "rrr_map_hd",
#         "rrr_map_hd_boundary"
#     shell:
#         "tar -xf esscolor.tar.gz; cd esscolor; gunzip meta.txt.gz; gunzip frequency_sorted.gz; cd"
   