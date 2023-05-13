configfile: "config.yaml"
SAMPLES = config["SAMPLES"]
EXTENSION=config["EXTENSION"]

import os

def dump_list(lst, ext, filename):
    with open(filename, 'w') as f:
        for line in lst:
            f.write(f"{line}{ext}\n")

def dump_list_with_pref(lst, ext, filename, pref):
    with open(filename, 'w') as f:
        for line in lst:
            f.write(f"{pref}{line}{ext}\n")


def write_to_file(lst, ext, listname, pre="" ):
    f = open(listname,"w")
    f.write((ext+"\n"+pre).join(lst)+ext)
    f.close()
    return pre+(ext+"\n").join(lst)+ext

def concat_list_kmers(lst):
    return ".kmers.gz ".join(lst)+".kmers.gz"

def generate_list_of_filenames(lst, ext):
    return  "esscs/"+(ext+"\nesscs/").join(lst)+ext


def num_colors(allsets):
    return len(allsets)
def read_val(filename):
    with open(filename) as f:
        lines = f.readlines()
        return lines[0].rstrip()


def get_ext_folder(EXTENSION):
    #"fungi/mers31/s5"
    replaced_path = EXTENSION.replace(".", "" )
    return os.path.abspath("../../")+"/"+replaced_path

def get_ext_folder_level0(EXTENSION):
    #"fungi/mers31/s5"
    replaced_path = EXTENSION.replace(".", "" )
    return os.path.abspath("../")+"/"+replaced_path    

def get_filesize(fname):
    file = open(fname)
    file.seek(0, os.SEEK_END)
    filebytes = file.tell()
    file.close()
    return int(filebytes)

def size_stat():
    coless_kb=(get_filesize("mega.essc")+get_filesize("rrrbv_1_delta.sdsl")+get_filesize("rrrbv_1_skip.sdsl")+get_filesize("rrrbv_1.sdsl"+get_filesize("rrr_bv_mapping.sdsl")))/1024.0
    print("megaessc size (KB)=", get_filesize("mega.essc")/1024.0)
    print("compressed bit size (KB)=", (get_filesize("rrr_bv_mapping.sdsl")+get_filesize("rrrbv.sdsl"))/1024.0)

    tot_indi=0
    for s in SAMPLES:
        tot_indi += get_filesize("esscs/"+s+EXTENSION+".essc")
    tot_indi /= 1024.0
    print("col-ess size (KB)=", coless_kb)
    print("indi-ess size (KB)=", tot_indi)
    f = open("stat_size", "w")
    f.write(str(coless_kb)+" "+str(tot_indi))
    f.close()
    

if config['ess'] == 'tip':    
    rule all:
        input:
            "mega.essc",
            "mega.essd",
            "ess_boundary_bit.txt",
            "col_bitmatrix",
            "uniq_ms.txt",
            "stat_m",
            "stat_nkmer_ess",
            #"bb_map",
            "esscolor.tar.gz",
            "size_esscolor_mb_tip"
else:
    rule all:
        input:
            # "ggout.fa.mfc",
            # "ggout.colors.dat.gz",
            # "ggout.stats.log.gz"
            # #expand("esscs/{sample}"+EXTENSION+".essd",sample=SAMPLES),
            # # expand("{sample}.kmers.gz",sample=SAMPLES),
            "mega.essc",
            "mega.essd",
            # # "rrrbv_1_delta.sdsl",
            # # "rrrbv_1.sdsl",
            # # "rrrbv_1_skip.sdsl",
            # # "rrr_bv_mapping.sdsl",
            # # "stat_size",
            "ess_boundary_bit.txt",
            "validate1",
            "col_bitmatrix",
            "uniq_ms.txt",
            "stat_m",
            "stat_nkmer_ess",
            #"bb_map",
            "esscolor.tar.gz",
            "size_esscolor_mb"


if config['option'] == 'from_essc':
    rule essc_to_essd:
        input:
            "esscs/{sample}"+EXTENSION+".essc"
        output:
            "esscs/{sample}"+EXTENSION+".essd"
        benchmark:
            "benchmarks/{sample}.essc_to_essd.txt"
        shell:
            "essDecompress {input} -f; "
elif config['option'] == 'create_essc':
    rule get_essc:
        input: 
            expand(get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION,sample=SAMPLES)
        output:
            get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION+".essc"
        params:
            k=config["k"],
            fol=get_ext_folder_level0(EXTENSION)+"/",
            ext=EXTENSION,
            ab=config["ab"]
        benchmark:
            "benchmarks/{sample}.f_to_essc.txt"
        shell:
            "essCompress -k {params.k} -i {params.fol}/{wildcards.sample}{params.ext} -a {params.ab}; "
    rule mv_essc:
        input: 
            get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION+".essc"
        output:
            "esscs/{sample}"+EXTENSION+".essc",
        params:
            fol=get_ext_folder_level0(EXTENSION)+"/",
            ext=EXTENSION,
        benchmark:
            "benchmarks/{sample}.mv_essc.txt"
        shell:
            "mkdir -p esscs;  mv {params.fol}{wildcards.sample}{params.ext}.essc esscs/"

    rule essc_to_essd2:
        input:
            "esscs/{sample}"+EXTENSION+".essc"
        output:
            "esscs/{sample}"+EXTENSION+".essd"
        benchmark:
            "benchmarks/{sample}.essc_to_essd.txt"
        shell:
            "essDecompress {input} -f"
elif config['option'] == '':
    if (config['unitig'] != 'ggcat'):
        rule f_to_megaessc:
            input:
                expand(get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION,sample=SAMPLES)
            params:
                k=config["k"],
                m=config["mem"],
                l=dump_list_with_pref(SAMPLES, EXTENSION, "mega", get_ext_folder_level0(EXTENSION)+"/")
            benchmark:
                "benchmarks/mega_from_f.essc.txt"
            output:
                "mega.essc"
            shell:
                "essCompress -i mega -k{params.k}"


if (config['option'] == 'from_essc' or config['option'] == 'create_essc' or config['option'] == 'from_essd'):
    rule essd_to_kmc:
        input:
            "esscs/{sample}"+EXTENSION+".essd"
        params:
            k=config["k"],
            m=config["mem"]
        benchmark:
            "benchmarks/{sample}.essd_to_kmc.txt"
        output:
            #"{sample}.kmc.kmc_pre",
            #"{sample}.kmc.kmc_suf"
            temp("{sample}.kmc.kmc_pre"),
            temp("{sample}.kmc.kmc_suf")
        shell:
            #"mkdir -p kmc_tmp_dir_{wildcards.sample}; kmc -k{params.k} -m{params.m} -ci1 -fa {input} {wildcards.sample}.kmc kmc_tmp_dir_{wildcards.sample}/; rm -rf kmc_tmp_dir_{wildcards.sample}"
            "mkdir -p kmc_tmp_dir_{wildcards.sample}; kmc -k{params.k} -m{params.m} -ci1 -fa {input} {wildcards.sample}.kmc kmc_tmp_dir_{wildcards.sample}/; rm -rf kmc_tmp_dir_{wildcards.sample}"

    rule essd_to_megaessc:
        input:
            essd=expand("esscs/{sample}"+EXTENSION+".essd", sample=SAMPLES)  
        params:
            k=config["k"],
            m=config["mem"],
            l=generate_list_of_filenames(SAMPLES, EXTENSION+".essd")
        benchmark:
            "benchmarks/mega.essc.txt"
        output:
            temp("mega.essc")
        shell:
            "echo \"{params.l}\"> mega;essCompress -i mega -k{params.k}"
else:
    rule f_to_kmc:
        input: 
            expand(get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION,sample=SAMPLES)
        output:
            # "{sample}.kmc.kmc_pre",
            # "{sample}.kmc.kmc_suf",
            temp("{sample}.kmc.kmc_pre"),
            temp("{sample}.kmc.kmc_suf")
        params:
            k=config["k"],
            fol=get_ext_folder_level0(EXTENSION)+"/",
            ext=EXTENSION,
            ab=config["ab"],
            m=config["mem"]
        benchmark:
            "benchmarks/{sample}.f_to_kmc.txt"
        shell:
            #"mkdir -p kmc_tmp_dir_{wildcards.sample}; kmc -k{params.k} -m{params.m} -ci1 -fa {input} {wildcards.sample}.kmc kmc_tmp_dir_{wildcards.sample}/; rm -rf kmc_tmp_dir_{wildcards.sample}"
            "mkdir -p kmc_tmp_dir_{wildcards.sample}; kmc -k{params.k} -m{params.m} -ci{params.ab} -fm {params.fol}/{wildcards.sample}{params.ext} {wildcards.sample}.kmc kmc_tmp_dir_{wildcards.sample}/; rm -rf kmc_tmp_dir_{wildcards.sample};"
      
# rule essd_to_jf:
#     input:
#         "esscs/{sample}"+EXTENSION+".essd"
#     output:
#         "{sample}.kmers.gz"
#     params:
#         k=config["k"]
#     shell:
#         "jellyfish count -C -m {params.k} -s 100000 -o {wildcards.sample}.jf <(cat {input}); jellyfish dump -c {wildcards.sample}.jf | sort -k 1 > {wildcards.sample}.kmers; gzip {wildcards.sample}.kmers" 



              
if config['matrix_generator'] == 'jc':
    rule kmc_to_kmcascii:
        input:
            "{sample}.kmc.kmc_pre",
            "{sample}.kmc.kmc_suf"
        benchmark:
            "benchmarks/{sample}.kmc_to_kmcascii.txt"
        output:
            "{sample}.kmers.gz"
        shell:
            "kmc_tools transform {wildcards.sample}.kmc dump -s {wildcards.sample}.kmers ;   gzip {wildcards.sample}.kmers"
            #"kmc_tools transform {wildcards.sample}.kmc dump {wildcards.sample}.ukmers; sort -k 1 {wildcards.sample}.ukmers -o {wildcards.sample}.kmers;  gzip {wildcards.sample}.kmers"
    rule kmc_to_joincounts:
        input:
            kmers=expand("{sample}.kmers.gz", sample=SAMPLES)
        output:
            "jc_matrix.tsv"
        params:
            l=concat_list_kmers(SAMPLES),
            l2=dump_list(SAMPLES, EXTENSION, "meta.txt"),
        benchmark:
            "benchmarks/joincounts.txt"
        shell:
            "joinCounts {params.l} > jc_matrix_unfixed ; cat jc_matrix_unfixed | head -n -1 > {output}; rm jc_matrix_unfixed "
    rule kmerlist_to_joinedlist:
        input:
            "ess_kmer_id.txt",
            "jc_matrix.tsv"
        output:
            "col_bitmatrix"
        benchmark:
            "benchmarks/kmerlist_to_joinedlist.txt"
        shell:
            "bash ess_color_sorter_final.sh ess_kmer_id.txt jc_matrix.tsv"
    
    rule megaessd_to_kmerlist:
        input:
            "mega.essd"
        output:
            "ess_kmer_id.txt",
            "ess_boundary_bit.txt"
        params:
            k=config["k"]
        benchmark:
            "benchmarks/kmerlist.txt"
        shell:
            "ess_color_kmerlist -k {params.k} -u mega.essd"
if config['matrix_generator'] == 'genmatrix':
    rule kmc_to_genmatrix:
        input:
            expand("{sample}.kmc.kmc_pre", sample=SAMPLES),
            expand("{sample}.kmc.kmc_suf", sample=SAMPLES),
            "mega.essd",
        output:
            temp("col_bitmatrix"),
            temp("ess_boundary_bit.txt"),
            temp("stat_nkmer_ess")
        params:
            l=dump_list(SAMPLES, ".kmc", "list_kmc"),
            l2=dump_list(SAMPLES, EXTENSION, "meta.txt"),
            k=config["k"],
        benchmark:
            "benchmarks/genmatrix.txt"
        shell:
            "/usr/bin/time  -f \"%M %e %U %S\" --output-file=kb_sec_genmatrix ~/s/proj4/git/ESSColor/bin/genmatrix  -c list_kmc -o {output} -s -l mega.essd -k {params.k}"



if (config['unitig'] != 'ggcat'):
    rule megaessc_to_megaessd:
        input:
            "mega.essc"
        output:
            "mega.essd"
        benchmark:
            "benchmarks/mega.essd.txt"
        shell:
            "essDecompress mega.essc"







if config['matrix_generator'] == 'jc':
    rule stat_nkmer_ess:
        input:
            "ess_kmer_id.txt"
        output:
            "stat_nkmer_ess"
        benchmark:
            "benchmarks/stat_nkmer_ess.txt"
        shell:
            "cat ess_kmer_id.txt | wc -l > stat_nkmer_ess"  
    rule stat_nkmer_jc:
        input:
            "jc_matrix.tsv"
        output:
            "stat_nkmer_jc"
        benchmark:
            "benchmarks/stat_nkmer_jc.txt"
        shell:
            "cat jc_matrix.tsv | wc -l > stat_nkmer_jc"    
            #"cat jc_matrix.tsv | head -n -1 > jc_fixed ; mv jc_fixed jc_matrix.tsv; cat jc_matrix.tsv | wc -l > stat_nkmer_jc"    

    rule validate_jc_kmers:
        input:
            "stat_nkmer_ess",
            "stat_nkmer_jc"
        output:
            "validate1"
        shell:
            "bash ess_color_validate_kmer.sh stat_nkmer_ess stat_nkmer_jc"

if config['matrix_generator'] == 'genmatrix':
    rule stat_nkmer_genmatrix:
        input:
            "col_bitmatrix"
        output:
            temp("stat_nkmer_genmatrix")
        benchmark:
            "benchmarks/stat_nkmer_genmatrix.txt"
        shell:
            "cat col_bitmatrix | wc -l > stat_nkmer_genmatrix"    
       
    rule validate_genmatrix_kmers:
        input:
            "stat_nkmer_ess",
            "stat_nkmer_genmatrix"
        output:
            "validate1"
        shell:
            "cmp --silent stat_nkmer_ess stat_nkmer_genmatrix && echo '### SUCCESS: Files Are Identical! ###' > validate1 || echo \"files are different\""
            

rule stat_uniq_colclass:
    input:
        "validate1",
        "col_bitmatrix"
    output:
        temp("uniq_ms.txt"),
        temp("stat_m")
    benchmark:
        "benchmarks/stat_uniq_colclass.txt"
    shell:
        "cat col_bitmatrix | sort -k 1 -n | uniq > temp_uniq; mv temp_uniq uniq_ms.txt;  cat uniq_ms.txt | wc -l > stat_m"

rule compress:
    input:
        "ess_boundary_bit.txt",
        "validate1",
        "col_bitmatrix",
        "uniq_ms.txt",
        "stat_m",
        "stat_nkmer_ess"
    params:
        c=len(SAMPLES),
    benchmark:
        "benchmarks/compress.txt"
    output:
        #"bb_map",
        temp("rrr_main"),
        temp("rrr_local_table"),
        temp("rrr_map_hd"),
        temp("rrr_map_hd_boundary"),
        temp("frequency_sorted")
    shell:
        "/usr/bin/time  -f \"%M %e %U %S\" --output-file=kb_sec_colcompress essColorAuxMatrixCompress -i uniq_ms.txt -d col_bitmatrix -c {params.c} -m $(cat stat_m) -k $(cat stat_nkmer_ess) -s ess_boundary_bit.txt -x 16"

rule zip_compress:
    input:
        "rrr_main",
        "rrr_local_table",
        "rrr_map_hd",
        "rrr_map_hd_boundary",
        "frequency_sorted",
        "mega.essc",
        "meta.txt"
    output:
        temp("esscolor"),
        "esscolor.tar.gz"
    benchmark:
        "benchmarks/final_gzip.txt"
    shell: 
        "mkdir -p esscolor; gzip -v9 meta.txt; gzip -v9 frequency_sorted; mv frequency_sorted.gz rrr_main rrr_local_table rrr_map_hd rrr_map_hd_boundary mega.essc meta.txt.gz esscolor/; tar cf esscolor.tar esscolor/;  gzip -v9 esscolor.tar; "

    
if config['ess'] == 'tip':
    rule zip_compress_size_tip:
        input: 
            "esscolor.tar.gz"
        output:
            "size_esscolor_mb_tip"
        shell:
            "ls -l | grep esscolor.tar.gz | awk '{{print $5/1024.0/1024.0}}' >  size_esscolor_mb_tip;"
else:
    rule zip_compress_size:
        input: 
            "esscolor.tar.gz"
        output:
            "size_esscolor_mb"
        shell:
            "ls -l | grep esscolor.tar.gz | awk '{{print $5/1024.0/1024.0}}' >  size_esscolor_mb;"
            #" nkmer=$(cat stat_nkmer_ess); ls -l | grep ess_color.tar.gz | awk -v nk=$nkmer '{{print $5*8.0/$nk}}' > size_esscolor_bitskmer"  
# rule all_stat:
#     input: 
#         "rrrbv_1_delta.sdsl",
#         "rrrbv_1.sdsl",
#         "rrrbv_1_skip.sdsl",
#         "rrr_bv_mapping.sdsl"
#     output:
#         "stat_size"
#     run:
#         size_stat()

# rule f_to_ggcat_unitig:
#     input:
#         expand(get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION,sample=SAMPLES)
#     params:
#         k=config["k"],
#         m=config["mem"],
#         ab=config["ab"],
#         l=dump_list_with_pref(SAMPLES, EXTENSION, "list_fa", get_ext_folder_level0(EXTENSION)+"/")
#     benchmark:
#         "benchmarks/f_to_ggcat.txt"
#     output:
#         "ggout.fa.mfc",
#         "ggout.colors.dat.gz",
#         "ggout.stats.log.gz"
#         kmers.esstip 
#     shell:
#         "/usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_ggcat_build ggcat build -k {params.k} -j 8 -l list_fa -o gg_unitigs.fa -s{params.ab} -p -e; essAuxMFCompressC ggout.fa; rm ggout.fa; gzip ggout.colors.dat; gzip ggout.stats.log; "
#         #ls -ls gg* | awk \'{sum = sum + $6} END {print sum/1024.0/1024.0}\' > stat_mb_gg "
#         /usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_tip  essAuxCompress -k {params.k} -i gg_unitigs.fa -t 1; kmers.esstip 

if config['unitig'] == 'ggcat':
    rule f_to_ggcatess:
        input:
            expand(get_ext_folder_level0(EXTENSION)+"/{sample}"+EXTENSION,sample=SAMPLES)
        params:
            k=config["k"],
            m=config["mem"],
            ab=config["ab"],
            l=dump_list_with_pref(SAMPLES, EXTENSION, "list_fa", get_ext_folder_level0(EXTENSION)+"/")
        benchmark:
            "benchmarks/f_to_ggcatess.txt"
        output:
            temp("gg_unitigs.fa"),
            temp("gg_unitigs.stats.log"),
            temp("kb_sec_ggcat_build")
        shell:
            "/usr/bin/time  -f \"%M %e %U %S\" --output-file=kb_sec_ggcat_build ggcat build -k {params.k} -j 8 -l list_fa -o gg_unitigs.fa -s{params.ab} -p -e"
    
if config['ess'] == 'tip':
    rule ggcat_unitig_to_ess_tip:
        input:
            "gg_unitigs.fa" 
        params:
            k=config["k"],
            m=config["mem"],
        benchmark:
            "benchmarks/ggcat_unitig_to_ess_tip.txt"
        output:
            temp("mega.essc"),
            temp("mega.essd"),
            temp("kb_sec_essauxc_tip"),
            temp("kb_sec_essauxd_tip"),
            temp("kb_sec_mfc_ess_tip")
        shell:
            "/usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_essauxc_tip  essColorAuxKmersCompress -k {params.k} -i {input} -t 1; /usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_essauxd_tip  essColorAuxKmersDecompress -i kmers.esstip 1; mv kmers.esstip.spss mega.essd; /usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_mfc_ess_tip   essAuxMFCompressC kmers.esstip; mv kmers.esstip.mfc mega.essc"

else:
    rule ggcat_unitig_to_ess:
        input:
            "gg_unitigs.fa" 
        params:
            k=config["k"],
            m=config["mem"],
        benchmark:
            "benchmarks/ggcat_unitig_to_ess.txt"
        output:
            temp("mega.essc"),
            temp("mega.essd"),
            temp("kb_sec_essauxc"),
            temp("kb_sec_essauxd"),
            temp("kb_sec_mfc_ess")
        shell:
            "/usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_essauxc  essColorAuxKmersCompress -k {params.k} -i {input} -t 0; /usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_essauxd  essColorAuxKmersDecompress -i kmers.ess 1; mv kmers.ess.spss mega.essd; /usr/bin/time  -f \"%M\t%e\" --output-file=kb_sec_mfc_ess   essAuxMFCompressC kmers.ess; mv kmers.ess.mfc mega.essc"
