//
//  external.hpp
//  hla-mapper
//
//  Created by Erick Castelli on 20/02/20.
//  Copyright © 2020 GeMBio.Unesp. All rights reserved.
//

#ifndef external_h
#define external_h

#include <iostream>
#include <cstdio>
#include <map>
#include <unordered_map>

using namespace std;

extern string ostype;
extern string Program_name;
extern string Program_company;
extern string Program_version;
extern string Program_author;
extern string Program_date;
extern string Program_website;
extern string Program_update_url;


extern map <string,int> samtools_versions;
extern map <string,int> bwa_versions;
extern map <string,int> whats_versions;
extern map <string,int> blast_versions;
extern map <string,int> freebayes_versions;
extern map <string,int> star_versions;

extern unsigned long v_memory;

extern int screen_size;
extern string v_message;
extern vector<string> warnings;
extern string v_system_out;
extern vector<string> v_system_out_list;
extern string v_command;;
extern vector <string> selected_reads;
extern int motif_size;
extern int v_buffer;
extern string v_genes;
extern vector<string> v_gene_list;
extern int v_mm_max; //-n BWA ALN
extern int v_mm_open; //-o BWA ALN
extern int v_size;
extern double v_tolerance;
extern vector <string> top_allele_pairs;
extern map <double,string> allele_pair_score;

extern std::map <string,int> position_hg38;
extern std::map <string,string> chr_hg38;
extern std::map <string,string> chr_size_hg38;
extern std::map <string,int> ref_size;
extern std::map <string,int> gene_opt_start;
extern std::map <string,int> gene_opt_end;
extern map <string,string> gene_type;


extern double read_mean_size;
extern map <string,double> read_count_gene;

extern string v_threads;
extern string v_output;
extern string v_db;
extern string v_dbprepare;
extern string v_r0;
extern string v_r1;
extern string v_r2;
extern string v_bam;
extern string v_mapout;
extern string v_sample;
extern string homedir;
extern string v_bwa;
extern string v_samtools;
extern string v_bedtools;
//extern string v_gatk;
extern string v_whats;
extern string v_blast;
//extern string v_java;
extern string v_freebayes;
extern string v_star;
extern string v_freebayes;
extern string v_bcftools;
extern double v_mtrim_error;
extern string v_bed;
extern int v_rnaseq;
extern int v_lowmem;
extern int v_callint;
extern int v_forceindex;

extern int v_debug;
extern int v_verbose;
extern int v_bypass_version_check;


extern map <pair<string,string>, int> sequence_list_r1;
extern map <pair<string,string>, int> sequence_list_r2;
extern map <string, int> sequence_size_r1;
extern map <string, int> sequence_size_r2;

extern map <pair<string,string>,int> reads_mapped_to_allele_A;
extern map <pair<string,string>,int> reads_mapped_to_allele_B;
extern map <pair<string,string>,int> reads_possibly_mapped;
extern map <string,int> reads_mapped_nm;
extern map <string,string> reads_mapped_gene;


extern map <pair<string,string>, int> sequence_address;

extern int v_skiptyping;
extern int v_typeall;
extern int v_assembling;
extern int v_intergenic;
extern int v_quiet;
extern int v_useconfig;
extern int keep_sam;
extern int minselect;
extern int v_skip_unmapped;
extern string v_map_type;
extern int v_use_local;
extern int downsampling;
extern int v_type_after_map;

extern int usechr;
extern int v_exome;
extern int v_multiphasing;
extern string v_target;


#endif /* external_h */
