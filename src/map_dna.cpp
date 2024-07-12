//
//  map_dna.cpp
//  hla-mapper
//
//  Created by Erick Castelli on 20/02/20.
//  Copyright © 2020 GeMBio.Unesp. All rights reserved.
//
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <string>
#include <thread>
#include <unordered_map>
#include <cstring>
#include <cassert>
#include <future>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <mutex>

#include "map_dna.hpp"
#include "external.hpp"
#include "functions.hpp"
#include "ThreadPool.hpp"
#include "preselect.hpp"
#include "typing.hpp"

mutex mtx_hg38_select;
mutex mtx_selection_general;
mutex mtxg;

vector <string> selected_results_r1;
vector <string> selected_results_r2;
vector <string> selected_results_trim_r1;
vector <string> selected_results_trim_r2;
//unordered_map <string,int> motifs_master;


void read_sam_dna (int frame, string sam, string gene)
{
    ifstream input( sam.c_str() );
    for( std::string line; getline( input, line ); )
    {
        if (line.substr(0,1) == "@") {continue;}
        vector<string> sam_line;
        boost::split(sam_line,line,boost::is_any_of("\t"));
  
        int mm = 10000;
        if ((sam_line[2] == "*") || (sam_line[5] == "*")) {continue;}
        
        vector<string> nm_value;
        boost::split(nm_value,sam_line[12],boost::is_any_of(":"));
        mm = stoi(nm_value[2]);
        

        string cigar = splitcigar(sam_line[5]);
        vector<string> cigardata;
        boost::split(cigardata,cigar,boost::is_any_of(","));
        string c1 = cigardata[0];
        string c2 = cigardata[cigardata.size()-1];
        
        
        
        if ((c1 != "") && (mm != 10000))
        {
            if ((c1.substr(c1.size()-1,1) == "S") || (c1.substr(c1.size()-1,1) == "H"))
            {mm = mm + stoi(c1.substr(0,c1.size()));}
        }
        if ((c2 != "") && (mm != 10000))
        {
            if ((c2.substr(c2.size()-1,1) == "S") || (c2.substr(c2.size()-1,1) == "H"))
            {mm = mm + stoi(c2.substr(0,c2.size()));}
        }
        
        if ((sam_line[0].substr(sam_line[0].size()-2,2) == "/1") || (sam_line[0].substr(sam_line[0].size()-2,2) == "/2")){sam_line[0] = sam_line[0].substr(0,sam_line[0].size()-2);}

        
        pair <string,string> key = make_pair(gene,sam_line[0]);
        
        if (frame == 1) {
            if (sequence_list_r1.find(key) != sequence_list_r1.end())
            {
                if (sequence_list_r1[key] > mm + 1)
                {
                    sequence_list_r1[key] = mm + 1;
//                    sequence_size_r1[sam_line[0]] = sam_line[9].size();
                }
            }
            else {
                sequence_list_r1[key] = mm + 1;
 //               sequence_size_r1[sam_line[0]] = sam_line[9].size();
            }
        }
        
        if (frame == 2) {
            if (sequence_list_r2.find(key) != sequence_list_r2.end())
            {
                if (sequence_list_r2[key] > mm +1)
                {
                    sequence_list_r2[key] = mm +1;
   //                 sequence_size_r2[sam_line[0]] = sam_line[9].size();
                }
            }
            else {
                sequence_list_r2[key] = mm +1;
  //              sequence_size_r2[sam_line[0]] = sam_line[9].size();
            }
        }
        
    }
    input.close();
    removefile(sam,v_debug);
}








void main_dna_map ()
{

    int v_check = 0;
    if (v_r0 != "") {v_check = 1;}
    if ((v_r1 != "") && (v_r2 != "")) {v_check = 1;}
    if (v_bam != "") {v_check = 1;}
 
    if (((v_db == "") || (v_check == 0)) || (v_sample == ""))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::dna";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     hla-mapper dna r1=R1.gz r2=R2.gz sample=your_sample_name <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Usage:     hla-mapper dna bam=your.bam sample=your_sample_name <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "sample      sample name or identification", 1, v_quiet);
        screen_message (screen_size, 2, "r0          a single-ended fastq (fq, fastq, or gz, ignore r1/r2/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "r1          a paired-ended forward fastq (fq, fastq, or gz, ignore r0/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "r2          a paired-ended reverse fastq (fq, fastq, or gz, ignore r0/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "bam         a BAM/CRAM file (ignore r0/r1/r2)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "db          path to an hla-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "output      output folder", 1, v_quiet);
        screen_message (screen_size, 2, "threads     number of threads [default: " + v_threads + "]", 1, v_quiet);
 //       screen_message (screen_size, 2, "bed         BED file to override the one in the database", 1, v_quiet);
        screen_message (screen_size, 2, "buffer      number of sequences in buffer [default: " + to_string(v_buffer) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "error       threshold for nucleotide quality trimming [default: " + to_string(v_mtrim_error).substr(0,4) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "tolerance   fraction of mismatches allowed [Default: " + to_string(v_tolerance).substr(0,4) + ", from 0.01 to 0.10]", 1, v_quiet);
        screen_message (screen_size, 2, "downsample  downsampling for the adjustment procedure [default: " + to_string(downsampling) + "]", 1, v_quiet);

 /*
        if (v_bwa != "") {v_message = "[found at " + v_bwa + "]";}
        if (v_bwa == "") {v_message = "(!!! bwa not detected !!!)";}
        v_message = "bwa         path to BWA " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        
        
        if (v_samtools != "") {v_message = "[found at " + v_samtools + "]";}
        if (v_samtools == "") {v_message = "[!!!Samtools not detected!!!]";}
        v_message = "samtools    path to Samtools " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

   */
        screen_message (screen_size, 2, "", 1, v_quiet);
//        screen_message (screen_size, 2, "--use-local-db        also consider the local database of sequences", 1, v_quiet);
        screen_message (screen_size, 2, "--type            call the typing algorithm after mapping", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-unmapped   skip retrieving unmapped reads from BAM/CRAM", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-adjust     skip the adjustment procedure [not recommended]", 1, v_quiet);
//        screen_message (screen_size, 2, "--type-all            search for the best match (alleles) for all genes [very slow]", 1, v_quiet);
//        screen_message (screen_size, 2, "--assemble            perform assembling for typed genes (under development)", 1, v_quiet);
//        screen_message (screen_size, 2, "--intergenic          map intergenic sequences using BWA", 1, v_quiet);
//        screen_message (screen_size, 2, "--keep-original-sam   keep the unoptimezed SAM file", 1, v_quiet);
        screen_message (screen_size, 2, "--low-mem         force low memory mode for sequence selection", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet           quiet mode", 1, v_quiet);
//        screen_message (screen_size, 2, "--noconfig            ignore the pre-configuration file", 1, v_quiet);

        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    if (v_sample == "") {
        v_r1 = "";
        v_r0 = "";
        v_r2 = "";
        v_bam = "";
        warnings.push_back("You must indicate a sample name or id (sample=)");
        main_dna_map();
        return;
    }
    
    if ((v_r1 != "") && (v_r1 == v_r2)) {warnings.push_back("r1 and r2 must be different fastq files"); main_dna_map(); return;}

    if (v_bam != "") {
        ifstream file_bam(v_bam.c_str());
        if (!file_bam)
        {
            warnings.push_back ("Error accessing the .bam file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            v_r2 = "";
            v_bam = "";
            main_dna_map();
            return;
        }
        file_bam.close();
    }

    if (v_r1 != "") {
        ifstream file_r1(v_r1.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Error accessing the r1 file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            main_dna_map();
            return;
        }
        file_r1.close();
    }

    if (v_r0 != "") {
        ifstream file_r1(v_r0.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Error accessing the r0 file.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            main_dna_map();
            return;
        }
        file_r1.close();
    }
    
    
    if (v_r2 != "") {
        ifstream file_r2(v_r2.c_str());
        if (!file_r2)
        {
            warnings.push_back ("Error accessing the r2 file.");
            v_sample = "";
            v_r2 = "";
            main_dna_map();
            return;
        }
        file_r2.close();
    }
    
    
    if (((ends_with(v_r1,".fq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fq")) && (ends_with(v_r1,".gz"))))
    {
        
        warnings.push_back ("hla-mapper cannot deal with .fastq and .gz at the same time.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    if (((ends_with(v_r1,".fastq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fastq")) && (ends_with(v_r1,".gz"))))
    {
        
        warnings.push_back ("hla-mapper cannot deal with .fastq and .gz at the same time.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    if (((v_r1 != "") && (v_r2 == "")) || ((v_r1 == "") && (v_r2 != "")))
    {
        warnings.push_back ("You must indicate r1 and r2 when using paired-end sequencing data.");
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    // checking database
    string v_db_info = v_db + "/db_dna.info";
    boost::replace_all(v_db_info, "\\ ", " ");

    if (! fileExists(v_db_info))
    {
        v_message = "Error accessing database " + v_db_info;
        cout << v_message << endl;
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }


    ifstream file_db(v_db_info.c_str());
    int db_version_ok = 0;
    for( std::string line; getline( file_db, line ); )
    {
        if (line == "hla-mapper:4.3+") {
            db_version_ok = 1;
            continue;
        }
        vector<string> db_data;
        boost::split(db_data,line,boost::is_any_of(":"));
        if (db_data[0] == "GENE")
        {
            v_genes = v_genes + db_data[1] + ",";
            chr_hg38[db_data[1]] = db_data[2];
            position_hg38[db_data[1]] = stoi(db_data[3]);
            chr_size_hg38[db_data[1]] = db_data[4];
            gene_opt_start[db_data[1]] = stoi(db_data[5]);
            gene_opt_end[db_data[1]] = stoi(db_data[6]);
            gene_type[db_data[1]] = db_data[7];
        }

 /*
        if (db_data[0] == "INTERGENIC")
        {
            if (db_data[1] == "FALSE") {v_intergenic = 0;}
        }
        if (db_data[0] == "TYPEALL")
        {
            if (db_data[1] == "FALSE") {v_typeall = 0;}
        }
        if (db_data[0] == "ASSEMBLE")
        {
            if (db_data[1] == "FALSE") {v_assembling = 0;}
        }
*/
        if (db_data[0] == "SIZE")
        {
            v_size = stoi(db_data[1]);
        }
    }
    file_db.close();
    
    
    if (db_version_ok == 0) {
        v_message = "The database is not compatible with this version of hla-mapper.";
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_dna_map();
        return;
    }
    
    
    if (v_bam != "")
    {
        if (v_bed != "")
        {
            ifstream file_bed(v_bed.c_str());
            if (!file_bed)
            {
                warnings.push_back ("Error accessing the BED file.");
                v_sample = "";
                v_r1 = "";
                v_bam = "";
                v_r2 = "";
                main_dna_map();
                return;
            }
            file_bed.close();
        }
    }
    
    
    
    
       v_rnaseq = 0;
       v_callint = 1;
       
    /*
       if (v_intergenic == 1)
       {
           string file = v_db + "/reference/dna/hg38/hg38_non_alt.fa.bwt";
           if (! fileExists(file))
           {
               warnings.push_back ("You need to run preparedb before using --intergenic");
               v_sample = "";
               v_r1 = "";v_r2 = "";
               v_bam = "";
               v_r0 = "";
               main_dna_map();
               return;
           }
       }
    */
    
    
    // criando output folder
    if (v_output == "")
    {
        if (v_bam != "") {v_output = findfilepath(v_bam) + "/hla-mapper/";}
        if (v_r0 != "") {v_output = findfilepath(v_r0) + "/hla-mapper/";}
        if (v_r1 != "") {v_output = findfilepath(v_r1) + "/hla-mapper/";}
    }
 
    
 
    
    v_command = "mkdir " + v_output;
    v_system_out = GetStdoutFromCommand(v_command);
    v_command = "mkdir " + v_output + "/log/";
    v_system_out = GetStdoutFromCommand(v_command);
    
    string v_sample_sub = v_sample.substr(0, v_sample.size()-1);
    
    
    string v_general_log = v_output + v_sample_sub + ".hla-mapper.log";
    ofstream general_log;
    general_log.open (v_general_log.c_str());
    
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = Program_name + "::dna";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log << v_message << endl;
 
    v_message = "Version " + Program_version;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log << v_message << endl;
   

    v_message = "Number of threads: " + v_threads;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log << v_message << endl;
        
    
    v_message = "Sample: " + v_sample_sub;
    general_log << v_message << endl;
    v_message = "Database: " + v_db;
    general_log << v_message << endl;
    v_message = "Output: " + v_output;
    general_log << v_message << endl;
    
    v_message = "Files (R0): " + v_r0;
    general_log << v_message << endl;
    v_message = "Files (R1): " + v_r1;
    general_log << v_message << endl;
    v_message = "Files (R2): " + v_r2;
    general_log << v_message << endl;
    v_message = "Files (BAM/CRAM): " + v_bam;
    general_log << v_message << endl;

        
    v_message = "Trimming error: " + to_string(v_mtrim_error);
    general_log << v_message << endl;
    v_message = "Tolerance: " + to_string(v_tolerance);
    general_log << v_message << endl;
    v_message = "Minimum read size: " + to_string(v_size);
    general_log << v_message << endl;

    if (v_bam != "")
    {
        if (v_skip_unmapped == 1) {
            v_message = "Retrieving unmapped: no";
        }
        if (v_skip_unmapped == 0) {
            v_message = "Retrieving unmapped: yes";
        }
        general_log << v_message << endl;
    }

    if (v_skiptyping == 1) {
        v_message = "Performing adjustments: no";
    }
    if (v_skiptyping == 0) {
        v_message = "Performing adjustments: yes";
    }
    general_log << v_message << endl;
    
    /*
    if (v_intergenic == 1) {v_message = "Map intergenic reads: yes";}
    if (v_intergenic == 0) {v_message = "Map intergenic reads: no";}
    general_log << v_message << endl;
*/

    
    main_preselect();
    if (v_callint == 2) {return;}
    
 
    v_message = "Type: " + v_map_type;
    general_log << v_message << endl;
    
    

    
    screen_message (screen_size, 0, "Sorting sequences ...", 2, v_quiet);

    int loop = stoi(v_threads);
    boost::algorithm::erase_all(v_genes, " ");
    boost::split(v_gene_list,v_genes,boost::is_any_of(","));

    string selectcleanR1 = v_output + v_sample + "selected.trim.R1.fastq";
    string selectcleanR2 = v_output + v_sample + "selected.trim.R2.fastq";
    if (v_map_type == "single") {selectcleanR1 = v_output + v_sample + "selected.trim.R0.fastq";}
    
    selected_results_r1.clear();
    selected_results_r2.clear();
    selected_results_trim_r1.clear();
    selected_results_trim_r2.clear();
    
    map <string,int> r1_trim_id_map;
    unordered_map <string,string> r1_trim_seq_map;
    unordered_map <string,string> r2_trim_seq_map;
    unordered_map <string,string> r1_trim_qual_map;
    unordered_map <string,string> r2_trim_qual_map;

    ifstream r1sort (selectcleanR1.c_str());
    ifstream r2sort (selectcleanR2.c_str());
    
    for( std::string line; getline( r1sort, line ); )
    {
        string idr1 = line;
        getline( r1sort, line );
        string seqr1 = line;
        getline( r1sort, line );
        string infor1 = line;
        getline( r1sort, line );
        string qualr1 = line;

        getline( r2sort, line );
        string idr2 = line;
        getline( r2sort, line );
        string seqr2 = line;
        getline( r2sort, line );
        string infor2 = line;
        getline( r2sort, line );
        string qualr2 = line;

        r1_trim_id_map[idr1] = 1;
        r1_trim_seq_map[idr1] = seqr1;
        r2_trim_seq_map[idr1] = seqr2;
        r1_trim_qual_map[idr1] = qualr1;
        r2_trim_qual_map[idr1] = qualr2;
    }
    r1sort.close();
    r2sort.close();
    
    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i ==  "") {continue;}
        screen_message (screen_size, 0, "Sorting sequences for " + *i + " ...", 2, v_quiet);
        
        unordered_map <string,int> motifs;

        string motifdbfile = v_db + "/motif/dna/loci/" + *i + ".txt";
        boost::replace_all(motifdbfile, "\\ ", " ");

        ifstream motifdb (motifdbfile.c_str());
        for( std::string line; getline( motifdb, line ); )
        {
            if (line == "") {continue;}
            motifs[line] = 1;
        }
        motifdb.close();
        
        if (v_use_local == 1)
        {
            string motifdbfile = "'" + v_db + "/local/dna/" + *i + "/" + *i + ".txt'";
            ifstream motifdb (motifdbfile.c_str());
            for( std::string line; getline( motifdb, line ); )
            {
                if (line == "") {continue;}
                motifs[line] = 1;
            }
            motifdb.close();
        }
        vector <string> select_data_sorted_r1;
        vector <string> select_data_sorted_r2;

        if (v_map_type == "paired") {
 
            ThreadPool pool(loop);
            std::vector< std::future<int> > results;
              
            for( auto & item : r1_trim_id_map)
              {
                  string id = item.first;
                  string r1seq = r1_trim_seq_map[id];
                  string r2seq = r2_trim_seq_map[id];
                  string r1qual = r1_trim_qual_map[id];
                  string r2qual = r2_trim_qual_map[id];

                  results.emplace_back(
                                         pool.enqueue([id, r1seq, r2seq, r1qual, r2qual, &motifs, &select_data_sorted_r1, &select_data_sorted_r2]
                          {
                          if (r1seq.size() < motif_size) {return 1;}
                          if (r2seq.size() < motif_size) {return 1;}
                          for (int c = 0; c < (r1seq.size() - motif_size); c++)
                          {
                              string sub = r1seq.substr(c,motif_size);
                              if ( motifs.find(sub) != motifs.end() ) {
                                  mtx_selection_general.lock();
                                  select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
                                  select_data_sorted_r2.push_back(id + "\n" + r2seq + "\n+\n" + r2qual);
                                  mtx_selection_general.unlock();
                                  return 1;
                              }
                            }
                              return 1;
                          })
                          );
              }
              for(auto && result: results){result.get();} // waiting for all threads
            
            ofstream OUTselectR1;
            string selectR1 = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
            OUTselectR1.open (selectR1.c_str());

            ofstream OUTselectR2;
            string selectR2 = v_output + v_sample + "sorted_" + *i + "_R2.fastq";
            OUTselectR2.open (selectR2.c_str());
            
            for (int a = 0; a < select_data_sorted_r1.size(); a++)
            {
                  OUTselectR1 << select_data_sorted_r1[a] << endl;
                  OUTselectR2 << select_data_sorted_r2[a] << endl;
            }
            OUTselectR1.close();
            OUTselectR2.close();
            select_data_sorted_r1.clear();
            select_data_sorted_r2.clear();

            screen_message (screen_size, 0, "Sorting sequences for " + *i + " ... done", 2, v_quiet);
            
        }
        
 
           if (v_map_type == "single") {
    
               ThreadPool pool(loop);
               std::vector< std::future<int> > results;
                 
               for( auto & item : r1_trim_id_map)
                 {
                     string id = item.first;
                     string r1seq = r1_trim_seq_map[id];
                     string r1qual = r1_trim_qual_map[id];
 
                     results.emplace_back(
                                            pool.enqueue([id, r1seq, r1qual, &motifs, &select_data_sorted_r1]
                             {
                             if (r1seq.size() < motif_size) {return 1;}
                             for (int c = 0; c < (r1seq.size() - motif_size); c++)
                             {
                                 string sub = r1seq.substr(c,motif_size);
                                 if ( motifs.find(sub) != motifs.end() ) {
                                     mtx_selection_general.lock();
                                     select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
                                     mtx_selection_general.unlock();
                                     return 1;
                                 }
                               }
                                 return 1;
                             })
                             );
                 }
                 for(auto && result: results){result.get();} // waiting for all threads
               
               ofstream OUTselectR1;
               string selectR1 = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
               OUTselectR1.open (selectR1.c_str());

               for (int a = 0; a < select_data_sorted_r1.size(); a++)
               {
                     OUTselectR1 << select_data_sorted_r1[a] << endl;
               }
               OUTselectR1.close();
               select_data_sorted_r1.clear();
               select_data_sorted_r2.clear();

               screen_message (screen_size, 0, "Sorting sequences for " + *i + " ... done", 2, v_quiet);
               
           }
        
    }

    screen_message (screen_size, 0, "Sorting sequences: done", 1, v_quiet);
    
    
    
    

    
    
    
 
    
    
    
    
    
    
    
    
    
    
    screen_message (screen_size, 0, "Scoring sequences ...", 2, v_quiet);
    
    boost::algorithm::erase_all(v_genes, " ");
    boost::split(v_gene_list,v_genes,boost::is_any_of(","));

    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i ==  "") {continue;}
        screen_message (screen_size, 0, "Scoring sequences for " + *i + " ...", 2, v_quiet);
        
       string v_ref = "'" + v_db + "/mapper/dna/" + *i + "/" + *i + ".fas' ";
       boost::replace_all(v_ref, "\\ ", " ");
       string v_sam1 = v_output + *i + ".mapper.1.sam";
       string v_sam2 = v_output + *i + ".mapper.2.sam";
           
       string v_r0_sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
       string v_r1_sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
       string v_r2_sort = v_output + v_sample + "sorted_" + *i + "_R2.fastq";
        
       string v_sai1 = v_output + "tmp1.sai ";
       string v_sai2 = v_output + "tmp2.sai ";


       string v_log = v_output + "/log/" + *i + "_mapper_aln.1.log";
       v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r1_sort + " > " + v_sai1 + " 2>" + v_log;
       
        if (v_map_type == "single") {v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r0_sort + " > " + v_sai1 + " 2>" + v_log;}
       system (v_command.c_str());
       
       v_log = v_output + "/log/" + *i + "_mapper_sam.1.log";
       v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r1_sort + " > " + v_sam1 + " 2>" + v_log;
       
       if (v_map_type == "single") {v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r0_sort + " > " + v_sam1 + " 2>" + v_log;}
       system (v_command.c_str());
       v_command = " rm " + v_sai1;
       system (v_command.c_str());

        if (v_map_type == "paired") {
            v_log = v_output + "/log/" + *i + "_mapper_aln.2.log";
            v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r2_sort + " > " + v_sai2 + " 2>" + v_log;
            system (v_command.c_str());
            v_log = v_output + "/log/" + *i + "_mapper_sam.2.log";
            v_command = v_bwa + " samse " + v_ref + v_sai2 + v_r2_sort + " > " + v_sam2 + " 2>" + v_log;
            system (v_command.c_str());
            v_command = " rm " + v_sai2;
            system (v_command.c_str());
        }
        
        thread r1 (read_sam_dna, 1, v_sam1, *i);
        thread r2 (read_sam_dna, 2, v_sam2, *i);
        r1.join();
        r2.join();
        

        
        
        if (v_use_local == 1) {
            v_ref = "'" + v_db + "/local/dna/" + *i + "/" + *i + ".fas' ";
            boost::replace_all(v_ref, "\\ ", " ");

            v_sam1 = v_output + *i + ".local.1.sam";
            v_sam2 = v_output + *i + ".local.2.sam";
                
            
            v_log = v_output + "/log/" + *i + "_local_aln.1.log";
            v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r1_sort + " > " + v_sai1 + " 2>" + v_log;
            if (v_map_type == "single") {v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r0_sort + " > " + v_sai1 + " 2>" + v_log;}
            system (v_command.c_str());
            v_log = v_output + "/log/" + *i + "_local_sam.1.log";
            v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r1_sort + " > " + v_sam1 + " 2>" + v_log;
            if (v_map_type == "single") {v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r0_sort + " > " + v_sam1 + " 2>" + v_log;}
            system (v_command.c_str());
            v_command = " rm " + v_sai1;
            system (v_command.c_str());
            
            if (v_map_type == "paired") {
                v_log = v_output + "/log/" + *i + "_local_aln.2.log";
                v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r2_sort + " > " + v_sai2 + " 2>" + v_log;
                system (v_command.c_str());
                v_log = v_output + "/log/" + *i + "_local_sam.2.log";
                v_command = v_bwa + " samse " + v_ref + v_sai2 + v_r2_sort + " > " + v_sam2 + " 2>" + v_log;
                system (v_command.c_str());
                v_command = " rm " + v_sai2;
                system (v_command.c_str());
            }
            
            thread r1local (read_sam_dna, 1, v_sam1, *i);
            thread r2local (read_sam_dna, 2, v_sam2, *i);
            r1local.join();
            r2local.join();
        }
        
    }
    
  
    v_message = "Scoring reads considering other sequences ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

     
    string v_ref = "'" + v_db + "/others/dna/hla_gen_others.fasta' ";
    boost::replace_all(v_ref, "\\ ", " ");

    string v_sam1 = v_output + "others.1.sam";
    string v_sam2 = v_output + "others.2.sam";
    string v_r0_sort = v_output + v_sample + "selected.trim.R0.fastq";
    string v_r1_sort = v_output + v_sample + "selected.trim.R1.fastq";
    string v_r2_sort = v_output + v_sample + "selected.trim.R2.fastq";
    string v_sai1 = v_output + "tmp1.sai ";
    string v_sai2 = v_output + "tmp2.sai ";

     
    string v_log = v_output + "/log/others_aln.1.log";
    v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r1_sort + " > " + v_sai1 + " 2>" + v_log;
    if (v_map_type == "single") {v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r0_sort + " > " + v_sai1 + " 2>" + v_log;}
    system (v_command.c_str());
    v_log = v_output + "/log/others_sam.1.log";
    v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r1_sort + " > " + v_sam1 + " 2>" + v_log;
    if (v_map_type == "single") {v_command = v_bwa + " samse " + v_ref + v_sai1 + v_r0_sort + " > " + v_sam1 + " 2>" + v_log;}
    system (v_command.c_str());
    v_command = " rm " + v_sai1;
    system (v_command.c_str());

    if (v_map_type == "paired"){
        v_log = v_output + "/log/others_aln.2.log";
        v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_r2_sort + " > " + v_sai2 + " 2>" + v_log;
        system (v_command.c_str());
        v_log = v_output + "/log/others_sam.2.log";
        v_command = v_bwa + " samse " + v_ref + v_sai2 + v_r2_sort + " > " + v_sam2 + " 2>" + v_log;
        system (v_command.c_str());
        v_command = " rm " + v_sai2;
        system (v_command.c_str());
    }
    thread r1 (read_sam_dna, 1, v_sam1, "others");
    thread r2 (read_sam_dna, 2, v_sam2, "others");
    r1.join();
    r2.join();
     
    v_message = "Scoring reads: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    
    
    

    

   
    
    
    
    
    
    
    
    
    
    
    
    
    //Addressing sequences
     v_message = "Comparing scores ...";
     screen_message (screen_size, 0, v_message, 2, v_quiet);
     
     vector<string> v_target_list_sub;
     boost::algorithm::erase_all(v_genes, " ");
     string tmp = v_genes + ",others";
     boost::split(v_target_list_sub,tmp,boost::is_any_of(","));
     
         
    float progress = 0;
    float totalreads = sequence_size_r1.size();
 
    std::map <string, string> address_list_r1;
    ThreadPool poolscore(stoi(v_threads));
    std::vector< std::future<int> > resultscore;
    int readcount = 0;
    vector <string> table;
    
    if (v_map_type == "paired") {
        for(map<string,int>::iterator it = sequence_size_r1.begin(); it != sequence_size_r1.end(); ++it)
        {
            string id = it->first;
            pair <string,string> key = make_pair ("others",id);
            int sub1 = 0;
            int sub2 = 0;
            if (sequence_list_r1.find(key) != sequence_list_r1.end()){sub1 = sequence_list_r1[key];}
            if (sequence_list_r2.find(key) != sequence_list_r2.end()){sub2 = sequence_list_r2[key];}
            if ((sub1 != 0) && (sub2 == 0)){sequence_list_r2[key] = sub1;}
            if ((sub2 != 0) && (sub1 == 0)){sequence_list_r1[key] = sub2;}
        }
    }
    
     for(map<string,int>::iterator it = sequence_size_r1.begin(); it != sequence_size_r1.end(); ++it)
     {
         string id = it->first;
         
         resultscore.emplace_back(
         poolscore.enqueue([readcount, id, v_target_list_sub, &address_list_r1, &table, &progress, totalreads]
           {
         
             map <string,int> nm_sum;
             string table_str = id;
         
             int min_nm_sum = 100000;
             string address_to = "";
             int count = 0;
         
             for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
             {
                 if (*g == "") {continue;}
             
                 nm_sum[*g] = 10000;
                 pair <string,string> key = make_pair (*g,id);
             
                 string scores = "-";
             
                 if (v_map_type == "paired"){
                     
                     if ((sequence_list_r1.find(key) != sequence_list_r1.end()) && (sequence_list_r2.find(key) != sequence_list_r2.end())){
                         if (((sequence_list_r1[key] != 10000) && (sequence_list_r2[key] != 10000)) && ((sequence_list_r1[key] != 0) && (sequence_list_r2[key] != 0)))
                         {
                             int v_max_r1 = 0;
                             int v_max_r2 = 0;
                             v_max_r1 = int((sequence_size_r1[id] * v_tolerance));
                             v_max_r2 = int((sequence_size_r2[id] * v_tolerance));
                             
//                           if (*g == "others"){
//                           if ((sequence_list_r1[key] == 1) && (sequence_list_r2[key] > v_max_r2))
//                             {mtxg.lock();sequence_list_r2[key] = 1;mtxg.unlock();}
//                             if ((sequence_list_r2[key] == 1) && (sequence_list_r1[key] > v_max_r1))
//                             {mtxg.lock();sequence_list_r1[key] = 1;mtxg.unlock();}
//                            }
                             
                             if (((sequence_list_r1[key]-1) <= v_max_r1) && ((sequence_list_r2[key]-1) <= v_max_r2))
                             {
                                 scores = to_string(sequence_list_r1[key]) + "," + to_string(sequence_list_r2[key]);
                                 if (min_nm_sum > (sequence_list_r1[key] + sequence_list_r2[key])) {min_nm_sum = (sequence_list_r1[key] + sequence_list_r2[key]);count++;}
                                 nm_sum[*g] = sequence_list_r1[key] + sequence_list_r2[key];
                             }
                         }
                     }
                 }
                 
                 if (v_map_type == "single"){
                     if (sequence_list_r1.find(key) != sequence_list_r1.end()){
                         if ((sequence_list_r1[key] != 10000) && (sequence_list_r1[key] != 0))
                         {
                             int v_max_r1 = 0;
                             v_max_r1 = int((sequence_size_r1[id] * v_tolerance));
                             if ((sequence_list_r1[key]-1) <= v_max_r1)
                             {
                                 scores = to_string(sequence_list_r1[key]);
                                 if (min_nm_sum > sequence_list_r1[key] ) {min_nm_sum = (sequence_list_r1[key]);count++;}
                                 nm_sum[*g] = sequence_list_r1[key];
                             }
                         }
                    }
                 }
                
                 table_str.append ("\t" + scores);
         }
         
             if (count >= 1)
             {
                 int count_min_nm_sum = 0;
                 for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
                 {
                     if (nm_sum[*g] == min_nm_sum)
                     {
                         count_min_nm_sum++;
                         address_to = address_to + "," + *g;
                         pair <string,string> key = make_pair (*g,id);
                         mtxg.lock();
                         sequence_address[key] = 1;
                         mtxg.unlock();
                     }
                 }
                 
                 if (address_to != "") {
                     mtxg.lock();
                     address_list_r1[id] = address_to;
                     mtxg.unlock();
                 }
                 
             }
            
             
             string subaddress = address_to;
             if (subaddress == "") {subaddress = ",-";}
             subaddress = subaddress.substr(1);
             
             table_str.append("\t" + subaddress);
             table_str.append("\n");
             mtxg.lock();
             table.push_back(table_str);
             progress++;
             float ratio = ((progress / totalreads) * 100);
             v_message = "Comparing scores ... " + to_string(ratio) + " % ";
             screen_message (screen_size, 0, v_message, 2, v_quiet);
             mtxg.unlock();
             return 1;
           })
         );
     }
    for(auto && result: resultscore){result.get();} // waiting for all threads
    sequence_list_r1.clear();
    sequence_list_r2.clear();
    sequence_size_r1.clear();
    sequence_size_r2.clear();
    
    
    string v_add = v_output + v_sample + "addressing_table.txt";
    ofstream ADD;
    ADD.open (v_add.c_str());
        
    ADD << "Read";
    for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
    {
       if (*g == "") {continue;}
       ADD << "\t" << *g;
    }
    ADD << "\tTarget" << endl;
    for (auto & item : table)
    {
        ADD << item;
    }
    ADD.close();
     
     v_message = "Comparing scores: done";
     screen_message (screen_size, 0, v_message, 1, v_quiet);
     
     
  
    
 
    


    
    
    
    
    
    
    
    
    
    // Creating fastq after filtering
    unordered_map <string,string> original_data_r1;
    unordered_map <string,string> original_data_r2;
    unordered_map <string,string> original_data_r1_trim;
    unordered_map <string,string> original_data_r2_trim;
    vector <int> list_of_read_sizes;
    
    v_message = "Making fastq files ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

    string r1_original = v_output + v_sample + "selected.R1.fastq";
    string r2_original = v_output + v_sample + "selected.R2.fastq";
    string r1_original_trim = v_output + v_sample + "selected.trim.R1.fastq";
    string r2_original_trim = v_output + v_sample + "selected.trim.R2.fastq";
    if (v_map_type == "single") {r1_original = v_output + v_sample + "selected.R0.fastq";}
    if (v_map_type == "single") {r1_original_trim = v_output + v_sample + "selected.trim.R0.fastq";}

  
  

    
    ifstream R1( r1_original.c_str());
    ifstream R2( r2_original.c_str());
    for( std::string line; getline( R1, line ); )
    {
        string id = line;
        getline( R1, line );
        string seq = line;
        getline( R1, line );
        string info = line;
        getline( R1, line );
        string qual = line;
        list_of_read_sizes.push_back(seq.size());
        vector<string> read_id;
        boost::split(read_id,id,boost::is_any_of(" "));
//        if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}

        
        original_data_r1[read_id[0]] = id + "\n" + seq + "\n" + info + "\n" + qual;
 
        getline( R2, line );
        id = line;
        getline( R2, line );
        seq = line;
        getline( R2, line );
        info = line;
        getline( R2, line );
        qual = line;
        original_data_r2[read_id[0]] = id + "\n" + seq + "\n" + info + "\n" + qual;
    }
    R1.close();
    R2.close();
    
    int total_size = 0;
    int count_reads = 0;
    for (auto &&item: list_of_read_sizes)
    {
        total_size = total_size + item;
        count_reads++;
    }
    if (count_reads > 0) {read_mean_size = (total_size / count_reads);}
    if (count_reads == 0) {read_mean_size = 0;}
    list_of_read_sizes.clear();

    
    ifstream R1t( r1_original_trim.c_str());
    ifstream R2t( r2_original_trim.c_str());
    for( std::string line; getline( R1t, line ); )
    {
        string id = line;
        getline( R1t, line );
        string seq = line;
        getline( R1t, line );
        string info = line;
        getline( R1t, line );
        string qual = line;
        vector<string> read_id;
        boost::split(read_id,id,boost::is_any_of(" "));
//        if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}

        original_data_r1_trim[read_id[0]] = seq + "\n" + info + "\n" + qual;
        getline( R2t, line );
        id = line;
        getline( R2t, line );
        seq = line;
        getline( R2t, line );
        info = line;
        getline( R2t, line );
        qual = line;
        original_data_r2_trim[read_id[0]] = seq + "\n" + info + "\n" + qual;
     }
     R1t.close();
     R2t.close();
    


    if (v_map_type == "paired") {
       for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
       {
           if (*i == "") {continue;}
           string gene = *i;
           
           v_message = "Making fastq files for " + *i + " ...";
           screen_message (screen_size, 0, v_message, 2, v_quiet);
           
           string sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";
           
           map <string,int> ids;
           ifstream idsrec( sort.c_str());
           for( std::string line; getline( idsrec, line ); )
           {
               string id = line;
               getline( R1t, line );
               string seq = line;
               getline( R1t, line );
               string info = line;
               getline( R1t, line );
               string qual = line;
               
               vector <string> read_id;
               boost::split(read_id,id,boost::is_any_of(" "));
//               if ((read_id[0].substr(read_id[0].size()-2,2) == "/1") || (read_id[0].substr(read_id[0].size()-2,2) == "/2")){read_id[0] = read_id[0].substr(0,read_id[0].size()-2);}
               ids[read_id[0]] = 1;
           }
           idsrec.close();
  
           int size = gene_opt_end[gene] - gene_opt_start[gene];
           int downS = (downsampling * size) / read_mean_size;
           int downS_counter = 0;
           
           string v_out_R1 = v_output + v_sample + *i + "_R1.fastq";
           ofstream OUT_R1;
           OUT_R1.open (v_out_R1.c_str());
           
           string v_out_R2 = v_output + v_sample + *i + "_R2.fastq";
           ofstream OUT_R2;
           OUT_R2.open (v_out_R2.c_str());

           string v_out_clean_R1 = v_output + v_sample + *i + ".trim.R1.fastq";
           ofstream OUT_clean_R1;
           OUT_clean_R1.open (v_out_clean_R1.c_str());

           string v_out_clean_R2 = v_output + v_sample + *i + ".trim.R2.fastq";
           ofstream OUT_clean_R2;
           OUT_clean_R2.open (v_out_clean_R2.c_str());

           
           for (auto & id : ids)
           {
               
               vector <string> data_split_r1;
               boost::split(data_split_r1,original_data_r1[id.first],boost::is_any_of("\n"));

               pair <string,string> key = make_pair (gene,id.first.substr(1));
               if (sequence_address.find(key) == sequence_address.end()) {continue;}
   
               vector <string> data_split_r1_trim;
               boost::split(data_split_r1_trim,original_data_r1_trim[id.first],boost::is_any_of("\n"));
   
                vector <string> data_split_r2;
                boost::split(data_split_r2,original_data_r2[id.first],boost::is_any_of("\n"));
                
                vector <string> data_split_r2_trim;
                boost::split(data_split_r2_trim,original_data_r2_trim[id.first],boost::is_any_of("\n"));
               
               

               OUT_R1 << data_split_r1[0] << endl;
               OUT_R1 << data_split_r1[1] << endl;
               OUT_R1 << data_split_r1[2] << endl;
               OUT_R1 << data_split_r1[3] << endl;
               OUT_R2 << data_split_r2[0] << endl;
               OUT_R2 << data_split_r2[1] << endl;
               OUT_R2 << data_split_r2[2] << endl;
               OUT_R2 << data_split_r2[3] << endl;

               if (downS_counter < downS) {
                           OUT_clean_R1 << id.first + " 1:trim" << endl;
                           OUT_clean_R1 << data_split_r1_trim[0] << endl;
                           OUT_clean_R1 << data_split_r1_trim[1] << endl;
                           OUT_clean_R1 << data_split_r1_trim[2] << endl;
                           OUT_clean_R2 << id.first + " 2:trim" << endl;
                           OUT_clean_R2 << data_split_r2_trim[0] << endl;
                           OUT_clean_R2 << data_split_r2_trim[1] << endl;
                           OUT_clean_R2 << data_split_r2_trim[2] << endl;
               }
                downS_counter++;
                continue;
           }
           OUT_R1.close();
           OUT_R2.close();
           OUT_clean_R1.close();
           OUT_clean_R2.close();
       }
    }

    
    if (v_map_type == "single") {
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i == "") {continue;}
            string gene = *i;
            
            v_message = "Making fastq files for " + *i + " ...";
            screen_message (screen_size, 0, v_message, 2, v_quiet);

            string sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";
            
            map <string,int> ids;
            ifstream idsrec( sort.c_str());
            for( std::string line; getline( idsrec, line ); )
            {
                string id = line;
                getline( R1t, line );
                string seq = line;
                getline( R1t, line );
                string info = line;
                getline( R1t, line );
                string qual = line;
                vector <string> read_id;
                boost::split(read_id,id,boost::is_any_of(" "));
                ids[read_id[0]] = 1;
            }
            idsrec.close();
            
            
            
            

            int size = gene_opt_end[gene] - gene_opt_start[gene];
            int downS = (downsampling * size) / read_mean_size;
            int downS_counter = 0;

            
            string v_out_R1 = v_output + v_sample + *i + "_R0.fastq";
            ofstream OUT_R1;
            OUT_R1.open (v_out_R1.c_str());
            
            string v_out_singlet = v_output + v_sample + *i + ".trim.R0.fastq";
            ofstream OUT_singlet;
            OUT_singlet.open (v_out_singlet.c_str());
            
            for (auto & id : ids)
            {
                vector <string> data_split_r1;
                boost::split(data_split_r1,original_data_r1[id.first],boost::is_any_of("\n"));
                
                pair <string,string> key = make_pair (gene,id.first.substr(1));
                if (sequence_address.find(key) == sequence_address.end()) {continue;}
    
                vector <string> data_split_r1_trim;
                boost::split(data_split_r1_trim,original_data_r1_trim[id.first],boost::is_any_of("\n"));
  
 
                OUT_R1 << data_split_r1[0] << endl;
                OUT_R1 << data_split_r1[1] << endl;
                OUT_R1 << data_split_r1[2] << endl;
                OUT_R1 << data_split_r1[3] << endl;

                if (downS_counter < downS) {
                    OUT_singlet << id.first << " 1:trim" << endl;
                    OUT_singlet << data_split_r1_trim[0] << endl;
                    OUT_singlet << data_split_r1_trim[1] << endl;
                    OUT_singlet << data_split_r1_trim[2] << endl;
                }
                downS_counter++;
                continue;
            }
            OUT_R1.close();
            OUT_singlet.close();
        }
    }
    original_data_r2.clear();
    original_data_r1.clear();
    original_data_r1_trim.clear();
    original_data_r2_trim.clear();
    
   
    v_message = "Making fastq files: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);


    

 
    
    
    
   
     
    
    
    
    
    
    //Typing
    
    
    int mincov = 0;
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i == "KIR2DL4") {mincov = 5;}
        if (*i == "HLA-G") {mincov = 10;}
    }

    map <string,string> typing_result_alleleA;
    map <string,string> typing_result_alleleB;

    if (v_skiptyping == 0) {
        v_message = "Performing adjustments ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
        
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i == "") {continue;}
            if (v_typeall == 0) {if (gene_type[*i] != "TYPE") {continue;}}
            
                string gene = *i;
                v_message = "Performing adjustments ... " + gene + " ... ";
                screen_message (screen_size, 0, v_message, 2, v_quiet);
           
                string fq1 = "";
                string fq2 = "";
                if (v_map_type == "paired")
                {
                    fq1 = v_output + v_sample + *i + ".trim.R1.fastq";
                    fq2 = v_output + v_sample + *i + ".trim.R2.fastq";
                }
                if (v_map_type == "single")
                {
                    fq1 = v_output + v_sample + *i + ".trim.R0.fastq";
                    fq2 = "";
                }
            
            
                string result = typing_dna(gene, fq1, fq2, 20, 5, v_map_type, mincov);
                if (v_debug == 1) {cout << result << endl;}
                
                mtx_selection_general.lock();
                if (result.substr(0,4) == "fail")
                {
                    typing_result_alleleA[gene] = "";
                    typing_result_alleleA[gene] = "";
                    mtx_selection_general.unlock();
                    continue;
                }
                vector <string> results;
                boost::split(results,result,boost::is_any_of(","));
                vector <string> alleles;
                boost::split(alleles,results[0],boost::is_any_of("\t"));
                typing_result_alleleA[gene] = alleles[0];
                typing_result_alleleB[gene] = alleles[1];
                mtx_selection_general.unlock();
        }
        
        general_log << endl;
/*
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i == "") {continue;}

            if (typing_result_alleleA[*i] != "") {
                general_log << "Possible genotype for " << *i << ": " << typing_result_alleleA[*i] << ", " << typing_result_alleleB[*i] << endl;
            }
            else {
                if ((gene_type[*i] == "TYPE") || (v_typeall == 1))
                {general_log << "Possible genotype for " << *i << ": failed" << endl;}
                if (gene_type[*i] != "TYPE")
                {general_log << "Possible genotype for " << *i << ": disabled" << endl;}
            }

        }
        general_log << endl;
 */
        v_message = "Performing adjustments: done";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        
    }
 
    
    for (auto & item : reads_possibly_mapped)
    {
        string gene = item.first.first;
        string read = item.first.second;
        int nm = item.second;
        
        if (reads_mapped_nm.find(read) == reads_mapped_nm.end())
        {
            reads_mapped_nm[read] = nm;
            reads_mapped_gene[read] = gene;
        }
        else
        {
            if (reads_mapped_nm[read] == nm)
            {
                reads_mapped_gene[read] = reads_mapped_gene[read] + ";" + gene;
            }
            if (reads_mapped_nm[read] > nm)
            {
                reads_mapped_nm[read] = nm;
                reads_mapped_gene[read] = gene;
            }
        }
    }
    
    

    


    
    
    
    
    
    

    
    

    
    
    // Final mapping
    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        
        if (*i == "") {continue;}
        v_message = "Mapping " + *i + " ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        
        

        /*
         string v_reference = "'" + v_db + "/reference/dna/loci/" + *i + ".fas'";
         boost::replace_all(v_reference, "\\ ", " ");
         string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
         string fastq_assembly = v_output + v_sample + *i + ".assemble.fastq";
         string v_sam_assembly = v_output + v_sample + *i + ".assemble.sam";;
         string v_sam_assembly_correct = v_output + v_sample + *i + ".assemble.corrected.sam";;
         
         
         v_log = v_output + "log/" + v_sample + *i + ".assemble.log";
         v_command = v_bwa + " mem -t " + v_threads + " -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + fastq_assembly + " > " + v_sam_assembly + " 2>" + v_log;
         system (v_command.c_str());


         ifstream sama (v_sam_assembly.c_str());
         ofstream OUT1;
         OUT1.open (v_sam_assembly_correct.c_str());
         int v_pos_correct = position_hg38[*i] - 1;
         
         for( std::string line; getline( sama, line ); )
         {
             if (line == "") {continue;}
             vector<string> sam_line;
             boost::split(sam_line,line,boost::is_any_of("\t"));
             
             
             if (sam_line[0] == "@SQ")
             {
                 OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
                 OUT1 << "@CO\thla-mapper " << Program_version << ", human genome version hg38 " << endl;
                 OUT1 << "@CO\tartificial haplotypes" << endl;
                 continue;
             }
             
             if ((sam_line[0].substr(0,1) == "@") || (sam_line[2] == "*"))
             {
                 OUT1 << line << endl;
                 continue;
             }
             
             string seq = sam_line[9];
             if (seq.length() < v_size) {continue;}

             sam_line[2] = chr_hg38[*i];
             int correct = stoi(sam_line[3]) + v_pos_correct;
             sam_line[3] = to_string(correct);
             OUT1 << sam_line[0];
             for (int helper = 1; helper < sam_line.size(); helper++)
             {OUT1 << "\t" << sam_line[helper];}
             OUT1 << endl;
         }
         OUT1.close();
        sama.close();
        
         string v_bam_assembly = v_output + v_sample + *i + "_assemble.bam";;
         v_log = v_output + "log/" + v_sample + *i + ".assemble.sort.log";
         v_command = v_samtools + " sort -@ " + v_threads + " " + v_sam_assembly_correct + " > " + v_bam_assembly + " 2>" + v_log;
         system (v_command.c_str());
         v_command = v_samtools + " index " + v_bam_assembly + " 2>" + v_log;
         system (v_command.c_str());

         removefile(v_sam_assembly_correct,v_debug);
         removefile(v_sam_assembly,v_debug);

        */
        
            string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
            string v_reference = "'" + v_db + "/reference/dna/loci/" + *i + ".fas'";
            boost::replace_all(v_reference, "\\ ", " ");

            string in0 = v_output + v_sample + *i + "_R0.fastq";
            string in1 = v_output + v_sample + *i + "_R1.fastq";
            string in2 = v_output + v_sample + *i + "_R2.fastq";
            v_sam1 = v_output + v_sample + *i + ".tmp.sam";
            v_sam2 = v_output + v_sample + *i + ".unique.sam";
            string v_sam3 = v_output + v_sample + *i + ".adjusted.sam";
            v_log = v_output + "log/" + v_sample + *i + ".map.log";
             
             
            if (v_map_type == "paired") {
                v_command = v_bwa + " mem -t " + v_threads + " -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in1 + " " + in2 + " > " + v_sam1 + " 2>" + v_log;
                if (*i == "HLA-DRB1")
                {
                    v_command = v_bwa + " mem -t " + v_threads + " -k 15 -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in1 + " " + in2 + " > " + v_sam1 + " 2>" + v_log;
                }
                
            }
            if (v_map_type == "single") {
                v_command = v_bwa + " mem -t " + v_threads + " -B 2 -O 3,3 -L 3,3 -R '@RG\\tID:" + v_tag_sub + "\\tLB:" + v_tag_sub + "\\tSM:" + v_tag_sub + "\\tPL:illumina\\tPU:" + v_tag_sub + "' " + v_reference + " " + in0 + " > " + v_sam1 + " 2>" + v_log;
            }
            system (v_command.c_str());
             

            ifstream reference (v_reference.c_str());
            string refseq = "";
            for( std::string line; getline( reference, line ); )
            {
                if (line.substr(0,1) == ">") {continue;}
                if (line == "") {continue;}
                refseq = refseq + line;
            }
            reference.close();
            ref_size[*i] = refseq.length();
            refseq = "";
             
             
            fstream sam (v_sam1.c_str());
        
            ofstream OUT1;
            OUT1.open (v_sam2.c_str());
            ofstream OUT2;
            OUT2.open (v_sam3.c_str());
        
             
            int count_reads = 0;
            int v_pos_correct = position_hg38[*i] - 1;
            unordered_map <string,int> adjusted_for_secondary;
            unordered_map <string,int> adjusted_for_primary;
             
            for( std::string line; getline( sam, line ); )
            {
                if (line == "") {continue;}
                vector<string> sam_line;
                boost::split(sam_line,line,boost::is_any_of("\t"));
                 
                if (sam_line[0] == "@SQ")
                {
                    if (usechr == 0) {
                        OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
                    }

                    if (usechr == 1) {
                        size_t found = chr_hg38[*i].find("chr");
                        if (found==std::string::npos){OUT1 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
                        else {OUT1 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
                    }

                    OUT1 << "@CO\thla-mapper " << Program_version << ", human genome version hg38 " << endl;
                    
                    if (usechr == 0) {
                        OUT2 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
                    }
                    if (usechr == 1) {
                        size_t found = chr_hg38[*i].find("chr");
                        if (found==std::string::npos){OUT2 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
                        else {OUT2 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;}
                    }
                    OUT2 << "@CO\thla-mapper " << Program_version << ", human genome version hg38 " << endl;
                    continue;
                }
                 
                if ((sam_line[0].substr(0,1) == "@") || (sam_line[2] == "*"))
                {
                    OUT1 << line << endl;
                    OUT2 << line << endl;
                    continue;
                }


                string seq = sam_line[9];
                if (seq.length() < v_size) {continue;}

                sam_line[2] = chr_hg38[*i];
                if (usechr == 1)
                {
                    size_t found = chr_hg38[*i].find("chr");
                    if (found==std::string::npos){sam_line[2] = "chr" + chr_hg38[*i];}
                }
             
                string tmp = address_list_r1[sam_line[0]] + ",";
                vector<string> address_multi_hits;
                boost::split(address_multi_hits,tmp,boost::is_any_of(","));
                 
 
                if (address_multi_hits.size() < 4)
                {
                    OUT1 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
                    OUT1 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
                    OUT1 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
                    OUT1 << (stoi(sam_line[7]) + v_pos_correct);
                    for (int helper = 8; helper < sam_line.size(); helper++)
                    {OUT1 << "\t" << sam_line[helper];}
                    OUT1 << endl;
                    OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
                    OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
                    OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
                    OUT2 << (stoi(sam_line[7]) + v_pos_correct);
                    for (int helper = 8; helper < sam_line.size(); helper++)
                    {OUT2 << "\t" << sam_line[helper];}
                    OUT2 << endl;
                    count_reads++;
                    continue;
                }
        
                
                if (address_multi_hits.size() >= 4)
                {
                    int adjust = 0;

                    pair <string,string> key = make_pair(*i,sam_line[0]);
                    if (reads_mapped_gene.find(sam_line[0]) != reads_mapped_gene.end())
                    {
                        if (reads_mapped_gene[sam_line[0]] == *i)
                        {
                            adjust = 1;
                        }
                    }
 
                    if (adjust == 1) {
                            OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
                            OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
                            OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
                            OUT2 << (stoi(sam_line[7]) + v_pos_correct);
                            for (int helper = 8; helper < sam_line.size(); helper++)
                            {OUT2 << "\t" << sam_line[helper];}
                            OUT2 << endl;
                            count_reads++;
                            adjusted_for_primary[sam_line[0]] = 1;
                    }
                         
                        // for paired
                        if (sam_line[1] == "163") {sam_line[1] = "419";}
                        if (sam_line[1] == "83") {sam_line[1] = "339";}
                        if (sam_line[1] == "99") {sam_line[1] = "355";}
                        if (sam_line[1] == "147") {sam_line[1] = "403";}
                        if (sam_line[1] == "81") {sam_line[1] = "337";}
                        if (sam_line[1] == "167") {sam_line[1] = "423";}
                        if (sam_line[1] == "145") {sam_line[1] = "401";}
                        if (sam_line[1] == "97") {sam_line[1] = "353";}
                        if (sam_line[1] == "185") {sam_line[1] = "441";}
                        if (sam_line[1] == "73") {sam_line[1] = "329";}
                        if (sam_line[1] == "113") {sam_line[1] = "369";}
                        if (sam_line[1] == "117") {sam_line[1] = "373";}
                        if (sam_line[1] == "121") {sam_line[1] = "377";}
                        if (sam_line[1] == "133") {sam_line[1] = "389";}
                        if (sam_line[1] == "137") {sam_line[1] = "393";}
                        if (sam_line[1] == "177") {sam_line[1] = "433";}
                        if (sam_line[1] == "181") {sam_line[1] = "437";}
                        if (sam_line[1] == "161") {sam_line[1] = "417";}

                        // for single
                        if (sam_line[1] == "16") {sam_line[1] = "272";}
                        if (sam_line[1] == "0") {sam_line[1] = "256";}

                         
                        adjusted_for_secondary[sam_line[0]] = 1;
                         
                        OUT1 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
                        OUT1 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
                        OUT1 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
                        OUT1 << (stoi(sam_line[7]) + v_pos_correct);
                        for (int helper = 8; helper < sam_line.size(); helper++)
                        {OUT1 << "\t" << sam_line[helper];}
                        OUT1 << endl;
                         
                        if (adjust == 0) {
                               OUT2 << sam_line[0] << "\t" << sam_line[1] << "\t" << sam_line[2] << "\t";
                               OUT2 << (stoi(sam_line[3]) + v_pos_correct) << "\t";
                               OUT2 << sam_line[4] << "\t" << sam_line[5] << "\t" << sam_line[6] << "\t";
                               OUT2 << (stoi(sam_line[7]) + v_pos_correct);
                               for (int helper = 8; helper < sam_line.size(); helper++)
                               {OUT2 << "\t" << sam_line[helper];}
                               OUT2 << endl;
                           }
                     
                }
            }
            sam.close();
            OUT1.close();
            OUT2.close();

            string log = v_output + v_sample + *i + ".log";
            ofstream log_mapping;
            log_mapping.open (log.c_str());

            log_mapping << endl << "List of reads marked as secondary" << endl;
            for (auto & item : adjusted_for_secondary)
            {
                 log_mapping << item.first << endl;
            }
            log_mapping << endl << endl;
            
            if (typing_result_alleleA[*i] != "") {
                log_mapping << "Possible genotype for " << *i << ": " << typing_result_alleleA[*i] << ", " << typing_result_alleleB[*i] << endl;
            }
        
            log_mapping << endl;
            log_mapping << "List of reads adjusted for primary after pseudo-typing:" << endl;
            for (auto & item : adjusted_for_primary)
            {
                 log_mapping << item.first << endl;
            }
            log_mapping << endl;
             
             
            log_mapping << "Number of mapped reads (only primary mappings): " << to_string(count_reads) << endl;
             

//            double cov = ((double(count_reads) * read_mean_size) / ref_size[*i]);
        
//           log_mapping << "Depth: " << to_string(cov) << endl;
//            general_log << *i << " depth: " << to_string(cov) << endl;
            read_count_gene[*i] = double(count_reads);
             
            log_mapping.close();
        
    }
    v_message = "Mapping: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    
    
    
    
    reads_mapped_nm.clear();
    reads_mapped_gene.clear();
    read_count_gene.clear();
    reads_possibly_mapped.clear();
    sequence_list_r1.clear();
    sequence_list_r2.clear();
    sequence_size_r1.clear();
    sequence_size_r2.clear();
    sequence_address.clear();
    read_count_gene.clear();
    
    
    v_message = "Merging files ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);
    
    map <string,int> sq_data;
    map <string,int> rg_data;
    unordered_map <string,int> optimized_reads;
    vector <string> adjusted_data;
    vector <string> unique_data;
    vector <string> hg38_data;
    
    
    v_message = "Merging files ... reading new alignments ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

    thread t1 ([&sq_data,&optimized_reads,&adjusted_data,&rg_data]{
    
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".adjusted.sam";
            ifstream input( file.c_str() );
            for( std::string line; getline( input, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,3) == "@SQ") {
 //                   if (v_bam == ""){
                        mtxg.lock();sq_data[line] = 1;mtxg.unlock();
  //                  }
                    continue;
                }
                if (line.substr(0,3) == "@RG") {mtxg.lock();rg_data[line] = 1;mtxg.unlock();continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                
                if (usechr == 1)
                {
                    if (data[2] == "6") {data[2] = "chr6";}
                    if (data[2] == "19") {data[2] = "chr19";}
                    line = data[0];
                    for (int a = 1; a < data.size(); a++)
                    {
                        line.append("\t" + data[a]);
                    }
                }
                
                mtxg.lock();
                optimized_reads[data[0]] = 1;
                mtxg.unlock();
                adjusted_data.push_back(line);
            }
            input.close();
        }
    });
    
    
    thread t2 ([&sq_data,&optimized_reads,&unique_data,&rg_data]{

        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".unique.sam";
            ifstream input( file.c_str() );
            for( std::string line; getline( input, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,3) == "@SQ")
                {
 //                   if (v_bam == ""){
                        mtxg.lock();sq_data[line] = 1;mtxg.unlock();
   //                 }
                    continue;
                }
                if (line.substr(0,3) == "@RG") {mtxg.lock();rg_data[line] = 1;mtxg.unlock();continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                
                if (usechr == 1)
                {
                    if (data[2] == "6") {data[2] = "chr6";}
                    if (data[2] == "19") {data[2] = "chr19";}
                    line = data[0];
                    for (int a = 1; a < data.size(); a++)
                    {
                        line.append("\t" + data[a]);
                    }
                }
                
                mtxg.lock();
                optimized_reads[data[0]] = 1;
                mtxg.unlock();
                unique_data.push_back(line);
            }
            input.close();
        }
     });
    
    t1.join();
    t2.join();
    
    
 
    
    
    if (v_bam != "")
    {
        general_log << endl << "Reads forced to MQ=0 after optimization:" << endl;
        v_message = "Merging files ... reading original alignments ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        string hg38_select = v_output + v_sample + "hg38.sam";

        string v_bed = "'" + v_db + "/bed/target_dna.bed'";
        boost::replace_all(v_bed, "\\ ", " ");

        v_command = v_samtools + " view -h -@ " + v_threads + " -ML " + v_bed + " " + v_bam;
        v_system_out = GetStdoutFromCommand(v_command);
        
        vector <string> hg38_full;
        boost::split(hg38_full,v_system_out,boost::is_any_of("\n"));
        
        for( auto &line : hg38_full)
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == "[") {continue;}
            if (line.substr(0,3) == "@SQ") {
                //sq_data[line] = 1;
                continue;
            }
            if (line.substr(0,3) == "@RG") {continue;}
            if (line.substr(0,1) == "@") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            if ((data[2] == "*") || (data[1] == "*")) {continue;}
            if ((data[3] == "*") || (data[3] == "0")) {continue;}
            
            
            if (optimized_reads.find(data[0]) == optimized_reads.end())
            {
                
                for (int a = 11; a < data.size(); a++)
                {
                    if (data[a].substr(0,5) == "RG:Z:") {data[a] = "RG:Z:" + v_sample_sub;break;}
                }
                
                int pos = stoi(data[3]);
                int helper = 0;
                for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
                {
                    if (*i ==  "") {continue;}
                    if ((pos >= gene_opt_start[*i]) && (pos <= gene_opt_end[*i]))
                    {
                        helper = 1;
                        break;
                    }
                }
                if (helper == 1){data[4] = "0"; general_log << data[0] << endl;}
                
                string news = data[0];
                for (int a = 1; a < data.size(); a++)
                {news.append("\t" + data[a]);}
                hg38_data.push_back(news);
                continue;
            }
        }
     }
    
    
    
    
    
    v_message = "Merging files ... writing merged files ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

    string merged_sam_adjusted = v_output + v_sample_sub + ".adjusted.sam";
    ofstream sam_adjusted;
    sam_adjusted.open (merged_sam_adjusted .c_str());
    sam_adjusted <<  "@CO\thla-mapper " << Program_version << ", human genome version hg38" << endl;
    
    for (auto &item : sq_data) {sam_adjusted << item.first << endl;}
    for (auto &item : rg_data) {sam_adjusted << item.first << endl;}
    for (auto &item : adjusted_data) {sam_adjusted << item << endl;}
    for (auto &item : hg38_data) {sam_adjusted << item << endl;}
    sam_adjusted.close();

    string v_bam_out = v_output + v_sample_sub + ".adjusted.bam";
    v_log = v_output + "/log/sort.log";
    v_command = v_samtools + " sort -@ " + v_threads + " -m 1g " + merged_sam_adjusted  + " > " +  v_bam_out + " 2>" + v_log;
    system (v_command.c_str());
    v_command = v_samtools + " index " + v_bam_out + " 2>" + v_log;
    system (v_command.c_str());

    
    string merged_sam_unique = v_output + v_sample_sub + ".unique.sam";
    ofstream sam_unique;
    sam_unique.open (merged_sam_unique .c_str());
    sam_unique <<  "@CO\thla-mapper " << Program_version << ", human genome version hg38" << endl;
    for (auto &item : sq_data) {sam_unique << item.first << endl;}
    for (auto &item : rg_data) {sam_unique << item.first << endl;}
    for (auto &item : unique_data) {sam_unique << item << endl;}
    for (auto &item : hg38_data) {sam_unique << item << endl;}
    sam_unique.close();
    
    v_bam_out = v_output + v_sample_sub + ".unique.bam";
    v_log = v_output + "/log/sort.log";
    v_command = v_samtools + " sort -@ " + v_threads + " -m 1g " + merged_sam_unique  + " > " +  v_bam_out + " 2>" + v_log;
    system (v_command.c_str());
    v_command = v_samtools + " index " + v_bam_out + " 2>" + v_log;
    system (v_command.c_str());
    
    v_message = "Merging files: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

      
  
    
    
    v_message = "Cleaning temporary files ... ";
    screen_message (screen_size, 0, v_message, 2, v_quiet);
    
    removefile(v_output + v_sample + "hg38.selected.sam",v_debug);
    if (v_debug == 1) {keep_sam = 1;}
    removefile(v_output + v_sample + "hg38.sam",keep_sam);
    removefile(v_output + v_sample + "hg38.bam",v_debug);
    removefile(v_output + v_sample + "hg38.unmapped.sam",v_debug);
    removefile(v_output + v_sample + "hg38.unmapped.bam",v_debug);

    removefile(v_output + v_sample_sub + ".adjusted.sam",v_debug);
    removefile(v_output + v_sample_sub + ".unique.sam",v_debug);
    removefile(v_output + v_sample + "R0.fastq",v_debug);
    removefile(v_output + v_sample + "R1.fastq",v_debug);
    removefile(v_output + v_sample + "R2.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R0.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R1.fastq",v_debug);
    removefile(v_output + v_sample + "selected.trim.R2.fastq",v_debug);
    removefile(v_output + v_sample + "mapped_r0.tmp.fq",v_debug);
    removefile(v_output + v_sample + "mapped_r1.tmp.fq",v_debug);
    removefile(v_output + v_sample + "mapped_r2.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r0.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r1.tmp.fq",v_debug);
    removefile(v_output + v_sample + "unmapped_r2.tmp.fq",v_debug);

    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i ==  "") {continue;}
        removefile(v_output + v_sample + *i + ".trim.R1.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R0.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R2.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".adjusted.sam",v_debug);
        removefile(v_output + v_sample + *i + ".unique.sam",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R0.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R1.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R2.fastq",v_debug);

 //       string file = v_output + v_sample + *i + "_assemble.bam";
 //       if (filesize(file.c_str()) == 0) {removefile(file,v_debug);}

        removefile(v_output + v_sample + *i + ".tmp.sam",v_debug);
 
    }
    
    /*
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
       {
           if (*i ==  "") {continue;}
           string region = chr_hg38[*i] + ":" + to_string(gene_opt_start[*i]) + "-" + to_string(gene_opt_end[*i]);
           string tags = "-g 65 -g 69 -g 73 -g 77 -g 81 -g 83 -g 97 -g 99 -g 113 -g 117 -g 121 -g 129 -g 133 -g 137 -g 141 -g 145 -g 147 -g 161 -g 163 -g 177 -g 181 -g 185 -g 196";
           
           string v_command = v_samtools + " depth -a -q 20 -r " + region + " " + tags + " " + v_output + v_sample_sub + ".adjusted.bam";
           v_system_out = GetStdoutFromCommand(v_command);
           
           vector <string> data;
           boost::split(data,v_system_out,boost::is_any_of("\n"));
           float sum = 0;
           for (auto & item : data)
           {
               if (item == "") {continue;}
               vector <string> values;
               boost::split(values,item,boost::is_any_of("\t"));
               sum = sum + stoi(values[2]);
           }
           float depth = 0;
           if (data.size() > 0) {depth = sum / data.size();}
           general_log << *i << "\tdepth\t" << depth << endl;
       }
    */
    
    
    general_log.close();
    int quiet_mem = v_quiet;
    if (v_type_after_map == 1) {
        v_message = "Calling the typing algorithm ...";
        screen_message (screen_size, 0, v_message, 0, v_quiet);
        v_mapout =   v_output;
        v_quiet = 1;
        main_typing();
        v_quiet = quiet_mem;
        v_message = "Calling the typing algorithm: done";
        screen_message (screen_size, 0, v_message, 1, v_quiet);

    }
    
}
    
      



       


    





    
    
    
    
    
    
