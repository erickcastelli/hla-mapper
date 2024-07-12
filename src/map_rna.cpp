
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

#include "map_rna.hpp"
#include "external.hpp"
#include "functions.hpp"
#include "ThreadPool.hpp"
#include "preselect.hpp"

mutex mtx_hg38_select_rna;
mutex mtx_selection_general_rna;
mutex mtxg_rna;
unordered_map <string,int> adjusted_for_secondary_inside_genes;

vector <string> selected_results_r1_rna;
vector <string> selected_results_r2_rna;
vector <string> selected_results_trim_r1_rna;
vector <string> selected_results_trim_r2_rna;
//unordered_map <string,int> motifs_master_rna;



void read_sam_rna (string sam, string gene)
{
    ifstream input( sam.c_str() );
    for( std::string line; getline( input, line ); )
    {
        if (line.substr(0,1) == "@") {continue;}
        vector<string> sam_line;
        boost::split(sam_line,line,boost::is_any_of("\t"));
  
        vector<string> read_id_split;
        boost::split(read_id_split,sam_line[0],boost::is_any_of("$"));
        string read_id = read_id_split[0];
        pair <string,string> key = make_pair(gene,read_id);
        
        int mm = 0;

        if (sam_line[2] == "*") {
            
            if (sequence_list_r1.find(key) != sequence_list_r1.end())
            {
                sequence_list_r1[key] = sequence_list_r1[key] + sam_line[9].size();
            }
            else {
                sequence_list_r1[key] = sam_line[9].size();
            }
            
        }

        else {
            
            vector<string> nm_value;
            boost::split(nm_value,sam_line[12],boost::is_any_of(":"));
            mm = stoi(nm_value[2]);
            
            string cigar = splitcigar(sam_line[5]);
            vector<string> cigardata;
            boost::split(cigardata,cigar,boost::is_any_of(","));
            string c1 = cigardata[0];
            string c2 = cigardata[cigardata.size()-1];
            
            if ((c1 != "") && (mm != sam_line[9].size()))
            {
                if ((c1.substr(c1.size()-1,1) == "S") || (c1.substr(c1.size()-1,1) == "H"))
                {mm = mm + stoi(c1.substr(0,c1.size()));}
            }
            if ((c2 != "") && (mm != sam_line[9].size()))
            {
                if ((c2.substr(c2.size()-1,1) == "S") || (c2.substr(c2.size()-1,1) == "H"))
                {mm = mm + stoi(c2.substr(0,c2.size()));}
            }
            
            if (sequence_list_r1.find(key) != sequence_list_r1.end())
            {
                sequence_list_r1[key] = sequence_list_r1[key] + mm;
            }
            else {
                sequence_list_r1[key] = mm;
            }
        }
    }
    input.close();
    removefile(sam,v_debug);
}





void main_rna_map ()
{

    int v_check = 0;
    if (v_r0 != "") {v_check = 1;}
    if ((v_r1 != "") || (v_r2 != "")) {v_check = 1;}
    if (v_bam != "") {v_check = 1;}
    if (v_sample != "") {v_check = 1;}

    v_tolerance = 0.08;

    if ((v_db == "") || (v_check == 0))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::rna";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        v_message = "Usage:     hla-mapper rna r1=R1.gz r2=R2.gz sample=sample_name db=path <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "           hla-mapper rna bam=sample.bam sample=sample_name db=path <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "sample      sample name or identification", 1, v_quiet);
        screen_message (screen_size, 2, "db          path to an hla-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "r0          a single-ended fastq (fq, fastq, or gz, ignore r1/r2/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "r1          a paired-ended forward fastq (fq, fastq, or gz, ignore r0/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "r2          a paired-ended reverse fastq (fq, fastq, or gz, ignore r0/bam)", 1, v_quiet);
        screen_message (screen_size, 2, "bam         a BAM/CRAM file (ignore r0/r1/r2)", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "output      output folder", 1, v_quiet);
        screen_message (screen_size, 2, "threads     number of threads [default: " + v_threads + "]", 1, v_quiet);
//        screen_message (screen_size, 2, "bed         BED file to override the one in the database", 1, v_quiet);
        screen_message (screen_size, 2, "buffer      number of sequences in buffer [default: " + to_string(v_buffer) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "error       threshold for nucleotide quality trimming [default: " + to_string(v_mtrim_error).substr(0,4) + "]", 1, v_quiet);
        screen_message (screen_size, 2, "tolerance   fraction of mismatches allowed [Default: " + to_string(v_tolerance).substr(0,4) + ", from 0.01 to 0.10]", 1, v_quiet);
 //       screen_message (screen_size, 2, "downsample  downsampling for the typing procedure [default: " + to_string(downsampling) + "]", 1, v_quiet);

                        
        if (v_bwa != "") {v_message = "[found at " + v_bwa + "]";}
        if (v_bwa == "") {v_message = "(!!! bwa not detected !!!)";}
        v_message = "bwa         path to BWA " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        if (v_star != "") {v_message = "[found at " + v_star + "]";}
        if (v_star == "") {v_message = "(!!! STAR not detected !!!)";}
        v_message = "star        path to star " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);
        
        
        if (v_samtools != "") {v_message = "[found at " + v_samtools + "]";}
        if (v_samtools == "") {v_message = "[!!!Samtools not detected!!!]";}
        v_message = "samtools    path to Samtools " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        
        
        screen_message (screen_size, 2, "", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-unmapped       skip retrieving unmapped reads from BAM/CRAM", 1, v_quiet);
        screen_message (screen_size, 2, "--skip-adjust         skip the typing/adjust procedure [not recommended]", 1, v_quiet);
//        screen_message (screen_size, 2, "--low-mem             force low memory mode for sequence preselection", 1, v_quiet);
        screen_message (screen_size, 2, "--quiet               quiet mode", 1, v_quiet);
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
        main_rna_map();
        return;
    }
    
    if ((v_r1 != "") && (v_r1 == v_r2)) {warnings.push_back("r1 and r2 must be different fastq files"); main_rna_map(); return;}

    if (v_bam != "") {
        ifstream file_bam(v_bam.c_str());
        if (!file_bam)
        {
            warnings.push_back ("Could not access the .bam file or file does not exist.");
            v_sample = "";
            v_r1 = "";
            v_r2 = "";
            v_bam = "";
            main_rna_map();
            return;
        }
        file_bam.close();
    }

    if (v_r1 != "") {
        ifstream file_r1(v_r1.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Could not access the r1 file or file does not exist.");
            v_sample = "";
            v_r1 = "";
            main_rna_map();
            return;
        }
        file_r1.close();
    }

    if (v_r0 != "") {
        ifstream file_r1(v_r0.c_str());
        if (!file_r1)
        {
            warnings.push_back ("Could not access the r0 file or file does not exist.");
            v_sample = "";
            v_r1 = "";
            v_r0 = "";
            main_rna_map();
            return;
        }
        file_r1.close();
    }
    
    
    if (v_r2 != "") {
        ifstream file_r2(v_r2.c_str());
        if (!file_r2)
        {
            warnings.push_back ("Could not access the r2 file or file does not exist.");
            v_sample = "";
            v_r2 = "";
            main_rna_map();
            return;
        }
        file_r2.close();
    }
    
    
    if (((ends_with(v_r1,".fq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fq")) && (ends_with(v_r1,".gz"))))
    {
        
        warnings.push_back ("hla-mapper cannot deal with .fastq and .gz at the same time.");
        v_sample = "";
        v_r1 = "";
        main_rna_map();
        return;
    }
    
    if (((ends_with(v_r1,".fastq")) && (ends_with(v_r2,".gz"))) || ((ends_with(v_r2,".fastq")) && (ends_with(v_r1,".gz"))))
    {
        
        warnings.push_back ("hla-mapper cannot deal with .fastq and .gz at the same time.");
        v_sample = "";
        v_r1 = "";
        main_rna_map();
        return;
    }
    
    if (((v_r1 != "") && (v_r2 == "")) || ((v_r1 == "") && (v_r2 != "")))
    {
        warnings.push_back ("You must indicate r1 and r2 when using paired-end sequencing data.");
        v_sample = "";
        v_r1 = "";
        main_rna_map();
        return;
    }
    
    // checking database
    string v_db_info = v_db + "/db_rna.info";
    boost::replace_all(v_db_info, "\\ ", " ");

    if (! fileExists(v_db_info))
    {
        v_message = "Could not access database " + v_db_info;
        cout << v_message << endl;
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_rna_map();
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
        main_rna_map();
        return;
    }
    
    
    if (v_bam != "")
    {
        if (v_bed != "")
        {
            ifstream file_bed(v_bed.c_str());
            if (!file_bed)
            {
                warnings.push_back ("Could not access the BED file or file does not exist.");
                v_sample = "";
                v_r1 = "";
                v_bam = "";
                v_r2 = "";
                main_rna_map();
                return;
            }
            file_bed.close();
        }
    }
    
    
    
    
   v_rnaseq = 1;
   v_callint = 1;
    
    // criando output folder
    if (v_output == "")
    {
        if (v_bam != "") {v_output = findfilepath(v_bam) + "hla-mapper/";}
        if (v_r0 != "") {v_output = findfilepath(v_r0) + "hla-mapper/";}
        if (v_r1 != "") {v_output = findfilepath(v_r1) + "hla-mapper/";}
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
    v_message = Program_name + "::rna";
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


    v_rnaseq = 1;
    main_preselect();
    if (v_callint == 2) {return;}
 
    
    
    
    
    
    
    
    screen_message (screen_size, 0, "Sorting sequences ...", 2, v_quiet);

    int loop = stoi(v_threads);
    boost::algorithm::erase_all(v_genes, " ");
    boost::split(v_gene_list,v_genes,boost::is_any_of(","));

    string selectcleanR1 = v_output + v_sample + "selected.trim.R1.fastq";
    string selectcleanR2 = v_output + v_sample + "selected.trim.R2.fastq";
    if (v_map_type == "single") {selectcleanR1 = v_output + v_sample + "selected.trim.R0.fastq";}
    
    selected_results_r1_rna.clear();
    selected_results_r2_rna.clear();
    selected_results_trim_r1_rna.clear();
    selected_results_trim_r2_rna.clear();
    
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

        string motifdbfile = v_db + "/motif/rna/loci/" + *i + ".txt";
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
            string motifdbfile = v_db + "/local/rna/" + *i + "/" + *i + ".txt";
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
                                  mtx_selection_general_rna.lock();
                                  select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
                                  select_data_sorted_r2.push_back(id + "\n" + r2seq + "\n+\n" + r2qual);
                                  mtx_selection_general_rna.unlock();
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
                                     mtx_selection_general_rna.lock();
                                     select_data_sorted_r1.push_back(id + "\n" + r1seq + "\n+\n" + r1qual);
                                     mtx_selection_general_rna.unlock();
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
    
    

    
    
    
    
    
    
    
//Detecting splicing sites
    ThreadPool poolsplice(stoi(v_threads));
    std::vector< std::future<int> > results_splice;
    int count = 0;
    
    v_message = "Detecting splicing sites ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);


    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i == "") {continue;}
        string gene = *i;
        string type = v_map_type;
        
        count++;

        results_splice.emplace_back(
            poolsplice.enqueue([count, gene, v_sample_sub, type]
           {
            string v_reference = "'" + v_db + "/reference/rna/loci/" + gene + "/'";
            boost::replace_all(v_reference, "\\ ", " ");
 
            string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
            string in0 = v_output + v_sample + "sorted_" + gene + "_R0.fastq";
            string in1 = v_output + v_sample + "sorted_" + gene + "_R1.fastq";
            string in2 = v_output + v_sample + "sorted_" + gene + "_R2.fastq";
            string in_spliced = v_output + v_sample + "sorted_spliced_" + gene + ".fastq";
            string v_log = v_output + "log/" + v_sample + gene + ".splice.log";
            
            string mnmax = "70";
            string readLmax = "1";
            string mismatchN = "999";
            string seed = "60";
            
            if (((gene == "HLA-G") || (gene == "HLA-E")) || (gene == "HLA-F"))
            {
                mnmax = "30";
                readLmax = "0.1";
                mismatchN = "999";
                seed = "50";
            }

            if (((gene == "HLA-DMA") || (gene == "HLA-DMB")) || (gene == "HLA-DOA"))
            {
                mnmax = "30";
                readLmax = "0.1";
                mismatchN = "10";
                seed = "50";
            }
     
            if (((gene == "HLA-DOB") || (gene == "MICA")) || (gene == "MICB"))
            {
                mnmax = "30";
                readLmax = "0.1";
                mismatchN = "10";
                seed = "50";
            }
    
            if (((gene == "TAP1") || (gene == "TAP2")) || (gene == "HLA-H"))
            {
                mnmax = "50";
                readLmax = "0.1";
                mismatchN = "999";
                seed = "50";
            }

            if ((gene == "HLA-DRA") || (gene == "HLA-DPA1"))
            {
                mnmax = "50";
                readLmax = "0.1";
                mismatchN = "999";
                seed = "50";
            }
            
            
            string rg_string = "ID:" + v_sample_sub + " LB:" + v_sample_sub + " SM:" + v_sample_sub;
            string cmd = "";
            
            if (type == "paired") {
                cmd = v_star + " --runThreadN 1 --genomeDir " + v_reference + " --readFilesIn " + in1 + " " + in2 + " --outFileNamePrefix " + v_output + v_sample + gene + "_splice_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " >" + v_log;
            }
            if (type == "single") {
                cmd = v_star + " --runThreadN 1 --genomeDir " + v_reference + " --readFilesIn " + in0 + " --outFileNamePrefix " + v_output + v_sample + gene + "_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " >" + v_log;
            }
            string out = GetStdoutFromCommand(cmd);

            string v_sam = v_output + v_sample + gene + "_splice_Aligned.out.sam";
            fstream sam (v_sam.c_str());

            ofstream sorted_spliced;
            sorted_spliced.open (in_spliced.c_str());
        
            for( std::string line; getline( sam, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector<string> sam_line;
                boost::split(sam_line,line,boost::is_any_of("\t"));
                string cigar = sam_line[5];
                if (cigar.find('N') != std::string::npos)
                {
                    
                    string split = splitNcigar(cigar);
                    vector<string> cigar_N_splited;
                    boost::split(cigar_N_splited,split,boost::is_any_of(","));
                    int start = 0;
                    int block = 1;
                    for (auto item : cigar_N_splited)
                    {
                        int read_size = 0;
                        string subcigar = splitcigar(item);
                        vector<string> cigar_N_splited_sub;
                        boost::split(cigar_N_splited_sub,subcigar,boost::is_any_of(","));
                        for (auto sub : cigar_N_splited_sub)
                        {
                            if (sub == "") {continue;}
                            string type = sub.substr(sub.length() - 1,1);
                            string size = sub.substr(0, sub.length() - 1);
                            if (type == "N") {continue;}
                            if (type == "D") {continue;}
                            read_size = read_size + stoi(size);
                        }
                        string subread = sam_line[9].substr(start,read_size);
                        string subqual = sam_line[10].substr(start,read_size);
                        
                        sorted_spliced << "@" << sam_line[0] << "$" << sam_line[1] << "$" << block  << endl;
                        sorted_spliced << subread << endl;
                        sorted_spliced << "+" << endl;
                        sorted_spliced << subqual << endl;
                        
                        start = read_size;
                        block++;
                    }
                    
                    
                }
                else {
                    sorted_spliced << "@" << sam_line[0] << "$" << sam_line[1] << "$1"  << endl;
                    sorted_spliced << sam_line[9] << endl;
                    sorted_spliced << "+" << endl;
                    sorted_spliced << sam_line[10] << endl;
                }
            }
            sam.close();
            sorted_spliced.close();
 
            return 1;
          })
        );
    }
    for(auto && result: results_splice){result.get();} // waiting for all threads
    v_message = "Detecting splicing sites: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    screen_message (screen_size, 0, "Scoring sequences ...", 2, v_quiet);
    
    boost::algorithm::erase_all(v_genes, " ");
    boost::split(v_gene_list,v_genes,boost::is_any_of(","));

    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {

        if (*i ==  "") {continue;}
        
        screen_message (screen_size, 0, "Scoring sequences for " + *i + " ...", 2, v_quiet);
         
        string v_ref = "'" + v_db + "/mapper/rna/" + *i + "/" + *i + ".fas' ";
        boost::replace_all(v_ref, "\\ ", " ");

        string v_sam = v_output + *i + ".mapper.sam";
        string v_sort = v_output + v_sample + "sorted_spliced_" + *i + ".fastq";
        string v_sai = v_output + "tmp.sai ";
        string v_log = v_output + "/log/" + *i + "_mapper_aln.log";

        v_command = v_bwa + " aln -o " + to_string(v_mm_open) + " -t " + v_threads + " -n " + to_string(v_mm_max) + " " + v_ref + v_sort + " > " + v_sai + " 2>" + v_log;
        system (v_command.c_str());
        
        v_log = v_output + "/log/" + *i + "_mapper_sam.log";
        v_command = v_bwa + " samse " + v_ref + v_sai + v_sort + " > " + v_sam + " 2>" + v_log;
        system (v_command.c_str());
        
        v_command = " rm " + v_sai;
        system (v_command.c_str());

        read_sam_rna(v_sam, *i);
    }
    
    v_message = "Scoring reads considering other sequences ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);

     
    string v_ref = "'" + v_db + "/others/rna/hla_gen_others.fasta' ";
    boost::replace_all(v_ref, "\\ ", " ");

    string v_sam1 = v_output + "others.sam";
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
    read_sam_rna(v_sam1, "others");
    read_sam_rna(v_sam2, "others");

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
    std::map <string, string> others_list;
    ThreadPool poolscore(stoi(v_threads));
    std::vector< std::future<int> > resultscore;
    int readcount = 0;
    vector <string> table;
    
    
     for(map<string,int>::iterator it = sequence_size_r1.begin(); it != sequence_size_r1.end(); ++it)
     {
         string id = it->first;
         
         resultscore.emplace_back(
         poolscore.enqueue([readcount, id, v_target_list_sub, &address_list_r1, &table, &progress, totalreads]
           {
         
             map <string,int> nm_sum;
             string table_str = id;
             
             int min_score_targets = 10;
             int min_score_others = 30;
         
             int invalid_read = 1000;
             string address_to = "";
             int count = 0;
         
             int min_score = invalid_read;
             for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
             {
                 if (*g == "") {continue;}
             
                 nm_sum[*g] = invalid_read;
                 pair <string,string> key = make_pair (*g,id);
             
                 string scores = "-";
             
                 if (sequence_list_r1.find(key) == sequence_list_r1.end()) {nm_sum[*g] = invalid_read;}
                 
                 if (sequence_list_r1.find(key) != sequence_list_r1.end()){
                     nm_sum[*g] = sequence_list_r1[key];
                 }
                 scores = to_string(nm_sum[*g]);
                 if (scores == to_string(invalid_read)) {scores = "-";}
                 table_str.append ("\t" + scores);
                 if (scores != "-") {count++;}
                 if (nm_sum[*g] < min_score) {min_score = nm_sum[*g];}
             }
  
             
             if (count >= 1)
             {
                 int count_min_nm_sum = 0;
                 for( std::vector<string>::const_iterator g = v_target_list_sub.begin(); g != v_target_list_sub.end(); ++g)
                 {
                     if (*g == "") {continue;}
                     
                     if ((*g != "others") && (nm_sum[*g] == min_score))
                     {
                         float max = (sequence_size_r1[id]*2) * v_tolerance;
//                         if (min_score < min_score_targets)
                         if (min_score <= max)
                         {
                             count_min_nm_sum++;
                             address_to = address_to + "," + *g;
                             pair <string,string> key = make_pair (*g,id);
                             mtxg_rna.lock();
                             sequence_address[key] = 1;
                             mtxg_rna.unlock();
                         }
                     }
                     if ((*g == "others") && (nm_sum[*g] == min_score))
                     {
                         if (min_score < min_score_others)
                         {
                             count_min_nm_sum++;
                             address_to = address_to + "," + *g;
                             pair <string,string> key = make_pair (*g,id);
                             mtxg_rna.lock();
                             sequence_address[key] = 1;
                             mtxg_rna.unlock();
                         }
                     }
                     
                 }
                 
                 if (address_to != "") {
                     mtxg_rna.lock();
                     address_list_r1[id] = address_to.substr(1);
                     mtxg_rna.unlock();
                 }
                 
             }
             string sub = "";
             if (address_to == "") {sub = "-";}
             else {sub = address_to.substr(1);}
             table_str.append("\t" + sub);
             table_str.append("\n");
 
          
          mtxg_rna.lock();
             table.push_back(table_str);
             progress++;
             float ratio = ((progress / totalreads) * 100);
             v_message = "Comparing scores ... " + to_string(ratio) + " % ";
             screen_message (screen_size, 0, v_message, 2, v_quiet);
             mtxg_rna.unlock();
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
    
    
    for (auto & item : address_list_r1)
    {
       if (item.second.find("others") != string::npos) {
           others_list[item.first] = 1;
       }
    }
    
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
           
           string sort = "";
           if (v_map_type == "paired") {sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";}
           if (v_map_type == "single") {sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";}
           
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

            string sort = "";
            if (v_map_type == "paired") {sort = v_output + v_sample + "sorted_" + *i + "_R1.fastq";}
            if (v_map_type == "single") {sort = v_output + v_sample + "sorted_" + *i + "_R0.fastq";}
            
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
                original_data_r1_trim[read_id[0]] = seq + "\n" + info + "\n" + qual;
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
                    fq1 = v_output + v_sample + *i + "_R1.fastq";
                    fq2 = v_output + v_sample + *i + "_R2.fastq";
                }
            
                if (v_map_type == "single")
                {
                    fq1 = v_output + v_sample + *i + "_R0.fastq";
                }
            
                string result = typing_rna(gene, fq1, fq2);
                if (v_debug == 1) {cout << result << endl;}
                mtx_selection_general_rna.lock();
                
                if (result == "")
                {
                typing_result_alleleA[gene] = "";
                typing_result_alleleA[gene] = "";
                mtx_selection_general_rna.unlock();
                continue;
                }
            
                if (result.substr(0,4) == "fail")
                {
                    typing_result_alleleA[gene] = "";
                    typing_result_alleleA[gene] = "";
                    mtx_selection_general_rna.unlock();
                    continue;
                }
            
                vector <string> results;
                boost::split(results,result,boost::is_any_of(","));
                vector <string> alleles;
                boost::split(alleles,results[0],boost::is_any_of("\t"));
                typing_result_alleleA[gene] = alleles[0];
                typing_result_alleleB[gene] = alleles[1];
                mtx_selection_general_rna.unlock();
        }
        
        general_log << endl;
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i == "") {continue;}
 
            if (typing_result_alleleA[*i] != "") {
                general_log << "List of sequences used for alignment adjustment " << *i << ": " << typing_result_alleleA[*i] << ", " << typing_result_alleleB[*i] << endl;
            }
            else {
                if ((gene_type[*i] == "TYPE") || (v_typeall == 1))
                {general_log << "Adjustment" << *i << ": failed" << endl;}
                if (gene_type[*i] != "TYPE")
                {general_log << "Adjustment" << *i << ": disabled" << endl;}
            }

        }
        general_log << endl;
        
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
    

    
/*
 // Mapping
     
     
     for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
     {
         
         if (*i == "") {continue;}
         v_message = "Mapping " + *i + " ...";
         screen_message (screen_size, 0, v_message, 2, v_quiet);
         
             string v_reference = "'" + v_db + "/reference/rna/loci/" + *i + "/'";
             boost::replace_all(v_reference, "\\ ", " ");

             string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
             string in0 = v_output + v_sample + *i + "_R0.fastq";
             string in1 = v_output + v_sample + *i + "_R1.fastq";
             string in2 = v_output + v_sample + *i + "_R2.fastq";
             v_sam1 = v_output + v_sample + *i + ".tmp.sam";
             v_sam2 = v_output + v_sample + *i + ".unique.sam";
             string v_sam3 = v_output + v_sample + *i + ".adjusted.sam";
             v_log = v_output + "log/" + v_sample + *i + ".map.log";
             
             int v_pos_correct = position_hg38[*i] - 1;
             
             string mnmax = "70";
             string readLmax = "1";
             string mismatchN = "999";
             string seed = "60";


             if (((*i == "HLA-G") || (*i == "HLA-E")) || (*i == "HLA-F"))
             {
                 mnmax = "10";
                 readLmax = "0.1";
                 mismatchN = "999";
                 seed = "50";
             }

             if (((*i == "HLA-DMA") || (*i == "HLA-DMB")) || (*i == "HLA-DOA"))
             {
                 mnmax = "20";
                 readLmax = "0.1";
                 mismatchN = "10";
                 seed = "50";
             }
          
             if (((*i == "HLA-DOB") || (*i == "MICA")) || (*i == "MICB"))
             {
                 mnmax = "20";
                 readLmax = "0.1";
                 mismatchN = "10";
                 seed = "50";
             }
         
             if (((*i == "TAP1") || (*i == "TAP2")) || (*i == "HLA-H"))
             {
                 mnmax = "50";
                 readLmax = "0.1";
                 mismatchN = "999";
                 seed = "50";
             }
   
             if ((*i == "HLA-DRA") || (*i == "HLA-DPA1"))
             {
                 mnmax = "50";
                 readLmax = "0.1";
                 mismatchN = "999";
                 seed = "50";
             }
         
             string rg_string = "ID:" + v_sample_sub + " LB:" + v_sample_sub + " SM:" + v_sample_sub;
         
         if (v_map_type == "paired") {
                 v_command = v_star + " --runThreadN " + v_threads + " --genomeDir " + v_reference + " --readFilesIn " + in1 + " " + in2 + " --outFileNamePrefix " + v_output + v_sample + *i + "_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " --peOverlapNbasesMin 1 >" + v_log;
             }
         

         if (v_map_type == "single") {
                 v_command = v_star + " --runThreadN " + v_threads + " --genomeDir " + v_reference + " --readFilesIn " + in0 + " --outFileNamePrefix " + v_output + v_sample + *i + "_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " --peOverlapNbasesMin 1 >" + v_log;
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
             
             v_sam1 = v_output + v_sample + *i + "_Aligned.out.sam";
              
             fstream sam (v_sam1.c_str());
             ofstream OUT1;
             OUT1.open (v_sam2.c_str());
             ofstream OUT2;
             OUT2.open (v_sam3.c_str());
         
              
             int count_reads = 0;
             v_pos_correct = position_hg38[*i] - 1;
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
                         OUT1 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
                     }
                     OUT1 << "@CO\thla-mapper " << Program_version << ", human genome version hg38 " << endl;
                     
                     if (usechr == 0) {
                         OUT2 << "@SQ\tSN:" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
                     }
                     if (usechr == 1) {
                         OUT2 << "@SQ\tSN:" << "chr" << chr_hg38[*i] << "\tLN:" << chr_size_hg38[*i] << endl;
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
                 if (usechr == 1) {sam_line[2] = "chr" + chr_hg38[*i];}
                 
              
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

                         if (sam_line[1] == "163") {sam_line[1] = "419";}
                         if (sam_line[1] == "83") {sam_line[1] = "339";}
                         if (sam_line[1] == "147") {sam_line[1] = "403";}
                         if (sam_line[1] == "99") {sam_line[1] = "355";}
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

                     
                         sam_line[4] = "3";
                          
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
             log_mapping << endl;
              
              
             log_mapping << "List of reads adjusted for primary alignment after finding the best allele combination" << endl;
             for (auto & item : adjusted_for_primary)
             {
                  log_mapping << item.first << endl;
             }
             log_mapping << endl;
              
             read_count_gene[*i] = double(count_reads);
             log_mapping.close();
         
     }
     v_message = "Mapping: done";
     screen_message (screen_size, 0, v_message, 1, v_quiet);
 */
    
    
    
    
    


    
    
// Mapping

    ThreadPool poolmap(stoi(v_threads));
    std::vector< std::future<int> > resultsmap;
    int count_map = 0;
    v_message = "Mapping reads to gene-specific references ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);


    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        
        if (*i == "") {continue;}
        string gene = *i;
        string type = v_map_type;
        
        count_map++;
        resultsmap.emplace_back(
                               poolmap.enqueue([count_map, gene, v_sample_sub, type, &address_list_r1]
                {
        
                string v_reference = "'" + v_db + "/reference/rna/loci/" + gene + "/'";
                boost::replace_all(v_reference, "\\ ", " ");

                string v_tag_sub = v_sample.substr(0, v_sample.size()-1);
                string in0 = v_output + v_sample + gene + "_R0.fastq";
                string in1 = v_output + v_sample + gene + "_R1.fastq";
                string in2 = v_output + v_sample + gene + "_R2.fastq";
                string v_sam1 = v_output + v_sample + gene + ".tmp.sam";
                string v_sam2 = v_output + v_sample + gene + ".unique.sam";
                string v_sam3 = v_output + v_sample + gene + ".adjusted.sam";
                string v_log = v_output + "log/" + v_sample + gene + ".map.log";
            
                int v_pos_correct = position_hg38[gene] - 1;
            
                string mnmax = "70";
                string readLmax = "1";
                string mismatchN = "999";
                string seed = "60";


                if (((gene == "HLA-G") || (gene == "HLA-E")) || (gene == "HLA-F"))
                {
                    mnmax = "30";
                    readLmax = "0.1";
                    mismatchN = "999";
                    seed = "50";
                }

                if (((gene == "HLA-DMA") || (gene == "HLA-DMB")) || (gene == "HLA-DOA"))
                {
                    mnmax = "30";
                    readLmax = "0.1";
                    mismatchN = "10";
                    seed = "50";
                }
         
                if (((gene == "HLA-DOB") || (gene == "MICA")) || (gene == "MICB"))
                {
                    mnmax = "30";
                    readLmax = "0.1";
                    mismatchN = "10";
                    seed = "50";
                }
        
                if (((gene == "TAP1") || (gene == "TAP2")) || (gene == "HLA-H"))
                {
                    mnmax = "50";
                    readLmax = "0.1";
                    mismatchN = "999";
                    seed = "50";
                }
  
                if ((gene == "HLA-DRA") || (gene == "HLA-DPA1"))
                {
                    mnmax = "50";
                    readLmax = "0.1";
                    mismatchN = "999";
                    seed = "50";
                }
        
                string rg_string = "ID:" + v_sample_sub + " LB:" + v_sample_sub + " SM:" + v_sample_sub;
                string cmd = "";
                if (type == "paired") {
                    cmd = v_star + " --runThreadN 1 --genomeDir " + v_reference + " --readFilesIn " + in1 + " " + in2 + " --outFileNamePrefix " + v_output + v_sample + gene + "_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " >" + v_log;
                }
        

                if (type == "single") {
                    cmd = v_star + " --runThreadN 1 --genomeDir " + v_reference + " --readFilesIn " + in0 + " --outFileNamePrefix " + v_output + v_sample + gene + "_ --seedPerWindowNmax " + seed + " --outFilterMismatchNmax " + mismatchN + " --outFilterMultimapNmax " + mnmax + " --outFilterMismatchNoverReadLmax " + readLmax + " --outSAMunmapped Within --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --outSAMattrRGline " + rg_string + " >" + v_log;
                    }
            GetStdoutFromCommand(cmd);


            ifstream reference (v_reference.c_str());
            string refseq = "";
            for( std::string line; getline( reference, line ); )
            {
                if (line.substr(0,1) == ">") {continue;}
                if (line == "") {continue;}
                refseq = refseq + line;
            }
            reference.close();
            mtxg_rna.lock();
            ref_size[gene] = refseq.length();
            mtxg_rna.unlock();
            refseq = "";
            
            v_sam1 = v_output + v_sample + gene + "_Aligned.out.sam";
             
            fstream sam (v_sam1.c_str());
            ofstream OUT1;
            OUT1.open (v_sam2.c_str());
            ofstream OUT2;
            OUT2.open (v_sam3.c_str());
        
            int count_reads = 0;
            v_pos_correct = position_hg38[gene] - 1;
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
                        OUT1 << "@SQ\tSN:" << chr_hg38[gene] << "\tLN:" << chr_size_hg38[gene] << endl;
                    }
                    if (usechr == 1) {
                        OUT1 << "@SQ\tSN:" << "chr" << chr_hg38[gene] << "\tLN:" << chr_size_hg38[gene] << endl;
                    }
                    OUT1 << "@CO\thla-mapper " << Program_version << ", human genome version hg38 " << endl;
                    
                    if (usechr == 0) {
                        OUT2 << "@SQ\tSN:" << chr_hg38[gene] << "\tLN:" << chr_size_hg38[gene] << endl;
                    }
                    if (usechr == 1) {
                        OUT2 << "@SQ\tSN:" << "chr" << chr_hg38[gene] << "\tLN:" << chr_size_hg38[gene] << endl;
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

                sam_line[2] = chr_hg38[gene];
                if (usechr == 1) {sam_line[2] = "chr" + chr_hg38[gene];}
                
             
                string tmp = address_list_r1[sam_line[0]] + ",";
                vector<string> address_multi_hits;
                boost::split(address_multi_hits,tmp,boost::is_any_of(","));
                
                if (address_multi_hits.size() < 3)
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

                
                
//                cout << reads_mapped_nm["SAMPLE_0001.622.A*68.01.01.04"] << endl;
//                cout << reads_mapped_gene["SAMPLE_0001.622.A*68.01.01.04"] << endl;

                 
                if (address_multi_hits.size() >= 3)
                {

                    int adjust = 0;
                    
                    pair <string,string> key = make_pair(gene,sam_line[0]);
                    if (reads_mapped_gene.find(sam_line[0]) != reads_mapped_gene.end())
                    {
                        if (reads_mapped_gene[sam_line[0]] == gene)
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

                        if (sam_line[1] == "163") {sam_line[1] = "419";}
                        if (sam_line[1] == "83") {sam_line[1] = "339";}
                        if (sam_line[1] == "147") {sam_line[1] = "403";}
                        if (sam_line[1] == "99") {sam_line[1] = "355";}
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

                    
                        sam_line[4] = "3";
                         
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

            string log = v_output + v_sample + gene + ".log";
            ofstream log_mapping;
            log_mapping.open (log.c_str());

            log_mapping << endl << "List of reads marked as secondary" << endl;
            for (auto & item : adjusted_for_secondary)
            {
                 log_mapping << item.first << endl;
            }
            log_mapping << endl;
             
             
            log_mapping << "List of reads adjusted for primary" << endl;
            for (auto & item : adjusted_for_primary)
            {
                 log_mapping << item.first << endl;
            }
            log_mapping << endl;
             
            read_count_gene[gene] = double(count_reads);
            log_mapping.close();
            
            return 1;
        })
        );
    }
    for(auto && result: resultsmap){result.get();} // waiting for all threads
    v_message = "Mapping reads to gene-specific references: done";
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
    
    unordered_map <string,int> optimized_reads;
    vector <string> adjusted_data;
    vector <string> unique_data;
    vector <string> hg38_data;
    
    
    v_message = "Merging files ... reading new alignments ...";
    screen_message (screen_size, 0, v_message, 2, v_quiet);


    if (v_bam == "")
    {
        vector <string> adjusted_data;
        vector <string> unique_data;
        vector <string> head;
        string v_tag_sub = v_sample.substr(0, v_sample.size()-1);

        
        head.push_back("@HD\tVN:1.4\n");
        head.push_back("@SQ\tSN:6\tLN:170805979\n");
        head.push_back("@CO\thla-mapper " + Program_version + ", human genome version hg38\n");
        head.push_back("@RG\tID:" + v_tag_sub + "\tLB:" + v_tag_sub + "\tSM:" + v_tag_sub + "\tPL:illumina\tPU:" + v_tag_sub + "\n");
        
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".adjusted.sam";
            ifstream input( file.c_str() );
             for( std::string line; getline( input, line ); )
             {
                 if (line == "") {continue;}
                 if (line.substr(0,1) == "@") {continue;}
                 adjusted_data.push_back(line);
             }
             input.close();

            file = v_output + v_sample + *i + ".unique.sam";
            ifstream input2( file.c_str() );
             for( std::string line; getline( input2, line ); )
             {
                 if (line == "") {continue;}
                 if (line.substr(0,1) == "@") {continue;}
                 adjusted_data.push_back(line);
             }
             input2.close();
        }
        

        string out1 = v_output + "/" + v_sample_sub + ".adjusted.tmp.sam";
        ofstream outmerged;
        outmerged.open (out1.c_str());
        for (auto & item : head)
        {outmerged << item;}
        for (auto & item : adjusted_data)
        {outmerged << item << endl;}
        outmerged.close();
        
        string out2 = v_output + "/" + v_sample_sub + ".unique.tmp.sam";
        ofstream outmergedunique;
        outmergedunique.open (out2.c_str());
        for (auto & item : head)
        {outmergedunique << item;}
        for (auto & item : adjusted_data)
        {outmergedunique << item << endl;}
        outmergedunique.close();
        

        string outsort = v_output + "/" + v_sample_sub + ".adjusted.bam";
        string log = v_output + "/log/merge.log";
        v_command = v_samtools + " sort " + out1 + " > " + outsort +  " 2>" + log;
        system (v_command.c_str());

        v_command = v_samtools + " index " + outsort;
        system (v_command.c_str());

        removefile(out1, v_debug);
        
        outsort = v_output + "/" + v_sample_sub + ".unique.bam";
        v_command = v_samtools + " sort " + out2 + " > " + outsort  +  " 2>" + log;;
        system (v_command.c_str());

        v_command = v_samtools + " index " + outsort;
        system (v_command.c_str());

        removefile(out2, v_debug);
    }
    
    
    if (v_bam != "")
    {
        
        for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
        {
            if (*i ==  "") {continue;}
            string file = v_output + v_sample + *i + ".adjusted.sam";
            ifstream input( file.c_str() );
            for( std::string line; getline( input, line ); )
            {
                if (line == "") {continue;}
                if (line.substr(0,1) == "@") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                if (data[2] == "*") {continue;}
                if (data[3] == "0") {continue;}
                
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
                
                optimized_reads[data[0]] = 1;
                adjusted_data.push_back(line);
            }
            input.close();
        }
        
        ofstream secondary_warning;
        string secondary_file = v_output + v_sample + "reads_overlapping_genes_but_optimization_failed.txt";
        secondary_warning.open (secondary_file.c_str());
        
        
        v_message = "Merging files ... reading original alignments ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);
        string hg38_not_target = v_output + v_sample + "hg38_not_target.sam";
        string hg38_target = v_output + v_sample + "hg38_target.sam";
        string hg38_unmapped = v_output + v_sample + "hg38_unmapped.sam";
        string v_log = v_output + "/log/hg38.view.log";
        string target = v_db + "/bed/target_rna.bed";
        string nottarget = v_db + "/bed/not_target_rna.bed";
        
        v_command = v_samtools + " view -@ " + v_threads + " -h -ML " + nottarget + " " + v_bam + " > " + hg38_not_target + " 2> " + v_log;
        system(v_command.c_str());

        v_command = v_samtools + " view -@ " + v_threads + " -h -ML " + target + " " + v_bam + " > " + hg38_target + " 2> " + v_log;
        system(v_command.c_str());

        v_command = v_samtools + " view -@ " + v_threads + " -h -f 4 " + v_bam + " > " + hg38_unmapped + " 2> " + v_log;
        system(v_command.c_str());

        
        ifstream input( hg38_target.c_str() );
        ofstream appendoriginal;
        appendoriginal.open(hg38_not_target, std::ios_base::app);
        
        for( std::string line; getline( input, line ); )
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == "@") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            if (optimized_reads.find(data[0]) != optimized_reads.end()) {continue;}
            if (optimized_reads.find(data[0]) == optimized_reads.end())
            {
                int secondary = 0;
                int others_detect = 0;
                
                if (others_list.find(data[0]) != others_list.end()) {others_detect = 1;}

                if (others_detect == 0) {
                    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
                    {
                        if (*i ==  "") {continue;}
                        if (((stoi(data[3]) > gene_opt_start[*i]) && (stoi(data[3]) < gene_opt_end[*i])) && (data[2] == chr_hg38[*i]))
                        {
                            secondary = 1;
                            break;
                        }
                    }
                }
                
                if ((secondary == 0) && (others_detect == 0)) { appendoriginal << line << endl; continue;}
  
                
                if ((secondary == 1) && (others_detect == 0))
                {
                    secondary_warning << data[0] << endl;

                    if (data[1] == "163") {data[1] = "419";}
                    if (data[1] == "83") {data[1] = "339";}
                    if (data[1] == "99") {data[1] = "355";}
                    if (data[1] == "147") {data[1] = "403";}
                    if (data[1] == "81") {data[1] = "337";}
                    if (data[1] == "167") {data[1] = "423";}
                    if (data[1] == "145") {data[1] = "401";}
                    if (data[1] == "97") {data[1] = "353";}
                    if (data[1] == "185") {data[1] = "441";}
                    if (data[1] == "73") {data[1] = "329";}
                    if (data[1] == "113") {data[1] = "369";}
                    if (data[1] == "117") {data[1] = "373";}
                    if (data[1] == "121") {data[1] = "377";}
                    if (data[1] == "133") {data[1] = "389";}
                    if (data[1] == "137") {data[1] = "393";}
                    if (data[1] == "177") {data[1] = "433";}
                    if (data[1] == "181") {data[1] = "437";}
                    
 //                   data[4] = "3";

                    adjusted_for_secondary_inside_genes[data[0]] = 1;
                    
                    string newline = data[0];
                    for (int helper = 1; helper < data.size(); helper++)
                    {newline = newline + "\t" + data[helper];}
                    appendoriginal << newline << endl;
                    continue;
                }
            }
        }
        input.close();
        
        
        
       ifstream unmap( hg38_unmapped.c_str() );
       for( std::string line; getline( unmap, line ); )
       {
           if (line == "") {continue;}
           if (line.substr(0,1) == "@") {continue;}
           vector <string> data;
           boost::split(data,line,boost::is_any_of("\t"));
           if (optimized_reads.find(data[0]) != optimized_reads.end()) {continue;}
           if (optimized_reads.find(data[0]) == optimized_reads.end())
           {
                   appendoriginal << line << endl;
                   continue;
           }
       }
       unmap.close();
        
        
        for (auto &item : adjusted_data) {appendoriginal << item << endl;}
        appendoriginal.close();

        float mem = ((v_memory / 3) * 2) / stoi(v_threads);
        if (mem > 3) {mem = 3;}
        int memint = int(mem);
        secondary_warning.close();
        v_message = "Merging files ... sorting and indexing BAM ...";
        screen_message (screen_size, 0, v_message, 2, v_quiet);

        
        string v_bam_out = v_output + v_sample_sub + ".adjusted.bam";
        v_log = v_output + "/log/sort.log";
        v_command = v_samtools + " sort -@ " + v_threads + " -m " + to_string(memint) + "g " + hg38_not_target  + " > " +  v_bam_out + " 2>" + v_log;
        system (v_command.c_str());
        v_command = v_samtools + " index " + v_bam_out + " 2>" + v_log;
        system (v_command.c_str());

     }
    v_message = "Merging files: done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

      
  
    
    
    v_message = "Cleaning temporary files ... ";
    screen_message (screen_size, 0, v_message, 2, v_quiet);
    
    removefile(v_output + v_sample + "hg38.selected.sam",v_debug);
    if (v_debug == 1) {keep_sam = 1;}
    removefile(v_output + v_sample + "hg38.sam",keep_sam);
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

    
    string hg38_not_target = v_output + v_sample + "hg38_not_target.sam";
    string hg38_target = v_output + v_sample + "hg38_target.sam";
    string hg38_unmapped_1 = v_output + v_sample + "hg38_unmapped.sam";
    string hg38_unmapped_2 = v_output + v_sample + "hg38.unmapped.sam";

    v_command = "rm -rf " + v_output + "/rna_type/";
    system (v_command.c_str());

    
    removefile(hg38_target,v_debug);
    removefile(hg38_not_target,v_debug);
    removefile(hg38_unmapped_1,v_debug);
    removefile(hg38_unmapped_2,v_debug);

    
    if (v_map_type == "paired") {
        v_command = "gzip -f " + v_output + v_sample + "selected.R1.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
        v_command = "gzip -f " + v_output + v_sample + "selected.R2.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
    }
    if (v_map_type == "single") {
        v_command = "gzip -f " + v_output + v_sample + "selected.R0.fastq";
        v_system_out = GetStdoutFromCommand(v_command);
    }
    
    
    for( std::vector<string>::const_iterator i = v_gene_list.begin(); i != v_gene_list.end(); ++i)
    {
        if (*i ==  "") {continue;}
        
        removefile(v_output + v_sample + *i + ".trim.R1.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R0.fastq",v_debug);
        removefile(v_output + v_sample + *i + ".trim.R2.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R0.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R1.fastq",v_debug);
        removefile(v_output + v_sample + "sorted_" + *i + "_R2.fastq",v_debug);
        removefile(v_output + v_sample + *i + "_Aligned.out.sam",v_debug);
        removefile(v_output + v_sample + *i + "_Log.out",v_debug);
        removefile(v_output + v_sample + *i + ".adjusted.sam",v_debug);
        removefile(v_output + v_sample + *i + ".unique.sam",v_debug);
        removefile(v_output + v_sample + *i + "_splice_Aligned.out.sam",v_debug);
        removefile(v_output + v_sample + "sorted_spliced_" + *i + ".fastq",v_debug);
        
        if (v_map_type == "paired") {
            v_command = "gzip -f " + v_output + v_sample + *i + "_R1.fastq";
            v_system_out = GetStdoutFromCommand(v_command);
            v_command = "gzip -f " + v_output + v_sample + *i + "_R2.fastq";
            v_system_out = GetStdoutFromCommand(v_command);
        }
        if (v_map_type == "single") {
            v_command = "gzip -f " + v_output + v_sample + *i + "_R0.fastq";
            v_system_out = GetStdoutFromCommand(v_command);
        }
            
        
    }
  
    if (v_debug != 1)
    {
        v_command = "rm " + v_output + "/*.out";
        system (v_command.c_str());
        v_command = "rm " + v_output + "/*.tab";
        system (v_command.c_str());
    }

    v_message = "Cleaning temporary files ... done";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    general_log.close();
    
    
    
}
    
      



       


    





    
    
    
    
    
    
