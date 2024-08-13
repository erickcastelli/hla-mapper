
//  hla-mapper
//
//  Created by Erick Castelli
//  Copyright © 2022 GeMBio.Unesp. All rights reserved.
//


#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <boost/algorithm/string.hpp>
#include <thread>
#include <pwd.h>

#include <sys/ioctl.h> //for screen_size
#include <stdio.h> //for screen_size
#include <unistd.h> //for screen_size

#include "external.hpp"
#include "preselect.hpp"
#include "functions.hpp"
#include "map_dna.hpp"
#include "typing.hpp"
#include "setup.hpp"
#include "map_rna.hpp"


using namespace std;

auto clock_start = std::chrono::steady_clock::now();

// max threads
int v_concurentThreadsSupported = std::thread::hardware_concurrency() / 2;
string v_threads = std::to_string(v_concurentThreadsSupported);


unsigned long v_memory = getTotalSystemMemory() / long(1024) / long(1024) / long(1024);

//Program identification
string Program_name = "hla-mapper";
string Program_company = "GeMBio/Unesp";
string Program_version = "4.5.1";
string Program_author = "Erick C. Castelli";
string Program_date = "Aug 13th 2024";
string Program_website = "www.castelli-lab.net/apps/hla-mapper";


map <string,int> samtools_versions;
map <string,int> bwa_versions;
map <string,int> whats_versions;
map <string,int> blast_versions;
map <string,int> freebayes_versions;
map <string,int> star_versions;



// Internal variables
string ostype = "";
int screen_size = 80;
string v_message = "";
string v_system_out = "";
vector<string> warnings;
string configfile = "";
vector<string> v_system_out_list;
string v_command;
vector <string> selected_reads;
int motif_size = 20;
int v_buffer = 1000000;
string v_genes = ""; // list of genes avaliable
vector<string> v_gene_list;
int v_mm_max = 20; //-n BWA ALN
int v_mm_open = 1; //-o BWA ALN
int v_size = 50;
double v_tolerance = 0.05; // tolerance
double v_mtrim_error = 0.08f; // error limit (mtrim)
int keep_sam = 0;
int minselect = 30;
string v_bed = "";
string v_map_type = "";
int downsampling = 30;
int v_rnaseq = 0;
int v_lowmem = 0;
int v_callint = 0;
int v_forceindex = 0;
int v_exome = 0;
int v_multiphasing = 1;
string v_target = "";

map <string,int> position_hg38;
map <string,string> chr_hg38;
map <string,string> chr_size_hg38;
map <string,int> gene_opt_start;
map <string,int> gene_opt_end;
map <string,string> gene_type;
map <string,int> ref_size;

map <pair<string,string>,int> reads_mapped_to_allele_A;
map <pair<string,string>,int> reads_mapped_to_allele_B;
map <pair<string,string>,int> reads_possibly_mapped;
map <string,int> reads_mapped_nm;
map <string,string> reads_mapped_gene;

double read_mean_size;
map <string,double> read_count_gene;

map <pair<string,string>, int> sequence_list_r1;
map <pair<string,string>, int> sequence_list_r2;
map <string, int> sequence_size_r1;
map <string, int> sequence_size_r2;

map <pair<string,string>, int> sequence_address;

int v_skiptyping = 0;
int v_typeall = 0;
int v_assembling = 0;
int v_intergenic = 0;
int v_quiet = 0;
int v_useconfig = 0;
int v_skip_unmapped = 0;
int v_use_local = 0;
int v_type_after_map = 0;


//Main parameters
string v_output = "";
string v_db = "";
string v_dbprepare = "";
string v_r0 = "";
string v_r1 = "";
string v_r2 = "";
string v_bam = "";
string v_mapout = "";
string v_sample = "";
string homedir = "";
int usechr = 0;


// Only for developers
int v_debug = 0; // defines the debug mode
int v_verbose = 0; // defines the verbose mode
int v_bypass_version_check = 0; //bypass bwa and samtools version


string v_bwa = "";
string v_samtools = "";
string v_whats = "";
string v_blast = "";
string v_star = "STAR";
string v_freebayes = "freebayes";
string v_bcftools = "bcftools";


void main_help(void)
{
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);

    v_message = "Program:  " + Program_name;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:  " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "Contact:  Erick C. Castelli <erick.castelli@unesp.br>", 1, v_quiet);

    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "Usage:    hla-mapper <command> [options]", 1, v_quiet);

    screen_message (screen_size, 0, "", 1, v_quiet);
    
    screen_message (screen_size, 0, "Commands: setup       configure hla-mapper", 1, v_quiet);
    screen_message (screen_size, 0, "          dna         map/align DNA sequences (WGS, WES, Amplicons)", 1, v_quiet);
//    screen_message (screen_size, 0, "          rna         (beta) map/align RNA sequences (RNASeq)", 1, v_quiet);
    screen_message (screen_size, 0, "          select      preselect sequences", 1, v_quiet);
//    screen_message (screen_size, 0, "          type        (beta) typing algorithm", 1, v_quiet);
 
    screen_message (screen_size, 0, "", 1, v_quiet);
    return;
}




void load_config (void)
{
    if (v_useconfig == 0) {
        ifstream config(configfile.c_str());
        if (config)
        {
            for( std::string line; getline( config, line ); )
            {
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                vector<string> v_set;
                boost::split(v_set,line,boost::is_any_of("="));
                if (v_set[0] == "bwa") {v_bwa = v_set[1];}
                if (v_set[0] == "samtools") {v_samtools = v_set[1];}
                if (v_set[0] == "star") {v_star = v_set[1];}
                if (v_set[0] == "db") {v_db = v_set[1] + "/";}
                if (v_set[0] == "whatshap") {v_whats = v_set[1];}
                if (v_set[0] == "blast") {v_blast = v_set[1];}
                if (v_set[0] == "freebayes") {v_freebayes = v_set[1];}
                if (v_set[0] == "star") {v_star = v_set[1];}
            }
        }
        config.close();
    }
    return;
}




int main(int argc, const char * argv[]) {

    
    
    #if defined (__linux__)
        ostype = "linux";
    #endif
    #if defined (__APPLE__)
        ostype = "mac";
    #endif
    
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    //screen_size = w.ws_col;
    
    
    if (stoi(v_threads) < 1) {v_threads = "2";}
    if (stoi(v_threads) > 10) {v_threads = "10";}

    
    
    
    
    samtools_versions["1.12"] = 1;
    samtools_versions["1.13"] = 1;
    samtools_versions["1.14"] = 1;
    samtools_versions["1.19"] = 1;
    samtools_versions["1.19.2"] = 1;
    samtools_versions["1.20"] = 1;

    bwa_versions["0.7.17"] = 1;
    bwa_versions["0.7.16"] = 1;
    
    blast_versions["2.12.0+"] = 1;
    blast_versions["2.11.0+"] = 1;
    blast_versions["2.15.0+"] = 1;

    whats_versions["2.3"] = 1;
    whats_versions["2.2"] = 1;

    freebayes_versions["v1.3.6"] = 1;

    star_versions["2.7.10a"] = 1;

    
    
        
    
    
    homedir = getpwuid(getuid())->pw_dir;
    configfile = homedir + "/.hla-mapper";
    
    
    if (fileExists(configfile))
    {
        load_config();
    }
    else
    {
        cout << endl << endl;
        cout << "Aparently this is the first time you use hla-mapper." << endl;
        cout << "Starting configuration ..." << endl;
        main_setup();
        return 0;
    }
    
    
    if (! fileExists(v_bwa)) {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0; }
    if (! fileExists(v_samtools)) {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0;}
    if (! fileExists(v_db)) {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0;}

    if (v_blast != "DISABLED") {
        if (! fileExists(v_blast))
        {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0;}
    }
            
    if (v_whats != "DISABLED") {
        if (! fileExists(v_whats)) {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0;}
    }
    
    if (v_freebayes != "DISABLED") {
        if (! fileExists(v_freebayes)) {cout << endl << "You need to run 'hla-mapper setup' again." << endl << endl; main_setup(); return 0;}
    }

    
    int command_ok = 1;
    
    if (argc > 2)
    {
       int a;
       for (a = 2; a < argc; a++)
       {
           string str = argv[a];
           
           if (str.find("--quiet") != string::npos)
           {
               v_quiet= 1;
               continue;
           }

           else if (str.find("--exome") != string::npos)
           {
               v_exome= 1;
               continue;
           }
           
           else if (str.find("r0=") != string::npos)
           {
               v_r0 = str.substr(3);
               if (! fileExists(v_r0))
               {
                   v_message = "Invalid r0: " + v_r0;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str.find("r1=") != string::npos)
           {
               v_r1 = str.substr(3);
               if (! fileExists(v_r1))
               {
                   v_message = "Invalid r1: " + v_r1;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str.find("r2=") != string::npos)
           {
               v_r2 = str.substr(3);
               if (! fileExists(v_r2))
               {
                   v_message = "Invalid r2: " + v_r2;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str.find("bam=") != string::npos)
           {
               v_bam = str.substr(4);
               if (! fileExists(v_bam))
               {
                   v_message = "Invalid bam: " + v_bam;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str.find("map=") != string::npos)
           {
               v_mapout = str.substr(4);
               if (! fileExists(v_mapout))
               {
                   v_message = "Invalid map folder: " + v_mapout;
                   v_mapout = "";
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str.find("bed=") != string::npos)
           {
               v_bed = str.substr(4);
               if (! fileExists(v_bed))
               {
                   v_message = "Invalid bed: " + v_bed;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           else if (str.find("threads=") != string::npos)
           {
               v_threads = str.substr(8);
               if (std::stoi(v_threads) < 1) {v_threads = 1;}
               continue;
           }

           else if (str.find("db=") != string::npos)
           {
               v_db = str.substr(3);
               char last_ch = v_db.back();
               char ch = '/';
               if (last_ch != ch)
               {
                   v_db = v_db + "/";
               }
               v_dbprepare = v_db;
               
               if (! fileExists(v_db))
               {
                   v_message = "Invalid db: " + v_db;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               
               continue;
           }
           
           else if (str.find("buffer=") != string::npos)
           {
               v_buffer = stoi(str.substr(7));
               if (v_buffer < 10000) {v_buffer = 10000;}
               continue;
           }

           else if (str.find("downsample=") != string::npos)
           {
               downsampling = stoi(str.substr(11));
               if (downsampling < 5) {downsampling = 5;}
               continue;
           }
           
           else if (str.find("output=") != string::npos)
           {
               v_output = str.substr(7);
               char last_ch = v_output.back();
               char ch = '/';
               if (last_ch != ch)
               {
                   v_output = v_output + "/";
               }
               continue;
           }
           
           else if (str.find("sample=") != string::npos)
           {
               v_sample = str.substr(7);
               if (v_sample != "") {v_sample = v_sample + "_";}
               continue;
           }
           
           else if (str.find("--noconfig") != string::npos)
           {
               v_useconfig = 1;
               continue;
           }
           else if (str.find("--nomultiphasing") != string::npos)
           {
               v_multiphasing = 0;
               continue;
           }
           

           else if (str.find("--type") != string::npos)
           {
               v_type_after_map = 1;
               continue;
           }
           
           else if (str.find("--rna") != string::npos)
           {
               v_rnaseq = 1;
               continue;
           }

           else if (str.find("--low-mem") != string::npos)
           {
               v_lowmem = 1;
               continue;
           }
           
           else if (str.find("--keep-original-sam") != string::npos)
           {
               keep_sam = 1;
               continue;
           }
           
           else if (str.find("samtools=") != string::npos)
           {
               v_samtools = str.substr(9);
               if (! fileExists(v_samtools))
               {
                   v_message = "Invalid samtools: " + v_samtools;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str.find("blast=") != string::npos)
           {
               v_blast = str.substr(6);
               if (! fileExists(v_blast))
               {
                   v_message = "Invalid blast: " + v_blast;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
/*
           else if (str.find("gatk=") != string::npos)
           {
               v_gatk = str.substr(5);
               if (! fileExists(v_gatk))
               {
                   v_message = "Invalid GATK: " + v_gatk;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
 
 */
           else if (str.find("whatshap=") != string::npos)
           {
               v_whats = str.substr(9);
               if (! fileExists(v_whats))
               {
                   v_message = "Invalid WhatsHap: " + v_whats;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }
           
           
           else if (str.find("bwa=") != string::npos)
           {
               v_bwa = str.substr(4);
               if (! fileExists(v_bwa))
               {
                   v_message = "Invalid bwa: " + v_bwa;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }


           else if (str.find("star=") != string::npos)
           {
               v_star = str.substr(5);
               if (! fileExists(v_star))
               {
                   v_message = "Invalid star: " + v_db;
                   warnings.push_back(v_message);
                   command_ok = 0;
               }
               continue;
           }

           else if (str.find("target=") != string::npos)
           {
               v_target = str.substr(7);
               continue;
           }
           
           else if (str.find("--debug") != string::npos)
           {
               v_debug = 1;
               continue;
           }

           
           else if (str.find("--verbose") != string::npos)
           {
               v_verbose = 1;
               continue;
           }
           
           else if (str.find("--bypass-version-check") != string::npos)
           {
               v_bypass_version_check = 1;
               continue;
           }
           
           /*
           else if (str.find("--forcereindexing") != string::npos)
           {
               v_forceindex = 1;
               continue;
           }
           

           else if (str.find("--intergenic") != string::npos)
           {
               v_intergenic = 1;
               continue;
           }
 
           else if (str.find("--fullmapping") != string::npos)
           {
               v_intergenic = 1;
               continue;
           }
*/

           else if (str.find("--skip-adjust") != string::npos)
           {
               v_skiptyping = 1;
               continue;
           }

 //          else if (str.find("--use-local-db") != string::npos)
 //          {
 //              v_use_local = 1;
 //              continue;
 //          }
           
           else if (str.find("--skip-unmapped") != string::npos)
           {
               v_skip_unmapped = 1;
               continue;
           }
           
/*
           else if (str.find("--type-all") != string::npos)
            {
                v_typeall = 1;
                continue;
            }

           else if (str.find("--assemble") != string::npos)
           {
               v_assembling = 1;
               continue;
           }
 */
           else if (str.find("size=") != string::npos)
           {
               v_size = stoi(str.substr(5));
               if (v_size < 30) {v_size = 30;}
               continue;
           }
           
           else if (str.find("tolerance=") != string::npos)
           {
               v_tolerance = stod(str.substr(10));
               continue;
           }

           else if (str.find("error=") != string::npos)
           {
               v_mtrim_error = stod(str.substr(6));
               continue;
           }
           
           else
           {
               v_message = "Unknown option: " + str;
               warnings.push_back(v_message);
               command_ok = 0;
           }
       }
    }


   
    
    if (argc == 1)
    {
        main_help();
 
        if (warnings.size() > 0)
        {
            v_message = "Warning:  " + warnings[0];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        if (warnings.size() > 1)
        {
            for(int a = 1; a < warnings.size(); a++)
                v_message = "          " + warnings[a];
                screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        screen_message (screen_size, 0, "", 1, v_quiet);
  
        return (0);
    }
    
    
    
    if (command_ok == 0) {v_r1 = "";v_r0 = "";v_bam = "";}
    command_ok = 0;
    
    
    if (strcmp(argv[1],"dna") == 0)
    {
            
        main_dna_map();
        command_ok = 1;
    }
    
 
    if (strcmp(argv[1],"rna") == 0)
    {
        main_rna_map();
        command_ok = 1;
    }
 
    if (strcmp(argv[1],"select") == 0)
    {
        main_preselect();
        command_ok = 1;
    }
    
    if (strcmp(argv[1],"type") == 0)
    {
        main_typing();
        command_ok = 1;
    }

    if (strcmp(argv[1],"setup") == 0)
    {
        main_setup();
        command_ok = 1;
    }
    
    if ((command_ok == 0))
    {
        string str = argv[1];
        v_message = "Unknown command: " + str;
        warnings.push_back(v_message);
        main_help();
        if (warnings.size() > 0)
        {
            screen_message (screen_size, 0, "", 1, v_quiet);
            v_message = "Warning:  " + warnings[0];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
        if (warnings.size() > 1)
        {
            for (int a = 1; a < warnings.size(); a++) {
                v_message = "          " + warnings[a];
                screen_message (screen_size, 0, v_message, 1, v_quiet);
            }
        }
        screen_message (screen_size, 0, "", 1, v_quiet);
        return (0);
    }

    
    
    
    auto clock_end = std::chrono::steady_clock::now();
    auto diff = clock_end - clock_start;
    v_message = "Elapsed time: " + to_string(((std::chrono::duration <double, std::milli> (diff).count())/1000)) + " s";
    warnings.push_back(v_message);
    
    
    if (warnings.size() > 0)
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Warning:  " + warnings[0];
        screen_message (screen_size, 0, v_message, 1, v_quiet);
    }
    
    if (warnings.size() > 1)
    {
        for (int a = 1; a < warnings.size(); a++) {
            v_message = "          " + warnings[a];
            screen_message (screen_size, 0, v_message, 1, v_quiet);
        }
    }
    screen_message (screen_size, 0, "", 1, v_quiet);

    return 0;
}
