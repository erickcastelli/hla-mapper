#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <string>
#include <thread>
#include <mutex>
#include <numeric>

#include "ThreadPool.hpp"
#include "functions.hpp"
#include "typing.hpp"
#include "external.hpp"

using namespace std;

string sample = "";
string reference = "";
int fraglimit = 750;
string bam = "";
string resources = "";
string sequence_master = "";
string chr = "";
float todo = 0;
//float done = 0;
int max_nm_limit = 6;
map <pair<string,string>,string> summary_genomic;
map <pair<string,string>,string> summary_cds;
string map_type = "";
string use_chr = "";
int limit_block_size = 30;


mutex mtx_type;






string downsample (string currentbam, int target, string currentout, string currentgene, string sample)
{
    
    string cmd = "samtools mpileup --output-QNAME " + currentbam;
    string data = GetStdoutFromCommand(cmd);
    string outsam = currentout + "/" + sample + "." + currentgene + ".downsampled.sam";
    string outbam = currentout + "/" + sample + "." + currentgene + ".downsampled.bam";

    map <string,int> used;

    vector <string> lines;
    boost::split(lines,data,boost::is_any_of("\n"));
    for (auto item : lines)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        int cov = stoi(fields[3]);
        vector <string> reads;
        boost::split(reads,fields[6],boost::is_any_of(","));
        
        if (cov <= target) {
            for (auto item : reads) {used[item] = 1;}
            continue;
        }
        
        if (cov > target) {
            int already_used = 0;
            vector <string> possibilities;
            for (auto item : reads)
            {
                if (used.find(item) != used.end()){already_used++;}
                else {possibilities.push_back(item);}
            }
            if (already_used >= target) {continue;}
            int number_to_select = target - already_used;
            if (possibilities.size() < number_to_select)
            {
                for (auto item : possibilities) {used[item] = 1;}
                continue;
            }
            

            for (int a = 1; a <= number_to_select; )
            {
                int random = rand() % possibilities.size();
                if (used.find(possibilities[random]) == used.end())
                {
                    used[possibilities[random]] = 1;
                    a++;
                    continue;
                }
                else {continue;}
            }
            continue;

        }

    }
    
    cmd = "samtools view -h " + currentbam;
    data = GetStdoutFromCommand(cmd);

    vector <string > newdata;
    lines.clear();
    boost::split(lines,data,boost::is_any_of("\n"));
    for (auto item : lines)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        if (item.substr(0,1) == "@") {newdata.push_back(item); continue;}
        vector <string> fields;
        boost::split(fields,item,boost::is_any_of("\t"));
        if (used.find(fields[0]) != used.end()) newdata.push_back(item); continue;
    }
    
    ofstream newdataout;
    newdataout.open(outsam.c_str());
    for (auto item : newdata)
    {
        newdataout << item << endl;
    }
    newdataout.close();
    
    cmd = "samtools sort " + outsam + " > " + outbam;
    GetStdoutFromCommand(cmd);
    
    cmd = "samtools index " + outbam;
    GetStdoutFromCommand(cmd);

    return outbam;
}


void load_reference () {
    ifstream input(reference );
    for( std::string line; getline( input, line ); )
    {
        if (line.substr(0,1) == ">"){continue;}
        if (line.substr(0,1) == ""){continue;}
        sequence_master += line;
    }
    input.close();
}


double randfrom(double min, double max)
{
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}

string combinations (int blocks)
{
  double limit = (pow(2,blocks) / 2);
  map <string,int> results;

  for (int count = 1; count < 10000000; count++)
  {
      string line = "1-h1";
      for (int a = 2; a <= blocks; a++)
      {
        int vetor = int(randfrom(1.0, 3.0));
        line.append(";" + to_string(a) + "-h" + to_string(vetor));
      }
      results[line] = 1;
      if (results.size() == limit) {break;}
  }
  string data;
  for (auto combination : results)
  {
    data.append(combination.first + "\n");
  }
  return data;
}


void generate_fasta (string vstart, string vend, string vcffile, string outfas)
{

    int $higher = 0;
    int $lower = 0;
    if (((vstart == "0") || (vend == "0")))
    {
        vector<string> line_data;
        ifstream input(vcffile);
        for( std::string line; getline( input, line ); )
        {
            if (line.substr(0,1) == "#"){continue;}
            if (line == ""){continue;}
            boost::split(line_data,line,boost::is_any_of("\t"));
            if ($higher == 0) {$higher = stoi(line_data[1]);}
            if ($lower == 0) {$lower = stoi(line_data[1]);}
            if (stoi(line_data[1]) >= $higher){$higher = stoi(line_data[1]);}
        }
        input.close();
    }
    if (vstart != "0") {$lower = stoi(vstart);}
    if (vend != "0") {$higher = stoi(vend);}

    string sequence = sequence_master.substr($lower-1, ($higher - $lower)+1);
    
    
    
    //loading vcf data
    vector<string> head;
    vector<string> subdata;
    vector<string> genotypes;
    vector<string> positions;
    std::map <pair<string,string>, string> snp_data;
    std::map <string,int> samples;
    std::map <pair<string,string>,string> alt_data;
    
    int snps = 0;
    int format_pos = 0;
    int ref_pos = 0;
    int alt_pos = 0;
    
    ifstream vcf(vcffile);
    for( std::string line; getline( vcf, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,2) == "##"){continue;}
        if (line.substr(0,4) == "#CHR")
        {
            boost::split(head,line,boost::is_any_of("\t"));
            int a = 0;
            for(vector<string>::iterator sample = head.begin();sample!=head.end();++sample)
            {
                if (*sample == "REF"){ref_pos = a;}
                if (*sample == "ALT"){alt_pos = a;}
                if (*sample == "FORMAT"){format_pos = a;break;}
                a++;
            }
            continue;
        }
        
        boost::split(subdata,line,boost::is_any_of("\t"));
        if ((stoi(subdata[1]) < stoi(vstart)) && (vstart != "0")) {continue;}
        if ((stoi(subdata[1]) > stoi(vend) ) && (vend != "0")) {continue;}
        
        positions.push_back (subdata[1]);
        snps++;
        
        vector <string> alternatives;
        boost::split(alternatives,subdata[alt_pos],boost::is_any_of(","));
        
        pair <string,string> key = make_pair(subdata[1],"0");
        alt_data[key] = subdata[ref_pos];
        
        int a = 1;
        for(vector<string>::iterator alt = alternatives.begin();alt!=alternatives.end();++alt)
        {
            pair <string,string> key = make_pair(subdata[1],to_string(a));
            alt_data[key] = *alt;
            a++;
        }
        
        for (a = format_pos + 1; a < subdata.size(); a++)
        {
            boost::split(genotypes,subdata[a],boost::is_any_of(":"));
            pair <string,string> key = make_pair(head[a],subdata[1]);
            snp_data[key] = genotypes[0];
            
        }
    }
    vcf.close();

    
    for (int a = format_pos + 1; a < subdata.size(); a++)
    {
        samples[head[a]] = 1;
    }
    
    ofstream myfile;
    myfile.open (outfas);
    
    if (snps > 0) {
        for (int a = format_pos + 1; a < subdata.size(); a++)
        {
            myfile << ">" << head[a] << "_h1\n";
            string draft = sequence;
            int correct = 0;
            for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
            {
                pair <string,string> key = make_pair(head[a],*pos);
                string gen = snp_data[key];
                vector<string> alleles;
                boost::split(alleles,gen,boost::is_any_of("|"));
                
                if (alleles.size() == 2)
                {
                    if (alleles[0] == "0") {continue;}
                    else
                    {
                        pair <string,string> key1 = make_pair(*pos,"0");
                        pair <string,string> key2 = make_pair(*pos,alleles[0]);
                        
                        if (alt_data[key2] == "*") {continue;}
                        
                        int b = stoi(*pos)+correct-$lower;
                        
                        draft.replace(b,alt_data[key1].size(),alt_data[key2]);
                        correct = correct + (alt_data[key2].size() - alt_data[key1].size());
                    }
                }
            }
            myfile << draft << endl;
            
            
            myfile << ">" << head[a] << "_h2\n";
            draft = sequence;
            correct = 0;
            for(vector<string>::iterator pos = positions.begin();pos!=positions.end();++pos)
            {
                pair <string,string> key = make_pair(head[a],*pos);
                string gen = snp_data[key];
                vector<string> alleles;
                boost::split(alleles,gen,boost::is_any_of("|"));
                
                if (alleles.size() == 2)
                {
                    
                    if (alleles[1] == "0") {continue;}
                    else
                    {
                        pair <string,string> key1 = make_pair(*pos,"0");
                        pair <string,string> key2 = make_pair(*pos,alleles[1]);
                        
                        if (alt_data[key2] == "*") {continue;}
                        
                        
                        int b = stoi(*pos)+correct-$lower;
                        
                        draft.replace(b,alt_data[key1].size(),alt_data[key2]);
                        correct = correct + (alt_data[key2].size() - alt_data[key1].size());
                    }
                }
            }
            myfile << draft << endl;
        }
    }
    
    
    if (snps == 0) {
        string samplename = head[format_pos + 1];
        myfile << ">" << samplename << "_h1\n";
        string draft = sequence;
        myfile << draft << endl;
            
        myfile << ">" << samplename << "_h2\n";
        myfile << draft << endl;
    }
    
    
    myfile.close();
}




/*
string run_gatk (string currentgene, string currentout, string gvcf, string vcf, string bed, string sample, string currentbam)
{
    
    string log = currentout + "/log/" + sample + "." + currentgene + ".gatk.log";
    string cmd = v_java + " -Xmx4g -jar " + v_gatk + " -T UnifiedGenotyper -rf NotPrimaryAlignment -rf DuplicateRead -glm BOTH -dt NONE -R " + reference + " -I " + currentbam + " -o " + vcf + " -L " + bed + " > " + log + " 2>" + log;
    GetStdoutFromCommand(cmd);
    if (filesize(vcf.c_str()) == 0) {return "fail";}

    ifstream vcfin;
    ofstream intervals;
    string intervalsfile = vcf + ".intervals";
    intervals.open(intervalsfile.c_str());
    vcfin.open(vcf.c_str());
    int do_realign = 0;
    for( std::string line; getline( vcfin, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        if ((data[3].size() > 1) || (data[4].size() > 1))
        {
            if (data[4].find(",") != string::npos) {continue;}
            int start = stoi(data[1]) - 10;
            int end = stoi(data[1]) + 10;
            intervals << data[0] << ":" << start << "-" << end << endl;
            do_realign = 1;
        }
    }
    vcfin.close();
    intervals.close();
    
    if (do_realign == 1)
    {
        string bamout = currentbam;
        boost::replace_all(bamout, ".bam", ".realigned.bam");
        
        string cmd = v_java + " -Xmx4g -jar " + v_gatk + " -T IndelRealigner -R " + reference + " -I " + currentbam + " -o " + bamout + " -targetIntervals " + intervalsfile + " >> " + log + " 2>>" + log;
        GetStdoutFromCommand(cmd);
       
        cmd = v_java + " -Xmx4g -jar " + v_gatk + " -T UnifiedGenotyper -rf NotPrimaryAlignment -rf DuplicateRead -glm BOTH -dt NONE -R " + reference + " -I " + bamout + " -o " + vcf + " -L " + bed + " >> " + log + " 2>>" + log;
        GetStdoutFromCommand(cmd);
    }
    return currentbam;
}
*/



string run_freebayes (string currentgene, string currentout, string vcf, string bed, string sample, string currentbam, string type)
{
    
    string log = currentout + "/log/" + sample + "." + currentgene + ".freebayes." + type + ".log";
    string cmd = v_freebayes + " -f " + reference + " -t " + bed + " " + currentbam + " > " + vcf + " 2>" + log;
    if (v_debug) {cout << cmd << endl << endl;}
    GetStdoutFromCommand(cmd);
    if (filesize(vcf.c_str()) == 0) {return "fail";}

    return currentbam;
}






string check_coverage (string currentgene, string bed, string currentbam)
{
    string cmd = v_samtools + " depth -a -b " + bed + " " + currentbam;
    string depth = GetStdoutFromCommand(cmd);
    vector <string> depth_data;
    boost::split(depth_data,depth,boost::is_any_of("\n"));
    float count = 0;
    float zero = 0;
    float sum = 0;
    for (auto item : depth_data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> data;
        boost::split(data,item,boost::is_any_of("\t"));
        if (data[2] == "0") {continue;}
        count++;
        sum = sum + stof(data[2]);
    }
    float threshold = 0.2 * (sum / count);
    if (threshold < 1) {threshold = 1;}
        
    
    count = 0;
    for (auto item : depth_data)
    {
        if (item == "") {continue;}
        if (item.substr(0,1) == "[") {continue;}
        vector <string> data;
        boost::split(data,item,boost::is_any_of("\t"));
        if (stof(data[2]) < threshold) {zero++;}
        count++;
    }
    float ratio = 1 - (zero / count);
    if ((ratio < 0.95) && (ratio >= 0.2)) {return "CDS";}
    if (ratio >= 0.95) {return "genomic";}
    return "fail";
}



void checking_vcf (string vcftmp, string vcf)
{
    ifstream vcftemp(vcftmp.c_str());
    ofstream vcfout;
    vcfout.open (vcf.c_str());
 
    float quality = 0;
    float count = 0;
    for( std::string line; getline( vcftemp, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        if (chr == "") {chr = data[0];}
        quality = quality + stof(data[5]);
        count++;
    }
    float limit_unknown = (quality / count) * 0.05;

    vcftemp.clear();
    vcftemp.seekg (0, ios::beg);

    for( std::string line; getline( vcftemp, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {vcfout << line << endl; continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        if (chr == "") {chr = data[0];}
        data[7] = ".";
        data[8] = "GT:AD:DP:GQ";
        vector <string> fields;
        boost::split(fields,data[9],boost::is_any_of(":"));
        boost::replace_all(fields[0], "|", "/");
        if (stoi(fields[3]) < 20) {fields[0] = "./.";}
        if ((fields[1] == ".") || (fields[2] == ".")){fields[0] = "./.";}
        
        float cov = stof(fields[2]);
        vector <string> allele_cov;
        boost::split(allele_cov,fields[1],boost::is_any_of(","));
        vector <string> alleles;
        boost::split(alleles,fields[0],boost::is_any_of("/"));

        if (stof(data[5]) < limit_unknown) {alleles[0] = "."; alleles[1] = ".";}

        //heterozigosis
        if ((alleles[0] != alleles[1]) && (alleles[0] != "."))
        {
            if ((stof(allele_cov[stoi(alleles[0])]) / cov ) < 0.2)
            {
                alleles[0] = ".";
//                alleles[1] = ".";
            }
            if ((stof(allele_cov[stoi(alleles[1])]) / cov ) < 0.2)
            {
//                alleles[0] = ".";
                alleles[1] = ".";
            }
        }

        //homozygosis
        if ((alleles[0] == alleles[1]) && (alleles[0] != "."))
        {
            if (cov < 8) {alleles[0] = "."; alleles[1] = ".";}
            if ((stof(allele_cov[stoi(alleles[0])]) / cov ) < 0.95)
            {
 //               alleles[0] = ".";
                alleles[1] = ".";
            }
        }
        
        fields[0] = alleles[0] + "/" + alleles[1];

        string newline = data[0];
        for (int a = 1; a <= 8; a++){newline.append("\t" + data[a]);}
        newline.append("\t" + fields[0]);
        for (int a = 1; a <= 3; a++){newline.append(":" + fields[a]);}
        vcfout << newline << endl; continue;
    }
    vcftemp.close();
    vcfout.close();
}



void checking_vcf_freebayes (string vcftmp, string vcf, string type, string bed)
{
    
    
    map <int,int> valid_positions;
    map <int,int> valid_position_limit;
    map <int,int> valid_segments;
    //LOADING BED
    ifstream bedfile(bed.c_str());
    int count = 1;
    for( std::string line; getline( bedfile, line ); )
    {
        if (line == "") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        int start = stoi(fields[1]);
        int end = stoi(fields[2]);
        valid_segments[start] = end;
        
        for (int a = start; a <= end; a++)
        {
            valid_positions[a] = count;
            valid_position_limit[a] = end;
        }
        count++;
    }
    bedfile.close();
    
    
    ifstream vcftemp(vcftmp.c_str());
    ofstream vcfout;
    vcfout.open (vcf.c_str());
 
    float quality = 0;
    count = 0;
    for( std::string line; getline( vcftemp, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        if (chr == "") {chr = data[0];}
        quality = quality + stof(data[5]);
        count++;
    }
    float limit_unknown = (quality / count) * 0.05;

    vcftemp.clear();
    vcftemp.seekg (0, ios::beg);

    for( std::string line; getline( vcftemp, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {vcfout << line << endl; continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        if (chr == "") {chr = data[0];}
        data[7] = ".";
        data[8] = "GT:DP:AD:RO";
        vector <string> fields;
        boost::split(fields,data[9],boost::is_any_of(":"));
        boost::replace_all(fields[0], "|", "/");
   //     if (stoi(fields[3]) < 20) {fields[0] = "./.";}
        if ((fields[1] == ".") || (fields[2] == ".")){fields[0] = "./.";}
        
        float cov = stof(fields[1]);
        vector <string> allele_cov;
        boost::split(allele_cov,fields[2],boost::is_any_of(","));
        vector <string> alleles;
        boost::split(alleles,fields[0],boost::is_any_of("/"));

/*
        // low qual
        if (stof(data[5]) < limit_unknown) {alleles[0] = "."; alleles[1] = ".";}
*/
        //heterozigosis
        if ((alleles[0] != alleles[1]) && (alleles[0] != "."))
        {
            if ((stof(allele_cov[stoi(alleles[0])]) / cov ) < 0.2)
            {
                alleles[0] = ".";
            }
            if ((stof(allele_cov[stoi(alleles[1])]) / cov ) < 0.2)
            {
                alleles[1] = ".";
            }
        }
/*
        //homozygosis
        if ((alleles[0] == alleles[1]) && (alleles[0] != "."))
        {
            if (cov < 8) {alleles[0] = "."; alleles[1] = ".";}
            if (alleles[0] != ".") {
                if ((stof(allele_cov[stoi(alleles[0])]) / cov ) < 0.95)
                {
                    alleles[1] = ".";
                }
            }
        }
        fields[0] = alleles[0] + "/" + alleles[1];
*/
        string newline = data[0];
        for (int a = 1; a <= 8; a++){newline.append("\t" + data[a]);}
        newline.append("\t" + fields[0]);
        for (int a = 1; a <= 3; a++){newline.append(":" + fields[a]);}
        vcfout << newline << endl; continue;
    }
    vcftemp.close();
    vcfout.close();
    
    
    removefile(vcf + ".gz", 0);
    string cmd = "bgzip " + vcf;
    GetStdoutFromCommand(cmd);
    cmd = "tabix -f -p vcf " + vcf + ".gz";
    GetStdoutFromCommand(cmd);

    int block = 1;
    for (auto item : valid_segments)
    {
        string cmd = v_bcftools + " view " + vcf + ".gz " + chr + ":" + to_string(item.first) + "-" + to_string(item.second) + " > " + vcf + ".part" + to_string(block) + ".vcf";
        GetStdoutFromCommand(cmd);
        block++;
    }
    
    
    
}



void  run_whatshap (string vcf, string whats, string whatstmp, string currentout, string currentgene, string sample, float frag_mean, float frag_std, string currentbam, string bed)
{
    
    map <int,int> valid_positions;
    map <int,int> valid_position_limit;
    map <int,int> valid_segments;
    //LOADING BED
    ifstream bedfile(bed.c_str());
    int count = 1;
    for( std::string line; getline( bedfile, line ); )
    {
        if (line == "") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        int start = stoi(fields[1]);
        int end = stoi(fields[2]);
        valid_segments[start] = end;
        
        for (int a = start; a <= end; a++)
        {
            valid_positions[a] = count;
            valid_position_limit[a] = end;
        }
        count++;
    }
    bedfile.close();
    
    
    
    int block = 1;
    for (auto item : valid_segments)
    {
        string blockvcf = vcf + ".part" + to_string(block) + ".vcf";
        string whatsblock = vcf + ".part" + to_string(block) + ".whatshap.vcf";
        string log = currentout + "/log/" + sample + "." + currentgene + ".whatshap.cds.block" + to_string(block) + ".log";
        string cmd = v_whats + " phase --ignore-read-groups --tag PS --indel --reference " + reference + " -o " + whatsblock + " " + blockvcf + " " + currentbam + " > " + log;
        GetStdoutFromCommand(cmd);
        block++;
    }
    
    
    float limit = frag_mean + frag_std;
    
    vector <string> head;
    vector <string> snps;
    map <int,int> number_of_heterogygous_snps_in_each_block;
    
    block = 1;
    for (auto item : valid_segments)
    {
        number_of_heterogygous_snps_in_each_block[block] = 0;
        string whatsblock = vcf + ".part" + to_string(block) + ".whatshap.vcf";
        ifstream file;
        file.open(whatsblock.c_str());
        
        for( std::string line; getline( file, line ); )
        {
           if (line.substr(0,1) == "#") {
               continue;
           }
            if (line == "") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            vector <string> fields;
            boost::split(fields,data[9],boost::is_any_of(":"));
            vector <string> alleles;
            boost::split(alleles,fields[0],boost::is_any_of("|"));
            if (alleles.size() < 2) {boost::split(alleles,fields[0],boost::is_any_of("/"));}
            if (alleles[0] == alleles[1]) {continue;}
            if ((alleles[0] != ".") && (alleles[1] != ".")) {
                if (alleles[0] != alleles[1]) {
                    number_of_heterogygous_snps_in_each_block[block]++;
                }
            }
        }
        file.close();
        block++;
    }
    
    
    block = 1;
    for (auto item : valid_segments)
    {
        
        string whatsblock = vcf + ".part" + to_string(block) + ".whatshap.vcf";
        ifstream file;
        file.open(whatsblock.c_str());
        string last_ps = "";
        int last_pos = 0;
        
        for( std::string line; getline( file, line ); )
        {
           if (line.substr(0,1) == "#") {
               if (block == 1) {head.push_back(line);}
               continue;
           }
            
            if (line == "") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            vector <string> fields;
            boost::split(fields,data[9],boost::is_any_of(":"));
            vector <string> alleles;
            boost::split(alleles,fields[0],boost::is_any_of("|"));
            if (alleles.size() < 2) {boost::split(alleles,fields[0],boost::is_any_of("/"));}
            if (alleles.size() < 2) {snps.push_back(line);continue;}
            if (alleles[0] == alleles[1]) {snps.push_back(line);continue;}
            
            int pos = stoi(data[1]);
            if (last_pos == 0) {last_pos = pos;}
            if (fields.size() == 4) {fields.push_back(data[1]);}
            if (last_ps == "") {last_ps = fields[4];}
            
            if ((fields[4] == last_ps) && ((pos - last_pos) > limit))
            {
                int check = 0; //fase falsa , =1 fase verdadeira
                if (check == 0){last_ps = data[1];}
            }
            
//            cout << "Pos: " << pos << endl;
//            cout << "Block: " << block << endl;
//            cout << "N_block: " << number_of_heterogygous_snps_in_each_block[block] << endl;
//            cout << fields[0] << endl;
            
            if (number_of_heterogygous_snps_in_each_block[block] == 1)
            {
                //cout << fields[0] << " " << alleles[0] << " " << alleles[1] << " " << pos << endl;
                fields[0] = alleles[0] + "|" + alleles[1];
                last_ps = data[1];
            }
            data[8] = "GT:DP:AD:RO:PS";
            data[9] = fields[0] + ":" + fields[1] + ":" + fields[2] + ":" + fields[3] + ":" + last_ps;
            
            
            string newline = "";
            for (auto item: data)
            {
                newline.append("\t" + item);
            }
            snps.push_back(newline.substr(1));
            last_pos = pos;
        }
        file.close();
        block++;
    }
    
    
    ofstream whatsout;
    whatsout.open(whats.c_str());
    for (auto item : head){whatsout << item << endl;}
    for (auto item : snps){whatsout << item << endl;}
    whatsout.close();
}



void run_local_phasing (string whats, string phased, string currentout, string currentgene, string phaselog, float frag_mean, string currentbam)
{

    mtx_type.lock();
    map <int,string> phasedata;
    vector <string> head;
    map <int,int> unphased;
    map <int,string> hetero_phased;
    map <int,string> hetero_unphased;

    ofstream log;
    log.open (phaselog.c_str());
    
    ifstream vcf(whats.c_str());
    for( std::string line; getline( vcf, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "#") {head.push_back(line);continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));
        vector <string> fields;
        boost::split(fields,data[9],boost::is_any_of(":"));
        vector <string> alleles;
        boost::split(alleles,fields[0],boost::is_any_of("|"));
        if (alleles.size() == 1) {boost::split(alleles,fields[0],boost::is_any_of("/"));}
        string newline;
        
        if ((alleles[0] == ".") || (alleles[1] == ".")) {
            phasedata[stoi(data[1])] = line;
            continue;
        }
        
        
        if (alleles[0] == alleles[1]) {
            fields[0] = alleles[0] + "|" + alleles[1];
            newline = data[0];
            for (int a = 1; a <= 8; a++) {newline.append("\t" + data[a]);}
            newline.append("\t" + fields[0]);
            for (int a = 1; a <= 3; a++) {newline.append(":" + fields[a]);}
            phasedata[stoi(data[1])] = newline;
            continue;
        }

        if (alleles[0] != alleles[1]) {
            std::size_t found = fields[0].find("|");
            if (found!=std::string::npos) {
                hetero_phased[stoi(data[1])] = fields[4];
                phasedata[stoi(data[1])] = line;
                continue;
            }
            else {
            hetero_unphased[stoi(data[1])] = 1;
            phasedata[stoi(data[1])] = line;
            unphased[stoi(data[1])] = 1;
            continue;
            }
        }
    }
    vcf.close();

    
    int round = 1;
    for (auto position : unphased)
    {
        log << "Phasing multi-allelic variant at position " << position.first << endl;
        int start = 0;
        int end = 0;
        int target = position.first;
        
        for (int a = 1; a <= 100; a++)
        {
            int phased_count = 0;
            int unphased_count = 0;
            map <string,int> phase_set;
            for (int b = (target-1); b < (target + a); b++)
            {
                if (hetero_phased.find(b) != hetero_phased.end())
                {
                    phased_count++;
                    phase_set[hetero_phased[b]] = 1;
                }
                if (hetero_unphased.find(b) != hetero_unphased.end()) {unphased_count++;}
            }
            if ((phased_count > 0) && (unphased_count > 0))
            {
                if (phase_set.size() == 1) {
                    start = target;
                    end = target + a;
                    if ((end - start) >= 20) {break;}
                }
            }
        }
        
        if (start == 0)
        {
            for (int a = 1; a <= 100; a++)
            {
                int phased_count = 0;
                int unphased_count = 0;
                map <string,int> phase_set;
                for (int b = target+1; b > (target - a);)
                {
                    if (hetero_phased.find(b) != hetero_phased.end())
                    {
                        phased_count++;
                        phase_set[hetero_phased[b]] = 1;
                    }
                    if (hetero_unphased.find(b) != hetero_unphased.end()) {unphased_count++;}
                    b = b - 1;
                }
                if ((phased_count > 0) && (unphased_count > 0))
                {
                    if (phase_set.size() == 1) {
                        start = target - a;
                        end = target + 1;
                        if ((end - start) >= 20) {break;}
                    }
                }
            }
        }
        
        
        log << "Phasing multi-allelic variant at position " << position.first << ", region " << start << " to " << end << endl;

        if (start == 0)
        {
            log << "Phasing multi-allelic variant at position " << position.first << " : failed: no close phased heterozygotes" << endl;
            continue;
        }

        string ps = "";
        string testA = "";
        string testB = "";
        string cmd = "mkdir " + currentout + "/phasing/";
        GetStdoutFromCommand(cmd);
        cmd = "mkdir " + currentout + "/phasing/round." + to_string(round);
        GetStdoutFromCommand(cmd);
        string tmpvcfA = currentout + "/phasing/round." + to_string(round) + "/round." + to_string(round) + ".A.vcf";
        string tmpvcfB = currentout + "/phasing/round." + to_string(round) + "/round." + to_string(round) + ".B.vcf";
        ofstream outvcfA;
        outvcfA.open (tmpvcfA.c_str());
        for (auto item : head) {outvcfA << item << endl;}
        
        for (int a = start; a <= end; a++)
        {
            if (phasedata.find(a) != phasedata.end())
            {
                vector <string> data;
                boost::split(data,phasedata[a],boost::is_any_of("\t"));
                vector <string> fields;
                boost::split(fields,data[9],boost::is_any_of(":"));
                ps = fields[4];
            }
        }
        
        for (int a = start; a <= end; a++)
        {
            if (phasedata.find(a) != phasedata.end())
            {
                if (a == target)
                {
                    vector <string> data;
                    boost::split(data,phasedata[a],boost::is_any_of("\t"));
                    vector <string> fields;
                    boost::split(fields,data[9],boost::is_any_of(":"));
                    vector <string> alleles;
                    boost::split(alleles,fields[0],boost::is_any_of("/"));
                    fields[0] = alleles[0] + "|" + alleles[1];
                    
                    string newline = data[0];
                    data[8] = "GT:DP:AD:RO:PS";
                    for (int a = 1; a <= 8; a++){newline.append("\t" + data[a]);}
                    newline.append("\t" + fields[0]);
                    for (int a = 1; a <= 3; a++){newline.append(":" + fields[a]);}
                    outvcfA << newline << endl;
                    testA = newline + ":" + ps;
                    continue;
                }
                
                outvcfA << phasedata[a] << endl;
                
                if (ps == "") {
                    vector <string> data;
                    boost::split(data,phasedata[a],boost::is_any_of("\t"));
                    vector <string> fields;
                    boost::split(fields,data[9],boost::is_any_of(":"));
                    vector <string> alleles;
                    boost::split(alleles,fields[0],boost::is_any_of("|"));
                    if (alleles.size() != 2) {continue;}
                    if (alleles[0] != alleles[1]) {
                    ps = fields[4];
                    }
                }
            }
        }
        outvcfA.close();
        
        ofstream outvcfB;
        outvcfB.open (tmpvcfB.c_str());
        for (auto item : head) {outvcfB << item << endl;}
        for (int a = start; a <= end; a++)
        {
            if (phasedata.find(a) != phasedata.end())
            {
                if (a == target)
                {
                    vector <string> data;
                    boost::split(data,phasedata[a],boost::is_any_of("\t"));
                    vector <string> fields;
                    boost::split(fields,data[9],boost::is_any_of(":"));
                    vector <string> alleles;
                    boost::split(alleles,fields[0],boost::is_any_of("/"));
                    fields[0] = alleles[1] + "|" + alleles[0];
                    
                    string newline = data[0];
                    data[8] = "GT:DP:AD:RO:PS";
                    for (int a = 1; a <= 8; a++){newline.append("\t" + data[a]);}
                    newline.append("\t" + fields[0]);
                    for (int a = 1; a <= 3; a++){newline.append(":" + fields[a]);}
                    outvcfB << newline << endl;
                    testB = newline + ":" + ps;;
                    continue;
                }
                outvcfB << phasedata[a] << endl;
            }
        }
        outvcfB.close();
        
        string outfasA = currentout + "/phasing/round." + to_string(round) + "/round." + to_string(round) + ".A.fas";
        string outfasB = currentout + "/phasing/round." + to_string(round) + "/round." + to_string(round) + ".B.fas";
        generate_fasta(to_string(start), to_string(end), tmpvcfA, outfasA);
        generate_fasta(to_string(start), to_string(end), tmpvcfB, outfasB);
        
        ifstream fasA(outfasA);
        string line;
        getline( fasA, line );
        getline( fasA, line );
        string seqh1A = line;
        string seqh2A = line;
        fasA.close();
        ifstream fasB(outfasB);
        getline( fasB, line );
        getline( fasB, line );
        string seqh1B = line;
        string seqh2B = line;
        fasB.close();
        boost::to_upper(seqh1A);
        boost::to_upper(seqh2A);
        boost::to_upper(seqh1B);
        boost::to_upper(seqh2B);

        string seqh1Arev = reverse_and_complement(seqh1A);
        string seqh2Arev = reverse_and_complement(seqh2A);
        string seqh1Brev = reverse_and_complement(seqh1B);
        string seqh2Brev = reverse_and_complement(seqh2B);

        
        cmd = v_samtools + " view " + currentbam + " " + chr + ":" + to_string(start) + "-" + to_string(end);
        string sam = GetStdoutFromCommand(cmd);
        vector <string> samdata;
        boost::split(samdata,sam,boost::is_any_of("\n"));
        
        map <string,int> hits;
        hits[seqh1A] = 0;
        hits[seqh2A] = 0;
        hits[seqh1B] = 0;
        hits[seqh2B] = 0;
        
        for (auto item : samdata)
        {
            if (item == "") {continue;}
            if (item.substr(0,1) == "[") {continue;}
            vector <string> data;
            boost::split(data,item,boost::is_any_of("\t"));
            string seq = data[9];
            
            if (seq.find(seqh1A) !=std::string::npos) {hits[seqh1A]++;}
            if (seq.find(seqh1Arev) !=std::string::npos) {hits[seqh1A]++;}
            if (seq.find(seqh2A) !=std::string::npos) {hits[seqh2A]++;}
            if (seq.find(seqh2Arev) !=std::string::npos) {hits[seqh2A]++;}
            if (seq.find(seqh1B) !=std::string::npos) {hits[seqh1B]++;}
            if (seq.find(seqh1Brev) !=std::string::npos) {hits[seqh1B]++;}
            if (seq.find(seqh2B) !=std::string::npos) {hits[seqh2B]++;}
            if (seq.find(seqh2Brev) !=std::string::npos) {hits[seqh2B]++;}
        }

        int groupA = 0;
        if ((hits[seqh1A] > 0) && (hits[seqh2A] > 0))
        { groupA = hits[seqh1A] + hits[seqh2A];}
        
        int groupB = 0;
        if ((hits[seqh1B] > 0) && (hits[seqh2B] > 0))
        {groupB = hits[seqh1B] + hits[seqh2B];}
        
        if (groupA > groupB)
        {
            phasedata[target] = testA;
            hetero_phased[target] = ps;
            hetero_unphased.erase(target);
            log << "Phasing multi-allelic variant at position " << position.first << " : done" << endl;
        }

        if (groupB > groupA)
        {
            phasedata[target] = testB;
            hetero_phased[target] = ps;
            hetero_unphased.erase(target);
            log << "Phasing multi-allelic variant at position " << position.first << " : done" << endl;
        }

        if (groupB == groupA)
        {
            log << "Phasing multi-allelic variant at position " << position.first << " : failed" << endl;
        }
        round++;
    }
    
    log.close();
    
    
    ofstream myfile;
    myfile.open (phased);
    for (auto item : head){ myfile << item << endl;}
    for (auto item : phasedata) {myfile << item.second << endl;}
    myfile.close();
    mtx_type.unlock();
}







void plot_graf(string sample, string currentout, string currentgene, string type, string alleleA, string alleleB, string currentbam)
{
    string bestAseq = "";
    string bestBseq = "";
    string refbwa = "";
    if (type == "cds") {refbwa = resources + "/fasta/" + currentgene + "/" + currentgene + "_nuc.bwa.fas";}
    if (type == "genomic") {refbwa = resources + "/fasta/" + currentgene + "/" + currentgene + "_gen.bwa.fas";}
    
    
    ifstream ref(refbwa.c_str());
    for( std::string line; getline( ref, line ); )
    {
        string id = line.substr(1);
        getline(ref, line );
        string seq = line;
        if (id == alleleA) {
            bestAseq = seq;
        }
        if (id == alleleB) {
            bestBseq = seq;
        }
    }
    ref.close();
    
    string bestAconverted = alleleA;
    boost::replace_all(bestAconverted, "*", "_");
    boost::replace_all(bestAconverted, ":", "_");
    string bestBconverted = alleleB;
    boost::replace_all(bestBconverted, "*", "_");
    boost::replace_all(bestBconverted, ":", "_");

    string cmd = "mkdir " + currentout + "/references/";
    GetStdoutFromCommand(cmd);
    
    map <string,string> references;
    string outrefA = currentout + "/references/" + bestAconverted + ".reference." + type + ".fas";
    ofstream outA;
    outA.open(outrefA);
    outA << ">" << alleleA << endl;
    outA << bestAseq << endl;
    references[alleleA] = bestAseq;
    outA.close();
    string outrefB = currentout + "/references/" + bestBconverted + ".reference." + type + ".fas";
    ofstream outB;
    outB.open(outrefB);
    outB << ">" << alleleB << endl;
    outB << bestBseq << endl;
    references[alleleB] = bestBseq;
    outB.close();
 
    cmd = v_bwa + " index " + outrefA;
    GetStdoutFromCommand(cmd);
    cmd = v_bwa + " index " + outrefB;
    GetStdoutFromCommand(cmd);
 
    
    cmd  = v_samtools + " view " + currentbam;
    string readdata = GetStdoutFromCommand(cmd);
    

    vector <string> data;
    boost::split(data,readdata,boost::is_any_of("\n"));
    map <string,int> validreads;
    for (auto line : data)
    {
        if (line.substr(0,1) == "[") {continue;}
        if (line == "") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        validreads["@"+ fields[0]] = 1;
    }
    
    string r1 = "";
    string r2 = "";
    
 
    if (map_type == "paired") {
        
        string subgene = currentgene;
        if (subgene == "HLA-DMA") {subgene = "HLA-DM";}
        if (subgene == "HLA-DMB") {subgene = "HLA-DM";}
        if (subgene == "TAP1") {subgene = "TAP";}
        if (subgene == "TAP2") {subgene = "TAP";}
        if (subgene == "HLA-DPA1") {subgene = "HLA-DPAB1";}
        if (subgene == "HLA-DPB1") {subgene = "HLA-DPAB1";}

       

        r1 = v_mapout + "/" + sample + "_" + subgene + "_R1.fastq";
        r2 = v_mapout + "/" + sample + "_" + subgene + "_R2.fastq";
        
        string r1new = currentout + "/" + sample + "." + currentgene + ".downsampled.R1.fastq";
        string r2new = currentout + "/" + sample + "." + currentgene + ".downsampled.R2.fastq";
        ifstream r1original;
        ofstream r1newout;
        r1original.open(r1.c_str());
        r1newout.open(r1new.c_str());
        for( std::string line; getline( r1original, line ); )
        {
            string id = line;
            vector <string> line_id;
            boost::split(line_id,id,boost::is_any_of(" "));
            id = line_id[0];
            getline(r1original, line);
            string seq = line;
            getline(r1original, line);
            string info = line;
            getline(r1original, line);
            string qual = line;
            if (validreads.find(id) != validreads.end())
            {
 
                r1newout << id << endl << seq << endl << info << endl << qual << endl;
            }
        }
        r1newout.close();
        r1original.close();

        ifstream r2original;
        ofstream r2newout;
        r2original.open(r2.c_str());
        r2newout.open(r2new.c_str());
        for( std::string line; getline( r2original, line ); )
        {
            string id = line;
            vector <string> line_id;
            boost::split(line_id,id,boost::is_any_of(" "));
            id = line_id[0];
            getline(r2original, line);
            string seq = line;
            getline(r2original, line);
            string info = line;
            getline(r2original, line);
            string qual = line;
            if (validreads.find(id) != validreads.end())
            {
                r2newout << id << endl << seq << endl << info << endl << qual << endl;
            }
        }
        r2newout.close();
        r2original.close();
        r1 = r1new ;
        r2 = r2new ;
    }
    
    
    
    
    if (map_type == "single") {
        r1 = v_mapout + "/" + sample + "_" + currentgene + "_R0.fastq";
        string r1new = currentout + "/" + sample + "." + currentgene + ".downsampled.R1.fastq";
        ifstream r1original;
        ofstream r1newout;
        r1original.open(r1.c_str());
        r1newout.open(r1new.c_str());
        for( std::string line; getline( r1original, line ); )
        {
            string id = line;
            vector <string> line_id;
            boost::split(line_id,id,boost::is_any_of(" "));
            id = line_id[0];
            getline(r1original, line);
            string seq = line;
            getline(r1original, line);
            string info = line;
            getline(r1original, line);
            string qual = line;
            if (validreads.find(id) != validreads.end())
            {
                r1newout << id << endl << seq << endl << info << endl << qual << endl;
            }
        }
        r1newout.close();
        r1original.close();
        r1 = r1new ;
    }
    

    string samAdata;
    string samBdata;
    string samA = currentout + "/" + sample + "." + currentgene + ".bwa.alelleA.sam";
    string samB = currentout + "/" + sample + "." + currentgene + ".bwa.alelleB.sam";
    
    if (map_type == "paired") {
        string sam = currentout + "/" + sample + "." + currentgene + ".bwa.alelleA.sam";
        cmd = v_bwa + " mem -B 2 -O 3,3 -L 2,2 -U 5 -o " + samA + " " + outrefA + " " + r1 + " " + r2;
        string out = GetStdoutFromCommand(cmd);
        cmd = v_bwa + " mem -B 2 -O 3,3 -L 2,2 -U 5 -o " + samB + " " + outrefB + " " + r1 + " " + r2;
        if (v_debug) {cout << cmd << endl << endl;}
        out = GetStdoutFromCommand(cmd);
    }
 
    if (map_type == "single") {
        cmd = v_bwa + " mem -B 2 -O 3,3 -L 2,2 -U 5 -o " + samA + " " + outrefA + " " + r1;
        string out = GetStdoutFromCommand(cmd);
        cmd = v_bwa + " mem -B 2 -O 3,3 -L 2,2 -U 5 -o " + samB + " " +  outrefB + " " + r1;
        if (v_debug) {cout << cmd << endl << endl;}
        out = GetStdoutFromCommand(cmd);
    }
    
    string headA = "";
    string headB = "";
    map <string,int> reads_A;
    map <string,int> reads_B;
    map <string,string> readdb_A;
    map <string,string> readdb_B;
    map <string,int> reads;
    data.clear();


    
    ifstream samAfile(samA.c_str());
    for( std::string line; getline( samAfile, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "[") {continue;}
        if (line.substr(0,1) == "@") {headA.append(line + "\n"); continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        if (fields[2] == "*") {continue;}
        if (fields[3] == "0") {continue;}
        int flag = stoi(fields[1]);
        if (flag > 200) {continue;}
        vector <string> nm;
        boost::split(nm,fields[11],boost::is_any_of(":"));
        if (stoi(nm[2]) > max_nm_limit) {continue;}
        reads_A[fields[0]] = reads_A[fields[0]] + stoi(nm[2]);
        readdb_A[fields[0]].append(line + "\n");
        reads[fields[0]] = 1;
    }
    samAfile.close();
    
    ifstream samBfile(samB.c_str());
    for( std::string line; getline( samBfile, line ); )
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "[") {continue;}
        if (line.substr(0,1) == "@") {headB.append(line + "\n"); continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        if (fields[2] == "*") {continue;}
        if (fields[3] == "0") {continue;}
        int flag = stoi(fields[1]);
        if (flag > 200) {continue;}
        vector <string> nm;
        boost::split(nm,fields[11],boost::is_any_of(":"));
        if (stoi(nm[2]) > max_nm_limit) {continue;}
        reads_B[fields[0]] = reads_B[fields[0]] + stoi(nm[2]);
        readdb_B[fields[0]].append(line + "\n");
        reads[fields[0]] = 1;
    }
    samBfile.close();
   
    
    map <string,int> doublehits;
    string outsamA = currentout + "/" + sample + "." + currentgene + "." + bestAconverted + "." + type + ".sam";
    string outsamB = currentout + "/" + sample + "." + currentgene + "." + bestBconverted + "." + type + ".sam";
    if (bestAconverted != bestBconverted) {
        ofstream OUTA;
        ofstream OUTB;
        OUTA.open (outsamA);
        OUTB.open (outsamB);
        OUTA << headA;
        OUTB << headB;
        for (auto read : reads)
        {
            string item = read.first;
            int scoreA = 1000;
            int scoreB = 1000;
            if (reads_A.find(item) != reads_A.end()) {scoreA = reads_A[item];}
            if (reads_B.find(item) != reads_B.end()) {scoreB = reads_B[item];}
            if (scoreA < scoreB)
            {
                OUTA << readdb_A[item];
            }
            if (scoreB < scoreA)
            {
                OUTB << readdb_B[item];
            }
            if (scoreB == scoreA)
            {
                OUTA << readdb_A[item];
                OUTB << readdb_B[item];
                doublehits[item] = 1;
            }
        }
        OUTA.close();
        OUTB.close();
    }
    if (bestAconverted == bestBconverted) {
        string outsamA = currentout + "/" + sample + "." + currentgene + "." + bestAconverted + "." + type + ".sam";
        ofstream OUTA;
        OUTA.open (outsamA);
        OUTA << headA;
        for (auto read : reads)
        {
            string item = read.first;
            OUTA << readdb_A[item];
            doublehits[item] = 1;
        }
        OUTA.close();
    }
    

    
    string bamA = currentout + "/" + sample + "." + currentgene + "." + bestAconverted + "." + type + ".bam";
    string bamB = currentout + "/" + sample + "." + currentgene + "." + bestBconverted + "." + type + ".bam";
    cmd = v_samtools + " sort " + outsamA + " > " + bamA;
    GetStdoutFromCommand(cmd);
    cmd = v_samtools + " index " + bamA;
    GetStdoutFromCommand(cmd);
    cmd = v_samtools + " sort " + outsamB + " > " + bamB;
    GetStdoutFromCommand(cmd);
    cmd = v_samtools + " index " + bamB;
    GetStdoutFromCommand(cmd);
    cmd = v_samtools + " mpileup -f " + outrefA + " -A -aa --output-QNAME " + bamA;
    string pileA = GetStdoutFromCommand(cmd);
    cmd = v_samtools + " mpileup -f " + outrefB + " -A -aa --output-QNAME " + bamB;
    string pileB = GetStdoutFromCommand(cmd);

    
    
    map <string,string> alleles;
    alleles[alleleA] = bestAconverted;
    alleles[alleleB] = bestBconverted;
    
    string outRs = currentout + "/" + sample + "." + currentgene + "." + type + ".R";
    string png = currentout + "/" + sample + "." + currentgene + "." + type + ".png";
    
    std::ofstream outfile;
    outfile.open(outRs);
    outfile << endl << endl;

    outfile << "if (!require(\"ggplot2\")) install.packages(\"ggplot2\", repos = \"http://cran.us.r-project.org\")" << endl;
    outfile << "if (!require(\"patchwork\")) install.packages(\"patchwork\", repos = \"http://cran.us.r-project.org\")" << endl;
    outfile << "library(ggplot2)" << endl;
    outfile << "library(patchwork)" << endl;

    outfile << "setwd(\"" << currentout << "\")" << endl;
    outfile << "png(file=\"" << png << "\",width=1200, height=700)" << endl;

    for (auto allele : alleles)
    {
        string alleleconv = allele.second;
        string ref = references[allele.first];
        string otherpredicted = "";
        if (allele.first == alleleA) {otherpredicted = alleleB;}
        if (allele.first == alleleB) {otherpredicted = alleleA;}
        string pile = "";
        if (allele.first == alleleA) {pile = pileA;}
        if (allele.first == alleleB) {pile = pileB;}
        
        string cmd = "mkdir " + currentout + "/metrics/";
        GetStdoutFromCommand(cmd);
        
        string outcov = currentout + "/metrics/" + sample + "." + currentgene + "." + alleleconv + "." + type + ".cov.txt";
        string outcovadj = currentout + "/metrics/" + sample + "." + currentgene + "." + alleleconv + "." + type + ".covadj.txt";
        string outerror = currentout + "/metrics/" + sample + "." + currentgene + "." + alleleconv + "." + type + ".error.txt";
        string outredundance = currentout + "/metrics/" + sample + "." + currentgene + "." + alleleconv + "." + type + ".redundance.txt";
 
        ofstream OUT;
        OUT.open(outcov.c_str());
        ofstream OUTE;
        OUTE.open(outerror.c_str());
        ofstream OUTR;
        OUTR.open(outredundance.c_str());
        ofstream OUTA;
        OUTA.open(outcovadj.c_str());
        
        vector <string> piledata;
        boost::split(piledata,pile,boost::is_any_of("\n"));
        float mean = 0;
        float valid = 0;
        
        map <int,int> covdata;
        
        for (auto line : piledata)
        {
            if (line == "") {continue;}
            if (line.substr(0,1) == "[") {continue;}
            vector <string> data;
            boost::split(data,line,boost::is_any_of("\t"));
            int pos = stoi(data[1]);
            if (ref.substr((pos - 1),1) == "N") {
                OUT << data[0] << "\t" << pos << "\t" << "0" << endl;
                OUTE << data[0] << "\t" << pos << "\t" << "0" << endl;
                OUTR << data[0] << "\t" << pos << "\t" << "0" << endl;
                OUTA << data[0] << "\t" << pos << "\t" << "0" << endl;
                continue;
            }
            if (ref.substr((pos - 1),1) == "a") {continue;}
            if (ref.substr((pos - 1),1) == "t") {continue;}
            if (ref.substr((pos - 1),1) == "c") {continue;}
            if (ref.substr((pos - 1),1) == "g") {continue;}
            string aln = data[4];
            int total = stoi(data[3]);
            
            vector <string> readline;
            boost::split(readline,data[6],boost::is_any_of(","));
            float redundance_level_support = 0;
            for (auto read : readline) {
                if (doublehits.find(read) != doublehits.end()) {redundance_level_support++;}
            }
            float redundance_level = redundance_level_support / float(readline.size());
            
            int cov = 0;
            for (int a = 0; a <= aln.size(); a++)
            {
                if ((aln.substr(a,1) == ".") || (aln.substr(a,1) == ",")) {cov++;}
            }
            int error = total - cov;
            
            float cov_adj = int(float(cov) * (1 -(redundance_level * 0.5)));
            covdata[pos] = cov_adj;
            
            mean = mean + cov;
            valid++;
            OUT << data[0] << "\t" << pos << "\t" << cov << endl;
            OUTE << data[0] << "\t" << pos << "\t" << error << endl;
            OUTR << data[0] << "\t" << pos << "\t" << redundance_level << endl;
            OUTA << data[0] << "\t" << pos << "\t" << cov_adj << endl;
        }
        OUT.close();
        OUTA.close();
        OUTE.close();
        OUTR.close();
        
        float meanfinal = mean / valid;
        
        
        outfile << "cov_" + alleleconv + " <- read.table(\"" << outcov << "\",h=F)" << endl;
        outfile << "covadj_" + alleleconv + " <- read.table(\"" << outcovadj << "\",h=F)" << endl;
        outfile << "x_" + alleleconv + " = cov_" + alleleconv + "$V2" << endl;
        outfile << "y1_" + alleleconv + " = cov_" + alleleconv + "$V3" << endl;
        outfile << "y2_" + alleleconv + " = covadj_" + alleleconv + "$V3" << endl;
        
        outfile << "p_" + alleleconv + " <- ggplot(cov_" + alleleconv + ") + geom_area(aes(x=x_" + alleleconv + ", y=y1_" + alleleconv + "), fill = \"darkblue\", alpha = 0.85) + geom_area(aes(x=x_" + alleleconv + ", y=y2_" + alleleconv + "), fill = \"lightsteelblue1\", alpha = 0.85, linetype = 1, size =0.5 ,colour=\"black\") +";
        
        string title  = sample + ": " + allele.first + "       (other prediction for this sample: " + otherpredicted + ")\")";
        
        outfile << "ggtitle(\"" + title + " + xlab(\"Position\") + ylab(\"Read depth\") + ";
        outfile << "geom_hline(yintercept=" << meanfinal << ", linetype=\"dashed\", color = \"red\")";
        outfile << " + labs(caption = \"hla-mapper   red line: mean read depth   light blue: uniquely mapped reads    dark blue: multi-map reads\")";
        
        outfile << endl;
//        outfile << "p_" + alleleconv << endl;
        
        string line = "";
        string subline = "";
        for (auto pos : covdata)
        {
            float ratio = float(pos.second) / meanfinal;
            string item = to_string(pos.first) + ", " + to_string(pos.second) + "\n";
            if (ratio < 0.05) {line.append (item);}
            if (pos.second == 0) {subline.append(","+ to_string(pos.first));}
        }
        
        mtx_type.lock();
        pair <string,string> key = make_pair(currentgene, allele.first);
        if (subline != "") {subline = subline.substr(1);}
        if (type == "genomic") {summary_genomic[key] = subline;}
        if (type == "cds") {summary_cds[key] = subline;}
        mtx_type.unlock();
        
        if (line != "")
        {
            string resultout = currentout + "/" + sample + "." + currentgene + "." + type + ".results.txt";
            ofstream outres;
            outres.open(resultout, std::ios_base::app);
            outres << endl << allele.first << endl << "List of positions with coverage < 5% of the average depth:" << endl;
            outres << "Position,Depth" << endl;
            outres << line << endl;
            outres.close();
        }
    }
    
    string list  = "";
    for (auto allele : alleles)
    {
        string alleleconv = allele.second;
        list.append("/p_" + alleleconv);
    }
    outfile << "figure = " << list.substr(1) << endl;
    outfile << "figure" << endl;
    outfile << "dev.off()" << endl;
    outfile.close();
    cmd = "Rscript " + outRs + " > null";
    GetStdoutFromCommand(cmd);

    
}



void clean_files(string currentout, string currentgene, string sample)
{
    if (v_debug == 1) {return;}
    string cmd = "rm " + currentout + "/" + sample + "*.tmp.*";
    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "*.sam";
    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "." + currentgene + ".bam";
    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "." + currentgene + ".bam.bai";
    GetStdoutFromCommand(cmd);

    cmd = "rm " + currentout +  "/" + sample + "." + currentgene + ".alleleA.sam";
    GetStdoutFromCommand(cmd);

    cmd = "rm " + currentout +  "/" + sample + "." + currentgene + ".alleleB.sam";
    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "*part*";

    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "*whatshap*";
    GetStdoutFromCommand(cmd);
    cmd = "rm " + currentout +  "/" + sample + "*.fas.blast.txt";
    GetStdoutFromCommand(cmd);

    cmd = "rm " + currentout +  "/" + sample + "*.fastq";
    GetStdoutFromCommand(cmd);
    
}


std::string base_name(std::string const & path)
{
  return path.substr(path.find_last_of("/\\") + 1);
}


string check_fragment_size(string currentbam, string bed)
{
    vector <double> v;
    string cmd = v_samtools + " view " + currentbam + " -L " + bed;
    string sam = GetStdoutFromCommand(cmd);
    vector <string> data;
    boost::split(data,sam,boost::is_any_of("\n"));
    for (auto line : data)
    {
        if (line == "") {continue;}
        if (line.substr(0,1) == "[") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));

        if (map_type == "single"){v.push_back(fields[9].size());}
        if (map_type == "paired"){
            if (stoi(fields[8]) > 0) {v.push_back(stoi(fields[8]));}
        }
    }
    
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   std::bind2nd(std::minus<double>(), mean));
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size());
    return to_string(mean) + "\t" + to_string(stdev);
}



string do_typing_cds (string sample, string output, string gene, string passed_bam)
{
    
    if (v_debug) {mtx_type.lock(); cout << gene << ": enter typing" << endl << endl; mtx_type.unlock();}
    
    string currentout = output;
    string currentgene = gene;
    string currentbam = output + "/" + sample + "." + gene + ".bam";
    string bed = "";
    if (use_chr == "yes") {bed = resources + "/bed/" + currentgene + ".chr.cds.bed";}
    if (use_chr == "no") {bed = resources + "/bed/" + currentgene + ".nochr.cds.bed";}
    string vcftmp = currentout + "/" + sample + "." + currentgene + ".tmp.cds.vcf";
    string vcf = currentout + "/" + sample + "." + currentgene + ".cds.vcf";
    string whatstmp= currentout + "/" + sample + "." + currentgene + ".whatshap.cds.tmp.vcf";
    string whats = currentout + "/" + sample + "." + currentgene + ".whatshap.cds.vcf";
    string phased = currentout + "/" + sample + "." + currentgene + ".phased.cds.vcf";
    string phaselog = currentout + "/log/" + sample + "." + currentgene + ".phased.cds.log";

    

    float frag_mean = 0;
    float frag_std = 0;
    vector <string> frag;
    string fragsize = check_fragment_size(currentbam, bed);
    if (v_debug) {mtx_type.lock(); cout << gene << ": " << fragsize << endl << endl; mtx_type.unlock();}
    boost::split(frag,fragsize,boost::is_any_of("\t"));
    frag_mean = stof(frag[0]);
    frag_std = stof(frag[1]);
    
    if (v_debug) {mtx_type.lock(); cout << gene << ": frag done" << endl << endl; mtx_type.unlock();}

    
    //FREEBAYES
    string definebam = run_freebayes(currentgene,currentout, vcftmp, bed, sample, currentbam, "cds");
    currentbam = definebam;
    if ((! fileExists(vcftmp.c_str())) || (filesize(vcftmp.c_str()) == 0)){return "fail";}
    
    if (v_debug) {mtx_type.lock(); cout << gene << ": freebayes done" << endl << endl; mtx_type.unlock();}
    
    map <int,int> valid_positions;
    map <int,int> valid_position_limit;
    map <int,int> valid_segments;
    //LOADING BED
    ifstream bedfile(bed.c_str());
    int count = 1;
    for( std::string line; getline( bedfile, line ); )
    {
        if (line == "") {continue;}
        vector <string> fields;
        boost::split(fields,line,boost::is_any_of("\t"));
        int start = stoi(fields[1]);
        int end = stoi(fields[2]);
        valid_segments[start] = end;
        
        for (int a = start; a <= end; a++)
        {
            valid_positions[a] = count;
            valid_position_limit[a] = end;
        }
        count++;
    }
    bedfile.close();
    
    if (v_debug) {mtx_type.lock(); cout << gene << ": valid positions done" << endl << endl; mtx_type.unlock();}
    
    
    //CHECKING VCF
    checking_vcf_freebayes (vcftmp,vcf, "cds", bed);
    
    if (v_debug) {mtx_type.lock(); cout << gene << ": check vcf done" << endl << endl; mtx_type.unlock();}

    // WHATSHAP
    run_whatshap (vcf,whats,whatstmp,currentout,currentgene,sample,frag_mean,frag_std,currentbam, bed);

    if (v_debug) {mtx_type.lock(); cout << gene << ": whats done" << endl << endl; mtx_type.unlock();}

 
    //LOCAl PHASING
    if ((! fileExists(whats.c_str())) || (filesize(whats.c_str()) == 0)){return "fail";}
    
    
    if (v_multiphasing == 1) {
        run_local_phasing (whats,phased,currentout,currentgene,phaselog,frag_mean,currentbam);

        if (v_debug) {mtx_type.lock(); cout << gene << ": local phasing done" << endl << endl; mtx_type.unlock();}
    }
    else {   phased = whats;}
    
    //BLOCK_DEFINITION
    
    vector <string> snp_data;
    ifstream vcfphased(phased.c_str());
    int last_exon = 1;
    for( std::string line; getline( vcfphased, line ); )
    {
        if (line.substr(0,1) == "#") {continue;}
        if (line == "") {continue;}
        snp_data.push_back(line);
    }
    vcfphased.close();

    if (v_debug) {mtx_type.lock(); cout << gene << ": block definition done" << endl << endl; mtx_type.unlock();}
  
    
    map <int,int> break_positions;
    for (auto item : valid_segments)
    {
        int start = item.first;
        int end = item.second;
        string last_ps = "";
        for (auto snp : snp_data)
        {
            vector <string> data;
            boost::split(data,snp,boost::is_any_of("\t"));
            if (stoi(data[1]) < start) {continue;}
            if (stoi(data[1]) > end) {continue;}
            vector <string> fields;
            boost::split(fields,data[9],boost::is_any_of(":"));
            vector <string> alleles;
            boost::split(alleles,fields[0],boost::is_any_of("|"));
            if (alleles.size() < 2)
            {
                int allele_size = data[3].size();
                if(data[4].size() > allele_size){allele_size = data[4].size();}
                for (int b = stoi(data[1]); b < stoi(data[1]) + allele_size; b++)
                {
                    break_positions[b] = 1;
                }
                // arrumar aqui
                continue;
            }
            if (alleles[0] == alleles[1]){continue;}
            if (alleles[0] != alleles[1])
            {
                if (last_ps == "") {last_ps = fields[4];}
            }
            if ((fields[4] != last_ps) && (alleles[0] != alleles[1]))
            {
                int allele_size = data[3].size();
                if(data[4].size() > allele_size){allele_size = data[4].size();}
                for (int b = stoi(data[1]); b < stoi(data[1]) + allele_size; b++)
                {
                    break_positions[b] = 1;
                }
                //break_positions[stoi(data[1])] = 1;
                continue;
                
            }
            last_ps = fields[4];
        }
    }



 
 map <int,string> blocks;
 int block = 1;
 int start = 0;
 int end = 0;
 int end_limit = 0;
    
 for (auto position : valid_positions)
 {
     int pos = position.first;
     if (end_limit == 0) {end_limit = valid_position_limit[pos];}
     
     if (break_positions.find(pos) == break_positions.end())
     {
         if (start == 0) {start = pos;}
         end = pos;
         if (end == end_limit)
         {
             end = end_limit;
             end_limit = 0;
             if ((end - start) > limit_block_size) {blocks[block] = to_string(start) + "," + to_string(end);}
             start = 0;
             end = 0;
             block++;
             continue;
         }
         if ((end - start) > limit_block_size) {blocks[block] = to_string(start) + "," + to_string(end);}
     }
     
     if (break_positions.find(pos) != break_positions.end())
     {
         if (start == 0) {start = pos;}
         if ((end - start) > limit_block_size) {
             blocks[block] = to_string(start) + "," + to_string(end);
             block++;
         }
         start = 0;
         end = 0;
     }
 }

    

 for (auto item : blocks)
 {
     string fas = currentout + "/" + sample + "." + currentgene + ".cds.block_" + to_string(item.first) + ".fas";
     vector <string> limits;
     boost::split(limits,item.second,boost::is_any_of(","));
     generate_fasta(limits[0], limits[1], phased, fas);
 }
    if (v_debug) {mtx_type.lock(); cout << gene << ": generate fastq done" << endl << endl; mtx_type.unlock();}
  


    map <string,int> allele_size;
    string currentdb = resources + "/fasta/" + currentgene + "/" + currentgene + "_nuc.blast.fas";
    ifstream db (currentdb.c_str());
    for( std::string line; getline( db, line ); )
    {
        string id = line.substr(1);
        vector <string> ids;
        boost::split(ids,id,boost::is_any_of("&"));
        getline( db, line );
        string seq = line;
        allele_size[ids[0]] = allele_size[ids[0]] + seq.size();
    }
    db.close();
    
    
    
 string resultout = currentout + "/" + sample + "." + currentgene + ".cds.results.txt";
 ofstream myfile;
 myfile.open (resultout);
 
 map <string,string> comp_h1;
 map <string,string> comp_h2;
 map <pair<string,string>,float> comp;

 for (auto item : blocks)
 {
     string fas = currentout + "/" + sample + "." + currentgene + ".cds.block_" + to_string(item.first) + ".fas";
     
     ifstream fasfile(fas.c_str());
     map <string,float> sequence_size;
     for( std::string line; getline( fasfile, line ); )
     {
         string id = line.substr(1);
         getline( fasfile, line );
         sequence_size[id] = float(line.length());

         getline( fasfile, line );
         id = line.substr(1);
         getline( fasfile, line );
         sequence_size[id] = float(line.length());
     }
     fasfile.close();

     
     
     
     string currentdb = resources + "/fasta/" + currentgene + "/" + currentgene + "_nuc.blast.fas";
     string cmd = v_blast + " -query " + fas + " -db " + currentdb + " -outfmt 7 -max_target_seqs 20000";
     string blast = GetStdoutFromCommand(cmd);
     
     ofstream blastout;
     string blastfile = fas + ".blast.txt";
     blastout.open(blastfile.c_str());
     blastout << blast;
     blastout.close();
}
 
    currentdb = resources + "/fasta/" + currentgene + "/" + currentgene + "_nuc.blast.fas";
    string outfolder = currentout;
    string cmd = "perl " + resources + "/scripts/check_blast_data_cds.pl -g " + currentgene + " -b " + outfolder + " -r " + currentdb;
    string bestguess = GetStdoutFromCommand(cmd);
    if (bestguess == "") {return "fail";}
    
    
    //string fas = currentout + "/" + sample + "." + currentgene + ".cds.block_" + to_string(item.first) + ".fas";
    
    
    vector <string> bestguess_fields;
    boost::split(bestguess_fields,bestguess,boost::is_any_of("\t"));

    string bestA = bestguess_fields[0];
    string bestB = bestguess_fields[1];
    string schemeA = bestguess_fields[2];
    string schemeB = bestguess_fields[3];

 //    myfile << "Number of mismatches for the best guess: " << score << endl;
     myfile << "Combination of the fasta blocks (sequence 1): " << schemeA << endl;
     myfile << "Combination of the fasta blocks (sequence 2): " << schemeB  << endl;
     myfile << "Best guess for alleles: " << bestA << "\t" << bestB << endl;

    plot_graf(sample,currentout,currentgene,"cds",bestA,bestB,currentbam);

 string return_text = bestA + "," + bestB + "," + to_string(frag_mean) + "," + to_string(frag_std) + "," + currentbam;
 return return_text;
}








string check_chr(string bam)
{
    string use_chr = "no";
    string cmd = v_samtools + " view -H " + bam;
    string head = GetStdoutFromCommand(cmd);
    vector <string> head_data;
    boost::split(head_data,head,boost::is_any_of("\n"));
    for (auto line : head_data)
    {
        if (line == "") {continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of("\t"));

        if (data[0] == "@SQ")
        {
            vector <string> fields;
            boost::split(fields,data[1],boost::is_any_of(":"));
            if (fields[0] == "SN")
            {
                if (fields[1].find("chr") != string::npos) {use_chr = "yes"; break;}
            }
        }
    }
    return use_chr;
}




void main_typing ()
{
    int v_check = 0;
    if (v_mapout != "") {v_check = 1;}
    int disabled = 0;
    if (v_freebayes == "DISABLED") {disabled = 1; v_check = 0;}
    if (v_blast == "DISABLED") {disabled = 1; v_check = 0;}
    if (v_whats == "DISABLED") {disabled = 1; v_check = 0;}

    if ((v_db == "") || (v_check == 0))
    {
        screen_message (screen_size, 0, "", 1, v_quiet);
        v_message = "Program:   " + Program_name + "::type";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        v_message = "Version:   " + Program_version + ", " + Program_date;
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        if (disabled == 1)
        {
            v_message = "ATTENTION: " + Program_name + "::type is disabled. Please rerun the setup.";
            screen_message (screen_size, 0, v_message, 1, v_quiet);
            screen_message (screen_size, 0, "", 1, v_quiet);
        }
        
        v_message = "Usage:     hla-mapper type map=hla_mapper_output_folder db=path <options>";
        screen_message (screen_size, 0, v_message, 1, v_quiet);
        screen_message (screen_size, 0, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Mandatory options:", 1, v_quiet);
        screen_message (screen_size, 2, "map         output folder from hla-mapper dna", 1, v_quiet);
//        screen_message (screen_size, 2, "db          path to an hla-mapper database", 1, v_quiet);
        screen_message (screen_size, 2, "", 1, v_quiet);
        
        screen_message (screen_size, 0, "Other options:", 1, v_quiet);
        screen_message (screen_size, 2, "threads     number of threads [default: " + v_threads + "]", 1, v_quiet);
  
        /*
        if (v_bwa != "") {v_message = "[found at " + v_bwa + "]";}
        if (v_bwa == "") {v_message = "(!!! bwa not detected !!!)";}
        v_message = "bwa         path to BWA " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);
       
        if (v_samtools != "") {v_message = "[found at " + v_samtools + "]";}
        if (v_samtools == "") {v_message = "[!!!Samtools not detected!!!]";}
        v_message = "samtools    path to Samtools " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        
        if (v_gatk != "") {v_message = "[found at " + v_gatk + "]";}
        if (v_gatk == "") {v_message = "(!!! gatk not detected !!!)";}
        v_message = "gatk        path to GATK 3.8 " + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        if (v_whats != "") {v_message = "[found at " + v_whats + "]";}
        if (v_whats == "") {v_message = "(!!! WhatsHap not detected !!!)";}
        v_message = "whatshap    path to WhatsHap" + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);

        if (v_blast != "") {v_message = "[found at " + v_blast + "]";}
        if (v_blast== "") {v_message = "(!!! blast not detected !!!)";}
        v_message = "blast       path to Blast" + v_message;
        screen_message (screen_size, 2, v_message, 1, v_quiet);
  */
        screen_message (screen_size, 2, "--quiet           quiet mode", 1, v_quiet);
        screen_message (screen_size, 2, "--nomultiphasing  disable phasing multi-allelic variants", 1, v_quiet);
 //       screen_message (screen_size, 2, "--exome     Set parameters for exome sequencing", 1, v_quiet);

        screen_message (screen_size, 0, "", 1, v_quiet);
        return;
    }
    
    
    string cmd = "ls " + v_mapout + "/*.hla-mapper.log";
    string file = GetStdoutFromCommand(cmd);
    boost::replace_all(file, "\n", "");
    
    
    ifstream logfile;
    logfile.open(file.c_str());
    int check_alg = 0;
    int check_ver = 0;
    string sample = "";
    
    
    for( std::string line; getline( logfile, line ); )
    {
        if (line == "hla-mapper::dna") {check_alg = 1; continue;}
        vector <string> data;
        boost::split(data,line,boost::is_any_of(":"));
        if (line.find("Version 4") != std::string::npos) {check_ver = 1; continue;}
        if (data[0] == "Sample") {sample = data[1].substr(1); continue;}
        if (data[0] == "Type") {map_type = data[1].substr(1); break;}
    }
    logfile.close();
    
    
    
    if (map_type == "") {map_type = "paired";} //assuming-paired for map jobs before 4.1.2
    
//    if (map_type == "") {warnings.push_back("hla-mapper typing is only compatible hla-mapper 4 mapping jobs.");}
    if (check_alg == 0) {warnings.push_back("hla-mapper typing is only compatible with algorithm DNA and hla-mapper 4.");};
    if (check_ver == 0) {warnings.push_back("hla-mapper typing is only compatible hla-mapper 4 mapping jobs.");};
    if (sample == "") {warnings.push_back("Was this mapping folder generated by hla-mapper dna?");};
    
    if (((check_alg == 0) || (check_ver == 0)) || (sample == "")) {return;}
    if (map_type == "") {return;}
    
    string v_db_info = v_db + "/db_dna.info";
    boost::replace_all(v_db_info, "\\ ", " ");

    if (! fileExists(v_db_info))
    {
        v_message = "Could not access database " + v_db_info;
        cout << v_message << endl;
        warnings.push_back (v_message);
        v_sample = "";
        v_r1 = "";
        main_typing();
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
    
    if (db_version_ok == 0) {cout << "Incompatible database version." << endl; return;}
    

    
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = Program_name + "::type";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
 
    v_message = "Version " + Program_version;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
   
    v_message = "Number of threads: " + v_threads;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
        
    v_message = "Sample: " + sample;
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    v_message = "Running ... ";
    screen_message (screen_size, 0, v_message, 1, v_quiet);


//    string msg = "Progress...";
//    screen_message (screen_size, 0, msg, 0, v_quiet);


    bam = v_mapout + "/" + sample + ".adjusted.bam";
    
    
    //*************************
    resources = v_db + "/type/resources/";
    
    use_chr = check_chr(bam);
    if (use_chr == "no") {
        reference = resources + "/genome/6.fasta";
    }
    if (use_chr == "yes") {
        reference = resources + "/genome/chr6.fasta";
    }

    cmd = "mkdir " + v_mapout + "/typing";
    GetStdoutFromCommand(cmd);



    string typing_type = "";
    string sequence_master = "";
    string chr = "";

    load_reference();

    boost::algorithm::erase_all(v_genes, " ");
    boost::split(v_gene_list,v_genes,boost::is_any_of(","));
    
    todo = 0;
    for (auto gene : v_gene_list)
    {
        if (gene == "") {continue;}
        if (! fileExists(resources + "/fasta/" + gene)) {continue;}
        todo++;
    }
    
//    ThreadPool pool(stoi(v_threads));
    ThreadPool pool(1);
    std::vector< std::future<int> > results;
    
    v_gene_list.clear();
    v_gene_list.push_back("HLA-A");
    v_gene_list.push_back("HLA-B");
    v_gene_list.push_back("HLA-C");
    v_gene_list.push_back("HLA-E");
    v_gene_list.push_back("HLA-F");
    v_gene_list.push_back("HLA-G");
//    v_gene_list.push_back("HLA-DOA");
//    v_gene_list.push_back("HLA-DOB");
//    v_gene_list.push_back("MICA");
//    v_gene_list.push_back("MICB");
    v_gene_list.push_back("HLA-DRB1");
//    v_gene_list.push_back("HLA-DRB5");
//    v_gene_list.push_back("HLA-DRA");
//    v_gene_list.push_back("HLA-DMA");
//    v_gene_list.push_back("HLA-DMB");
//    v_gene_list.push_back("TAP1");
//    v_gene_list.push_back("TAP2");
    v_gene_list.push_back("HLA-DPA1");
    v_gene_list.push_back("HLA-DPB1");
    v_gene_list.push_back("HLA-DQA1");
    v_gene_list.push_back("HLA-DQB1");

        
    
    float todo = v_gene_list.size();

    for (auto gene : v_gene_list)
    {
        if (gene == "") {continue;}
        if (v_target != "") {
            if (gene != v_target){continue;}
        }

        if (! fileExists(resources + "/fasta/" + gene)) {continue;}

        mtx_type.lock();
        string currentout = v_mapout + "/typing/" + gene + "/";
        cmd = "mkdir " + currentout;
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();

        mtx_type.lock();
        cmd = "mkdir " + currentout + "/log";
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();

        mtx_type.lock();
        cmd = "mkdir " + currentout + "/log";
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();
        
        mtx_type.lock();
        cmd = "mkdir " + currentout + "/metrics";
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();

        mtx_type.lock();
        cmd = "mkdir " + currentout + "/phasing";
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();

        mtx_type.lock();
        cmd = "mkdir " + currentout + "/references";
        GetStdoutFromCommand(cmd);
        mtx_type.unlock();
        
 //       results.emplace_back(
 //                              pool.enqueue([sample, currentout, gene, todo]
 //               {
 
            string currentgene = gene;

            string bed = "";
            if (use_chr == "yes") {bed = resources + "/bed/" + currentgene + ".chr.genomic.bed";}
            if (use_chr == "no") {bed = resources + "/bed/" + currentgene + ".nochr.genomic.bed";}
 
            
            string currentbam = currentout + "/" + sample + "." + gene + ".bam";
            string cmd = v_samtools + " view -hb -L " + bed + " " + bam + " > " + currentbam;
            GetStdoutFromCommand(cmd);
            cmd = v_samtools + " index " + currentbam;
            GetStdoutFromCommand(cmd);
 
            //string out = downsample (currentbam, 75, currentout, currentgene, sample);
            //currentbam = out;
            string typing_type = check_coverage(currentgene, bed, currentbam);
//            if (typing_type == "fail") {return 1;}
            if (typing_type == "fail") {continue;}
            if ((v_exome == 1) && (typing_type != "fail")) {typing_type = "CDS";}
            if ((typing_type == "CDS") || (typing_type == "genomic")) {
                
                do_typing_cds(sample,currentout,gene,bam);
            }
            
            
                string genomic = "";
                mtx_type.lock();
                clean_files(currentout,gene,sample);
                string msg = "  ... " + currentgene + " done";
                screen_message (screen_size, 0, msg, 1, v_quiet);
                mtx_type.unlock();
//            return 1;
//        })
//        );
    }
//    for(auto && result: results){result.get();} // waiting for all threads
    
/*
    if (summary_genomic.size() >= 1) {
        ofstream summary;
        string sum = v_mapout + "/typing/" + sample + ".genomic.summary.txt";
        summary.open(sum.c_str());
        summary << "GENE\tALLELE\tERRORS" << endl;
        for (auto item : summary_genomic)
        {
            string gene = item.first.first;
            string allele = item.first.second;
            string error = item.second;
            summary << gene << "\t" << allele << "\t" << error << endl;
        }
        summary.close();
    }
*/
    if (summary_cds.size() >= 1) {
        ofstream summary;
        string sum = v_mapout + "/typing/" + sample + ".CDS.summary.txt";
        summary.open(sum.c_str());
        summary << "GENE\tALLELE\tERRORS" << endl;
        for (auto item : summary_cds)
        {
            string gene = item.first.first;
            string allele = item.first.second;
            string error = item.second;
            summary << gene << "\t" << allele << "\t" << error << endl;
        }
        summary.close();
    }
    
    screen_message (screen_size, 0, "Task completed.", 1, v_quiet);
 
}
