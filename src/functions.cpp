

#include "functions.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
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
#include <mutex>
#include <sys/stat.h>
#include <math.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <future>
#include <unistd.h>

#include "external.hpp"
#include "ThreadPool.hpp"


using namespace std;
mutex mtx;


void screen_message (int size, int left, string message, int enter, int quiet)
{
    if(quiet == 1) {return;}
    if ((message.length()+left) > size) {message = message.substr(0,(size-left-1));}
    cout << "\r";
    for (int a = 0; a < left; a++){cout << " ";}
    cout << message;
    for (int a = 0; a < (((size-left)-message.length())); a++){cout << " ";}
    std::cout.flush();
    if (enter == 1) {cout << endl;}
    return;
}

bool ends_with(const std::string &filename, const std::string &ext)
{
    return ext.length() <= filename.length() &&
    std::equal(ext.rbegin(), ext.rend(), filename.rbegin());
}

string GetStdoutFromCommand(string cmd) {
    
    string data;
    FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
    cmd.append(" 2>&1");
    
    stream = popen(cmd.c_str(), "r");
    if (stream) {
        while (!feof(stream))
            if (fgets(buffer, max_buffer, stream) != NULL) data.append(buffer);
        pclose(stream);
    }
    return data;
}


string findfilepath (string v_file)
{
    string v_filename = v_file.substr(0,(v_file.find_last_of("/"))+1);
    return (v_filename);
}


int phred (char v_char) {
    return int(v_char) - 33;
}







string mtrim (string seq, string qual)
{
    int start = -1;
    int end = -1;
    int low = 0;
    int high = 0;
    
    
    
    for (int $a = 0; $a < qual.length(); $a++) {
        double chr_phred = double(phred(qual.at($a)));
        //        double chr_phred = phred(qual.at($a));
        double P = pow(double(10),(chr_phred / double(-10)));
        
        if ((P <= v_mtrim_error) && (start == -1)) {
            start = $a;
            end = $a;
            continue;
        }
        
        if ((P > v_mtrim_error) && (start >= 0)) {
            end = $a - 1;
            if ((end - start) > (high - low)) {
                low = start;
                high = end;
            }
            start = -1;
            end = -1;
            continue;
        }
    }
    
    if (start >= 0) {
        end = qual.length();
        if ((end - start) > (high - low)) {
            low = start;
            high = end;
        }
    }
    
    if (low >= high)
    {
        seq = "";
        qual = "";
    }
    
    if (low < high)
    {
        seq = seq.substr(low,(high-low+1));
        qual = string(seq.size(),'A');
    }
    
    string v_limits = "";
    if (seq.size() >= v_size) {v_limits = seq + "\n+\n" + qual;}

    return (v_limits);
}






double filesize(const char *filename)
{
    FILE *f = fopen(filename,"rb");  /* open the file in read only */
    
    long size = 0;
    if (fseek(f,0,SEEK_END)==0) /* seek was successful */
        size = ftell(f);
    fclose(f);
    return size;
}



string splitcigar (string cigar)
{
    boost::replace_all(cigar, "S", "S,");
    boost::replace_all(cigar, "H", "H,");
    boost::replace_all(cigar, "M", "M,");
    boost::replace_all(cigar, "D", "D,");
    boost::replace_all(cigar, "I", "I,");
    boost::replace_all(cigar, "N", "N,");
    boost::replace_all(cigar, "P", "P,");
    boost::replace_all(cigar, "X", "X,");
    return cigar;
}

string splitNcigar (string cigar)
{
    boost::replace_all(cigar, "N", "N,");
    return cigar;
}



bool fileExists(const std::string& filename)
{
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1)
    {
        return true;
    }
    return false;
}

void removefile (string v_file, int v_debug)
{
    if (v_debug == 1) {return;}
    if (fileExists(v_file)) {
        const int result = remove(v_file.c_str());
    }
}


int is_num (string str)
{
    if (str == "0") {return 1;}
    if (str == "1") {return 1;}
    if (str == "2") {return 1;}
    if (str == "3") {return 1;}
    if (str == "4") {return 1;}
    if (str == "5") {return 1;}
    if (str == "6") {return 1;}
    if (str == "7") {return 1;}
    if (str == "8") {return 1;}
    if (str == "9") {return 1;}
    return 0;
}

string decompose_mpileup (string cigar)
{
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '$'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '^'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), ']'), cigar.end());
    cigar.erase(std::remove(cigar.begin(), cigar.end(), '~'), cigar.end());
//    cigar.erase(std::remove(cigar.begin(), cigar.end(), 'N'), cigar.end());
//    cigar.erase(std::remove(cigar.begin(), cigar.end(), 'n'), cigar.end());

    string new_cigar = "";
    for (int a = 0; a < cigar.size(); a++)
    {
        string sub = cigar.substr(a,1);
        
        string next;
        if (a < cigar.size()) {next = cigar.substr(a+1,1);}
        
        if ((next == "+") || (next == "-"))
        {
            int start_search = a+2;
            int end_search = a+2;
            for (end_search = start_search; end_search < cigar.size();end_search++) {if (is_num(cigar.substr(end_search,1)) == 0) {break;}}

            int get = end_search - start_search;
            int size_n = stoi(cigar.substr(start_search, get));
            string size_te = cigar.substr(start_search, get);
            

            string value = "";
            value.append(cigar.substr(a,1));
            value.append(next);
            value.append(size_te);
            a = a + 2 + size_te.size();
            value.append(cigar.substr(a,size_n));
            a = a + size_n - 1;
            new_cigar.append(value + ";");
            continue;
        }
        
        if (sub == ".") {new_cigar.append(sub + ";");continue;}
        if (sub == ",") {new_cigar.append(sub + ";");continue;}
        if (sub == "A") {new_cigar.append(sub + ";");continue;}
        if (sub == "T") {new_cigar.append(sub + ";");continue;}
        if (sub == "C") {new_cigar.append(sub + ";");continue;}
        if (sub == "G") {new_cigar.append(sub + ";");continue;}
        if (sub == "a") {new_cigar.append(sub + ";");continue;}
        if (sub == "t") {new_cigar.append(sub + ";");continue;}
        if (sub == "c") {new_cigar.append(sub + ";");continue;}
        if (sub == "g") {new_cigar.append(sub + ";");continue;}
        if (sub == "N") {new_cigar.append(sub + ";");continue;}
        if (sub == "n") {new_cigar.append(sub + ";");continue;}
        if (sub == "*") {new_cigar.append(sub + ";");continue;}
    }
    std::transform(new_cigar.begin(), new_cigar.end(),new_cigar.begin(), ::toupper);
    return new_cigar;
}


size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2)
{
    // from https://rosettacode.org/wiki/Levenshtein_distance#C.2B.2B
    
    const size_t m(s1.size());
    const size_t n(s2.size());
 
    if( m==0 ) return n;
    if( n==0 ) return m;
 
    size_t *costs = new size_t[n + 1];
 
    for( size_t k=0; k<=n; k++ ) costs[k] = k;
 
    size_t i = 0;
    for ( std::string::const_iterator it1 = s1.begin(); it1 != s1.end(); ++it1, ++i )
    {
        costs[0] = i+1;
        size_t corner = i;
 
        size_t j = 0;
        for ( std::string::const_iterator it2 = s2.begin(); it2 != s2.end(); ++it2, ++j )
        {
            size_t upper = costs[j+1];
            if( *it1 == *it2 )
            {
          costs[j+1] = corner;
      }
            else
      {
        size_t t(upper<corner?upper:corner);
                costs[j+1] = (costs[j]<t?costs[j]:t)+1;
      }
 
            corner = upper;
        }
    }
 
    size_t result = costs[n];
    delete [] costs;
 
    return result;
}



string reverse_and_complement (string seq)
{
    string rev = "";
    for (int a = seq.size(); a >= 0;a--)
    {
        string chr = seq.substr(a,1);
        if (chr == "T") {rev.append("A");continue;}
        if (chr == "A") {rev.append("T");continue;}
        if (chr == "G") {rev.append("C");continue;}
        if (chr == "C") {rev.append("G");continue;}
        if (chr == "N") {rev.append("N");continue;}
    }
    return rev;
}









string typing_dna (string gene, string fastq1, string fastq2, int maxselect, int maxerror, string maptype, float mincov)
{

        if (filesize(fastq1.c_str()) < 25000) {return "fail:low cov low size";}

    
        int motif_size = 20;
        unordered_map <string,int> motifs;

        ifstream list;
        string reffile = v_db + "/typing/dna/" + gene + "/reference/" + gene + ".fa";
        boost::replace_all(reffile, "\\ ", " ");

        if (! fileExists(reffile.c_str())) {gene_type[gene] = "."; return "disabled";}
        string vlog = v_output + "/log/typing_select.log";
        string vout = v_output + "/" + v_sample + gene + ".typing_search.sam";
        string command = "";
        if (maptype == "paired") {command = v_bwa + " mem -a -t " + v_threads + " '" + reffile + "' " + fastq1 + " " + fastq2 + " > " + vout + " 2>" + vlog;}
        if (maptype == "single") {command = v_bwa + " mem -a -t " + v_threads + " '" + reffile + "' " + fastq1 + " > " + vout + " 2>" + vlog;}
        system(command.c_str());

    
        map <string,float> reference;
        unordered_map <string,float> counterhits;
        unordered_map <string,float> error;
        vector <int> readsizes;
        
        ifstream samsearch;
        samsearch.open (vout.c_str());
        for( std::string item; getline( samsearch, item ); )
        {
                if (item == "") {continue;}
                if (item.substr(0,1) == "[") {continue;}

                if (item.substr(0,3) == "@SQ") {
                                vector <string> samdata;
                                boost::split(samdata,item,boost::is_any_of("\t"));
                                string target = samdata[1].substr(3);
                                string size = samdata[2].substr(3);
                                reference[target] = stoi(size) - 120;
                                continue;
                }
                if (item.substr(0,1) == "@") {continue;}
                
                vector <string> samdata;
                boost::split(samdata,item,boost::is_any_of("\t"));

                if (samdata[2] == "*") {continue;}
                if (samdata[9] != "*") {
                        readsizes.push_back(samdata[9].size());
                        for (int a = 0; a < (samdata[9].size() - motif_size); a++)
                        {
                                string sub = samdata[9].substr(a,motif_size);
                                motifs[sub]++;
                        }
                }
                
                if (samdata.size() > 10) {
                        if (samdata[11].substr(0,2) == "NM") {
                                int nm = stoi(samdata[11].substr(5));
                                if (nm > maxerror) {continue;}
                                counterhits[samdata[2]]++;
                                error[samdata[2]] = error[samdata[2]] + nm;
                        }
                }
        }
        samsearch.close();
        removefile(vout,v_debug);
        

        int sumsizes = 0;
        for (auto & item : readsizes) {sumsizes = sumsizes + item;}
        float meanreadsize = float(sumsizes) / float(readsizes.size());
        
        float highercov = 0;
        unordered_map <string,int> banned;
        for (auto & item : reference)
        {
                string allele = item.first;
                float allelecov = (float(counterhits[allele]) * meanreadsize) / float(item.second);
                if (highercov < allelecov) {highercov = allelecov;}
        }
        
        
        if (highercov < mincov) {return "fail:low cov after BWA";}
        
        for (auto & item : reference)
        {
                string allele = item.first;
                float allelecov = (float(counterhits[allele]) * meanreadsize) / float(item.second);
                if (allelecov < (0.7 * highercov)) {banned[allele] = 1; continue;}
        }


        string motif_file = v_db + "/typing/dna/" + gene + "/motifs/" + gene + ".txt";
        boost::replace_all(motif_file, "\\ ", " ");

        ifstream motif;
        motif.open(motif_file.c_str());
        unordered_map <string,int> accepted;
        for( std::string line; getline( motif, line ); )
        {
                vector <string> item;
                boost::split(item,line,boost::is_any_of("\t"));
                string allele = item[0];
                
                if (banned.find(allele) != banned.end()) {continue;}
                
                for (int a = 1; a < item.size(); a++)
                {
                        string sub = item[a];
                        if (motifs.find(sub) != motifs.end()) {
                                if (motifs[sub] >= (0.2 * highercov)) {
                                        accepted[allele] = 1;
                                        break;
                                }
                        }
                }
        }
        motif.close();
        
        


        ThreadPool select_pool(stoi(v_threads));
        std::vector< std::future<int> > select_step_results;
        map <float,string> scores;
        for (auto & target : counterhits)
        {
                string current = target.first;
                select_step_results.emplace_back(
                                        select_pool.enqueue([&scores,accepted,current,&counterhits,&error,&reference]
                {
                        if (accepted.find(current) == accepted.end()) {return 1;}
                        float cov = (counterhits[current]*100) / reference[current];
                        float noise = error[current] / counterhits[current];
                        float score = cov / noise;
                        mtx.lock();
                        scores[score].append("," + current);
                        mtx.unlock();
                        return 1;
                })
                );
        }
        for(auto && result: select_step_results){result.get();}
        

        
        int count = 0;
        vector <string> selected;
        int best_score_selection = 0;
        for (auto it = scores.rbegin(); it != scores.rend(); it++) {
                                if (best_score_selection == 0) {best_score_selection = it->first;}
                                vector <string> item;
                                boost::split(item,it->second,boost::is_any_of(","));
                                for (auto & allele : item)
                                {
                                        if (allele == "") {continue;}
                                        selected.push_back(allele);
                                        count++;
                                }
                                if (count > maxselect) {break;}
        }

        if (v_debug == 1) {
                for (auto & allele : selected)
                {
                        cout << endl << gene << " allele -> " << allele;
                }
        }
        
        if (selected.size() == 0) {return "fail: no preselect";}



        unordered_map <string,string> head;
        unordered_map <string,int> reads;
        map <pair<string,string>,string> mapdata;
        map <pair<string,string>,int> nmdata;
        
        ThreadPool loading_pool(stoi(v_threads));
        std::vector< std::future<int> > loading_results;
        
        for (auto & allele : selected)
        {
                string current = allele;
                loading_results.emplace_back(
                                        loading_pool.enqueue([current,&mapdata,&nmdata,&head,gene,fastq1,fastq2,&reads,maxerror,maptype]
                {
                        string file = current;
                        replace( file.begin(), file.end(), ':', '_');
                        replace( file.begin(), file.end(), '*', '_');
                        
                        string ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa'";
                        boost::replace_all(ref, "\\ ", " ");

                        string cmd = "";
                        if (maptype == "paired"){cmd = v_bwa + " mem -v 1 -a " + ref + " " + fastq1 + " " + fastq2;}
                        if (maptype == "single"){cmd = v_bwa + " mem -v 1 -a " + ref + " " + fastq1;}
                        string out = GetStdoutFromCommand(cmd);
                        
                        vector <string> lines;
                        boost::split(lines,out,boost::is_any_of("\n"));
                        for (auto & item: lines)
                        {
                                if (item == "") {continue;}
                                if (item.substr(0,1) == "[") {continue;}
                                if (item.substr(1,2) == "SQ")
                                {
                                        vector <string> data;
                                        boost::split(data,item,boost::is_any_of("\t"));
                                        string allele = data[1].substr(3);
                                        mtx.lock();
                                        head[allele] = item;
                                        mtx.unlock();
                                        continue;
                                }
                                if (item.substr(0,1) == "@") {continue;}
                                vector <string> data;
                                boost::split(data,item,boost::is_any_of("\t"));
                                if (data[2] == "*") {continue;}
                                if (data.size() >= 11)
                                {
                                        int nm = stoi(data[11].substr(5));
                                        pair <string,string> key = make_pair(data[2],data[0]);
                                        mtx.lock();
                                        mapdata[key].append("\n" + item);
                                        nmdata[key] = nmdata[key] + nm;
                                        reads[data[0]] = 1;
                                        mtx.unlock();
                                }
                        }
                        return 1;
                })
                );
        }
        for(auto && result: loading_results){result.get();}
        if (reads.size() == 0) {return "fail: no loaded data";}




        unordered_map <string,int> used;
        int counter = 0;
        
        string best_hit = "";
        float best_score = 0;
        
        ThreadPool pool(stoi(v_threads));
        std::vector< std::future<int> > results;

        for (auto & alleleA : selected)
        {
                string currentA = alleleA;
                
                for (auto & alleleB : selected)
                {
                        string currentB = alleleB;
                        if (used.find(currentB) != used.end()) {continue;}
                        results.emplace_back(
                                                pool.enqueue([counter, gene, &head, currentA, currentB, &reads, &nmdata, &mapdata, &reference, maxerror,&best_score, &best_hit]
                        {
                                        string samA = "";
                                        string samB = "";
                                        unordered_map <string,int> readsA;
                                        unordered_map <string,int> readsB;
                                        unordered_map <string,int> reads_general;

                                        samA.append(head[currentA]);
                                        samB.append(head[currentB]);
                                
                                        for (auto & readid : reads)
                                        {
                                                int nmA = 1000;
                                                int nmB = 1000;
                                                pair <string,string> keyA = make_pair(currentA,readid.first);
                                                pair <string,string> keyB = make_pair(currentB,readid.first);
                                                if (nmdata.find(keyA) != nmdata.end()) {nmA = nmdata[keyA];}
                                                if (nmdata.find(keyB) != nmdata.end()) {nmB = nmdata[keyB];}
                                                if ((nmA == 1000) && (nmB == 1000)) {continue;}
                                                if ((nmA > 1) && (nmB == 1000)){continue;}
                                                if ((nmB > 1) && (nmA == 1000)){continue;}
                                                if ((nmA > maxerror) && (nmB > maxerror)){continue;}

                                                if (nmA <= nmB) {samA.append(mapdata[keyA]);readsA[readid.first] = 1;}
                                                if (nmA >= nmB) {samB.append(mapdata[keyB]);readsB[readid.first] = 1;}
                                                reads_general[readid.first] = 1;
                                        }
                                        

                                        
                                        string fileA = currentA;
                                        replace( fileA.begin(), fileA.end(), ':', '_');
                                        replace( fileA.begin(), fileA.end(), '*', '_');
                                        string fileB = currentB;
                                        replace( fileB.begin(), fileB.end(), ':', '_');
                                        replace( fileB.begin(), fileB.end(), '*', '_');
                                     
                                        string ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + fileA + ".fa'";
                                        boost::replace_all(ref, "\\ ", " ");

                                        
                                        string cmd = "";
                                        string pileA = "";
                                        string outfile = v_output + "/" + v_sample + "_pair_" + fileA + "_" + fileB + ".1.sam";
                                        mtx.lock();
                                        ofstream out;
                                        out.open (outfile.c_str());
                                        out << samA << endl;
                                        out.close();
                                        mtx.unlock();
                                        cmd = v_samtools + " sort " + outfile + " | " + v_samtools + " mpileup --verbosity 0 -a --reference " + ref + " -";
                                        pileA = GetStdoutFromCommand(cmd);
                                        mtx.lock(); removefile(outfile,0); mtx.unlock();
                                        
                                        string pileB = "";
                                        if (currentA != currentB) {
                                                ref = "'" + v_db + "/typing/dna/" + gene + "/alleles/" + fileB + ".fa'";
                                                boost::replace_all(ref, "\\ ", " ");

                                                cmd = "";
                                                outfile = v_output + "/" + v_sample + "_pair_" + fileA + "_" + fileB + ".2.sam";
                                                mtx.lock();
                                                out.open (outfile.c_str());
                                                out << samB << endl;
                                                out.close();
                                                mtx.unlock();
                                                cmd = v_samtools + " sort " + outfile + " | " + v_samtools + " mpileup --verbosity 0 -a --reference " + ref + " -";
                                                pileB = GetStdoutFromCommand(cmd);
                                                mtx.lock(); removefile(outfile,0); mtx.unlock();
                                        }
                                
                                        vector <string> pileAdata;
                                        vector <string> pileBdata;
                                        boost::split(pileAdata,pileA,boost::is_any_of("\n"));
                                        boost::split(pileBdata,pileB,boost::is_any_of("\n"));
                                        
                                        float covA_sum;
                                        float covA_count;
                                        float errorA_sum = 1;
                                        float errorA_count = 1;
                                        
                                        for (auto item : pileAdata)
                                        {
                                                if (item == "") {continue;}
                                                if (item.substr(0,1) == "[") {continue;}
                                                vector <string> data;
                                                boost::split(data,item,boost::is_any_of("\t"));
                                                if (((data[2] == "N") || (data[2] == "0")) || (data[2] == "*")) {continue;}
                                                float cov = stof(data[3]);
                                                if (cov == 0) {continue;}
                                                string cigar = data[4];
                                                float ref = 0;
                                                for (int a = 0; a < cigar.size(); a++)
                                                {
                                                    if ((cigar.substr(a,1) == ".") || (cigar.substr(a,1) == ",")) {ref++;}
                                                }
                                                
                                                float indel = 0;
                                                for (int a = 0; a < cigar.size(); a++)
                                                {
                                                        if ((cigar.substr(a,1) == "+") || (cigar.substr(a,1) == "-")) {indel++;}
                                                }
                                                errorA_sum = errorA_sum + ((cov - ref) + indel);
                                                if (((cov - ref) + indel) > 0) {errorA_count++;}
                                        }
                                        

                                        float covB_sum;
                                        float covB_count;
                                        float errorB_sum = 1;
                                        float errorB_count = 1;
                                        
                                        if (currentA != currentB) {
                                            for (auto item : pileBdata)
                                            {
                                                    if (item == "") {continue;}
                                                    if (item.substr(0,1) == "[") {continue;}
                                                    vector <string> data;
                                                    boost::split(data,item,boost::is_any_of("\t"));
                                                    if (((data[2] == "N") || (data[2] == "0")) || (data[2] == "*")) {continue;}
                                                    float cov = stof(data[3]);
                                                    if (cov == 0) {continue;}
                                                    string cigar = data[4];
                                                    float ref = 0;
                                                    for (int a = 0; a < cigar.size(); a++)
                                                    {
                                                        if ((cigar.substr(a,1) == ".") || (cigar.substr(a,1) == ",")) {ref++;}
                                                    }
                                                    
                                                    float indel = 0;
                                                    for (int a = 0; a < cigar.size(); a++)
                                                    {
                                                            if ((cigar.substr(a,1) == "+") || (cigar.substr(a,1) == "-")) {indel++;}
                                                    }
                                                    errorB_sum = errorB_sum + ((cov - ref) + indel);
                                                    if (((cov - ref) + indel) > 0) {errorB_count++;}
                                            }
                                        }
                                        
                                        
                                        if (currentA == currentB) {
                                            errorB_sum = errorA_sum;
                                            errorB_count = errorA_count;
                                        }
                                        
                                        float errorA = errorA_sum / errorA_count;
                                        float errorB = errorB_sum / errorB_count;
                                        float errorsum = errorA + errorB;

                                        float size = (max(reference[currentA], reference[currentB])) / 1000;
                                        float sizekb = ((float)((int)(size * 100))) / 100;
                                        float abund = float(reads_general.size()) / sizekb;
                                        
                                        float score = abund / errorsum;
                                            
                                        
                                        mtx.lock();
                                        if (v_debug == 1) {
                                            cout << currentA << "\t" << currentB << "\t" << errorA << "\t" << errorB << "\t" << errorsum << "\t" << abund << "\t" << score << endl;
                                        }
                                        if (best_score == score) {best_hit.append("," + currentA + "\t" + currentB);}
                                        if (best_score < score) {best_score = score; best_hit = currentA + "\t" + currentB;}
                                        mtx.unlock();
                                        return counter;
                        })
                        );
                                counter++;
                                used[currentA] = 1;
                }
        }
        for(auto && result: results){result.get();}

        if (best_hit == "") {return "fail: no best hit";}
        
        vector <string> typed;
         if (best_hit != "") {
            vector <string> subbest;
            boost::split(subbest,best_hit,boost::is_any_of(","));
            best_hit = subbest[0];
            boost::split(typed,best_hit,boost::is_any_of("\t"));
         }
        
        
        
        /*
        if ((v_map_type == "paired") && (v_assembling == 1)){
                    string config_assemble_A_file = v_output + v_sample + gene + ".assemble.config";
                    ofstream config_assemble_A;
                    config_assemble_A.open (config_assemble_A_file.c_str());
                   
                    config_assemble_A << "project = " << gene << "_assemble\n";
                    config_assemble_A << "job = genome,mapping,accurate\n";
                    config_assemble_A << "parameters = -NW:cdrn=no -GE:not=4 -SK:mmhr=3 -CO:mroir=yes -ED:mace=no -MI:ef1=yes -NW:cmrnl=no -NW:cac=warn ";
                    config_assemble_A << " SOLEXA_SETTINGS -CL:c3ppmsl=12:c3ppmea=3 -AL:egpl=low -CO:msr=no\n";
                    config_assemble_A << "\n";
                    config_assemble_A << "readgroup = gaps500\n";
                    config_assemble_A << "is_reference\n";

                    string fileA = typed[0];
                    replace( fileA.begin(), fileA.end(), ':', '_');
                    replace( fileA.begin(), fileA.end(), '*', '_');
                    string refA = v_db + "/typing/dna/" + gene + "/alleles/" + fileA + ".fa";

                    string fileB = typed[1];
                    replace( fileB.begin(), fileB.end(), ':', '_');
                    replace( fileB.begin(), fileB.end(), '*', '_');
                    string refB = v_db + "/typing/dna/" + gene + "/alleles/" + fileB + ".fa";
                    
                    string vrefout = v_output + "/" + v_sample + gene + ".assemble.ref.fa";
                    if (typed[0] != typed[1])
                    {
                        string cmd = "cat " + refA + " " + refB + " > " + vrefout;
                        GetStdoutFromCommand(cmd);
                        
                    }
                    else
                    {
                        string cmd = "cat " + refA + " > " + vrefout;
                        GetStdoutFromCommand(cmd);
                    }
                    
                    config_assemble_A << "data = " << vrefout << endl;
                    config_assemble_A << "technology = text\n";
                    config_assemble_A << "\n";
                    config_assemble_A << "readgroup = Samples\n";

                    string fq = " " + v_output + v_sample + gene + ".trim.R" + "*.fastq";
                    config_assemble_A << "data = " << fq << endl;
                    config_assemble_A << "technology = solexa\n";
                    config_assemble_A << "template_size = 300 1000\n";
                    config_assemble_A << "segment_placement = ---> <--- infoonly\n";
                    config_assemble_A << "segment_naming = solexa\n";
                    config_assemble_A.close();
                    

                    v_command = "mkdir " + v_output + gene + "_assemble/";
                    GetStdoutFromCommand(v_command);
                    string log = v_output + "/log/" + gene + ".mira.log";
                    v_command = v_mira + " --thread=" + v_threads + " --cwd=" + v_output + "/" + gene + "_assemble/ " + config_assemble_A_file + " > " + log + " 2>" + log;
                    system(v_command.c_str());
                    removefile(config_assemble_A_file,v_debug);
                    removefile(vrefout,v_debug);
                    
                    
                    map <string,string> assemble_result;
                    string assembled = v_output + gene + "_assemble/" + gene + "_assemble_assembly/" + gene + "_assemble_d_results/";
                    assembled = assembled + gene + "_assemble_out_AllStrains.padded.fasta";
         
                    ifstream assembly ( assembled.c_str() );
                    string id = "";
                    for( std::string line; getline( assembly, line ); )
                    {
                        if (line.substr(0,1) == ">") {id = line.substr(1,line.size()-4); continue;}
                        for (int a = 0; a < line.size(); a++)
                        {
                            if ((line.substr(a,1) == "T") || (line.substr(a,1) == "t")) {assemble_result[id].append("T");continue;}
                            if ((line.substr(a,1) == "A") || (line.substr(a,1) == "a")) {assemble_result[id].append("A");continue;}
                            if ((line.substr(a,1) == "C") || (line.substr(a,1) == "c")) {assemble_result[id].append("C");continue;}
                            if ((line.substr(a,1) == "G") || (line.substr(a,1) == "g")) {assemble_result[id].append("G");continue;}
                            assemble_result[id].append("N");
                        }
                    }
                    assembly.close();
                    
                    
                    string fq_assemble = v_output + v_sample + gene + ".assemble.fastq";
                    string fas_assemble = v_output + v_sample + gene + ".assemble.fasta";
                    ofstream OUTFQ;
                    OUTFQ.open (fq_assemble.c_str());
                    ofstream OUTFA;
                    OUTFA.open (fas_assemble.c_str());
                    
                    
                    
                    int allele = 1;
                    for (auto & item : assemble_result)
                    {
                        string seq = item.second;
                        vector <string> segments;
                        boost::split(segments,seq,boost::is_any_of("N"));
                        int segnumber = 1;
                        for (auto & seg : segments)
                        {
                            if (seg.size() < 50) {continue;}
                            OUTFA << ">allele_" << allele << "-segment_" << segnumber << endl;
                            OUTFA << seg << endl;
                            for (int a = 1; a <= 10; a++)
                            {
                                OUTFQ << "@allele_" << allele << ":segment_" << segnumber << ":copy_" << a << endl;
                                OUTFQ << seg << endl;
                                OUTFQ << "+" << endl;
                                OUTFQ << string(seg.size(),'A') << endl;
                            }
                            segnumber++;
                        }
                        allele++;
                    }
                    OUTFQ.close();
                    OUTFA.close();
                }
        */
        
        if (best_hit != "")
        {
            string fq;
            if (v_map_type == "single") {fq = " " + v_output + v_sample + gene + "_R0.fastq";}
            if (v_map_type == "paired") {fq = " " + v_output + v_sample + gene + "_R1.fastq " + v_output + v_sample + gene + "_R2.fastq";}
            
            vector <string> typed;
            boost::split(typed,best_hit,boost::is_any_of("\t"));
            string file = typed[0];
            replace( file.begin(), file.end(), ':', '_');
            replace( file.begin(), file.end(), '*', '_');
            string ref = " '" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa' ";
            boost::replace_all(ref, "\\ ", " ");

            v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
            string out = GetStdoutFromCommand(v_command);
            
            
            unordered_map <string,int> nm_allele_A;
            unordered_map <string,int> nm_allele_B;
            
            vector <string> samdata;
            boost::split(samdata,out,boost::is_any_of("\n"));
            for (auto & item : samdata)
            {
                if (item == "") {continue;}
                if (item.substr(0,1) == "[") {continue;}
                vector <string> data;
                boost::split(data,item,boost::is_any_of("\t"));
                    
                if (data[2] == "*") {continue;}
                if (data.size() >= 11)
                {
                    int nm = stoi(data[11].substr(5));
                    nm_allele_A[data[0]] = nm_allele_A[data[0]] + nm;
                }
                

                
            }
            for (auto &item : nm_allele_A)
            {
                if (item.second <= 2)
                {
                    pair <string,string> key = make_pair(gene,item.first);
                    reads_possibly_mapped[key] = item.second;
                }
            }
 
            if (typed[0] != typed[1]) {
                file = typed[1];
                replace( file.begin(), file.end(), ':', '_');
                replace( file.begin(), file.end(), '*', '_');
                ref = " '" + v_db + "/typing/dna/" + gene + "/alleles/" + file + ".fa' ";
                boost::replace_all(ref, "\\ ", " ");

                v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
                    
                out = GetStdoutFromCommand(v_command);
                
                samdata.clear();
                boost::split(samdata,out,boost::is_any_of("\n"));
                for (auto & item : samdata)
                {
                    if (item == "") {continue;}
                    if (item.substr(0,1) == "[") {continue;}
                            
                    vector <string> data;
                    boost::split(data,item,boost::is_any_of("\t"));
                    
                            
                    if (data[2] == "*") {continue;}
                    if (data.size() >= 11)
                    {
                        int nm = stoi(data[11].substr(5));
                        nm_allele_B[data[0]] = nm_allele_B[data[0]] + nm;
                    }
                }
                    
                for (auto &item : nm_allele_B)
                {
                    if (item.second <= 2 )
                    {
                        pair <string,string> key = make_pair(gene,item.first);
                        if (reads_possibly_mapped.find(key) != reads_possibly_mapped.end())
                        {
                            if (reads_possibly_mapped[key] > item.second)
                            {
                                reads_possibly_mapped[key] = item.second;
                            }
                        }
                        else {reads_possibly_mapped[key] = item.second;}
                    }
                }
            }
        }
        

        if (best_hit != "") {return best_hit;}
        return "fail: no best hit";
}









string typing_rna (string gene, string fastq1, string fastq2)
{
    string cmd = "mkdir " + v_output + "/rna_type/";
    string out = GetStdoutFromCommand(cmd);
    string v_out = v_output + "/rna_type/";
    
    string perl = "'" + v_db + "/scripts/type_rna.pl'";
    boost::replace_all(perl, "\\ ", " ");
    
    string ref = "'" + v_db + "/typing/rna/type_db/'";
    boost::replace_all(ref, "\\ ", " ");

    cmd = "perl " + perl + " -d " + ref + " -g " +  gene + " -o " +  v_out + " -a " + fastq1 + " -b " + fastq2 + " -t " + v_threads;
    out = GetStdoutFromCommand(cmd);

    
    
    if (v_debug == 1)
    {
        cout << endl << out << endl;
    }
    
    vector <string> type_return;
    boost::split(type_return,out,boost::is_any_of("\n"));

    string alleleA = "";
    string alleleB = "";
    
    int help = 0;
    for (auto & item : type_return)
    {
        if (item == "") {continue;}
        if (item.substr(0,4) == "Best") {help = 1; continue;}
        if (help == 1)
        {
            vector <string> type_data;
            boost::split(type_data,item,boost::is_any_of("\t"));
            alleleA = type_data[0];
            alleleB = type_data[1];
            break;
        }

    }
    

    if ((alleleA == "") || (alleleB == ""))
    {
        v_command = "rm -rf " + v_out;
        system (v_command.c_str());
        return "fail";
    }
    
    string best_hit = alleleA + "\t" + alleleB;
    if (best_hit != "")
    {
        string fq;
        if (v_map_type == "single") {fq = " " + v_output + v_sample + gene + "_R0.fastq";}
        if (v_map_type == "paired") {fq = " " + v_output + v_sample + gene + "_R1.fastq " + v_output + v_sample + gene + "_R2.fastq";}
        
        vector <string> typed;
        boost::split(typed,best_hit,boost::is_any_of("\t"));
        string file = typed[0];
        replace( file.begin(), file.end(), ':', '_');
        replace( file.begin(), file.end(), '*', '_');
        string ref = v_out + "/" + gene + "/reftmp_" + file + ".fas";
        string log = v_out + "/" + gene + "/reftmp_" + file + ".log";
        v_command = v_bwa + " mem -t " + v_threads + " " + ref + " " + fq + " 2> " + log;
        out = GetStdoutFromCommand(v_command);
        

        unordered_map <string,int> nm_allele_A;
        unordered_map <string,int> nm_allele_B;
        
        vector <string> data;
        boost::split(data,out,boost::is_any_of("\n"));

 
        for( auto & line : data)
        {
                if (line == "") {continue;}
                if (line.substr(0,1) == "@") {continue;}
                if (line.substr(0,1) == "[") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                if (data[2] == "*") {continue;}
                if (data.size() >= 11)
                {
                    int nm = stoi(data[11].substr(5));
                    nm_allele_A[data[0]] = nm_allele_A[data[0]] + nm;
                }
        }
        
        for (auto &item : nm_allele_A)
        {
            if (item.second <= 2)
            {
                pair <string,string> key = make_pair(gene,item.first);
                reads_possibly_mapped[key] = item.second;
            }
        }

        
        
        
        if (typed[0] != typed[1]) {
            
            file = typed[1];
            replace( file.begin(), file.end(), ':', '_');
            replace( file.begin(), file.end(), '*', '_');
            string ref = v_out + "/" + gene + "/reftmp_" + file + ".fas";
            string log = v_out + "/" + gene + "/reftmp_" + file + ".log";

  
            v_command = v_bwa + " mem -t " + v_threads + " " + ref + " " + fq + " 2> " + log;
            out = GetStdoutFromCommand(v_command);


            vector <string> data;
            boost::split(data,out,boost::is_any_of("\n"));

            for( auto & line : data)
            {
                if (line == "") {continue;}
                if (line.substr(0,1) == "@") {continue;}
                if (line.substr(0,1) == "[") {continue;}
                vector <string> data;
                boost::split(data,line,boost::is_any_of("\t"));
                
                if (data[2] == "*") {continue;}
                if (data.size() >= 11)
                {
                    int nm = stoi(data[11].substr(5));
                    nm_allele_B[data[0]] = nm_allele_B[data[0]] + nm;
                }
            }

                
            for (auto &item : nm_allele_B)
            {
                if (item.second <= 1 )
                {
                    pair <string,string> key = make_pair(gene,item.first);
                    if (reads_possibly_mapped.find(key) != reads_possibly_mapped.end())
                    {
                        if (reads_possibly_mapped[key] > item.second)
                        {
                            reads_possibly_mapped[key] = item.second;
                        }
                    }
                    else {reads_possibly_mapped[key] = item.second;}
                }
            }
        }
        
    }
   


    v_command = "rm -rf " + v_out;
    system (v_command.c_str());
    
    
    if (best_hit != "")
    {
        string fq;
        if (v_map_type == "single") {fq = " " + v_output + v_sample + gene + "_R0.fastq";}
        if (v_map_type == "paired") {fq = " " + v_output + v_sample + gene + "_R1.fastq " + v_output + v_sample + gene + "_R2.fastq";}
        
        vector <string> typed;
        boost::split(typed,best_hit,boost::is_any_of("\t"));
        string file = typed[0];
        replace( file.begin(), file.end(), ':', '_');
        replace( file.begin(), file.end(), '*', '_');
        string ref = " '" + v_db + "/typing/rna/" + gene + "/alleles/" + file + ".fa' ";
        boost::replace_all(ref, "\\ ", " ");

        v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
        string out = GetStdoutFromCommand(v_command);
        
        
        unordered_map <string,int> nm_allele_A;
        unordered_map <string,int> nm_allele_B;
        
        vector <string> samdata;
        boost::split(samdata,out,boost::is_any_of("\n"));
        for (auto & item : samdata)
        {
            if (item == "") {continue;}
            if (item.substr(0,1) == "[") {continue;}
            vector <string> data;
            boost::split(data,item,boost::is_any_of("\t"));
                
            if (data[2] == "*") {continue;}
            if (data.size() >= 11)
            {
                int nm = stoi(data[11].substr(5));
                nm_allele_A[data[0]] = nm_allele_A[data[0]] + nm;
            }
            

            
        }
        for (auto &item : nm_allele_A)
        {
            if (item.second <= 2)
            {
                pair <string,string> key = make_pair(gene,item.first);
                reads_possibly_mapped[key] = item.second;
            }
        }


        if (typed[0] != typed[1]) {
            file = typed[1];
            replace( file.begin(), file.end(), ':', '_');
            replace( file.begin(), file.end(), '*', '_');
            ref = " '" + v_db + "/typing/rna/" + gene + "/alleles/" + file + ".fa' ";
            boost::replace_all(ref, "\\ ", " ");

            v_command = v_bwa + " mem -v 1 -t " + v_threads + ref + fq;
                
            out = GetStdoutFromCommand(v_command);
            
            
            samdata.clear();
            boost::split(samdata,out,boost::is_any_of("\n"));
            for (auto & item : samdata)
            {
                if (item == "") {continue;}
                if (item.substr(0,1) == "[") {continue;}
                        
                vector <string> data;
                boost::split(data,item,boost::is_any_of("\t"));
                
                        
                if (data[2] == "*") {continue;}
                if (data.size() >= 11)
                {
                    int nm = stoi(data[11].substr(5));
                    nm_allele_B[data[0]] = nm_allele_B[data[0]] + nm;
                }
            }
                
            for (auto &item : nm_allele_B)
            {
                if (item.second <= 2 )
                {
                    pair <string,string> key = make_pair(gene,item.first);
                    if (reads_possibly_mapped.find(key) != reads_possibly_mapped.end())
                    {
                        if (reads_possibly_mapped[key] > item.second)
                        {
                            reads_possibly_mapped[key] = item.second;
                        }
                    }
                    else {reads_possibly_mapped[key] = item.second;}
                }
            }
        }
    }
    
  
    if ((alleleA != "") && (alleleB != "")) {return best_hit;}
    
    return "fail";
}


unsigned long long getTotalSystemMemory()
{
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    return pages * page_size;
}
