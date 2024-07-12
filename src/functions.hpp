//
//  functions.hpp
//  hla-mapper
//
//  Created by Erick Castelli on 20/02/20.
//  Copyright © 2020 GeMBio.Unesp. All rights reserved.
//

#ifndef functions_hpp
#define functions_hpp

#include <stdio.h>
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

using namespace std;

void screen_message (int size, int left, string message, int enter, int quiet);
bool ends_with(const std::string &filename, const std::string &ext);
string GetStdoutFromCommand(string cmd);
string findfilepath (string v_file);
int phred (char v_char);
string mtrim (string seq, string qual);
double filesize(const char *filename);
string splitcigar (string cigar);
void removefile (string v_file, int v_debug);
bool fileExists(const std::string& filename);
int is_num (string str);
string decompose_mpileup (string cigar);
size_t uiLevenshteinDistance(const std::string &s1, const std::string &s2);
string reverse_and_complement (string seq);
string typing_dna (string gene, string fastq1, string fastq2, int maxselect, int maxerror, string maptype, float mincov);
string typing_rna (string gene, string fastq1, string fastq2);
unsigned long long getTotalSystemMemory();
string splitNcigar (string cigar);
#endif /* functions_hpp */
