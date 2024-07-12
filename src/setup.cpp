//
//  setup.cpp
//  hla-mapper
//
//  Created by Erick Castelli on 23/08/21.
//  Copyright © 2021 GeMBio.Unesp. All rights reserved.
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

#include "setup.hpp"
#include "external.hpp"
#include "functions.hpp"


void main_setup() {
    
    screen_message (screen_size, 0, "", 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
    v_message = "Program:   " + Program_name + "::setup";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "Version:   " + Program_version + ", " + Program_date;
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);


    v_message = "hla-mapper needs some third-party software to work properly.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    v_message = "This is the list of the third-party software and their tested versions.";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);

    
    v_message = "For hla-mapper dna:";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    string list = "  Samtools, version(s)";
    for (auto item : samtools_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);

    
    list = "  BWA, version(s)";
    for (auto item : bwa_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);


    
    v_message = "For hla-mapper rna (same as dna plus):";
    screen_message (screen_size, 0, v_message, 1, v_quiet);

    list = "  STAR, version(s)";
    for (auto item : star_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);


    v_message = "For hla-mapper type (same as dna plus):";
    screen_message (screen_size, 0, v_message, 1, v_quiet);
    
    list = "  Freebayes, version(s)";
    for (auto item : freebayes_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    
    
    
    list = "  WhatsHap, version(s)";
    for (auto item : whats_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    
    
    list = "  BLAST, version(s)";
    for (auto item : blast_versions)
    {
        list.append (" " + item.first);
    }
    screen_message (screen_size, 0, list, 1, v_quiet);
    screen_message (screen_size, 0, "", 1, v_quiet);
 
    
    
    cout << "Shall we proceed with the setup process? yes (y) or close/cancel (c): ";
    string res = "";
    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res != "y") && (res != "Y")) {return;}

    
    
    
    
    

    
// ********** SAMTOOLS  ************
    screen_message (screen_size, 0, "", 1, v_quiet);
    string program = "Samtools";
    string command = "which samtools";
    string samtools = "";
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for all hla-mapper algorithms." << endl;
    string item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    
samtools_start:
    if (item != "")
    {
        cout << program << " path (automatic detection): " << item << endl;
        cout << "use this one (y), or change (c), or quit (q): ";
        string res = "";
        cin >> res;
        res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
        int valid = 0;
        if ((res == "y") || (res == "Y")) {valid = 1; samtools = item;}
        if ((res == "c") || (res == "C")) {valid = 1; item = "";}
        if ((res == "q") || (res == "Q")) {return;}
        if (valid == 0) {goto samtools_start;}
    }

samtools:
    if (item == "")
    {
        while (! fileExists(item)) {
            cout << "Please write the " << program << " path, or drag the binary below (or q to quit): " << endl;
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            if ((item == "q") || (item == "Q")) {cout << "Quiting!" << endl; return;}
            if (! fileExists(item))
            {
                cout << "Can't find this file." << endl;
            }
            else {cout << program << " path: " << item << endl;}
        }
    }
    samtools = item;

    item = GetStdoutFromCommand(samtools + " --version");
    vector <string> data;
    boost::split(data,item,boost::is_any_of("\n"));
    vector <string> sub;
    boost::split(sub,data[0],boost::is_any_of(" "));
    if (samtools_versions.find(sub[1]) != samtools_versions.end())
        {
            cout << "This Samtools version (" << sub[1] << ") is compatible with hla-mapper" << endl;
        }
    else {
        cout << "This Samtools version has not been tested with hla-mapper\nContinue anyway (y), quit (q) or return (r)?";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1;item = ""; goto samtools;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto samtools;}
    }

    

    
    
    
    
    
    
// ********** BWA  ************
    program = "BWA";
    command = "which bwa";
    string bwa = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for all hla-mapper algorithms." << endl;
    item = GetStdoutFromCommand(command);
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
 
bwa_start:
    if (item != "")
        {
            cout << program << " path (automatic detection): " << item << endl;
            cout << "use this one (y) or change (c): ";
            string res = "";
            cin >> res;
            int valid = 0;
            if ((res == "y") || (res == "Y")) {valid = 1; bwa = item;}
            if ((res == "c") || (res == "C")) {valid = 1;item = "";}
            if (valid == 0) {goto bwa_start;}
        }
    
bwa:
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path, or drag the binary bellow (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << "Quiting!" << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Can't find this file." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    
    bwa = item;
    item = GetStdoutFromCommand(bwa);
    data.clear();
    boost::split(data,item,boost::is_any_of("\n"));
    sub.clear();
    boost::split(sub,data[2],boost::is_any_of(" "));
    boost::split(data,sub[1],boost::is_any_of("-"));
    if (bwa_versions.find(data[0]) != bwa_versions.end())
        {
            cout << "This BWA version (" << data[0] << ") is compatible with hla-mapper" << endl;
        }
    else {
        cout << "This BWA version has not been tested with hla-mapper\nContinue anyway (y), quit (q), or return (r)? ";
        cin >> item;
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        int valid = 0;
        if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
        if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto bwa;}
        if ((item == "y") || (item == "Y")) {valid = 1;}
        if (valid == 0) {goto bwa;}
    }

    

    
    
    
    
    
// ********** Database  ************
    
    item = "";
    program = "Database";
    string db = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for all hla-mapper algorithms." << endl;
    
    if (item == "")
        {
            while (! fileExists(item)) {
                cout << "Please write the " << program << " path or drag the folder bellow (or q to quit): " << endl;
                cin >> item;
                item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                if ((item == "q") || (item == "Q"))  {cout << "Quiting!" << endl; return;}
                if (! fileExists(item))
                    {
                        cout << "Can't find this folder." << endl;
                    }
                else {cout << program << " path: " << item << endl << endl;}
            }
        }
    db = item;
    
    

    
    
    
    
    
    
string version = "";
    
// ********** Freebayes  ************
    
    program = "Freebayes";
    command = "which freebayes";
    string freebayes = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for hla-mapper type." << endl;
    cout << "If " + program + " is not available, [hla-mapper type] will be disabled." <<endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {freebayes = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }

    if (((res == "y") || (res == "Y")) && (freebayes != "DISABLED"))
    {
        item = GetStdoutFromCommand(command);
    try
    {
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    }
        catch (const std::exception& e)
    { item = "";}
    
    freebayes_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "use this one (y), or change (c), or quit (q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "Y")) {valid = 1; freebayes = item;}
                if ((res == "c") || (res == "C")) {valid = 1;item = "";}
                if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }
                if (valid == 0) {goto freebayes_start;}
            }
        
    freebayes:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write the " << program << " path, or drag the binary bellow (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if ((item == "q") || (item == "Q")) {cout << "Quiting!" << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Can't find this file." << endl;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        freebayes = item;
        
        try
        {
            item = GetStdoutFromCommand(freebayes);
        }
            catch (const std::exception& e)
        { item = ""; freebayes = "";}
        
        data.clear();
        
        try
        {
            boost::split(data,item,boost::is_any_of("\n"));
        }
            catch (const std::exception& e)
        { data.clear();}
        
        
        string version = "";
        for (auto value : data)
        {
            vector <string> version_data;
            boost::split(version_data,value,boost::is_any_of(":"));
            if (version_data[0] == "version")
            {
                version = version_data[1];
                version.erase(remove(version.begin(), version.end(), ' '), version.end());
            }
        }
        
        if (freebayes_versions.find(version) != freebayes_versions.end())
            {
                cout << "This Freebayes version (" << version << ") is compatible with hla-mapper" << endl;
            }
        else {
            cout << "This Freebayes version has not been tested with hla-mapper\nContinue anyway (y), quit (q), or return (r)? ";
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            int valid = 0;
            if ((res == "q") || (res == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
            if ((res == "r") || (res == "R")) {valid = 1; item = ""; goto freebayes;}
            if ((res == "y") || (res == "Y")) {valid = 1;}
            if (valid == 0) {goto freebayes;}
        }
    }
    
    
    
    
/*
    
    
    // ********** STAR  ************
    
    program = "STAR";
    command = "which STAR";
    string star = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for hla-mapper rna." << endl;
    cout << "If " + program + " is not available, [hla-mapper rna] will be disabled." <<endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {star = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }

    if (((res == "y") || (res == "Y")) && (star != "DISABLED"))
    {
    
    
        item = GetStdoutFromCommand(command);
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
     
    star_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "use this one (y), or change (c), or quit(q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "Y"))  {valid = 1; star = item;}
                if ((res == "c") || (res == "C")) {valid = 1;item = "";}
                if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }
                if (valid == 0) {goto star_start;}
            }
        
    star:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write the " << program << " path, or drag the binary bellow (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if (item == "q") {cout << "Quiting!" << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Can't find this file." << endl;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        star = item;
        item = GetStdoutFromCommand(star);
        data.clear();
        boost::split(data,item,boost::is_any_of("\n"));
        version = "";
        for (auto value : data)
        {
            vector <string> version_data;
            boost::split(version_data,value,boost::is_any_of("="));
            if (version_data[0] == "STAR version")
            {
                vector <string> subversion;
                boost::split(subversion,version_data[1],boost::is_any_of("_"));
                version = subversion[0];
            }
        }
        
        if (star_versions.find(version) != star_versions.end())
            {
                cout << "This STAR version (" << version << ") is compatible with hla-mapper" << endl;
            }
        else {
            cout << "This STAR version has not been tested with hla-mapper\nContinue anyway (y), quit (q), or return (r)? ";
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            int valid = 0;
            if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
            if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto star;}
            if ((item == "y") || (item == "Y")) {valid = 1;}
            if (valid == 0) {goto star;}
        }
    }
    
    
    
 */
    
    
    
    
   
    
    
    
    
// ********** WHATSHAP  ************
    
    program = "WhatsHap";
    command = "which whatshap";
    string whats = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for hla-mapper type." << endl;
    cout << "If " + program + " is not available, [hla-mapper type] will be disabled." <<endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {whats = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }

    if (((res == "y") || (res == "Y")) && (whats != "DISABLED"))
    {
    
        item = GetStdoutFromCommand(command);
        try
        {
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
        }
            catch (const std::exception& e)
        { item = "";}

    whats_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "use this one (y), or change (c), or quit (q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "Y")) {valid = 1; whats = item;}
                if ((res == "c") || (res == "C")) {valid = 1; item = "";}
                if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }
                if (valid == 0) {goto whats_start;}
            }
        
    whats:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write the " << program << " path, or drag the binary bellow (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if ((item == "q") || (item == "Q")) {cout << "Quiting!" << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Can't find this file." << endl;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        whats = item;
        
        item = GetStdoutFromCommand(whats + " --version");
        data.clear();
        boost::split(data,item,boost::is_any_of("\n"));
        
        if (whats_versions.find(data[0]) != whats_versions.end())
            {
                cout << "This WhatsHap version (" << data[0] << ") is compatible with hla-mapper" << endl;
            }
        else {
            cout << "This WhatsHap version has not been tested with hla-mapper\nContinue anyway (y), quit (q), or return (r)? ";
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            int valid = 0;
            if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
            if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto whats;}
            if ((item == "y") || (item == "Y")) {valid = 1;}
            if (valid == 0) {goto whats;}
        }
    }

    

    
    
    
    
// ********** BLAST  ************
    
    program = "BLAST";
    command = "which blastn";
    string blast = "";
    cout << endl;
    cout << "**** Configuring " << program << " ***" << endl;
    cout << program << " is mandatory for hla-mapper type." << endl;
    cout << "If " + program + " is not available, [hla-mapper type] will be disabled." <<endl;
    cout << "Shall we configure " + program + "? Yes (y), ignore it (i), or quit (q): ";

    cin >> res;
    res.erase(std::remove(res.begin(), res.end(), '\n'), res.end());
    if ((res == "i") || (res == "i")) {blast = "DISABLED";}
    if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }

    if (((res == "y") || (res == "Y")) && (blast != "DISABLED"))
    {
        
        item = GetStdoutFromCommand(command);
        item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    blast_start:
        if (item != "")
            {
                cout << program << " path (automatic detection): " << item << endl;
                cout << "use this one (y) or change (c), or quit (q): ";
                string res = "";
                cin >> res;
                int valid = 0;
                if ((res == "y") || (res == "y")) {valid = 1; blast = item;}
                if ((res == "c") || (res == "C")) {valid = 1; item = "";}
                if ((res == "q") || (res == "Q")) {cout << "Quiting!" << endl; return; }
                if (valid == 0) {goto blast_start;}
            }
        
    blast:
        if (item == "")
            {
                while (! fileExists(item)) {
                    cout << "Please write " << program << " path or drag the binary bellow (or q to quit): " << endl;
                    cin >> item;
                    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
                    if ((item == "q") || (item == "Q")) {cout << "Quiting!" << endl; return;}
                    if (! fileExists(item))
                        {
                            cout << "Can't find this file." << endl;
                        }
                    else {cout << program << " path: " << item << endl << endl;}
                }
            }
        
        blast = item;
        item = GetStdoutFromCommand(blast + " -version");
        data.clear();
        boost::split(data,item,boost::is_any_of("\n"));
        sub.clear();
        boost::split(sub,data[0],boost::is_any_of(" "));
        

        if (blast_versions.find(sub[1]) != blast_versions.end())
            {
                cout << "This BLAST version (" << sub[1] << ") is compatible with hla-mapper" << endl;
            }
        else {
            cout << "This BLAST version has not been tested with hla-mapper\nContinue anyway (y), quit (q), or return (r)? ";
            cin >> item;
            item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
            int valid = 0;
            if ((item == "q") || (item == "Q")) {valid = 1; cout << "Quiting!" << endl; return; }
            if ((item == "r") || (item == "R")) {valid = 1; item = ""; goto blast;}
            if ((item == "y") || (item == "Y")) {valid = 1;}
            if (valid == 0) {goto blast;}
        }
    }


    
    
    
    cout << endl << endl;
    cout << "**** Configuration list ***" << endl;
    cout << "Database: " << db << endl;
    cout << "Samtools: " << samtools << endl;
    cout << "BWA:      " << bwa << endl;
    cout << "WhatsHap: " << whats << endl;
    cout << "Freebayes: " << freebayes << endl;
//    cout << "STAR: " << star << endl;
    cout << "BLAST:    " << blast << endl;
    cout << endl;
    
write:
    cout << "Proceed writing configuration file? (y or n) ";
    cin >> item;
    item.erase(std::remove(item.begin(), item.end(), '\n'), item.end());
    int valid = 0;
    if (item == "n") {valid = 1; cout << "Quiting!" << endl; return; }
    if (item == "y") {valid = 1;}
    if (valid == 0) {goto write;}
    
    
    string homedir = getpwuid(getuid())->pw_dir;
    string configfile = homedir + "/.hla-mapper";
    
    ofstream myfile;
    myfile.open (configfile);
    myfile << "db=" + db + "\n";
    myfile << "samtools=" + samtools + "\n";
    myfile << "bwa=" + bwa + "\n";
    myfile << "whatshap=" + whats + "\n";
    myfile << "freebayes=" + freebayes + "\n";
//    myfile << "star=" + star + "\n";
    myfile << "blast=" + blast + "\n";
    myfile.close();
    return;
}
