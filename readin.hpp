// *****************************************************************************
// *** readin.h                                                              ***
// *** Jun Liu 2014-02-16                                                    ***
// *** Reading the input file                                                ***
// *****************************************************************************

#ifndef READIN_H
#define READIN_H

#include <fstream>
#include <sstream>
#include <cctype>
#include <iostream>
#include <cstdlib>

namespace jun{

template<typename T>

int val_in_str(T &val, std::string str){
    //put the first number in string str to val, and return 1,
    //if no number return 0
    std::string temp;
    T re;
    int flag=0;
    std::istringstream input_str(str),save;
    while(!input_str.eof()){
        input_str>>temp;
        if(isdigit(temp[0])||temp[0]=='-'||temp[0]=='+'){
            save.str(temp);
            save>>re;
            flag=1;
        }
    }
    if(flag)
        val=re;
    return flag;
}

int is_check_in(std::string str, std::string check);
    //check if the word check is in check

template<typename T>

inline int read(T &val, std::string file_name, std::string parameter){
    //read parameter from file_name, put the value to val, and return 1,
    //if failed return 0
    T re;
    std::string str;
    int flag=0;
    std::ifstream file_in(file_name.c_str());
    if(!file_in){
        std::cerr<<"Can not open input file!"<<std::endl;
        std::exit(0);
    }
    while(!file_in.eof()){
        getline(file_in,str);
        if(is_check_in(str,parameter)){
            val_in_str(re,str);
            flag=1;
        }
    }
    if(flag)
        val=re;
    file_in.close();
    return flag;
}

int count_lines(std::string infile);

}

#endif
