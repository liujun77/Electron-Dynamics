// *****************************************************************************
// *** readin.cpp                                                            ***
// *** Jun Liu 2014-02-16                                                    ***
// *** Reading the input file                                                ***
// *****************************************************************************

#include <sstream>
#include <fstream>

using namespace std;

namespace jun{

int is_check_in(string str, string check){
    //check if the word check is in str
    istringstream input_str(str);
    string temp;
    while(!input_str.eof()){
        input_str>>temp;
        if(temp[0]=='#')
            return 0;
        if(temp==check)
            return 1;
    }
    return 0;
}

int count_lines(string filename){
    int i=0;
    ifstream infile(filename.c_str());
    string line;
    while(getline(infile,line)){
        ++i;
    }
    infile.close();
    return i;
}

}
