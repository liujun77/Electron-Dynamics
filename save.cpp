// *****************************************************************************
// *** save.cpp                                                              ***
// *** Jun Liu 2014-08-18                                                    ***
// *** Saving data to file                                                   ***
// *****************************************************************************

#include <fstream>
#include <string>
#include <cstdarg>

using namespace std;

namespace jun{

int save(string filename, int col, int row, ...){
    ofstream file;
    file.open(filename.c_str());
    int i,j;
    va_list arg_p;
    for(i=0;i<row;++i){
        va_start(arg_p,row);
        for(j=0;j<col;++j)
            file<<*(va_arg(arg_p,double*)+i)<<"\t";
        file<<endl;
        va_end(arg_p);
    }
    file.close();
    return 0;
}

}
