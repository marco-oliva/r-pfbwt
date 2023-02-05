//
//  version.hpp.in
//
//  Copyright 2023s Marco Oliva. All rights reserved.
//

#ifndef version_hpp
#define version_hpp

#include <string>
#include <iostream>

namespace rpfbwt
{

struct Version
{
    static std::string RPFBWT_GIT_BRANCH;
    static std::string RPFBWT_GIT_COMMIT_HASH;
    static int RPFBWT_MAJOR;
    static int RPFBWT_MINOR;
    static int RPFBWT_PATCH;
    
    static void print();
};

}

#endif //version_hpp