// Copyright (c) 2013, Adam Harrison*
// http://www.ualberta.ca/~apharris/
// All rights reserved.

// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

// -Redistributions of source code must retain the above copyright notice, the footnote below, this list of conditions and the following disclaimer.
// -Redistributions in binary form must reproduce the above copyright notice, the footnote below, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
// -Neither the name of the University of Alberta nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// *This work originated as part of a Ph.D. project under the supervision of Dr. Dileepan Joseph at the Electronic Imaging Lab, University of Alberta.


#ifndef LATTICEEXCEPTION_H
#define LATTICEEXCEPTION_H

#include <stdexcept>
#include <string>

namespace LibMIA
{

class LatticeException: public std::runtime_error
{
    public:
        LatticeException(): std::runtime_error("LatticeException") {}
        LatticeException(const std::string& __arg): std::runtime_error(__arg) {}
    protected:
    private:
};


class LatticeParameterException: public std::runtime_error
{
    public:
        LatticeParameterException(): std::runtime_error("Invalid Parameters to Lattice Operation") {}
        LatticeParameterException(const std::string& __arg): std::runtime_error(__arg) {}
    protected:
    private:
};

class RankDeficientException: public std::runtime_error
{
    public:
        RankDeficientException(): std::runtime_error("Invalid Parameters to Lattice Operation") {}
        RankDeficientException(const std::string& __arg): std::runtime_error(__arg) {}
    protected:
    private:
};

}
#endif // LATTICEEXCEPTION_H
