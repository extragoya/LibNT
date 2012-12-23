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
