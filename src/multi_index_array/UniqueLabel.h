#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/slot/slot.hpp>

#if !defined(UNIQUE_LABEL)
#define UNIQUE_LABEL
#define BOOST_PP_VALUE 1
#include BOOST_PP_ASSIGN_SLOT(1)
#undef BOOST_PP_VALUE
#else
#define BOOST_PP_VALUE BOOST_PP_INC(BOOST_PP_SLOT(1))
#include BOOST_PP_ASSIGN_SLOT(1)
#undef BOOST_PP_VALUE
#endif
