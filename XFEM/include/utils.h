#ifndef UTILS_XFEM_H
#define UTILS_XFEM_H

#include <iostream>
#include <boost/lexical_cast.hpp>
#include <ctime>

#ifndef __func__
#define __func__  __FUNCTION__
#endif

#define DEBUG ::std::cerr << __FILE__ << ":" << __LINE__ << " : function " << __func__ << ::std::endl
#define WARNING(a) ::std::cerr << "WARNING IN " << __FILE__ << ":" << __LINE__ << " : function " << __func__ << " : " << a << ::std::endl
#define SHOW(a) ::std::cout << #a << " :  " << (a) << ::std::endl
#define NOT_IMPLEMENTED_METHOD ::std::cerr << this->getClassName() << "::" << __func__  << " : not implemented !" << ::std::endl; DEBUG
#define NOT_IMPLEMENTED_FUNCTION ::std::cerr << __func__  << " : not implemented !" << ::std::endl ; DEBUG
#define LOG_FUNCALL serr << getClassName() << "::" << __func__  << " called !" << sendl
#define VARIABLENAME(obj) #obj

#endif // UTILS_H
