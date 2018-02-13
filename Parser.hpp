/*
  @Author : Kai-Hsin Wu
  @Date   : 2017/06/30
  @Version: v1.0
  @ for C++ version < 11
*/


#ifndef PARSER_HPP_INCLUDED
#define PARSER_HPP_INCLUDED
#include <map>
#include <string>
//#include <typeinfo>
//#include <typeindex>
#include <fstream>
#include <vector>
#include <cstring>
//#include <unordered_map>

class Parser{
private:
   enum {
	INT_type,
	FLT_type,
	UINT_type,
	ULNG_type,
	LNG_type,
	DBL_type,
	STR_type,
	LLNG_type
   };
   std::map<unsigned int,std::string> DictType;
public :
   Parser(){
	Init();
   }
   void Init(){
        DictType[INT_type] = "int";
        DictType[FLT_type] = "float";
        DictType[UINT_type] = "unsigned int";
        DictType[ULNG_type] = "unsigned long";
        DictType[LNG_type] = "long";
        DictType[DBL_type] = "double";
        DictType[STR_type] = "string";
	DictType[LLNG_type] = "long long";
   }

    struct rElem{
        void* pVal;
        std::string key;
        int pType;
        bool isRead;
    };

    std::map<std::string,unsigned int > Matcher;
    std::vector< rElem > Keys;
    
    //template<class T>
    void Bind(const std::string &key, int &Var);
    void Bind(const std::string &key, float &Var);
    void Bind(const std::string &key, unsigned int &Var);
    void Bind(const std::string &key, unsigned long &Var);
    void Bind(const std::string &key, long &Var);
    void Bind(const std::string &key, double &Var);
    void Bind(const std::string &key, std::string &Var);
    void Bind(const std::string &key, long long &Var);




    
    ///Overload for different type:
    int TryRead(const std::string &rawkey,const std::string &readin_VarString);
    void Parse(const std::string &rc_fname);
    void PrintBinds();
    void PrintVars();
    void Remove(const std::string &keys);
    void Check_All();


};



#endif // SSE_HPP_INCLUDED
