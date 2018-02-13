/*
  @Author : Kai-Hsin Wu
  @Date   : 2017/06/30
  @Version: v1.0
  @ for C++ version < 11
*/


#include "Parser.hpp"
#include <iostream>
#include <iomanip>
using namespace std;

string RemoveBlnk(string istr){

    if(istr.empty()){ istr.clear(); return istr;}
    while(istr.at(0)==' '){
        istr.erase(istr.begin(),istr.begin()+1);

        if(istr.empty()){ istr.clear(); return istr;}
    }

    while(istr.at(istr.size()-1)==' '){
        istr.erase(istr.end()-1,istr.end());
        if(istr.empty()){ istr.clear(); return istr;}
    }

    if(istr.find(" ") != std::string::npos){ istr.clear(); return istr;}

    return istr;

}


int Parser::TryRead(const string &rawkey,const std::string &readin_VarString){

    map<string,unsigned int>::iterator it = Matcher.find(RemoveBlnk(rawkey));
    if(it == Matcher.end()) return 0;
    unsigned int idx = it->second;

    string rhs = RemoveBlnk(readin_VarString);
    if(rhs.empty()) return 0;


    //checkType:
    if(Keys[idx].pType == INT_type){
	int* tmp = (int*)Keys[idx].pVal;
	*tmp = atol(rhs.c_str());
    }else if(Keys[idx].pType == UINT_type){
	unsigned int* tmp = (unsigned int*)Keys[idx].pVal;
	*tmp = atol(rhs.c_str());
    }else if(Keys[idx].pType == LNG_type){
	long* tmp = (long*)Keys[idx].pVal;
	*tmp = atol(rhs.c_str());
    }else if(Keys[idx].pType == ULNG_type){
	unsigned long* tmp = (unsigned long*)Keys[idx].pVal;
	*tmp = atol(rhs.c_str());
    }else if(Keys[idx].pType == FLT_type){
	float* tmp = (float*)Keys[idx].pVal;
	*tmp = atof(rhs.c_str());
    }else if(Keys[idx].pType == DBL_type){
	double* tmp = (double*)Keys[idx].pVal;
	*tmp = atof(rhs.c_str());
    }else if(Keys[idx].pType == STR_type){
	string* tmp = (string*)Keys[idx].pVal;
	*tmp = rhs;
    }else if(Keys[idx].pType == LLNG_type){
	long long* tmp = (long long*)Keys[idx].pVal;
	*tmp = atoll(rhs.c_str());
    }else{
	cout << "[ERROR][Parser::TryRead] invalid Type." << endl;
	exit(1);
    }
    Keys[idx].isRead = 1;
    return 1;



}

void Parser::Parse(const string &rc_fname){
    FILE * pFile;
    char line[256];//buffer
    char *tmp;
    string lhs,rhs;
    pFile = fopen(rc_fname.c_str() , "rb");

    ///chk exist
    if(pFile != NULL){
        ///loop line
        while(1){

            ///chk EOF & read line
            if(fgets(line,sizeof(line),pFile)==NULL) break;

            if(line[0]=='#') continue;

                ///get lhs(key)
                tmp = strtok(line,"=");
                lhs = string(tmp);
                if(tmp != NULL){
                    ///get rhs(val)
                    tmp = strtok(NULL,"=");
                    if(tmp ==NULL ) continue;
                    tmp = strtok(tmp,"\n");
                    if(tmp == NULL) continue;
                    rhs = string(tmp);

		    ///tryread:
		    TryRead(lhs,rhs);
		}
	}

    	fclose(pFile);
    }else{
	cout << "[ERROR][Parser::Parse] invalid .rc file" << endl;
	exit(1);
    }
}

void Parser::PrintBinds(){
    cout << "====================\n";
    cout << "Total element bind: " << Keys.size()<<endl;

    for(int i=0;i<Keys.size();i++){
    	cout << Keys[i].key << "\t isRead " << Keys[i].isRead  << "\t Type: " << DictType[Keys[i].pType]  << endl;
    }
    //cout << "====================\n";

}

void Parser::PrintVars(){
    cout << "===============================\n";
    cout << "Total element bind: " << Keys.size() << endl;
    cout << "---------------------------------" << endl;
    string stat;

    for(int idx=0;idx<Keys.size();idx++){

        if(Keys[idx].isRead==0){

	     cout << setw(8) << left << Keys[idx].key << " [** No Parsed **] " << endl;
	     continue;
	}

	//checkType:
	if(Keys[idx].pType == INT_type){
	    int* tmp = (int*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << " = " <<  *tmp << endl;
	}else if(Keys[idx].pType == UINT_type){
	    unsigned int* tmp = (unsigned int*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType == LNG_type){
	    long* tmp = (long*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType == ULNG_type){
	    unsigned long* tmp = (unsigned long*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType == FLT_type){
	    float* tmp = (float*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType ==  DBL_type){
	    double* tmp = (double*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType == STR_type){
	    string* tmp = (string*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else if(Keys[idx].pType == LLNG_type){
	    long long* tmp = (long long*)Keys[idx].pVal;
	    cout << setw(8) << left << Keys[idx].key << "\t = " <<  *tmp << endl;
	}else{
	    cout << "[ERROR][Parser::PrintVars] invalid Type." << endl;
	    exit(1);
	}

	//cout << Keys[i].key << "\t Type: " << Keys[i].pType.name() << endl;
    }
    cout << "=================================\n";

}

void Parser::Remove(const std::string &key){
     map<std::string,unsigned int>::iterator it;
     it = Matcher.find(key);
     if(it!=Matcher.end()){
	unsigned int idx = it->second;
	Keys.erase(Keys.begin()+idx);
	Matcher.erase(it);
     }

}
void Parser::Check_All(){
    int ff =0;
    for(int i=0;i<Keys.size();i++)
    {

	if(Keys[i].isRead==0){
	    cout << "[ERROR][Parser::Check_All] < "  << Keys[i].key << " > No read." << endl;
	    ff = 1;
	}
    }
    if(ff) exit(1);

}



void Parser::Bind(const std::string &key, int &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",INT_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
 
}
void Parser::Bind(const std::string &key, float &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",FLT_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, unsigned int &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",UINT_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, unsigned long &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",ULNG_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, long &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",LNG_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, double &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",DBL_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, std::string &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",STR_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
void Parser::Bind(const std::string &key, long long &Var){
    Matcher[key] = Keys.size();
    rElem tmp = {NULL,"",LLNG_type,0 } ;
    tmp.pVal = (void*)&Var;
    tmp.key = key;
    Keys.push_back(tmp);
}
