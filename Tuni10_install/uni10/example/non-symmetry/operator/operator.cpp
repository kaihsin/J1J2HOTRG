#include "operator.h"
#include <cmath> 
#include <assert.h> 
#include <exception>

using namespace std;

const uni10_int32 CURRENT_SUPPORTED_SPIN_DIM = 5;

void spin_check(uni10_float32 spin){

  if(!(spin > 0 && floor(2 * spin) == 2 * spin)){
      std::ostringstream err;
      err<<"The given spin is not half integer.";
      throw std::runtime_error(err.str());
  }

}

uni10::Matrix<uni10_double64> matSp(uni10_float32 spin){

  spin_check(spin);

  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }

  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 1,\
      0, 0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
      0.0, 0.0, 1.0,\
      0.0, 0.0, 0.0};
    return sqrt(2) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,     0.0,     0.0, sqrt(6), 0.0,\
      0.0,     0.0,     0.0,     0.0, 2.0,\
      0.0,     0.0,     0.0,     0.0, 0.0,};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matSm(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 0,\
      1, 0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 0.0, 0.0,\
      1.0, 0.0, 0.0,\
      0.0, 1.0, 0.0};
    return sqrt(2) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     0.0,     0.0,     0.0, 0.0,\
      2.0,     0.0,     0.0,     0.0, 0.0,\
      0.0, sqrt(6),     0.0,     0.0, 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,     0.0,     0.0,     2.0, 0.0,};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::Matrix<uni10_double64> matSx(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0,   0.5,\
      0.5, 0  };
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
      1.0, 0.0, 1.0,\
      0.0, 1.0, 0.0};
    return (1.0 / sqrt(2)) * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
      2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0, sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0, sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,     2.0, 0.0,};
    return 0.5 * uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
	
  return uni10::Matrix<uni10_double64>();

}
/*
uni10::Matrix<uni10_double64> matSy(uni10_float32 spin){

  spin_check(spin);
  int dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 1,\
     -1, 0};
    return std::complex<uni10_double64>(0,1./2.0)*uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
     -1.0, 0.0, 1.0,\
      0.0,-1.0, 0.0};
    return std::complex<uni10_double64>(0,1./sqrt(2.0))* uni10::Matrix(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
     -2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,-sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0,-sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,    -2.0, 0.0,};
    return std::complex<uni10_double64>(0,1./2.0)*uni10::Matrix(dim, dim, mat_elem);
  }

  return uni10::Matrix();

}
*/
uni10::Matrix<uni10_double64> matSz(uni10_float32 spin){

  spin_check(spin);
  uni10_int32 dim = spin * 2 + 1;

  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0.5,  0,\
      0,   -0.5  };
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      1.0, 0.0,  0.0,\
      0.0, 0.0,  0.0,\
      0.0, 0.0, -1.0};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      2.0, 0.0, 0.0, 0.0, 0.0,\
      0.0, 1.0, 0.0, 0.0, 0.0,\
      0.0, 0.0, 0.0, 0.0, 0.0,\
      0.0, 0.0, 0.0,-1.0, 0.0,\
      0.0, 0.0, 0.0, 0.0,-2.0,};
    return uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
  }

  return uni10::Matrix<uni10_double64>();

}

uni10::UniTensor<uni10_double64> periodicSx(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> Sx = matSx(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sx.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sx.col()));
  uni10::UniTensor<uni10_double64> perSx(bondI);
  perSx.putBlock(Sx);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSx = otimes(perSx, Id);

  uni10::UniTensor<uni10_double64> tmp = perSx;
  vector<uni10_int32> labels = perSx.label();
  vector<uni10_int32> per_labels(labels.size());
  uni10_int32 bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(uni10_int32 l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSx += permute(tmp, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }

  return perSx;

}

/*
uni10::UniTensor periodicSy(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix Sy = matSy(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sy.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sy.col()));
  uni10::UniTensor perSy(bondI);
  perSy.putBlock(Sy);
  uni10::UniTensor Id(uni10::CTYPE, bondI);
  Id.identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSy = otimes(perSy, Id);

  uni10::UniTensor tmp = perSy;
  vector<int> labels = perSy.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSy += tmp.permute(uni10::CTYPE, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }

  return perSy;

}
*/

uni10::UniTensor<uni10_double64> periodicSz(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> Sz = matSz(spin);  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, Sz.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, Sz.col()));
  uni10::UniTensor<uni10_double64> perSz(bondI);
  perSz.putBlock(Sz);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    perSz = otimes(perSz, Id);
	
  uni10::UniTensor<uni10_double64> tmp = perSz;
  vector<int> labels = perSz.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    perSz += permute(tmp, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }

  return perSz;

}

uni10::UniTensor<uni10_double64> periodicsqSz(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> sqSz = uni10::dot(matSz(spin), matSz(spin));  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSz.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSz.col()));
  uni10::UniTensor<uni10_double64> persqSz(bondI);
  persqSz.putBlock(sqSz);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSz = otimes(persqSz, Id);

  uni10::UniTensor<uni10_double64> tmp = persqSz;
  vector<int> labels = persqSz.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();

  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSz += permute(tmp, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }

  return persqSz;

}

uni10::UniTensor<uni10_double64> periodicsqSx(uni10_uint64 siteNum, uni10_float32 spin){

  uni10::Matrix<uni10_double64> sqSx = uni10::dot(matSx(spin), matSx(spin));  
  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSx.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSx.col()));
  uni10::UniTensor<uni10_double64> persqSx(bondI);
  persqSx.putBlock(sqSx);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.identity();

  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSx = otimes(persqSx, Id);

  uni10::UniTensor<uni10_double64> tmp = persqSx;
  vector<uni10_int32> labels = persqSx.label();
  vector<uni10_int32> per_labels(labels.size());
  uni10_int32 bondNum = per_labels.size();

  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSx += permute(tmp, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }

  return persqSx;

}

uni10::UniTensor<uni10_double64> periodicsqSy(uni10_uint64 siteNum, uni10_float32 spin){

  spin_check(spin);
  int dim = spin * 2 + 1;
  uni10::Matrix<uni10_double64> matSy_Rpart;
  uni10::Matrix<uni10_double64> sqSy;
  if(dim > CURRENT_SUPPORTED_SPIN_DIM){
    std::ostringstream err;
    err<<"Spin = "<<spin <<" is not yet supported.";
    throw std::runtime_error(err.str());
  }
  if(dim == 2){
    uni10_double64 mat_elem[] = {\
      0, 1,\
     -1, 0};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/4.0) * uni10::dot(matSy_Rpart, matSy_Rpart);
  }
  else if(dim == 3){
    uni10_double64 mat_elem[] = {\
      0.0, 1.0, 0.0,\
     -1.0, 0.0, 1.0,\
      0.0,-1.0, 0.0};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/2.0) * uni10::dot(matSy_Rpart, matSy_Rpart);
  }
  else if(dim == 5){
    uni10_double64 mat_elem[] = {\
      0.0,     2.0,     0.0,     0.0, 0.0,\
     -2.0,     0.0, sqrt(6),     0.0, 0.0,\
      0.0,-sqrt(6),     0.0, sqrt(6), 0.0,\
      0.0,     0.0,-sqrt(6),     0.0, 2.0,\
      0.0,     0.0,     0.0,    -2.0, 0.0,};
    matSy_Rpart = uni10::Matrix<uni10_double64>(dim, dim, mat_elem);
    sqSy = (-1.0/4.0) * uni10::dot(matSy_Rpart, matSy_Rpart);
  }

  std::vector<uni10::Bond> bondI;
  bondI.push_back(uni10::Bond(uni10::BD_IN, sqSy.row())); 
  bondI.push_back(uni10::Bond(uni10::BD_OUT, sqSy.col()));
  uni10::UniTensor<uni10_double64> persqSy(bondI);
  persqSy.putBlock(sqSy);
  uni10::UniTensor<uni10_double64> Id(bondI);
  Id.identity();
  for(uni10_uint64 i = 0; i < siteNum - 1; i++)
    persqSy = otimes(persqSy, Id);
  uni10::UniTensor<uni10_double64> tmp = persqSy;
  vector<int> labels = persqSy.label();
  vector<int> per_labels(labels.size());
  int bondNum = per_labels.size();
  for(uni10_uint64 h = 0; h < siteNum - 1; h++){
    for(int l = 0; l < bondNum/2; l++){
      per_labels[l] = labels[(l + 1) % (bondNum/2)];
      per_labels[(bondNum/2) + l] = labels[bondNum/2 + ((l + 1) % (bondNum/2))];
    }
    persqSy += permute(tmp, per_labels, bondNum / 2);
    tmp.setLabel(labels);
  }
  return persqSy;
}


