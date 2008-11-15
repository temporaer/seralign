/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

#define DIMACS_PRINT 0

#include "sdpa_io.h"
#include <vector>
#include <algorithm>

namespace sdpa {

// 2008/02/27  kazuhide nakata 
#if 0  // not use
void IO::read(FILE* fpData, int m,
	      int SDP_nBlock,int* SDP_blockStruct,
	      int SOCP_nBlock,int* SOCP_blockStruct,
	      int LP_nBlock, 
	      int nBlock, int* blockStruct,
	      BlockType* blockType, int* blockNumber,
	      InputData& inputData, bool isDataSparse)
{
  inputData.initialize_bVec(m);
  read(fpData,inputData.b);
  long position = ftell(fpData);
  // C,A must be accessed "twice".

  // count numbers of elements of C and A
  int* SDP_CNonZeroCount;
  NewArray(SDP_CNonZeroCount,int,SDP_nBlock);
  int* SDP_ANonZeroCount;
  NewArray(SDP_ANonZeroCount,int,m*SDP_nBlock);
  
  // count numbers of elements of C and A
  int* SOCP_CNonZeroCount;
  NewArray(SOCP_CNonZeroCount,int,SOCP_nBlock);
  int* SOCP_ANonZeroCount;
  NewArray(SOCP_ANonZeroCount,int,m*SOCP_nBlock);
  
  // count numbers of elements of C and A
  bool* LP_CNonZeroCount;
  NewArray(LP_CNonZeroCount,bool,LP_nBlock);
  bool* LP_ANonZeroCount;
  NewArray(LP_ANonZeroCount,bool,m*LP_nBlock);
  //   initialize C and A
  read(fpData,m,
       SDP_nBlock, SDP_blockStruct, SDP_CNonZeroCount, SDP_ANonZeroCount,
       SOCP_nBlock, SOCP_blockStruct,
       SOCP_CNonZeroCount, SOCP_ANonZeroCount,
       LP_nBlock, LP_CNonZeroCount, LP_ANonZeroCount,
       nBlock, blockStruct, blockType, blockNumber,
       isDataSparse);
  //   rMessage(" C and A count over");
  inputData.initialize_CMat(SDP_nBlock, SDP_blockStruct,
			    SDP_CNonZeroCount,
			    SOCP_nBlock,  SOCP_blockStruct,
			    SOCP_CNonZeroCount,
			    LP_nBlock, LP_CNonZeroCount);
  inputData.initialize_AMat(m,SDP_nBlock, SDP_blockStruct,
			    SDP_ANonZeroCount,
			    SOCP_nBlock,  SOCP_blockStruct,
			    SOCP_ANonZeroCount,
			    LP_nBlock, LP_ANonZeroCount);
  DeleteArray(SDP_CNonZeroCount);
  DeleteArray(SDP_ANonZeroCount);
  DeleteArray(SOCP_CNonZeroCount);
  DeleteArray(SOCP_ANonZeroCount);
  DeleteArray(SOCP_ANonZeroCount);
  DeleteArray(LP_CNonZeroCount);
  DeleteArray(LP_ANonZeroCount);
  
  //   rMessage(" C and A initialize over");
  read(fpData, inputData, m, 
       SDP_nBlock, SDP_blockStruct, 
       SOCP_nBlock, SOCP_blockStruct, 
       LP_nBlock, 
       nBlock, blockStruct, blockType, blockNumber,
       position, isDataSparse);
  //   rMessage(" C and A have been read");
}
#endif

void IO::read(FILE* fpData, FILE* fpout, int& m, char* str)
{
  while (true) {
    volatile int dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (str[0]=='*' || str[0]=='"') {
      fprintf(fpout,"%s",str);
    } else {
      sscanf(str,"%d",&m);
      break;
    }
  }
}

void IO::read(FILE* fpData, int& nBlock)
{
  fscanf(fpData,"%d",&nBlock);
}

void IO::read(FILE* fpData, int nBlock, int* blockStruct)
{
  for (int l=0; l<nBlock; ++l) {
    fscanf(fpData,"%*[^0-9+-]%d",&blockStruct[l]);
  }
}

void IO::read(FILE* fpData, Vector& b)
{
  for (int k=0; k<b.nDim; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&b.ele[k]);
  }
}

void IO::read(FILE* fpData, DenseLinearSpace& xMat,
	      Vector& yVec, DenseLinearSpace& zMat,
	      int SDP_nBlock, int* SDP_blockStruct,
	      int SOCP_nBlock, int* SOCP_blockStruct,
	      int LP_nBlock,
	      int nBlock, int* blockStruct,
	      BlockType* blockType, int* blockNumber,
	      bool inputSparse)
{
  // yVec is opposite sign
  for (int k=0; k<yVec.nDim; ++k) {
    double tmp;
    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
    yVec.ele[k] = -tmp;
    // rMessage("yVec.ele[" << k << "] = " << tmp);
  }
  
  if (inputSparse) {
    // sparse case , zMat , xMat in this order
    int i,j,l,target;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&target)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
      #if 0
      rMessage("target = " << target
	       << ": l " << l
	       << ": i " << i
	       << ": j " << j
	       << ": value " <<value);
      #endif

      if (blockType[l-1] == btSDP) {
	int l2 = blockNumber[l-1];
	if (target==1) {
	  zMat.setElement_SDP(l2,i-1,j-1,value);
	} else {
	  xMat.setElement_SDP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btSOCP) {
	rError("io:: current version does not support SOCP");
	int l2 = blockNumber[l-1];
	if (target==1) {
	  zMat.setElement_SOCP(l2,i-1,j-1,value);
	} else {
	  xMat.setElement_SOCP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btLP) {
	if (i != j){
	  rError("io:: LP part  3rd element != 4th element\n"
		 "column should be the same as row in LP part.");
	}
	#if 0
	rMessage("l = " << l
		 << ": blockNumber[l-1] = " << blockNumber[l-1]
		 << ": index = " << blockNumber[l-1]+i-1
		 << ": i = " << i);
	#endif
	if (target==1) {
	  zMat.setElement_LP(blockNumber[l-1]+i-1,value);
	} else {
	  xMat.setElement_LP(blockNumber[l-1]+i-1,value);
	}
      }
    } // end of 'while (true)'
  } else {
    // dense case , zMat , xMat in this order
    // for SDP
    for (int l=0; l<nBlock; ++l) {
      if (blockType[l] == btSDP) {
	int l2 = blockNumber[l];
	int size = blockStruct[l];
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      zMat.setElement_SDP(l2,i,j,tmp);
	    }
	  }
	}
      }
      else if (blockType[l] == btSOCP) {
	rError("io:: current version does not support SOCP");
      }
      else if (blockType[l] == btLP) {
	int size  = blockStruct[l];
	int index = blockNumber[l];
	for (int j=0; j<size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    zMat.setElement_LP(index,tmp);
	  }
	  index++;
	}
      }
    }
    
    for (int l=0; l<nBlock; ++l) {
      if (blockType[l] == btSDP) {
	int l2 = blockNumber[l];
	int size = blockStruct[l];
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      xMat.setElement_SDP(l2,i,j,tmp);
	    }
	  }
	}
      }
      else if (blockType[l] == btSOCP) {
	rError("io:: current version does not support SOCP");
      }
      else if (blockType[l] == btLP) {
	int size  = blockStruct[l];
	int index = blockNumber[l];
	for (int j=0; j<size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    xMat.setElement_LP(index,tmp);
	  }
	  index++;
	}
      }
    }
  } // end of 'if (inputSparse)'
}

// 2008/02/27 kazuhide nakata   
// not use
void IO::read(FILE* fpData, int m, 
	      int SDP_nBlock, int* SDP_blockStruct,
	      int* SDP_CNonZeroCount, int* SDP_ANonZeroCount,
	      int SOCP_nBlock, int* SOCP_blockStruct,
	      int* SOCP_CNonZeroCount, int* SOCP_ANonZeroCount,
	      int LP_nBlock,
	      bool* LP_CNonZeroCount, bool* LP_ANonZeroCount,
	      int nBlock, int* blockStruct,
	      BlockType* blockType, int* blockNumber,
	      bool isDataSparse)
{
  // only count the numbers of C,A[k]
  for (int l=0; l<SDP_nBlock; ++l) {
    SDP_CNonZeroCount[l] = 0;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<SDP_nBlock; ++l) {
      SDP_ANonZeroCount[k*SDP_nBlock + l] = 0;
    }
  }
  for (int l=0; l<SOCP_nBlock; ++l) {
    SOCP_CNonZeroCount[l] = 0;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<SOCP_nBlock; ++l) {
      SOCP_ANonZeroCount[k*SOCP_nBlock + l] = 0;
    }
  }
  for (int l=0; l<LP_nBlock; ++l) {
    LP_CNonZeroCount[l] = false;
  }
  for (int k=0; k<m; ++k) {
    for (int l=0; l<LP_nBlock; ++l) {
      LP_ANonZeroCount[k*LP_nBlock + l] = false;
    }
  }
  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }

      if (blockType[l-1] == btSDP ) {	// SDP part
	int l2 = blockNumber[l-1];
	if (k==0) {
	  SDP_CNonZeroCount[l2]++;
	} else {
	  SDP_ANonZeroCount[(k-1)*SDP_nBlock+l2]++;
	}
      } else if (blockType[l-1] == btSOCP) {	// SOCP part
	rError("io:: current version does not support SOCP");
	int l2 = blockNumber[l-1];;
	if (k==0) {
	  SOCP_CNonZeroCount[l2]++;
	} else {
	  SOCP_ANonZeroCount[(k-1)*SOCP_nBlock+l2]++;
	}
      } else if (blockType[l-1] == btLP) { // LP part
	int l2 =blockNumber[l-1];
	if (k==0) {
	  LP_CNonZeroCount[l2+i-1] = true;
	} else {
	  LP_ANonZeroCount[(k-1)*LP_nBlock+l2+i-1] = true;
	}
      } else {
	rError("io::read not valid blockType");
      }
    }// end of 'while (true)'

  } else { // isDataSparse == false

    // k==0
    for (int l2=0; l2<nBlock; ++l2){
      if (blockType[l2] == btSDP) { // SDP part
	int l = blockNumber[l2];
	int size = SDP_blockStruct[l];
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      SDP_CNonZeroCount[l]++;
	    }
	  }
	}
      } else if (blockType[l2] == btSOCP) { // SOCP part
	rError("io:: current version does not support SOCP");
      } else if (blockType[l2] == btLP) { // LP part
	for (int j=0; j<blockStruct[l2]; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    LP_CNonZeroCount[blockNumber[l2]+j]++;
	  }
	}
      } else {
	rError("io::read not valid blockType");
      }
    }

    for (int k=0; k<m; ++k) {
      // k>0
      for (int l2=0; l2<nBlock; ++l2){
	if (blockType[l2] == btSDP) { // SDP part
	  int l = blockNumber[l2];
	  int size = SDP_blockStruct[l];
	  for (int i=0; i<size; ++i) {
	    for (int j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		SDP_ANonZeroCount[k*SDP_nBlock+l]++;
	      }
	    }
	  }
	} else if (blockType[l2] == btSOCP) { // SOCP part
	  rError("io:: current version does not support SOCP");
	} else if (blockType[l2] == btLP) { // LP part
	  for (int j=0; j<blockStruct[l2]; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (tmp!=0.0) {
	      LP_ANonZeroCount[k*LP_nBlock+blockNumber[l2]+j] = true;
	    }
	  }
	} else {
	  rError("io::read not valid blockType");
	}
      }
    }

  } // end of 'if (isDataSparse)'

}


// 2008/02/27 kazuhide nakata   
// not use
void IO::read(FILE* fpData,  InputData& inputData, int m, 
	      int SDP_nBlock, int* SDP_blockStruct, 
	      int SOCP_nBlock, int* SOCP_blockStruct, 
	      int LP_nBlock, 
	      int nBlock, int* blockStruct,
	      BlockType* blockType, int* blockNumber,
	      long position, bool isDataSparse)

{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     

      if (blockType[l-1] == btSDP) {	// SDP part
	int l2 = blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SDP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btSOCP) {	// SOCP part
	rError("io:: current version does not support SOCP");
	int l2 = blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btLP) { // LP part
	if (i != j){
	  rError("io:: LP part  3rd elemtn != 4th elemnt");
	}
	if (k==0) {
	  inputData.C.setElement_LP(blockNumber[l-1]+i-1,-value);
	} else {
	  inputData.A[k-1].setElement_LP(blockNumber[l-1]+i-1,value);
	}
      } else {
	rError("io::read not valid blockType");
      }
    } 
  } else {  // dense

    // k==0
    for (int l2=0; l2<nBlock; ++l2){
      if (blockType[l2] == btSDP) { // SDP part
	int l = blockNumber[l2];
	int size = SDP_blockStruct[l];
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      inputData.C.setElement_SDP(l,i,j,-tmp);
	    }
	  }
	}
      } else if (blockType[l2] == btSOCP) { // SOCP part
	rError("io:: current version does not support SOCP");
      } else if (blockType[l2] == btLP) { // LP part
	for (int j=0; j<blockStruct[l2]; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    inputData.C.setElement_LP(blockNumber[l2]+j,-tmp);
	  }
	}
      } else {
	rError("io::read not valid blockType");
      }
    }

    // k > 0
    for (int k=0; k<m; ++k) {
	  
      for (int l2=0; l2<nBlock; ++l2){
	if (blockType[l2] == btSDP) { // SDP part
	  int l = blockNumber[l2];
	  int size = SDP_blockStruct[l];
	  for (int i=0; i<size; ++i) {
	    for (int j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		inputData.A[k].setElement_SDP(l,i,j,tmp);
	      }
	    }
	  }
	} else if (blockType[l2] == btSOCP) { // SOCP part
	  rError("io:: current version does not support SOCP");
	} else if (blockType[l2] == btLP) { // LP part
	  for (int j=0; j<blockStruct[l2]; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (tmp!=0.0) {
	      inputData.A[k].setElement_LP(blockNumber[l2]+j,tmp);
	    }
	  }
	} else {
	  rError("io::read not valid blockType");
	}
      }

    } // for k

  } // end of 'if (isDataSparse)'

}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
#if 1
void IO::read(FILE* fpData, int m,
	      int SDP_nBlock,
              int* SDP_blockStruct,
	      int SOCP_nBlock,
              int* SOCP_blockStruct,
	      int LP_nBlock, 
              int nBlock, 
              int* blockStruct,
              BlockType* blockType, 
              int* blockNumber,
	      InputData& inputData, bool isDataSparse)
{
  inputData.initialize_bVec(m);
  read(fpData,inputData.b);
  long position = ftell(fpData);

  // C,A must be accessed "double".

  //   initialize block struct of C and A
  setBlockStruct(fpData, inputData, m,
                 SDP_nBlock, 
                 SDP_blockStruct,
                 SOCP_nBlock, 
                 SOCP_blockStruct,
                 LP_nBlock,
                 nBlock, blockStruct, blockType, blockNumber,
                 position, isDataSparse);
  //   rMessage(" C and A initialize over");
    
  setElement(fpData, inputData, m, 
             SDP_nBlock, 
             SDP_blockStruct, 
             SOCP_nBlock, 
             SOCP_blockStruct, 
             LP_nBlock, 
             nBlock, blockStruct, blockType, blockNumber,
             position, isDataSparse);
  //   rMessage(" C and A have been read");
}
#endif

// 2008/02/27 kazuhide nakata   
// without LP_ANonZeroCount
void IO::setBlockStruct(FILE* fpData, InputData& inputData, int m,
                        int SDP_nBlock, 
                        int* SDP_blockStruct,
                        int SOCP_nBlock, 
                        int* SOCP_blockStruct,
                        int LP_nBlock,
                        int nBlock, int* blockStruct, 
                        BlockType* blockType, int* blockNumber,
                        long position, bool isDataSparse)
{
  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  vector<int>* SDP_index;
  NewArray(SDP_index,vector<int>,m+1);
  vector<int>* SOCP_index;
  NewArray(SOCP_index,vector<int>,m+1);
  vector<int>* LP_index;
  NewArray(LP_index,vector<int>,m+1);

  // for SDP
  int SDP_sp_nBlock;
  int* SDP_sp_index;
  int* SDP_sp_blockStruct;
  int* SDP_sp_NonZeroNumber;
  NewArray(SDP_sp_index,int,SDP_nBlock);
  NewArray(SDP_sp_blockStruct,int,SDP_nBlock);
  NewArray(SDP_sp_NonZeroNumber,int,SDP_nBlock);
  // for SOCP
  int SOCP_sp_nBlock;
  int* SOCP_sp_blockStruct;
  int* SOCP_sp_index;
  int* SOCP_sp_NonZeroNumber;
  // for LP
  int LP_sp_nBlock;
  int* LP_sp_index;
  NewArray(LP_sp_index,int,LP_nBlock);

  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
      
      if (blockType[l-1] == btSDP) {	// SDP part
        int l2 = blockNumber[l-1];
        SDP_index[k].push_back(l2);
      } else if (blockType[l-1] == btSOCP) {	// SOCP part
        rError("io:: current version does not support SOCP");
        int l2 = blockNumber[l-1];;
        SOCP_index[k].push_back(l2);
      } else if (blockType[l-1] == btLP) { // LP part
        if (i!=j){
          printf("invalid data file k:%d, l:%d, i:%d, j:%d, value:%lf\n"
                 ,k,l,i,j,value);
          rError("IO::initializeLinearSpace");
        }
        int l2 =blockNumber[l-1];
        LP_index[k].push_back(l2+i-1);
      } else {
        rError("io::read not valid blockType");
      }
    }// end of 'while (true)'
    
  } else { // isDataSparse == false
    
    // k==0
    for (int l2=0; l2<nBlock; ++l2){
      if (blockType[l2] == btSDP) { // SDP part
        int l = blockNumber[l2];
        int size = SDP_blockStruct[l];
        for (int i=0; i<size; ++i) {
          for (int j=0; j<size; ++j) {
            double tmp;
            fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
            if (i<=j && tmp!=0.0) {
              SDP_index[0].push_back(l);
            }
          }
        }
      } else if (blockType[l2] == btSOCP) { // SOCP part
        rError("io:: current version does not support SOCP");
      } else if (blockType[l2] == btLP) { // LP part
        for (int j=0; j<blockStruct[l2]; ++j) {
          double tmp;
          fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
          if (tmp!=0.0) {
              LP_index[0].push_back(blockNumber[l2]+j);
          }
        }
      } else {
        rError("io::read not valid blockType");
      }
    }
    
    for (int k=0; k<m; ++k) {
      // k>0
      for (int l2=0; l2<nBlock; ++l2){
        if (blockType[l2] == btSDP) { // SDP part
          int l = blockNumber[l2];
          int size = SDP_blockStruct[l];
          for (int i=0; i<size; ++i) {
            for (int j=0; j<size; ++j) {
              double tmp;
              fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
              if (i<=j && tmp!=0.0) {
                SDP_index[k+1].push_back(l);
              }
            }
          }
        } else if (blockType[l2] == btSOCP) { // SOCP part
          rError("io:: current version does not support SOCP");
        } else if (blockType[l2] == btLP) { // LP part
          for (int j=0; j<blockStruct[l2]; ++j) {
            double tmp;
            fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
            if (tmp!=0.0) {
              LP_index[k+1].push_back(blockNumber[l2]+j);
            }
          }
        } else {
          rError("io::read not valid blockType");
        }
      }
    }
    
  } // end of 'if (isDataSparse)'

  NewArray(inputData.A,SparseLinearSpace,m);
  for (int k=0 ; k<m+1; k++){
    sort(SDP_index[k].begin(),SDP_index[k].end());
    SDP_sp_nBlock = 0;
    int previous_index = -1;
    int index;
    for (unsigned int i=0; i<SDP_index[k].size(); i++){
      index = SDP_index[k][i];
      if (previous_index != index){
        SDP_sp_index[SDP_sp_nBlock] = index;
        SDP_sp_blockStruct[SDP_sp_nBlock] = SDP_blockStruct[index];
        SDP_sp_NonZeroNumber[SDP_sp_nBlock] = 1;
        previous_index = index;
        SDP_sp_nBlock++;
      } else {
        SDP_sp_NonZeroNumber[SDP_sp_nBlock-1]++;
      }
    }

    // dummy initialization to surpress compiler warning
    SOCP_sp_nBlock        = 0;
    SOCP_sp_blockStruct   = NULL;
    SOCP_sp_index         = NULL;
    SOCP_sp_NonZeroNumber = NULL;
    
    sort(LP_index[k].begin(),LP_index[k].end());
    LP_sp_nBlock=0;
    previous_index = -1;
    for (unsigned int i=0; i<LP_index[k].size(); i++){
      index = LP_index[k][i];
      if (previous_index != index){
        LP_sp_index[LP_sp_nBlock] = index;
        previous_index = index;
        LP_sp_nBlock++;
      }
    }

    if (k==0){
      inputData.C.initialize(SDP_sp_nBlock, 
                             SDP_sp_index,
                             SDP_sp_blockStruct, 
                             SDP_sp_NonZeroNumber,
                             SOCP_sp_nBlock, 
                             SOCP_sp_blockStruct, 
                             SOCP_sp_index,
                             SOCP_sp_NonZeroNumber,
                             LP_sp_nBlock, 
                             LP_sp_index);
    } else {
      inputData.A[k-1].initialize(SDP_sp_nBlock, 
                                  SDP_sp_index,
                                  SDP_sp_blockStruct, 
                                  SDP_sp_NonZeroNumber,
                                  SOCP_sp_nBlock, 
                                  SOCP_sp_blockStruct, 
                                  SOCP_sp_index,
                                  SOCP_sp_NonZeroNumber,
                                  LP_sp_nBlock, 
                                  LP_sp_index);
    }
  }

  DeleteArray(SDP_index);
  DeleteArray(SOCP_index);
  DeleteArray(LP_index);

  DeleteArray(SDP_sp_index);
  DeleteArray(SDP_sp_blockStruct);
  DeleteArray(SDP_sp_NonZeroNumber);
  DeleteArray(SDP_sp_NonZeroNumber);
#if 0
  DeleteArray(SOCP_sp_index);
  DeleteArray(SOCP_sp_blockStruct);
  DeleteArray(SOCP_sp_NonZeroNumber);
#endif
  DeleteArray(LP_sp_index);
}


// 2008/02/27 kazuhide nakata   
// without LP_ANonZeroCount
void IO::setElement(FILE* fpData, InputData& inputData, int m, 
                    int SDP_nBlock, int* SDP_blockStruct, 
                    int SOCP_nBlock, int* SOCP_blockStruct, 
                    int LP_nBlock, 
                    int nBlock, int* blockStruct, 
                    BlockType* blockType, int* blockNumber,
                    long position, bool isDataSparse)
{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  if (isDataSparse) {
    int i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     

      if (blockType[l-1] == btSDP) {	// SDP part
	int l2 = blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SDP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btSOCP) {	// SOCP part
	rError("io:: current version does not support SOCP");
	int l2 = blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
	}
      } else if (blockType[l-1] == btLP) { // LP part
	if (i != j){
	  rError("io:: LP part  3rd element != 4th element\n"
		 "column should be same as row in LP part.");
	}
	if (k==0) {
	  inputData.C.setElement_LP(blockNumber[l-1]+i-1,-value);
	} else {
	  inputData.A[k-1].setElement_LP(blockNumber[l-1]+i-1,value);
	}
      } else {
	rError("io::read not valid blockType");
      }
    } 
  } else {  // dense

    // k==0
    for (int l2=0; l2<nBlock; ++l2){
      if (blockType[l2] == btSDP) { // SDP part
	int l = blockNumber[l2];
	int size = SDP_blockStruct[l];
	for (int i=0; i<size; ++i) {
	  for (int j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      inputData.C.setElement_SDP(l,i,j,-tmp);
	    }
	  }
	}
      } else if (blockType[l2] == btSOCP) { // SOCP part
	rError("io:: current version does not support SOCP");
      } else if (blockType[l2] == btLP) { // LP part
	for (int j=0; j<blockStruct[l2]; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    inputData.C.setElement_LP(blockNumber[l2]+j,-tmp);
	  }
	}
      } else {
	rError("io::read not valid blockType");
      }
    }

    // k > 0
    for (int k=0; k<m; ++k) {
	  
      for (int l2=0; l2<nBlock; ++l2){
	if (blockType[l2] == btSDP) { // SDP part
	  int l = blockNumber[l2];
	  int size = SDP_blockStruct[l];
	  for (int i=0; i<size; ++i) {
	    for (int j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		inputData.A[k].setElement_SDP(l,i,j,tmp);
	      }
	    }
	  }
	} else if (blockType[l2] == btSOCP) { // SOCP part
	  rError("io:: current version does not support SOCP");
	} else if (blockType[l2] == btLP) { // LP part
	  for (int j=0; j<blockStruct[l2]; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (tmp!=0.0) {
	      inputData.A[k].setElement_LP(blockNumber[l2]+j,tmp);
	    }
	  }
	} else {
	  rError("io::read not valid blockType");
	}
      }

    } // for k

  } // end of 'if (isDataSparse)'

}

void IO::printHeader(FILE* fpout, FILE* Display)
{
  if (fpout) {
    fprintf(fpout,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
  }
  if (Display) {
    fprintf(Display,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
  }
}

void IO::printOneIteration(int pIteration,
			    AverageComplementarity& mu,
			    RatioInitResCurrentRes& theta,
			    SolveInfo& solveInfo,
			    StepLength& alpha,
			    DirectionParameter& beta,
			    FILE* fpout,
			    FILE* Display)
{
  #if REVERSE_PRIMAL_DUAL
  if (Display) {
    fprintf(Display,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.dual, theta.primal,
	    -solveInfo.objValDual,-solveInfo.objValPrimal,
	    alpha.dual, alpha.primal, beta.value);
  }
  if (fpout) {
    fprintf(fpout,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.dual, theta.primal,
	    -solveInfo.objValDual,-solveInfo.objValPrimal,
	    alpha.dual, alpha.primal, beta.value);
  }
  #else
  if (Display) {
    fprintf(Display,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.primal, theta.dual,
	    solveInfo.objValPrimal, solveInfo.objValDual,
	    alpha.primal, alpha.dual, beta.value);
  }
  if (fpout) {
    fprintf(fpout,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.primal, theta.dual,
	    solveInfo.objValPrimal, solveInfo.objValDual,
	    alpha.primal, alpha.dual, beta.value);
  }
  #endif
}

void IO::printLastInfo(int pIteration,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       SolveInfo& solveInfo,
		       StepLength& alpha,
		       DirectionParameter& beta,
		       Residuals& currentRes,
		       Phase & phase,
		       Solutions& currentPt,
		       double cputime,
		       InputData& inputData,
                       WorkVariables& work,
		       ComputeTime& com,
		       Parameter& param,
		       FILE* fpout,
		       FILE* Display,
		       bool printTime)
{
  int nDim = currentPt.nDim;

  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta, fpout, Display);

  double mean = (fabs(solveInfo.objValPrimal)
		 + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // double dominator;
  double relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }

  double gap    = mu.current*nDim;
  double digits = 1000; // 1000 means infinity in this case
  digits = -log10(fabs(PDgap/mean));


  #if DIMACS_PRINT
  double tmp = 0.0;
  double b1 = 0.0;
  for (int k=0; k<inputData.b.nDim; ++k) {
    b1 = max(b1, fabs(inputData.b.ele[k]));
  }
  double c1 = 0.0;
  for (int l=0; l<inputData.C.SDP_sp_nBlock; ++l) {
    SparseMatrix& Cl = inputData.C.SDP_sp_block[l];
    if (Cl.type == SparseMatrix::SPARSE) {
      for (int i=0; i<Cl.NonZeroCount; ++i) {
	c1 = max(c1, fabs(Cl.sp_ele[i]));
      }
    } else if (Cl.type == SparseMatrix::DENSE) {
      for (int i=0; i<Cl.nRow*Cl.nCol; ++i) {
	c1 = max(c1, fabs(Cl.de_ele[i]));
      }
    }
  }
  for (int l=0; l<inputData.C.SOCP_sp_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  for (int l=0; l<inputData.C.LP_sp_nBlock; ++l) {
    c1 = max(c1, fabs(inputData.C.LP_sp_block[l]));
  }
  double p_norm;
  Lal::let(tmp,'=',currentRes.primalVec,'.',currentRes.primalVec);
  p_norm = sqrt(tmp);
  double d_norm = 0.0;
  for (int l=0; l<currentRes.dualMat.SDP_nBlock; ++l) {
    Lal::let(tmp,'=',currentRes.dualMat.SDP_block[l],
	     '.',currentRes.dualMat.SDP_block[l]);
    d_norm += sqrt(tmp);
  }
  for (int l=0; l<currentRes.dualMat.SOCP_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  tmp = 0.0;
  for (int l=0; l<currentRes.dualMat.LP_nBlock; ++l) {
    tmp += currentRes.dualMat.LP_block[l]
      * currentRes.dualMat.LP_block[l];
  }
  d_norm += sqrt(tmp);
  double x_min =  Jal::getMinEigen(currentPt.xMat,work);
  double z_min =  Jal::getMinEigen(currentPt.zMat,work);
					
  //  printf("b1:%e\n",b1);
  //  printf("c1:%e\n",c1);
  //  printf("p_norm:%e\n",p_norm);
  //  printf("d_norm:%e\n",d_norm);
  //  printf("x_min:%e\n",x_min);
  //  printf("z_min:%e\n",z_min);
  
  double ctx = solveInfo.objValPrimal;
  double bty = solveInfo.objValDual;
  double xtz = 0.0;
  Lal::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  double err1 = p_norm / (1+b1);
  double err2 = max( 0.0, - x_min / (1+b1));
  double err3 = d_norm / (1+c1);
  double err4 = max( 0.0, - z_min / (1+c1));
  double err5 = (ctx - bty) / (1 + fabs(ctx) + fabs(bty));
  double err6 = xtz / (1 + fabs(ctx) + fabs(bty));
    
  #endif
  
  
  if (Display) {
    fprintf(Display, "\n");
    phase.display(Display);
    fprintf(Display, "   Iteration = %d\n",  pIteration);
    fprintf(Display, "          mu = ");
    fprintf(Display, param.infPrint, mu.current);
    fprintf(Display, "\n");
    fprintf(Display, "relative gap = ");
    fprintf(Display, param.infPrint, relgap);
    fprintf(Display, "\n");
    fprintf(Display, "        gap = ");
    fprintf(Display, param.infPrint, gap);
    fprintf(Display, "\n");
    fprintf(Display, "     digits = ");
    fprintf(Display, param.infPrint, digits);
    fprintf(Display, "\n");
    #if REVERSE_PRIMAL_DUAL
    fprintf(Display, "objValPrimal = ");
    fprintf(Display, param.infPrint, -solveInfo.objValDual);
    fprintf(Display, "\n");
    fprintf(Display, "objValDual   = ");
    fprintf(Display, param.infPrint, -solveInfo.objValPrimal);
    fprintf(Display, "\n");
    fprintf(Display, "p.feas.error   = ");
    fprintf(Display, param.infPrint, currentRes.normDualMat);
    fprintf(Display, "\n");
    fprintf(Display, "d.feas.error   = ");
    fprintf(Display, param.infPrint, currentRes.normPrimalVec);
    fprintf(Display, "\n");
    #else
    fprintf(Display, "objValPrimal = ");
    fprintf(Display, param.infPrint, solveInfo.objValPrimal);
    fprintf(Display, "\n");
    fprintf(Display, "objValDual   = ");
    fprintf(Display, param.infPrint, solveInfo.objValDual);
    fprintf(Display, "\n");
    fprintf(Display, "p.feas.error   = ");
    fprintf(Display, param.infPrint, currentRes.normPrimalVec);
    fprintf(Display, "\n");
    fprintf(Display, "d.feas.error   = ");
    fprintf(Display, param.infPrint, currentRes.normDualMat);
    fprintf(Display, "\n");
    #endif
    if (printTime == true) {
      fprintf(Display, "total time   = %.6f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(Display, "\n");
    fprintf(Display, "* DIMACS_ERRORS * \n");
    fprintf(Display, "err1 = ");
    fprintf(Display, param.infPrint, err1);
    fprintf(Display,"  [||Ax-b|| / (1+||b||_1)]\n");
    fprintf(Display, "err2 = ");
    fprintf(Display, param.infPrint, err2);
    fprintf(Display,"  [max(0, -lambda(x) / (1+||b||_1))]\n");
    fprintf(Display, "err3 = ");
    fprintf(Display, param.infPrint, err3);
    fprintf(Display,"  [||A^Ty + z - c || / (1+||c||_1)]\n");
    fprintf(Display, "err4 = ");
    fprintf(Display, param.infPrint, err4);
    fprintf(Display,"  [max(0, -lambda(z) / (1+||c||_1))]\n");
    fprintf(Display, "err5 = ");
    fprintf(Display, param.infPrint, err5);
    fprintf(Display,"  [(<c,x> - by) / (1 + |<c,x>| + |by|)]\n");
    fprintf(Display, "err6 = ");
    fprintf(Display, param.infPrint, err6);
    fprintf(Display,"  [<x,z> / (1 + |<c,x>| + |by|)]\n");
    fprintf(Display, "\n");
    #endif
  }
  if (fpout) {
    fprintf(fpout, "\n");
    phase.display(fpout);
    fprintf(fpout, "   Iteration = %d\n",  pIteration);
    fprintf(fpout, "          mu = ");
    fprintf(fpout, param.infPrint, mu.current);
    fprintf(fpout, "\n");
    fprintf(fpout, "relative gap = ");
    fprintf(fpout, param.infPrint, relgap);
    fprintf(fpout, "\n");
    fprintf(fpout, "        gap = ");
    fprintf(fpout, param.infPrint, gap);
    fprintf(fpout, "\n");
    fprintf(fpout, "     digits = ");
    fprintf(fpout, param.infPrint, digits);
    fprintf(fpout, "\n");
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout, "objValPrimal = ");
    fprintf(fpout, param.infPrint, -solveInfo.objValDual);
    fprintf(fpout, "\n");
    fprintf(fpout, "objValDual   = ");
    fprintf(fpout, param.infPrint, -solveInfo.objValPrimal);
    fprintf(fpout, "\n");
    fprintf(fpout, "p.feas.error   = ");
    fprintf(fpout, param.infPrint, currentRes.normDualMat);
    fprintf(fpout, "\n");
    fprintf(fpout, "d.feas.error   = ");
    fprintf(fpout, param.infPrint, currentRes.normPrimalVec);
    fprintf(fpout, "\n");
    #else
    fprintf(fpout, "objValPrimal = ");
    fprintf(fpout, param.infPrint, solveInfo.objValPrimal);
    fprintf(fpout, "\n");
    fprintf(fpout, "objValDual   = ");
    fprintf(fpout, param.infPrint, solveInfo.objValDual);
    fprintf(fpout, "\n");
    fprintf(fpout, "p.feas.error   = ");
    fprintf(fpout, param.infPrint, currentRes.normPrimalVec);
    fprintf(fpout, "\n");
    fprintf(fpout, "d.feas.error   = ");
    fprintf(fpout, param.infPrint, currentRes.normDualMat);
    fprintf(fpout, "\n");
    #endif
    if (printTime == true) {
      fprintf(fpout, "total time   = %.6f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(fpout, "\n");
    fprintf(fpout, "* DIMACS_ERRORS * \n");
    fprintf(fpout, "err1 = ");
    fprintf(fpout, param.infPrint, err1);
    fprintf(fpout,"  [||Ax-b|| / (1+||b||_1)]\n");
    fprintf(fpout, "err2 = ");
    fprintf(fpout, param.infPrint, err2);
    fprintf(fpout,"  [max(0, -lambda(x) / (1+||b||_1))]\n");
    fprintf(fpout, "err3 = ");
    fprintf(fpout, param.infPrint, err3);
    fprintf(fpout,"  [||A^Ty + z - c || / (1+||c||_1)]\n");
    fprintf(fpout, "err4 = ");
    fprintf(fpout, param.infPrint, err4);
    fprintf(fpout,"  [max(0, -lambda(z) / (1+||c||_1))]\n");
    fprintf(fpout, "err5 = ");
    fprintf(fpout, param.infPrint, err5);
    fprintf(fpout,"  [(<c,x> - by) / (1 + |<c,x>| + |by|)]\n");
    fprintf(fpout, "err6 = ");
    fprintf(fpout, param.infPrint, err6);
    fprintf(fpout,"  [<x,z> / (1 + |<c,x>| + |by|)]\n");
    fprintf(fpout, "\n");
    #endif

    fprintf(fpout, "\n\nParameters are\n");
    param.display(fpout);
    com.display(fpout);

    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,-1.0);
    fprintf(fpout,"xMat = \n");
    currentPt.zMat.display(fpout);
    fprintf(fpout,"yMat = \n");
    currentPt.xMat.display(fpout);
    #else
    currentPt.display(fpout);
    #endif
    #endif
  }
}


void IO::printLastInfo(int pIteration,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       SolveInfo& solveInfo,
		       StepLength& alpha,
		       DirectionParameter& beta,
		       Residuals& currentRes,
		       Phase & phase,
		       Solutions& currentPt,
		       double cputime,
		       int nBlock,
		       int* blockStruct,
		       BlockType* blockType,
		       int* blockNumber,
		       InputData& inputData,
                       WorkVariables& work,
		       ComputeTime& com,
		       Parameter& param,
		       FILE* fpout,
		       FILE* Display,
		       bool printTime)
{
  int nDim = currentPt.nDim;

  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta, fpout, Display);

  double mean = (fabs(solveInfo.objValPrimal)
		 + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // double dominator;
  double relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }

  double gap    = mu.current*nDim;
  double digits = 1000; // 1000 means infinity in this case
  digits = -log10(fabs(PDgap/mean));


  #if DIMACS_PRINT
  double tmp = 0.0;
  double b1 = 0.0;
  for (int k=0; k<inputData.b.nDim; ++k) {
    b1 = max(b1, fabs(inputData.b.ele[k]));
  }
  double c1 = 0.0;
  for (int l=0; l<inputData.C.SDP_sp_nBlock; ++l) {
    SparseMatrix& Cl = inputData.C.SDP_sp_block[l];
    if (Cl.type == SparseMatrix::SPARSE) {
      for (int i=0; i<Cl.NonZeroCount; ++i) {
	c1 = max(c1, fabs(Cl.sp_ele[i]));
      }
    } else if (Cl.type == SparseMatrix::DENSE) {
      for (int i=0; i<Cl.nRow*Cl.nCol; ++i) {
	c1 = max(c1, fabs(Cl.de_ele[i]));
      }
    }
  }
  for (int l=0; l<inputData.C.SOCP_sp_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  for (int l=0; l<inputData.C.LP_sp_nBlock; ++l) {
    c1 = max(c1, fabs(inputData.C.LP_sp_block[l]));
  }
  double p_norm;
  Lal::let(tmp,'=',currentRes.primalVec,'.',currentRes.primalVec);
  p_norm = sqrt(tmp);
  double d_norm = 0.0;
  for (int l=0; l<currentRes.dualMat.SDP_nBlock; ++l) {
    Lal::let(tmp,'=',currentRes.dualMat.SDP_block[l],'.',currentRes.dualMat.SDP_block[l]);
    d_norm += sqrt(tmp);
  }
  for (int l=0; l<currentRes.dualMat.SOCP_nBlock; ++l) {
    rError("io:: current version does not support SOCP");
  }
  tmp = 0.0;
  for (int l=0; l<currentRes.dualMat.LP_nBlock; ++l) {
    tmp += currentRes.dualMat.LP_block[l] * currentRes.dualMat.LP_block[l];
  }
  d_norm += sqrt(tmp);
  double x_min =  Jal::getMinEigen(currentPt.xMat,work);
  double z_min =  Jal::getMinEigen(currentPt.zMat,work);
					
  //  printf("b1:%e\n",b1);
  //  printf("c1:%e\n",c1);
  //  printf("p_norm:%e\n",p_norm);
  //  printf("d_norm:%e\n",d_norm);
  //  printf("x_min:%e\n",x_min);
  //  printf("z_min:%e\n",z_min);
  
  double ctx = solveInfo.objValPrimal;
  double bty = solveInfo.objValDual;
  double xtz = 0.0;
  Lal::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  double err1 = p_norm / (1+b1);
  double err2 = max( 0.0, - x_min / (1+b1));
  double err3 = d_norm / (1+c1);
  double err4 = max( 0.0, - z_min / (1+c1));
  double err5 = (ctx - bty) / (1 + fabs(ctx) + fabs(bty));
  double err6 = xtz / (1 + fabs(ctx) + fabs(bty));
    
  #endif

  if (Display) {
    fprintf(Display, "\n");
    phase.display(Display);
    fprintf(Display, "   Iteration = %d\n",  pIteration);
    fprintf(Display, "          mu = ");
    fprintf(Display, param.infPrint, mu.current);
    fprintf(Display, "\n");
    fprintf(Display, "relative gap = ");
    fprintf(Display, param.infPrint, relgap);
    fprintf(Display, "\n");
    fprintf(Display, "         gap = ");
    fprintf(Display, param.infPrint, gap);
    fprintf(Display, "\n");
    fprintf(Display, "      digits = ");
    fprintf(Display, param.infPrint, digits);
    fprintf(Display, "\n");
    #if REVERSE_PRIMAL_DUAL
    fprintf(Display, "objValPrimal = ");
    fprintf(Display, param.infPrint, -solveInfo.objValDual);
    fprintf(Display, "\n");
    fprintf(Display, "objValDual   = ");
    fprintf(Display, param.infPrint, -solveInfo.objValPrimal);
    fprintf(Display, "\n");
    fprintf(Display, "p.feas.error = ");
    fprintf(Display, param.infPrint, currentRes.normDualMat);
    fprintf(Display, "\n");
    fprintf(Display, "d.feas.error = ");
    fprintf(Display, param.infPrint, currentRes.normPrimalVec);
    fprintf(Display, "\n");
    #else
    fprintf(Display, "objValPrimal = ");
    fprintf(Display, param.infPrint, solveInfo.objValPrimal);
    fprintf(Display, "\n");
    fprintf(Display, "objValDual   = ");
    fprintf(Display, param.infPrint, solveInfo.objValDual);
    fprintf(Display, "\n");
    fprintf(Display, "p.feas.error = ");
    fprintf(Display, param.infPrint, currentRes.normPrimalVec);
    fprintf(Display, "\n");
    fprintf(Display, "d.feas.error = ");
    fprintf(Display, param.infPrint, currentRes.normDualMat);
    fprintf(Display, "\n");
    #endif
    if (printTime == true) {
      fprintf(Display, "total time   = %.6f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(Display, "\n");
    fprintf(Display, "* DIMACS_ERRORS * \n");
    fprintf(Display, "err1 = ");
    fprintf(Display, param.infPrint, err1);
    fprintf(Display,"  [||Ax-b|| / (1+||b||_1)]\n");
    fprintf(Display, "err2 = ");
    fprintf(Display, param.infPrint, err2);
    fprintf(Display,"  [max(0, -lambda(x) / (1+||b||_1))]\n");
    fprintf(Display, "err3 = ");
    fprintf(Display, param.infPrint, err3);
    fprintf(Display,"  [||A^Ty + z - c || / (1+||c||_1)]\n");
    fprintf(Display, "err4 = ");
    fprintf(Display, param.infPrint, err4);
    fprintf(Display,"  [max(0, -lambda(z) / (1+||c||_1))]\n");
    fprintf(Display, "err5 = ");
    fprintf(Display, param.infPrint, err5);
    fprintf(Display,"  [(<c,x> - by) / (1 + |<c,x>| + |by|)]\n");
    fprintf(Display, "err6 = ");
    fprintf(Display, param.infPrint, err6);
    fprintf(Display,"  [<x,z> / (1 + |<c,x>| + |by|)]\n");
    fprintf(Display, "\n");
    #endif
  }
  if (fpout) {
    fprintf(fpout, "\n");
    phase.display(fpout);
    fprintf(fpout, "   Iteration = %d\n",  pIteration);
    fprintf(fpout, "          mu = ");
    fprintf(fpout, param.infPrint, mu.current);
    fprintf(fpout, "\n");
    fprintf(fpout, "relative gap = ");
    fprintf(fpout, param.infPrint, relgap);
    fprintf(fpout, "\n");
    fprintf(fpout, "        gap  = ");
    fprintf(fpout, param.infPrint, gap);
    fprintf(fpout, "\n");
    fprintf(fpout, "     digits  = ");
    fprintf(fpout, param.infPrint, digits);
    fprintf(fpout, "\n");
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout, "objValPrimal = ");
    fprintf(fpout, param.infPrint, -solveInfo.objValDual);
    fprintf(fpout, "\n");
    fprintf(fpout, "objValDual   = ");
    fprintf(fpout, param.infPrint, -solveInfo.objValPrimal);
    fprintf(fpout, "\n");
    fprintf(fpout, "p.feas.error = ");
    fprintf(fpout, param.infPrint, currentRes.normDualMat);
    fprintf(fpout, "\n");
    fprintf(fpout, "d.feas.error = ");
    fprintf(fpout, param.infPrint, currentRes.normPrimalVec);
    fprintf(fpout, "\n");
    #else
    fprintf(fpout, "objValPrimal = ");
    fprintf(fpout, param.infPrint, solveInfo.objValPrimal);
    fprintf(fpout, "\n");
    fprintf(fpout, "objValDual   = ");
    fprintf(fpout, param.infPrint, solveInfo.objValDual);
    fprintf(fpout, "\n");
    fprintf(fpout, "p.feas.error = ");
    fprintf(fpout, param.infPrint, currentRes.normPrimalVec);
    fprintf(fpout, "\n");
    fprintf(fpout, "d.feas.error = ");
    fprintf(fpout, param.infPrint, currentRes.normDualMat);
    fprintf(fpout, "\n");
    #endif
    if (printTime == true) {
      fprintf(fpout, "total time   = %.6f\n",cputime);
    }
    #if DIMACS_PRINT
    fprintf(fpout, "\n");
    fprintf(fpout, "* DIMACS_ERRORS * \n");
    fprintf(fpout, "err1 = ");
    fprintf(fpout, param.infPrint, err1);
    fprintf(fpout,"  [||Ax-b|| / (1+||b||_1)]\n");
    fprintf(fpout, "err2 = ");
    fprintf(fpout, param.infPrint, err2);
    fprintf(fpout,"  [max(0, -lambda(x) / (1+||b||_1))]\n");
    fprintf(fpout, "err3 = ");
    fprintf(fpout, param.infPrint, err3);
    fprintf(fpout,"  [||A^Ty + z - c || / (1+||c||_1)]\n");
    fprintf(fpout, "err4 = ");
    fprintf(fpout, param.infPrint, err4);
    fprintf(fpout,"  [max(0, -lambda(z) / (1+||c||_1))]\n");
    fprintf(fpout, "err5 = ");
    fprintf(fpout, param.infPrint, err5);
    fprintf(fpout,"  [(<c,x> - by) / (1 + |<c,x>| + |by|)]\n");
    fprintf(fpout, "err6 = ");
    fprintf(fpout, param.infPrint, err6);
    fprintf(fpout,"  [<x,z> / (1 + |<c,x>| + |by|)]\n");
    fprintf(fpout, "\n");
    #endif

    fprintf(fpout, "\n\nParameters are\n");
    param.display(fpout,param.infPrint);
    com.display(fpout);

    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,-1.0,param.xPrint);
    fprintf(fpout,"xMat = \n");
    displayDenseLinarSpaceLast(currentPt.zMat,
			       nBlock, blockStruct,
			       blockType, blockNumber, fpout,
			       param.XPrint);
    fprintf(fpout,"yMat = \n");
    displayDenseLinarSpaceLast(currentPt.xMat,
			       nBlock, blockStruct,
			       blockType, blockNumber, fpout,
			       param.YPrint);
    #else
    fprintf(fpout,"xMat = \n");
    displayDenseLinarSpaceLast(currentPt.xMat,
			       nBlock, blockStruct,
			       blockType, blockNumber, fpout,
			       param.YPrint);
    fprintf(fpout,"yVec = \n");
    currentPt.yVec.display(fpout,1.0,param.xPrint);
    fprintf(fpout,"zMat = \n");
    displayDenseLinarSpaceLast(currentPt.zMat,
			       nBlock, blockStruct,
			       blockType, blockNumber, fpout,
			       param.XPrint);
    #endif
    #endif
  }
}

void IO::displayDenseLinarSpaceLast(DenseLinearSpace& aMat,
				    int nBlock,
				    int* blockStruct,
				    BlockType* blockType,
				    int* blockNumber,
				    FILE* fpout,
				    char* printFormat)
{
  if (fpout == NULL) {
    return;
  }
  if (strcmp(printFormat,NO_P_FORMAT) == 0) {
    fprintf(fpout,"%s\n",NO_P_FORMAT);
    return;
  }
  fprintf(fpout,"{\n");
  for (int i=0; i<nBlock; i++){
    if (blockType[i] == btSDP) {
      int l = blockNumber[i];
      aMat.SDP_block[l].display(fpout,printFormat);
    } else if (blockType[i] == btSOCP) {
      rError("io:: current version does not support SOCP");
      int l = blockNumber[i];
      aMat.SOCP_block[l].display(fpout,printFormat);
    } else if (blockType[i] == btLP) {
      fprintf(fpout,"{");
      for (int l=0; l<blockStruct[i]-1; ++l) {
	fprintf(fpout,printFormat,aMat.LP_block[blockNumber[i]+l]);
	fprintf(fpout,",");
      }
      if (blockStruct[i] > 0) {
	fprintf(fpout,printFormat,
		aMat.LP_block[blockNumber[i]+blockStruct[i]-1]);
	fprintf(fpout,"}\n");
      } else {
	fprintf(fpout,"  }\n");
      }
    } else {
      rError("io::displayDenseLinarSpaceLast not valid blockType");
    }
  }
  fprintf(fpout,"}\n");
}

} // end of namespace 'sdpa'
