
#ifndef __RELEASE_VERSION_DETECTOR__
#define __RELEASE_VERSION_DETECTOR__

#define MW_FIRST_API_VERSION 700
#define R2017b 700
#define R2018a 800
#define MW_LATEST_API_VERSION 800


#define MW_REL2VER(A) A 

#if defined(MX_COMPAT_32) || defined(MEX_DOUBLE_HANDLE)

#if defined(MATLAB_MEXCMD_RELEASE) || defined(MATLAB_MEXSRC_RELEASE)

/* Errors! Legacy knobs cannot be used with release-based hard knobs */

#if defined(MX_COMPAT_32) && defined(MATLAB_MEXCMD_RELEASE)
#error "MEX command option -R20XXx is incompatible with MX_COMPAT_32"
#endif

#if defined(MEX_DOUBLE_HANDLE) && defined(MATLAB_MEXCMD_RELEASE)
#error "MEX command option -R20XXx is incompatible with MEX_DOUBLE_HANDLE"
#endif

#if defined(MX_COMPAT_32) && defined(MATLAB_MEXSRC_RELEASE)
#error "Source code macro MATLAB_MEXSRC_RELEASE is incompatible with MX_COMPAT_32"
#endif

#if defined(MEX_DOUBLE_HANDLE) && defined(MATLAB_MEXSRC_RELEASE)
#error "Source code macro MATLAB_MEXSRC_RELEASE is incompatible with MEX_DOUBLE_HANDLE"
#endif

#else

/* Legacy knobs are defined  */

#define MATLAB_TARGET_API_VERSION MW_FIRST_API_VERSION

#endif

#else /* defined(MX_COMPAT_32) || defined(MEX_DOUBLE_HANDLE) */

/* No Legacy knobs. Check release-based tag */

#if defined(MATLAB_MEXCMD_RELEASE)
#define MW_MEXCMD_VERSION MW_REL2VER(MATLAB_MEXCMD_RELEASE)
#if MW_MEXCMD_VERSION < MW_FIRST_API_VERSION
#error invalid MATLAB_MEXCMD_RELEASE definition
#endif
#endif

#if defined(MATLAB_MEXSRC_RELEASE)
#define MW_MEXSRC_VERSION MW_REL2VER(MATLAB_MEXSRC_RELEASE)
#if MW_MEXSRC_VERSION < MW_FIRST_API_VERSION
#error invalid MATLAB_MEXSRC_RELEASE definition
#endif
#endif
      
#if defined(MATLAB_DEFAULT_RELEASE)
#define MW_DEFAULT_VERSION MW_REL2VER(MATLAB_DEFAULT_RELEASE)
#if MW_DEFAULT_VERSION < MW_FIRST_API_VERSION
#error invalid MATLAB_DEFAULT_RELEASE definition
#endif
#endif

#if defined(MATLAB_MEXCMD_RELEASE) && defined(MATLAB_MEXSRC_RELEASE)
#if MW_MEXCMD_VERSION != MW_MEXSRC_VERSION
#error "MEX command option -R20XXx is incompatible with MATLAB_MEXSRC_RELEASE"
#endif
#endif

#if defined(MATLAB_MEXCMD_RELEASE) || defined(MATLAB_MEXSRC_RELEASE)

/* Check whether MEXCMD and MEXSRC release tags are compatible */

#if defined(MATLAB_MEXCMD_RELEASE)
#define MATLAB_TARGET_API_VERSION MW_MEXCMD_VERSION
#else
#define MATLAB_TARGET_API_VERSION MW_MEXSRC_VERSION
#endif

#else /* defined(MATLAB_MEXCMD_RELEASE) || defined(MATLAB_MEXSRC_RELEASE) */

#if defined(MATLAB_DEFAULT_RELEASE)
#define MATLAB_TARGET_API_VERSION MW_DEFAULT_VERSION
#else

/* None of the input macros are defined. Use LATEST. */
#define MATLAB_TARGET_API_VERSION MW_LATEST_API_VERSION

#endif /* defined(MATLAB_DEFAULT_RELEASE) */

#endif /* defined(MATLAB_MEXCMD_RELEASE) || defined(MATLAB_MEXSRC_RELEASE) */

#endif /* defined(MX_COMPAT_32) || defined(MEX_DOUBLE_HANDLE) */

#if defined(TARGET_API_VERSION)
#if MATLAB_TARGET_API_VERSION != TARGET_API_VERSION
#error MATLAB_TARGET_API_VERSION != TARGET_API_VERSION
#endif
#else
#define TARGET_API_VERSION MATLAB_TARGET_API_VERSION
#endif

#endif /* __RELEASE_VERSION_DETECTOR__ */
#if defined(WITH_COMMENTS)
/*
 * fintrf.h	- MATLAB/FORTRAN interface header file. This file
 *		  contains the declaration of the pointer type needed
 *		  by the MATLAB/FORTRAN interface.
 *
 * Copyright 1984-2011 The MathWorks, Inc.
 * All Rights Reserved.
 */
#endif
#if defined(__LP64__) || defined(_M_AMD64) || defined(__amd64)
# define mwpointer integer*8
# define mwPointer integer*8
# define MWPOINTER INTEGER*8
#else
# define mwpointer integer*4
# define mwPointer integer*4
# define MWPOINTER INTEGER*4
#endif

#if defined(MX_COMPAT_32)
# define mwsize  integer*4
# define mwSize  integer*4
# define MWSIZE  INTEGER*4
# define mwindex integer*4
# define mwIndex integer*4
# define MWINDEX INTEGER*4
# define mwsignedindex integer*4
# define mwSignedIndex integer*4
# define MWSIGNEDINDEX INTEGER*4
#else
# define mwsize  mwpointer
# define mwSize  mwpointer
# define MWSIZE  MWPOINTER
# define mwindex mwpointer
# define mwIndex mwpointer
# define MWINDEX MWPOINTER
# define mwsignedindex mwpointer
# define mwSignedIndex mwpointer
# define MWSIGNEDINDEX MWPOINTER
#endif

#if defined(WITH_COMMENTS)
/*
 * Compatibility layer for MEX files using the 32-bit mxArray APIs
 */
#endif

#if defined(MX_COMPAT_32)

#define mxCalcSingleSubscript mxCalcSingleSubscript700
#define mxCreateCellArray mxCreateCellArray700
#define mxCreateCellMatrix mxCreateCellMatrix700
#define mxCreateCharArray mxCreateCharArray700
#define mxCreateCharMatrixFromStrings mxCreateCharMatrixFromStrs700
#define mxCreateDoubleMatrix mxCreateDoubleMatrix700
#define mxCreateNumericArray mxCreateNumericArray700
#define mxCreateNumericMatrix mxCreateNumericMatrix700
#define mxCreateSparse mxCreateSparse700
#define mxCreateStructArray mxCreateStructArray700
#define mxCreateStructMatrix mxCreateStructMatrix700
#define mxGetCell mxGetCell700
#define mxGetDimensions mxGetDimensions700
#define mxGetField mxGetField700
#define mxGetFieldByNumber mxGetFieldByNumber700
#define mxGetProperty mxGetProperty700
#define mxGetIr mxGetIr700
#define mxGetJc mxGetJc700
#define mxGetNumberOfDimensions mxGetNumberOfDimensions700
#define mxGetNzmax mxGetNzmax700
#define mxGetString mxGetString700
#define mxSetCell mxSetCell700
#define mxSetDimensions mxSetDimensions700
#define mxSetField mxSetField700
#define mxSetFieldByNumber mxSetFieldByNumber700
#define mxSetProperty mxSetProperty700
#define mxSetIr mxSetIr700
#define mxSetJc mxSetJc700
#define mxSetM mxSetM700
#define mxSetN mxSetN700
#define mxSetNzmax mxSetNzmax700
#define mxGetM mxGetM700
#define mxGetN mxGetN700
#define mxGetNumberOfElements mxGetNumberOfElements700
#define mxGetElementSize mxGetElementSize700
#define mxMalloc mxMalloc700
#define mxCalloc mxCalloc700
#define mxRealloc mxRealloc700
#define mxCopyReal4ToPtr mxCopyReal4ToPtr700
#define mxCopyPtrToReal4 mxCopyPtrToReal4700
#define mxCopyReal8ToPtr mxCopyReal8ToPtr700
#define mxCopyPtrToReal8 mxCopyPtrToReal8700
#define mxCopyCharacterToPtr mxCopyCharacterToPtr700
#define mxCopyPtrToCharacter mxCopyPtrToCharacter700
#define mxCopyInteger1ToPtr mxCopyInteger1ToPtr700
#define mxCopyPtrToInteger1 mxCopyPtrToInteger1700
#define mxCopyInteger2ToPtr mxCopyInteger2ToPtr700
#define mxCopyPtrToInteger2 mxCopyPtrToInteger2700
#define mxCopyInteger4ToPtr mxCopyInteger4ToPtr700
#define mxCopyPtrToInteger4 mxCopyPtrToInteger4700
#define mxCopyInteger8ToPtr mxCopyInteger8ToPtr700
#define mxCopyPtrToInteger8 mxCopyPtrToInteger8700
#define mxCopyPtrToPtrArray mxCopyPtrToPtrArray700
#define mxCopyComplex16ToPtr mxCopyComplex16ToPtr700
#define mxCopyPtrToComplex16 mxCopyPtrToComplex16700
#define mxCopyComplex8ToPtr mxCopyComplex8ToPtr700
#define mxCopyPtrToComplex8 mxCopyPtrToComplex8700
#define mxCopyMWIndexToPtr mxCopyMWIndexToPtr700
#define mxCopyPtrToMWIndex mxCopyPtrToMWIndex700
#define mxCreateCharMatrixFromStrs mxCreateCharMatrixFromStrs700
#define MXCALCSINGLESUBSCRIPT MXCALCSINGLESUBSCRIPT700
#define MXCREATECELLARRAY MXCREATECELLARRAY700
#define MXCREATECELLMATRIX MXCREATECELLMATRIX700
#define MXCREATECHARARRAY MXCREATECHARARRAY700
#define MXCREATECHARMATRIXFROMSTRINGS MXCREATECHARMATRIXFROMSTRS700
#define MXCREATEDOUBLEMATRIX MXCREATEDOUBLEMATRIX700
#define MXCREATENUMERICARRAY MXCREATENUMERICARRAY700
#define MXCREATENUMERICMATRIX MXCREATENUMERICMATRIX700
#define MXCREATESPARSE MXCREATESPARSE700
#define MXCREATESTRUCTARRAY MXCREATESTRUCTARRAY700
#define MXCREATESTRUCTMATRIX MXCREATESTRUCTMATRIX700
#define MXGETCELL MXGETCELL700
#define MXGETDIMENSIONS MXGETDIMENSIONS700
#define MXGETFIELD MXGETFIELD700
#define MXGETFIELDBYNUMBER MXGETFIELDBYNUMBER700
#define MXGETPROPERTY MXGETPROPERTY700
#define MXGETIR MXGETIR700
#define MXGETJC MXGETJC700
#define MXGETNUMBEROFDIMENSIONS MXGETNUMBEROFDIMENSIONS700
#define MXGETNZMAX MXGETNZMAX700
#define MXGETSTRING MXGETSTRING700
#define MXSETCELL MXSETCELL700
#define MXSETDIMENSIONS MXSETDIMENSIONS700
#define MXSETFIELD MXSETFIELD700
#define MXSETFIELDBYNUMBER MXSETFIELDBYNUMBER700
#define MXSETPROPERTY MXSETPROPERTY700
#define MXSETIR MXSETIR700
#define MXSETJC MXSETJC700
#define MXSETM MXSETM700
#define MXSETN MXSETN700
#define MXSETNZMAX MXSETNZMAX700
#define MXGETM MXGETM700
#define MXGETN MXGETN700
#define MXGETNUMBEROFELEMENTS MXGETNUMBEROFELEMENTS700
#define MXGETELEMENTSIZE MXGETELEMENTSIZE700
#define MXMALLOC MXMALLOC700
#define MXCALLOC MXCALLOC700
#define MXREALLOC MXREALLOC700
#define MXCOPYREAL4TOPTR MXCOPYREAL4TOPTR700
#define MXCOPYPTRTOREAL4 MXCOPYPTRTOREAL4700
#define MXCOPYREAL8TOPTR MXCOPYREAL8TOPTR700
#define MXCOPYPTRTOREAL8 MXCOPYPTRTOREAL8700
#define MXCOPYCHARACTERTOPTR MXCOPYCHARACTERTOPTR700
#define MXCOPYPTRTOCHARACTER MXCOPYPTRTOCHARACTER700
#define MXCOPYINTEGER1TOPTR MXCOPYINTEGER1TOPTR700
#define MXCOPYPTRTOINTEGER1 MXCOPYPTRTOINTEGER1700
#define MXCOPYINTEGER2TOPTR MXCOPYINTEGER2TOPTR700
#define MXCOPYPTRTOINTEGER2 MXCOPYPTRTOINTEGER2700
#define MXCOPYINTEGER4TOPTR MXCOPYINTEGER4TOPTR700
#define MXCOPYPTRTOINTEGER4 MXCOPYPTRTOINTEGER4700
#define MXCOPYINTEGER8TOPTR MXCOPYINTEGER8TOPTR700
#define MXCOPYPTRTOINTEGER8 MXCOPYPTRTOINTEGER8700
#define MXCOPYPTRTOPTRARRAY MXCOPYPTRTOPTRARRAY700
#define MXCOPYCOMPLEX16TOPTR MXCOPYCOMPLEX16TOPTR700
#define MXCOPYPTRTOCOMPLEX16 MXCOPYPTRTOCOMPLEX16700
#define MXCOPYCOMPLEX8TOPTR MXCOPYCOMPLEX8TOPTR700
#define MXCOPYPTRTOCOMPLEX8 MXCOPYPTRTOCOMPLEX8700
#define MXCOPYMWINDEXTOPTR MXCOPYMWINDEXTOPTR700
#define MXCOPYPTRTOMWINDEX MXCOPYPTRTOMWINDEX700
#define MXCREATECHARMATRIXFROMSTRS MXCREATECHARMATRIXFROMSTRS700
#define mxcalcsinglesubscript mxcalcsinglesubscript700
#define mxcreatecellarray mxcreatecellarray700
#define mxcreatecellmatrix mxcreatecellmatrix700
#define mxcreatechararray mxcreatechararray700
#define mxcreatecharmatrixfromstrings mxcreatecharmatrixfromstrs700
#define mxcreatedoublematrix mxcreatedoublematrix700
#define mxcreatenumericarray mxcreatenumericarray700
#define mxcreatenumericmatrix mxcreatenumericmatrix700
#define mxcreatesparse mxcreatesparse700
#define mxcreatestructarray mxcreatestructarray700
#define mxcreatestructmatrix mxcreatestructmatrix700
#define mxgetcell mxgetcell700
#define mxgetdimensions mxgetdimensions700
#define mxgetfield mxgetfield700
#define mxgetfieldbynumber mxgetfieldbynumber700
#define mxgetproperty mxgetproperty700
#define mxgetir mxgetir700
#define mxgetjc mxgetjc700
#define mxgetnumberofdimensions mxgetnumberofdimensions700
#define mxgetnzmax mxgetnzmax700
#define mxgetstring mxgetstring700
#define mxsetcell mxsetcell700
#define mxsetdimensions mxsetdimensions700
#define mxsetfield mxsetfield700
#define mxsetfieldbynumber mxsetfieldbynumber700
#define mxsetproperty mxsetproperty700
#define mxsetir mxsetir700
#define mxsetjc mxsetjc700
#define mxsetm mxsetm700
#define mxsetn mxsetn700
#define mxsetnzmax mxsetnzmax700
#define mxgetm mxgetm700
#define mxgetn mxgetn700
#define mxgetnumberofelements mxgetnumberofelements700
#define mxgetelementsize mxgetelementsize700
#define mxmalloc mxmalloc700
#define mxcalloc mxcalloc700
#define mxrealloc mxrealloc700
#define mxcopyreal4toptr mxcopyreal4toptr700
#define mxcopyptrtoreal4 mxcopyptrtoreal4700
#define mxcopyreal8toptr mxcopyreal8toptr700
#define mxcopyptrtoreal8 mxcopyptrtoreal8700
#define mxcopycharactertoptr mxcopycharactertoptr700
#define mxcopyptrtocharacter mxcopyptrtocharacter700
#define mxcopyinteger1toptr mxcopyinteger1toptr700
#define mxcopyptrtointeger1 mxcopyptrtointeger1700
#define mxcopyinteger2toptr mxcopyinteger2toptr700
#define mxcopyptrtointeger2 mxcopyptrtointeger2700
#define mxcopyinteger4toptr mxcopyinteger4toptr700
#define mxcopyptrtointeger4 mxcopyptrtointeger4700
#define mxcopyinteger8toptr mxcopyinteger8toptr700
#define mxcopyptrtointeger8 mxcopyptrtointeger8700
#define mxcopyptrtoptrarray mxcopyptrtoptrarray700
#define mxcopycomplex16toptr mxcopycomplex16toptr700
#define mxcopyptrtocomplex16 mxcopyptrtocomplex16700
#define mxcopycomplex8toptr mxcopycomplex8toptr700
#define mxcopyptrtocomplex8 mxcopyptrtocomplex8700
#define mxcopymwindextoptr mxcopymwindextoptr700
#define mxcopyptrtomwindex mxcopyptrtomwindex700
#define mxcreatecharmatrixfromstrs mxcreatecharmatrixfromstrs700

#else

#ifndef __linux
#define mxCalcSingleSubscript mxCalcSingleSubscript730
#define mxCreateCellArray mxCreateCellArray730
#define mxCreateCellMatrix mxCreateCellMatrix730
#define mxCreateCharArray mxCreateCharArray730
#define mxCreateCharMatrixFromStrings mxCreateCharMatrixFromStrs730
#define mxCreateDoubleMatrix mxCreateDoubleMatrix730
#define mxCreateNumericArray mxCreateNumericArray730
#define mxCreateNumericMatrix mxCreateNumericMatrix730
#define mxCreateSparse mxCreateSparse730
#define mxCreateStructArray mxCreateStructArray730
#define mxCreateStructMatrix mxCreateStructMatrix730
#define mxGetCell mxGetCell730
#define mxGetDimensions mxGetDimensions730
#define mxGetField mxGetField730
#define mxGetFieldByNumber mxGetFieldByNumber730
#define mxGetProperty mxGetProperty730
#define mxGetIr mxGetIr730
#define mxGetJc mxGetJc730
#define mxGetNumberOfDimensions mxGetNumberOfDimensions730
#define mxGetNzmax mxGetNzmax730
#define mxGetString mxGetString730
#define mxSetCell mxSetCell730
#define mxSetDimensions mxSetDimensions730
#define mxSetField mxSetField730
#define mxSetFieldByNumber mxSetFieldByNumber730
#define mxSetProperty mxSetProperty730
#define mxSetIr mxSetIr730
#define mxSetJc mxSetJc730
#define mxSetM mxSetM730
#define mxSetN mxSetN730
#define mxSetNzmax mxSetNzmax730
#define mxGetM mxGetM730
#define mxGetN mxGetN730
#define mxGetNumberOfElements mxGetNumberOfElements730
#define mxGetElementSize mxGetElementSize730
#define mxMalloc mxMalloc730
#define mxCalloc mxCalloc730
#define mxRealloc mxRealloc730
#define mxCopyReal4ToPtr mxCopyReal4ToPtr730
#define mxCopyPtrToReal4 mxCopyPtrToReal4730
#define mxCopyReal8ToPtr mxCopyReal8ToPtr730
#define mxCopyPtrToReal8 mxCopyPtrToReal8730
#define mxCopyCharacterToPtr mxCopyCharacterToPtr730
#define mxCopyPtrToCharacter mxCopyPtrToCharacter730
#define mxCopyInteger1ToPtr mxCopyInteger1ToPtr730
#define mxCopyPtrToInteger1 mxCopyPtrToInteger1730
#define mxCopyInteger2ToPtr mxCopyInteger2ToPtr730
#define mxCopyPtrToInteger2 mxCopyPtrToInteger2730
#define mxCopyInteger4ToPtr mxCopyInteger4ToPtr730
#define mxCopyPtrToInteger4 mxCopyPtrToInteger4730
#define mxCopyInteger8ToPtr mxCopyInteger8ToPtr730
#define mxCopyPtrToInteger8 mxCopyPtrToInteger8730
#define mxCopyPtrToPtrArray mxCopyPtrToPtrArray730
#define mxCopyComplex16ToPtr mxCopyComplex16ToPtr730
#define mxCopyPtrToComplex16 mxCopyPtrToComplex16730
#define mxCopyComplex8ToPtr mxCopyComplex8ToPtr730
#define mxCopyPtrToComplex8 mxCopyPtrToComplex8730
#define mxCopyMWIndexToPtr mxCopyMWIndexToPtr730
#define mxCopyPtrToMWIndex mxCopyPtrToMWIndex730
#define mxCreateCharMatrixFromStrs mxCreateCharMatrixFromStrs730
#define MXCALCSINGLESUBSCRIPT MXCALCSINGLESUBSCRIPT730
#define MXCREATECELLARRAY MXCREATECELLARRAY730
#define MXCREATECELLMATRIX MXCREATECELLMATRIX730
#define MXCREATECHARARRAY MXCREATECHARARRAY730
#define MXCREATECHARMATRIXFROMSTRINGS MXCREATECHARMATRIXFROMSTRS730
#define MXCREATEDOUBLEMATRIX MXCREATEDOUBLEMATRIX730
#define MXCREATENUMERICARRAY MXCREATENUMERICARRAY730
#define MXCREATENUMERICMATRIX MXCREATENUMERICMATRIX730
#define MXCREATESPARSE MXCREATESPARSE730
#define MXCREATESTRUCTARRAY MXCREATESTRUCTARRAY730
#define MXCREATESTRUCTMATRIX MXCREATESTRUCTMATRIX730
#define MXGETCELL MXGETCELL730
#define MXGETDIMENSIONS MXGETDIMENSIONS730
#define MXGETFIELD MXGETFIELD730
#define MXGETFIELDBYNUMBER MXGETFIELDBYNUMBER730
#define MXGETPROPERTY MXGETPROPERTY730
#define MXGETIR MXGETIR730
#define MXGETJC MXGETJC730
#define MXGETNUMBEROFDIMENSIONS MXGETNUMBEROFDIMENSIONS730
#define MXGETNZMAX MXGETNZMAX730
#define MXGETSTRING MXGETSTRING730
#define MXSETCELL MXSETCELL730
#define MXSETDIMENSIONS MXSETDIMENSIONS730
#define MXSETFIELD MXSETFIELD730
#define MXSETFIELDBYNUMBER MXSETFIELDBYNUMBER730
#define MXSETPROPERTY MXSETPROPERTY730
#define MXSETIR MXSETIR730
#define MXSETJC MXSETJC730
#define MXSETM MXSETM730
#define MXSETN MXSETN730
#define MXSETNZMAX MXSETNZMAX730
#define MXGETM MXGETM730
#define MXGETN MXGETN730
#define MXGETNUMBEROFELEMENTS MXGETNUMBEROFELEMENTS730
#define MXGETELEMENTSIZE MXGETELEMENTSIZE730
#define MXMALLOC MXMALLOC730
#define MXCALLOC MXCALLOC730
#define MXREALLOC MXREALLOC730
#define MXCOPYREAL4TOPTR MXCOPYREAL4TOPTR730
#define MXCOPYPTRTOREAL4 MXCOPYPTRTOREAL4730
#define MXCOPYREAL8TOPTR MXCOPYREAL8TOPTR730
#define MXCOPYPTRTOREAL8 MXCOPYPTRTOREAL8730
#define MXCOPYCHARACTERTOPTR MXCOPYCHARACTERTOPTR730
#define MXCOPYPTRTOCHARACTER MXCOPYPTRTOCHARACTER730
#define MXCOPYINTEGER1TOPTR MXCOPYINTEGER1TOPTR730
#define MXCOPYPTRTOINTEGER1 MXCOPYPTRTOINTEGER1730
#define MXCOPYINTEGER2TOPTR MXCOPYINTEGER2TOPTR730
#define MXCOPYPTRTOINTEGER2 MXCOPYPTRTOINTEGER2730
#define MXCOPYINTEGER4TOPTR MXCOPYINTEGER4TOPTR730
#define MXCOPYPTRTOINTEGER4 MXCOPYPTRTOINTEGER4730
#define MXCOPYINTEGER8TOPTR MXCOPYINTEGER8TOPTR730
#define MXCOPYPTRTOINTEGER8 MXCOPYPTRTOINTEGER8730
#define MXCOPYPTRTOPTRARRAY MXCOPYPTRTOPTRARRAY730
#define MXCOPYCOMPLEX16TOPTR MXCOPYCOMPLEX16TOPTR730
#define MXCOPYPTRTOCOMPLEX16 MXCOPYPTRTOCOMPLEX16730
#define MXCOPYCOMPLEX8TOPTR MXCOPYCOMPLEX8TOPTR730
#define MXCOPYPTRTOCOMPLEX8 MXCOPYPTRTOCOMPLEX8730
#define MXCOPYMWINDEXTOPTR MXCOPYMWINDEXTOPTR730
#define MXCOPYPTRTOMWINDEX MXCOPYPTRTOMWINDEX730
#define MXCREATECHARMATRIXFROMSTRS MXCREATECHARMATRIXFROMSTRS730
#define mxcalcsinglesubscript mxcalcsinglesubscript730
#define mxcreatecellarray mxcreatecellarray730
#define mxcreatecellmatrix mxcreatecellmatrix730
#define mxcreatechararray mxcreatechararray730
#define mxcreatecharmatrixfromstrings mxcreatecharmatrixfromstrs730
#define mxcreatedoublematrix mxcreatedoublematrix730
#define mxcreatenumericarray mxcreatenumericarray730
#define mxcreatenumericmatrix mxcreatenumericmatrix730
#define mxcreatesparse mxcreatesparse730
#define mxcreatestructarray mxcreatestructarray730
#define mxcreatestructmatrix mxcreatestructmatrix730
#define mxgetcell mxgetcell730
#define mxgetdimensions mxgetdimensions730
#define mxgetfield mxgetfield730
#define mxgetfieldbynumber mxgetfieldbynumber730
#define mxgetproperty mxgetproperty730
#define mxgetir mxgetir730
#define mxgetjc mxgetjc730
#define mxgetnumberofdimensions mxgetnumberofdimensions730
#define mxgetnzmax mxgetnzmax730
#define mxgetstring mxgetstring730
#define mxsetcell mxsetcell730
#define mxsetdimensions mxsetdimensions730
#define mxsetfield mxsetfield730
#define mxsetfieldbynumber mxsetfieldbynumber730
#define mxsetproperty mxsetproperty730
#define mxsetir mxsetir730
#define mxsetjc mxsetjc730
#define mxsetm mxsetm730
#define mxsetn mxsetn730
#define mxsetnzmax mxsetnzmax730
#define mxgetm mxgetm730
#define mxgetn mxgetn730
#define mxgetnumberofelements mxgetnumberofelements730
#define mxgetelementsize mxgetelementsize730
#define mxmalloc mxmalloc730
#define mxcalloc mxcalloc730
#define mxrealloc mxrealloc730
#define mxcopyreal4toptr mxcopyreal4toptr730
#define mxcopyptrtoreal4 mxcopyptrtoreal4730
#define mxcopyreal8toptr mxcopyreal8toptr730
#define mxcopyptrtoreal8 mxcopyptrtoreal8730
#define mxcopycharactertoptr mxcopycharactertoptr730
#define mxcopyptrtocharacter mxcopyptrtocharacter730
#define mxcopyinteger1toptr mxcopyinteger1toptr730
#define mxcopyptrtointeger1 mxcopyptrtointeger1730
#define mxcopyinteger2toptr mxcopyinteger2toptr730
#define mxcopyptrtointeger2 mxcopyptrtointeger2730
#define mxcopyinteger4toptr mxcopyinteger4toptr730
#define mxcopyptrtointeger4 mxcopyptrtointeger4730
#define mxcopyinteger8toptr mxcopyinteger8toptr730
#define mxcopyptrtointeger8 mxcopyptrtointeger8730
#define mxcopyptrtoptrarray mxcopyptrtoptrarray730
#define mxcopycomplex16toptr mxcopycomplex16toptr730
#define mxcopyptrtocomplex16 mxcopyptrtocomplex16730
#define mxcopycomplex8toptr mxcopycomplex8toptr730
#define mxcopyptrtocomplex8 mxcopyptrtocomplex8730
#define mxcopymwindextoptr mxcopymwindextoptr730
#define mxcopyptrtomwindex mxcopyptrtomwindex730
#define mxcreatecharmatrixfromstrs mxcreatecharmatrixfromstrs730
#endif

#endif
#ifndef __MX_API_VER_HPP__
#define __MX_API_VER_HPP__

/* Current MATRIX published API version */
#define MX_CURRENT_API_VER z'08000000'

/* Backward compatible current MATRIX published API version */
#define MX_API_VER MX_CURRENT_API_VER

/* Backward compatible MATRIX published API versions */
#define MX_LAST_32BIT_VER z'07000000'
#define MX_LAST_SEPARATE_COMPLEX_VER z'07300000'

/* Required MEX-file MATRIX published API version */
#if TARGET_API_VERSION == 700
#if defined(MX_COMPAT_32)
#define MX_TARGET_API_VER MX_LAST_32BIT_VER
#else
#define MX_TARGET_API_VER MX_LAST_SEPARATE_COMPLEX_VER
#endif
#else
#define MX_TARGET_API_VER MX_CURRENT_API_VER
#endif

/*
 * The following macros enable conditional compilation based on the
 * target published API. The macros can be used in a single source file
 * that is intended to be built against multiple matrix API versions.
 *
 * MX_HAS_64BIT_ARRAY_DIMS evaluates to a non-zero value if array
 * dimensions are 64 bits wide.
 *
 * MX_HAS_INTERLEAVED_COMPLEX evaluates to a non-zero value if complex
 * array data is interleaved.
 *
 */
#define MX_HAS_64BIT_ARRAY_DIMS MX_TARGET_API_VER > MX_LAST_32BIT_VER
#define MX_HAS_INTERLEAVED_COMPLEX MX_TARGET_API_VER > MX_LAST_SEPARATE_COMPLEX_VER

#endif /* __MX_API_VER_HPP__ */
