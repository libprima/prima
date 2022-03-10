
#ifndef MXHBB__RELEASE_VERSION_DETECTOR__
#define MXHBB__RELEASE_VERSION_DETECTOR__

#define MW_FIRST_API_VERSION 700
#define R2017b 700
#define R2018a 800
#define R2018b 800
#define R2019a 800
#define R2019b 800
#define R201aa 800
#define R201ab 800
#define R201ba 800
#define R201bb 800
#define R201ca 800
#define R201cb 800
#define R201da 800
#define R201db 800
#define R201ea 800
#define R201eb 800
#define R201fa 800
#define R201fb 800
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

/* None of the input macros are defined. Use LATEST */
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

#endif /* MXHBB__RELEASE_VERSION_DETECTOR__ */
#if defined(WITH_COMMENTS)
/*
 * fintrf.h	- MATLAB/FORTRAN interface header file. This file
 *		  contains the declaration of the pointer type needed
 *		  by the MATLAB/FORTRAN interface.
 *
 * Copyright 1984-2018 The MathWorks, Inc.
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

#if defined(TARGET_API_VERSION)
#if !(TARGET_API_VERSION == 700 || TARGET_API_VERSION == 800)
#error invalid TARGET_VERSION_API definition
#elif defined(MEX_DOUBLE_HANDLE) && TARGET_API_VERSION != 700
#error It is illegal to use MEX_DOUBLE_HANDLE with linear versioning
#elif defined(MX_COMPAT_32) && TARGET_API_VERSION != 700
#error It is illegal to use MX_COMPAT_32 with linear versioning
#endif
#endif

#if TARGET_API_VERSION == 800

#if defined(WITH_COMMENTS)
/*
 * Interleaved Complex (IC) mxArray APIs
 */
#endif

#define mxGetPi mxGetPiIsDeprecated
#define MXGETPI MXGETPIISDEPRECATED
#define mxgetpi mxgetpiisdeprecated
#define mxSetPi mxSetPiIsDeprecated
#define MXSETPI MXSETPIISDEPRECATED
#define mxsetpi mxsetpiisdeprecated
#define mxGetImageData mxGetImagDataIsDeprecated
#define MXGETIMAGEDATA MXGETIMAGDATAISDEPRECATED
#define mxgetimagedata mxgetimagdataisdeprecated
#define mxSetImageData mxSetImagDataIsDeprecated
#define MXSETIMAGEDATA MXSETIMAGDATAISDEPRECATED
#define mxsetimagedata mxsetimagdataisdeprecated
#define mxCreateFull mxCreateFullDeprecated
#define MXCREATEFULL MXCREATEFULLDEPRECATED
#define mxcreatefull mxcreatefulldeprecated

#define mxFree mxFree800
#define MXFREE MXFREE800
#define mxfree mxfree800
#define mxDestroyArray mxDestroyArray800
#define MXDESTROYARRAY MXDESTROYARRAY800
#define mxdestroyarray mxdestroyarray800
#define mxCreateString mxCreateString800
#define MXCREATESTRING MXCREATESTRING800
#define mxcreatestring mxcreatestring800
#define mxCreateDoubleScalar mxCreateDoubleScalar800
#define MXCREATEDOUBLESCALAR MXCREATEDOUBLESCALAR800
#define mxcreatedoublescalar mxcreatedoublescalar800
#define mxIsDouble mxIsDouble800
#define MXISDOUBLE MXISDOUBLE800
#define mxisdouble mxisdouble800
#define mxIsSingle mxIsSingle800
#define MXISSINGLE MXISSINGLE800
#define mxissingle mxissingle800
#define mxIsNumeric mxIsNumeric800
#define MXISNUMERIC MXISNUMERIC800
#define mxisnumeric mxisnumeric800
#define mxIsInt64 mxIsInt64800
#define MXISINT64 MXISINT64800
#define mxisint64 mxisint64800
#define mxIsUint64 mxIsUint64800
#define MXISUINT64 MXISUINT64800
#define mxisuint64 mxisuint64800
#define mxIsInt32 mxIsInt32800
#define MXISINT32 MXISINT32800
#define mxisint32 mxisint32800
#define mxIsUint32 mxIsUint32800
#define MXISUINT32 MXISUINT32800
#define mxisuint32 mxisuint32800
#define mxIsInt16 mxIsInt16800
#define MXISINT16 MXISINT16800
#define mxisint16 mxisint16800
#define mxIsUint16 mxIsUint16800
#define MXISUINT16 MXISUINT16800
#define mxisuint16 mxisuint16800
#define mxIsInt8 mxIsInt8800
#define MXISINT8 MXISINT8800
#define mxisint8 mxisint8800
#define mxIsUint8 mxIsUint8800
#define MXISUINT8 MXISUINT8800
#define mxisuint8 mxisuint8800
#define mxIsChar mxIsChar800
#define MXISCHAR MXISCHAR800
#define mxischar mxischar800
#define mxIsLogical mxIsLogical800
#define MXISLOGICAL MXISLOGICAL800
#define mxislogical mxislogical800
#define mxIsSparse mxIsSparse800
#define MXISSPARSE MXISSPARSE800
#define mxissparse mxissparse800
#define mxIsStruct mxIsStruct800
#define MXISSTRUCT MXISSTRUCT800
#define mxisstruct mxisstruct800
#define mxIsCell mxIsCell800
#define MXISCELL MXISCELL800
#define mxiscell mxiscell800
#define mxIsClass mxIsClass800
#define MXISCLASS MXISCLASS800
#define mxisclass mxisclass800
#define mxisfinite mxisfinite800
#define MXISFINITE MXISFINITE800
#define mxisfinite mxisfinite800
#define mxIsInf mxIsInf800
#define MXISINF MXISINF800
#define mxisinf mxisinf800
#define mxIsNaN mxIsNaN800
#define MXISNAN MXISNAN800
#define mxisnan mxisnan800
#define mxIsEmpty mxIsEmpty800
#define MXISEMPTY MXISEMPTY800
#define mxisempty mxisempty800
#define mxIsFromGlobalWS mxIsFromGlobalWS800
#define MXISFROMGLOBALWS MXISFROMGLOBALWS800
#define mxisfromglobalws mxisfromglobalws800
#define mxGetNaN mxGetNaN800
#define MXGETNAN MXGETNAN800
#define mxgetnan mxgetnan800
#define mxGetClassId mxGetClassId800
#define MXGETCLASSID MXGETCLASSID800
#define mxgetclassid mxgetclassid800
#define mxGetClassName mxGetClassName800
#define MXGETCLASSNAME MXGETCLASSNAME800
#define mxgetclassname mxgetclassname800
#define mxGetNumberOfFields mxGetNumberOfFields800
#define MXGETNUMBEROFFIELDS MXGETNUMBEROFFIELDS800
#define mxgetnumberoffields mxgetnumberoffields800
#define mxGetFieldNameByNumber mxGetFieldNameByNumber800
#define MXGETFIELDNAMEBYNUMBER MXGETFIELDNAMEBYNUMBER800
#define mxgetfieldnamebynumber mxgetfieldnamebynumber800
#define mxGetFieldNumber mxGetFieldNumber800
#define MXGETFIELDNUMBER MXGETFIELDNUMBER800
#define mxgetfieldnumber mxgetfieldnumber800
#define mxAddField mxAddField800
#define MXADDFIELD MXADDFIELD800
#define mxaddfield mxaddfield800
#define mxRemoveField mxRemoveField800
#define MXREMOVEFIELD MXREMOVEFIELD800
#define mxremovefield mxremovefield800
#define mxClassIDFromClassName mxClassIDFromClassName800
#define MXCLASSIDFROMCLASSNAME MXCLASSIDFROMCLASSNAME800
#define mxclassidfromclassname mxclassidfromclassname800
#define mxGetN mxgetn800
#define MXGETN MXGETN800
#define mxgetn mxgetn800
#define mxCreateCellArray mxCreateCellArray800
#define MXCREATECELLARRAY MXCREATECELLARRAY800
#define mxcreatecellarray mxcreatecellarray800
#define mxCreateCellMatrix mxCreateCellMatrix800
#define MXCREATECELLMATRIX MXCREATECELLMATRIX800
#define mxcreatecellmatrix mxcreatecellmatrix800
#define mxCreateCharArray mxCreateCharArray800
#define MXCREATECHARARRAY MXCREATECHARARRAY800
#define mxcreatechararray mxcreatechararray800
#define mxCreateDoubleMatrix mxCreateDoubleMatrix800
#define MXCREATEDOUBLEMATRIX MXCREATEDOUBLEMATRIX800
#define mxcreatedoublematrix mxcreatedoublematrix800
#define mxCreateNumericArray mxCreateNumericArray800
#define MXCREATENUMERICARRAY MXCREATENUMERICARRAY800
#define mxcreatenumericarray mxcreatenumericarray800
#define mxCreateNumericMatrix mxCreateNumericMatrix800
#define MXCREATENUMERICMATRIX MXCREATENUMERICMATRIX800
#define mxcreatenumericmatrix mxcreatenumericmatrix800
#define mxCalcSingleSubscript mxCalcSingleSubscript800
#define MXCALCSINGLESUBSCRIPT MXCALCSINGLESUBSCRIPT800
#define mxcalcsinglesubscript mxcalcsinglesubscript800
#define mxCreateCharMatrixFromStrings mxCreateCharMatrixFromStrings800
#define MXCREATECHARMATRIXFROMSTRINGS MXCREATECHARMATRIXFROMSTRINGS800
#define mxcreatecharmatrixfromstrings mxcreatecharmatrixfromstrings800
#define mxCreateSparse mxCreateSparse800
#define MXCREATESPARSE MXCREATESPARSE800
#define mxcreatesparse mxcreatesparse800
#define mxCreateStructArray mxCreateStructArray800
#define MXCREATESTRUCTARRAY MXCREATESTRUCTARRAY800
#define mxcreatestructarray mxcreatestructarray800
#define mxCreateStructMatrix mxCreateStructMatrix800
#define MXCREATESTRUCTMATRIX MXCREATESTRUCTMATRIX800
#define mxcreatestructmatrix mxcreatestructmatrix800
#define mxGetCell mxGetCell800
#define MXGETCELL MXGETCELL800
#define mxgetcell mxgetcell800
#define mxGetDimensions mxGetDimensions800
#define MXGETDIMENSIONS MXGETDIMENSIONS800
#define mxgetdimensions mxgetdimensions800
#define mxGetField mxGetField800
#define MXGETFIELD MXGETFIELD800
#define mxgetfield mxgetfield800
#define mxGetFieldByNumber mxGetFieldByNumber800
#define MXGETFIELDBYNUMBER MXGETFIELDBYNUMBER800
#define mxgetfieldbynumber mxgetfieldbynumber800
#define mxGetProperty mxGetProperty800
#define MXGETPROPERTY MXGETPROPERTY800
#define mxgetproperty mxgetproperty800
#define mxGetIr mxGetIr800
#define MXGETIR MXGETIR800
#define mxgetir mxgetir800
#define mxGetJc mxGetJc800
#define MXGETJC MXGETJC800
#define mxgetjc mxgetjc800
#define mxGetNumberOfDimensions mxGetNumberOfDimensions800
#define MXGETNUMBEROFDIMENSIONS MXGETNUMBEROFDIMENSIONS800
#define mxgetnumberofdimensions mxgetnumberofdimensions800
#define mxGetNzmax mxGetNzmax800
#define MXGETNZMAX MXGETNZMAX800
#define mxgetnzmax mxgetnzmax800
#define mxGetString mxGetString800
#define MXGETSTRING MXGETSTRING800
#define mxgetstring mxgetstring800
#define mxSetCell mxSetCell800
#define MXSETCELL MXSETCELL800
#define mxsetcell mxsetcell800
#define mxSetField mxSetField800
#define MXSETFIELD MXSETFIELD800
#define mxsetfield mxsetfield800
#define mxSetFieldByNumber mxSetFieldByNumber800
#define MXSETFIELDBYNUMBER MXSETFIELDBYNUMBER800
#define mxsetfieldbynumber mxsetfieldbynumber800
#define mxSetProperty mxSetProperty800
#define MXSETPROPERTY MXSETPROPERTY800
#define mxsetproperty mxsetproperty800
#define mxSetIr mxSetIr800
#define MXSETIR MXSETIR800
#define mxsetir mxsetir800
#define mxSetJc mxSetJc800
#define MXSETJC MXSETJC800
#define mxsetjc mxsetjc800
#define mxGetM mxgetm800
#define MXGETM MXGETM800
#define mxgetm mxgetm800
#define mxGetNumberOfElements mxGetNumberOfElements800
#define MXGETNUMBEROFELEMENTS MXGETNUMBEROFELEMENTS800
#define mxgetnumberofelements mxgetnumberofelements800
#define mxMalloc mxMalloc800
#define MXMALLOC MXMALLOC800
#define mxmalloc mxmalloc800
#define mxCalloc mxCalloc800
#define MXCALLOC MXCALLOC800
#define mxcalloc mxcalloc800
#define mxRealloc mxRealloc800
#define MXREALLOC MXREALLOC800
#define mxrealloc mxrealloc800
#define mxCreateCharMatrixFromStrs mxCreateCharMatrixFromStrs800
#define MXCREATECHARMATRIXFROMSTRS MXCREATECHARMATRIXFROMSTRS800
#define mxcreatecharmatrixfromstrs mxcreatecharmatrixfromstrs800
#define mxCopyReal4ToPtr mxCopyReal4ToPtr800
#define MXCOPYREAL4TOPTR MXCOPYREAL4TOPTR800
#define mxcopyreal4toptr mxcopyreal4toptr800
#define mxCopyPtrToReal4 mxCopyPtrToReal4800
#define MXCOPYPTRTOREAL4 MXCOPYPTRTOREAL4800
#define mxcopyptrtoreal4 mxcopyptrtoreal4800
#define mxCopyReal8ToPtr mxCopyReal8ToPtr800
#define MXCOPYREAL8TOPTR MXCOPYREAL8TOPTR800
#define mxcopyreal8toptr mxcopyreal8toptr800
#define mxCopyPtrToReal8 mxCopyPtrToReal8800
#define MXCOPYPTRTOREAL8 MXCOPYPTRTOREAL8800
#define mxcopyptrtoreal8 mxcopyptrtoreal8800
#define mxCopyCharacterToPtr mxCopyCharacterToPtr800
#define MXCOPYCHARACTERTOPTR MXCOPYCHARACTERTOPTR800
#define mxcopycharactertoptr mxcopycharactertoptr800
#define mxCopyPtrToCharacter mxCopyPtrToCharacter800
#define MXCOPYPTRTOCHARACTER MXCOPYPTRTOCHARACTER800
#define mxcopyptrtocharacter mxcopyptrtocharacter800
#define mxCopyInteger1ToPtr mxCopyInteger1ToPtr800
#define MXCOPYINTEGER1TOPTR MXCOPYINTEGER1TOPTR800
#define mxcopyinteger1toptr mxcopyinteger1toptr800
#define mxCopyPtrToInteger1 mxCopyPtrToInteger1800
#define MXCOPYPTRTOINTEGER1 MXCOPYPTRTOINTEGER1800
#define mxcopyptrtointeger1 mxcopyptrtointeger1800
#define mxCopyInteger2ToPtr mxCopyInteger2ToPtr800
#define MXCOPYINTEGER2TOPTR MXCOPYINTEGER2TOPTR800
#define mxcopyinteger2toptr mxcopyinteger2toptr800
#define mxCopyPtrToInteger2 mxCopyPtrToInteger2800
#define MXCOPYPTRTOINTEGER2 MXCOPYPTRTOINTEGER2800
#define mxcopyptrtointeger2 mxcopyptrtointeger2800
#define mxCopyInteger4ToPtr mxCopyInteger4ToPtr800
#define MXCOPYINTEGER4TOPTR MXCOPYINTEGER4TOPTR800
#define mxcopyinteger4toptr mxcopyinteger4toptr800
#define mxCopyPtrToInteger4 mxCopyPtrToInteger4800
#define MXCOPYPTRTOINTEGER4 MXCOPYPTRTOINTEGER4800
#define mxcopyptrtointeger4 mxcopyptrtointeger4800
#define mxCopyInteger8ToPtr mxCopyInteger8ToPtr800
#define MXCOPYINTEGER8TOPTR MXCOPYINTEGER8TOPTR800
#define mxcopyinteger8toptr mxcopyinteger8toptr800
#define mxCopyPtrToInteger8 mxCopyPtrToInteger8800
#define MXCOPYPTRTOINTEGER8 MXCOPYPTRTOINTEGER8800
#define mxcopyptrtointeger8 mxcopyptrtointeger8800
#define mxCopyPtrToPtrArray mxCopyPtrToPtrArray800
#define MXCOPYPTRTOPTRARRAY MXCOPYPTRTOPTRARRAY800
#define mxcopyptrtoptrarray mxcopyptrtoptrarray800
#define mxCopyMWIndexToPtr mxCopyMWIndexToPtr800
#define MXCOPYMWINDEXTOPTR MXCOPYMWINDEXTOPTR800
#define mxcopymwindextoptr mxcopymwindextoptr800
#define mxCopyPtrToMWIndex mxCopyPtrToMWIndex800
#define MXCOPYPTRTOMWINDEX MXCOPYPTRTOMWINDEX800
#define mxcopyptrtomwindex mxcopyptrtomwindex800
#define mxCopyPtrToComplex16 mxCopyPtrToComplex16800
#define MXCOPYPTRTOCOMPLEX16 MXCOPYPTRTOCOMPLEX16800
#define mxcopyptrtocomplex16 mxcopyptrtocomplex16800
#define mxCopyComplex16ToPtr mxCopyComplex16ToPtr800
#define MXCOPYCOMPLEX16TOPTR MXCOPYCOMPLEX16TOPTR800
#define mxcopycomplex16toptr mxcopycomplex16toptr800
#define mxCopyPtrToComplex8 mxCopyPtrToComplex8800
#define MXCOPYPTRTOCOMPLEX8 MXCOPYPTRTOCOMPLEX8800
#define mxcopyptrtocomplex8 mxcopyptrtocomplex8800
#define mxCopyComplex8ToPtr mxCopyComplex8ToPtr800
#define MXCOPYCOMPLEX8TOPTR MXCOPYCOMPLEX8TOPTR800
#define mxcopycomplex8toptr mxcopycomplex8toptr800
#define mxGetElementSize mxGetElementSize800
#define MXGETELEMENTSIZE MXGETELEMENTSIZE800
#define mxgetelementsize mxgetelementsize800
#define mxSetM mxSetM800
#define MXSETM MXSETM800
#define mxsetm mxsetm800
#define mxSetN mxSetN800
#define MXSETN MXSETN800
#define mxsetn mxsetn800
#define mxSetDimensions mxSetDimensions800
#define MXSETDIMENSIONS MXSETDIMENSIONS800
#define mxsetdimensions mxsetdimensions800
#define mxSetNzmax mxSetNzmax800
#define MXSETNZMAX MXSETNZMAX800
#define mxsetnzmax mxsetnzmax800
#define mxGetPr mxGetPr800
#define MXGETPR MXGETPR800
#define mxgetpr mxgetpr800
#define mxSetPr mxSetPr800
#define MXSETPR MXSETPR800
#define mxsetpr mxsetpr800
#define mxGetData mxGetData800
#define MXGETDATA MXGETDATA800
#define mxgetdata mxgetdata800
#define mxSetData mxSetData800
#define MXSETDATA MXSETDATA800
#define mxsetdata mxsetdata800
#define mxGetScalar mxGetScalar800
#define MXGETSCALAR MXGETSCALAR800
#define mxgetscalar mxgetscalar800
#define mxDuplicateArray mxDuplicateArray800
#define MXDUPLICATEARRAY MXDUPLICATEARRAY800
#define mxduplicatearray mxduplicatearray800
#define mxIsComplex mxIsComplex800
#define MXISCOMPLEX MXISCOMPLEX800
#define mxiscomplex mxiscomplex800
#define mxMakeArrayReal mxMakeArrayReal800
#define MXMAKEARRAYREAL MXMAKEARRAYREAL800
#define mxmakearrayreal mxmakearrayreal800
#define mxMakeArrayComplex mxMakeArrayComplex800
#define MXMAKEARRAYCOMPLEX MXMAKEARRAYCOMPLEX800
#define mxmakearraycomplex mxmakearraycomplex800
#define mxGetDoubles mxGetDoubles800
#define MXGETDOUBLES MXGETDOUBLES800
#define mxgetdoubles mxgetdoubles800
#define mxSetDoubles mxSetDoubles800
#define MXSETDOUBLES MXSETDOUBLES800
#define mxsetdoubles mxsetdoubles800
#define mxGetComplexDoubles mxGetComplexDoubles800
#define MXGETCOMPLEXDOUBLES MXGETCOMPLEXDOUBLES800
#define mxgetcomplexdoubles mxgetcomplexdoubles800
#define mxSetComplexDoubles mxSetComplexDoubles800
#define MXSETCOMPLEXDOUBLES MXSETCOMPLEXDOUBLES800
#define mxsetcomplexdoubles mxsetcomplexdoubles800
#define mxGetSingles mxGetSingles800
#define MXGETSINGLES MXGETSINGLES800
#define mxgetsingles mxgetsingles800
#define mxSetSingles mxSetSingles800
#define MXSETSINGLES MXSETSINGLES800
#define mxsetsingles mxsetsingles800
#define mxGetComplexSingles mxGetComplexSingles800
#define MXGETCOMPLEXSINGLES MXGETCOMPLEXSINGLES800
#define mxgetcomplexsingles mxgetcomplexsingles800
#define mxSetComplexSingles mxSetComplexSingles800
#define MXSETCOMPLEXSINGLES MXSETCOMPLEXSINGLES800
#define mxsetcomplexsingles mxsetcomplexsingles800
#define mxGetInt8s mxGetInt8s800
#define MXGETINT8S MXGETINT8S800
#define mxgetint8s mxgetint8s800
#define mxSetInt8s mxSetInt8s800
#define MXSETINT8S MXSETINT8S800
#define mxsetint8s mxsetint8s800
#define mxGetUint8s mxGetUint8s800
#define MXGETUINT8S MXGETUINT8S800
#define mxgetuint8s mxgetuint8s800
#define mxSetUint8s mxSetUint8s800
#define MXSETUINT8S MXSETUINT8S800
#define mxsetuint8s mxsetuint8s800
#define mxGetInt16s mxGetInt16s800
#define MXGETINT16S MXGETINT16S800
#define mxgetint16s mxgetint16s800
#define mxSetInt16s mxSetInt16s800
#define MXSETINT16S MXSETINT16S800
#define mxsetint16s mxsetint16s800
#define mxGetUint16s mxGetUint16s800
#define MXGETUINT16S MXGETUINT16S800
#define mxgetuint16s mxgetuint16s800
#define mxSetUint16s mxSetUint16s800
#define MXSETUINT16S MXSETUINT16S800
#define mxsetuint16s mxsetuint16s800
#define mxGetInt32s mxGetInt32s800
#define MXGETINT32S MXGETINT32S800
#define mxgetint32s mxgetint32s800
#define mxSetInt32s mxSetInt32s800
#define MXSETINT32S MXSETINT32S800
#define mxsetint32s mxsetint32s800
#define mxGetUint32s mxGetUint32s800
#define MXGETUINT32S MXGETUINT32S800
#define mxgetuint32s mxgetuint32s800
#define mxSetUint32s mxSetUint32s800
#define MXSETUINT32S MXSETUINT32S800
#define mxsetuint32s mxsetuint32s800
#define mxGetInt64s mxGetInt64s800
#define MXGETINT64S MXGETINT64S800
#define mxgetint64s mxgetint64s800
#define mxSetInt64s mxSetInt64s800
#define MXSETINT64S MXSETINT64S800
#define mxsetint64s mxsetint64s800
#define mxGetUint64s mxGetUint64s800
#define MXGETUINT64S MXGETUINT64S800
#define mxgetuint64s mxgetuint64s800
#define mxSetUint64s mxSetUint64s800
#define MXSETUINT64S MXSETUINT64S800
#define mxsetuint64s mxsetuint64s800

#if defined(WITH_COMMENTS)
/*
 * MEX module APIs
 */
#endif
#define mexPrintf mexPrintf800
#define mexprintf mexprintf800
#define MEXPRINTF MEXPRINTF800
#define mexErrMsgIdAndTxt mexErrMsgIdAndTxt800
#define mexerrmsgidandtxt mexerrmsgidandtxt800
#define MEXERRMSGIDANDTXT MEXERRMSGIDANDTXT800
#define mexWarnMsgIdAndTxt mexWarnMsgIdAndTxt800
#define mexwarnmsgidandtxt mexwarnmsgidandtxt800
#define MEXWARNMSGIDANDTXT MEXWARNMSGIDANDTXT800
#define mexErrMsgTxt mexErrMsgTxt800
#define mexerrmsgtxt mexerrmsgtxt800
#define MEXERRMSGTXT MEXERRMSGTXT800
#define mexWarnMsgTxt mexWarnMsgTxt800
#define mexwarnmsgtxt mexwarnmsgtxt800
#define MEXWARNMSGTXT MEXWARNMSGTXT800
#define mexIsLocked mexIsLocked800
#define mexislocked mexislocked800
#define MEXISLOCKED MEXISLOCKED800
#define mexLock mexLock800
#define mexlock mexlock800
#define MEXLOCK MEXLOCK800
#define mexUnlock mexUnlock800
#define mexunlock mexunlock800
#define MEXUNLOCK MEXUNLOCK800
#define mexMakeArrayPersistent mexMakeArrayPersistent800
#define mexmakearraypersistent mexmakearraypersistent800
#define MEXMAKEARRAYPERSISTENT MEXMAKEARRAYPERSISTENT800
#define mexMakeMemoryPersistent mexMakeMemoryPersistent800
#define mexmakememorypersistent mexmakememorypersistent800
#define MEXMAKEMEMORYPERSISTENT MEXMAKEMEMORYPERSISTENT800
#define mexIsGlobal mexIsGlobal800
#define mexisglobal mexisglobal800
#define MEXISGLOBAL MEXISGLOBAL800
#define mexFunctionName mexFunctionName800
#define mexfunctionname mexfunctionname800
#define MEXFUNCTIONNAME MEXFUNCTIONNAME800
#define mexAtExit mexAtExit800
#define mexatexit mexatexit800
#define MEXATEXIT MEXATEXIT800

#define mexCallMATLAB mexCallMATLAB800
#define mexcallmatlab mexcallmatlab800
#define MEXCALLMATLAB MEXCALLMATLAB800
#define mexCallMATLABWithTrap mexCallMATLABWithTrap800
#define mexcallmatlabwithtrap mexcallmatlabwithtrap800
#define MEXCALLMATLABWITHTRAP MEXCALLMATLABWITHTRAP800

#define mexEvalString mexEvalString800
#define mexevalstring mexevalstring800
#define MEXEVALSTRING MEXEVALSTRING800
#define mexEvalStringWithTrap mexEvalStringWithTrap800
#define mexevalstringwithtrap mexevalstringwithtrap800
#define MEXEVALSTRINGWITHTRAP MEXEVALSTRINGWITHTRAP800
#define mexGetVariable mexGetVariable800
#define mexgetvariable mexgetvariable800
#define MEXGETVARIABLE MEXGETVARIABLE800
#define mexGetVariablePtr mexGetVariablePtr800
#define mexgetvariableptr mexgetvariableptr800
#define MEXGETVARIABLEPTR MEXGETVARIABLEPTR800
#define mexPutVariable mexPutVariable800
#define mexputvariable mexputvariable800
#define MEXPUTVARIABLE MEXPUTVARIABLE800

#define mexGet mexGetIsDeprecated
#define mexget mexgetisdeprecated
#define MEXGET MEXGETISDEPRECATED
#define mexSet mexSetIsDeprecated
#define mexset mexsetisdeprecated
#define MEXSET MEXSETISDEPRECATED
#define mexSetTrapFlag mexSetTrapFlagIsDeprecated
#define mexsettrapflag mexsettrapflagisdeprecated
#define MEXSETTRAPFLAG MEXSETTRAPFLAGISDEPRECATED

#if defined(WITH_COMMENTS)
/*
 * MAT module APIs
 */
#endif
#define matOpen matOpen800
#define matopen matopen800
#define MATOPEN MATOPEN800
#define matClose matClose800
#define matclose matclose800
#define MATCLOSE MATCLOSE800
#define matGetVariable matGetVariable800
#define matGetVariable matGetVariable800
#define MATGETVARIABLE MATGETVARIABLE800
#define matGetNextVariable matGetNextVariable800
#define matGetNextVariable matGetNextVariable800
#define MATGETNEXTVARIABLE MATGETNEXTVARIABLE800
#define matGetVariableInfo matGetVariableInfo800
#define matGetVariableInfo matGetVariableInfo800
#define MATGETVARIABLEINFO MATGETVARIABLEINFO800
#define matGetDir matGetDir800
#define matGetDir matGetDir800
#define MATGETDIR MATGETDIR800
#define matGetErrno matGetErrno800
#define matGetErrno matGetErrno800
#define MATGETERRNO MATGETERRNO800
#define matPutVariable matPutVariable800
#define matPutVariable matPutVariable800
#define MATPUTVARIABLE MATPUTVARIABLE800
#define matPutVariableAsGlobal matPutVariableAsGlobal800
#define matPutVariableAsGlobal matPutVariableAsGlobal800
#define MATPUTVARIABLEASGLOBAL MATPUTVARIABLEASGLOBAL800

#elif defined(MX_COMPAT_32)

#if defined(WITH_COMMENTS)
/*
 * Compatibility layer for MEX files using the 32-bit mxArray APIs
 */
#endif

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

#define mexGet mexGetIsDeprecated
#define mexget mexgetisdeprecated
#define MEXGET MEXGETISDEPRECATED
#define mexSet mexSetIsDeprecated
#define mexset mexsetisdeprecated
#define MEXSET MEXSETISDEPRECATED
#define mexSetTrapFlag mexSetTrapFlagIsDeprecated
#define mexsettrapflag mexsettrapflagisdeprecated
#define MEXSETTRAPFLAG MEXSETTRAPFLAGISDEPRECATED

#endif
/* Current MATRIX published API version */
#define MX_CURRENT_API_VER 0x08000000
#define FORT_MX_CURRENT_API_VER z'08000000'

/* Backward compatible MATRIX published API versions */
#define MX_LAST_32BIT_VER 0x07000000
#define MX_LAST_SEPARATE_COMPLEX_VER 0x07300000
#define FORT_MX_LAST_32BIT_VER z'07000000'
#define FORT_MX_LAST_SEPARATE_COMPLEX_VER z'07300000'

/* Required MEX-file MATRIX published API version */
#if TARGET_API_VERSION == 700
#if defined(MX_COMPAT_32)
#define MX_TARGET_API_VER MX_LAST_32BIT_VER
#define FORT_MX_TARGET_API_VER FORT_MX_LAST_32BIT_VER
#else
#define MX_TARGET_API_VER MX_LAST_SEPARATE_COMPLEX_VER
#define FORT_MX_TARGET_API_VER FORT_MX_LAST_SEPARATE_COMPLEX_VER
#endif
#else
#define MX_TARGET_API_VER MX_CURRENT_API_VER
#define FORT_MX_TARGET_API_VER FORT_MX_CURRENT_API_VER
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

#ifdef __GFORTRAN__
#define MX_HAS_64BIT_ARRAY_DIMS MX_TARGET_API_VER > MX_LAST_32BIT_VER
#define MX_HAS_INTERLEAVED_COMPLEX MX_TARGET_API_VER > MX_LAST_SEPARATE_COMPLEX_VER
#else
#define MX_HAS_64BIT_ARRAY_DIMS FORT_MX_TARGET_API_VER > FORT_MX_LAST_32BIT_VER
#define MX_HAS_INTERLEAVED_COMPLEX FORT_MX_TARGET_API_VER > FORT_MX_LAST_SEPARATE_COMPLEX_VER
#endif

