#ifdef CHEMAPP
/* ----------------------------------------------------------------------
 *  System      : ChemApp
 * ----------------------------------------------------------------------
 *  Module      : cacint.h (ChemApp C/C++ interface)
 *
 * ----------------------------------------------------------------------
 *  File        : $RCSfile: cacint.h,v $ 	
 *  Revision    : $Revision: 1.28 $	
 *  Last Change : $Date: 2006/06/02 09:21:45 $ by $Author: sp $  
 *
 *  Language    : C
 * ----------------------------------------------------------------------
 *  Subject     : This file contains the C/C++ interface to ChemApp 
 *                (header file)
 * ----------------------------------------------------------------------
 */

/*
This file is Copyright (C) GTT-Technologies, Herzogenrath, Germany.
It may only be used together with GTT-Technologies´ ChemApp software.
Unauthorised duplication and distribution both in printed and online
form, also in parts, is prohibited.
*/

/* ChemApp DLL function prototypes */

#ifndef _cacint
#define _cacint

#ifdef UNIX
#define LI  long            	
#define LIP long*       	
#define LNT long        	
#define DB  double		
#define DBP double*		
#define CHP char*		
#define CMT extern int		
#define VDP void*		
#define ftnlen long		/* FORTRAN string length type */
#else
#define LI  long         	/* unsigned int		*/
#define LIP long*       	/* unsigned int*	*/
#define LNT long        	/* unsigned int		*/
#define DB  double		/* double		*/
#define DBP double*		/* double*		*/
#define CHP char*		/* char*		*/
#define CMT void __stdcall	/* void __stdcall	*/
#define VDP void*		/* void*		*/
#endif
                                        
/* Length of a TQ String */
#define TQSTRLEN 25
/* Macro for defining TQStrings */
#define TQSTRING(x) char x[TQSTRLEN]
/* TQ Error Message Buffer */
extern char TQERRMSG[3][80]; 

#ifdef __cplusplus
extern "C" {
#endif
int tqini (LIP NOERR);
int tqopen(CHP FILE, LI LUN, LIP NOERR);
int tqclos(LI LUN, LIP NOERR);
int tqgio (CHP OPTION, LIP IVAL, LIP NOERR);
int tqcio (CHP OPTION, LI IVAL, LIP NOERR);
int tqrfil(LIP NOERR);
int tqgsu (CHP OPTION, CHP UNIT, LIP NOERR);
int tqcsu (CHP OPTION, CHP UNIT, LIP NOERR);
int tqinsc(CHP NAME, LIP INDEXS, LIP NOERR);
int tqgnsc(LI INDEXS, CHP NAME, LIP NOERR);
int tqnosc(LIP NSCOM, LIP NOERR);
int tqstsc(LI INDEXS, DBP STOI, DBP WMASS, LIP NOERR);
int tqcsc (CHP NAME, LIP NOERR);
int tqinp (CHP NAME, LIP INDEXP, LIP NOERR);
int tqgnp (LI INDEXP, CHP NAME, LIP NOERR);
int tqnop (LIP NPHASE, LIP NOERR);
int tqinpc(CHP NAME, LI INDEXP, LIP INDEXC, LIP NOERR);
int tqgnpc(LI INDEXP, LI INDEXC, CHP NAME, LIP NOERR);
int tqnopc(LI INDEXP, LIP NPCON, LIP NOERR);
int tqstpc(LI INDEXP, LI INDEXC, DBP STOI, DBP WMASS, LIP NOERR);
int tqgsp (LI INDEXP, CHP OPTION, LIP NOERR);
int tqcsp (LI INDEXP, CHP OPTION, LIP NOERR);
int tqgspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR);
int tqcspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR);
int tqsetc(CHP OPTION, LI INDEXP, LI INDEX, DB VAL, LIP NUMCON, LIP NOERR);
int tqremc(LI NUMCON, LIP NOERR);
int tqsttp(CHP IDENTS, DBP VALS, LIP NOERR);
int tqstca(CHP IDENTS, LI INDEXP, LI INDEXC, DB VAL, LIP NOERR);
int tqstec(CHP OPTION, LI INDEXP, DB VAL, LIP NOERR);
int tqstrm(CHP IDENTS, LIP NOERR);
int tqce  (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
int tqcel (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
int tqclim(CHP OPTION, DB VAL, LIP NOERR);
int tqgetr(CHP OPTION, LI INDEXP, LI INDEX, DBP VAL, LIP NOERR);
int tqgdpc(CHP OPTION, LI INDEXP, LI INDEXC,DBP VAL, LIP NOERR);
int tqshow(LIP NOERR);
int tqerr (CHP MESS, LIP NOERR);

int tqcprt(LIP NOERR);
int tqvers(LIP NVERS, LIP NOERR);
int tqsize(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,LIP NJ,LIP NK,LIP NOERR);
int tqmodl(LI INDEXP, CHP NAME, LIP NOERR);
int tqstxp(CHP IDENTS,CHP OPTION, DBP VAL, LIP NOERR);
int tqlite(LIP LITE, LIP NOERR);
int tqrbin(LIP NOERR);
int tqmap(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT, LIP NOERR);
int tqmapl(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT, LIP NOERR);
int tqpcis(LI INDEXP, LI INDEXC, LIP ISPERM, LIP NOERR);
int tqopna(CHP FILE, LI LUN, LIP NOERR);
int tqopnb(CHP FILE, LI LUN, LIP NOERR);
int tqnosl(LI INDEXP, LIP NSUBL, LIP NOERR);
int tqnolc(LI INDEXP, LI INDEXL, LIP NSLCON, LIP NOERR);
int tqinlc(CHP NAME, LI INDEXP, LI INDEXL, LIP INDEXC, LIP NOERR);
int tqgnlc(LI INDEXP, LI INDEXL, LI INDEXC, CHP NAME, LIP NOERR);
int tqgtlc(LI INDEXP, LI INDEXL, LI INDEXC, DBP VAL, LIP NOERR);
/*
int tqgopn (CHP FILE,LI LUN,CHP FFORM,CHP FSTAT,CHP FACC,LI RECL,
	    LIP IOSTAT,LIP NOERR);
*/
int tqbond(LI INDEXP, LI INDEXA, LI INDEXB, LI INDEXC, LI INDEXD, 
	   DBP VAL, LIP NOERR);
int tqused(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,
	   LIP NJ,LIP NK,LIP NOERR);
int tqgtrh(LIP TFHVER,
	   CHP TFHNWP,
	   LIP TFHVNW,
	   CHP TFHNRP,
	   LIP TFHVNR,
	   LIP TFHDTC,
	   LIP TFHDTE,
	   CHP TFHID,
	   CHP TFHUSR,
	   CHP TFHREM,
	   LIP NOERR);
int tqopnt(CHP FILE, LI LUN, LIP NOERR);
int tqrcst(LIP NOERR);
int tqgtid(CHP ID, LIP NOERR);
int tqgtnm(CHP NAME, LIP NOERR);
int tqgtpi(CHP PID, LIP NOERR);
int tqwstr(CHP OPTION, CHP CTXT, LIP NOERR);
int tqgted(LIP EDMON, LIP EDYEAR, LIP NOERR);
int tqgthi(CHP HASPT, LIP HASPID, LIP NOERR);
int tqcen (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
int tqcenl(CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR);
int tqwasc(CHP FILE, LIP NOERR);
int tqcdat(LI I1, LI I2, LI I3, LI I4, LI I5, DB VAL, LIP NOERR);
int tqchar(LI INDEXP, LI INDEXC, DBP VAL, LIP NOERR);
int tqcnsc (LI INDEXS, CHP NAME, LIP NOERR);

#ifdef __cplusplus
};
#endif          /* __cplusplus        */
#endif 		/* _cacint */
#endif