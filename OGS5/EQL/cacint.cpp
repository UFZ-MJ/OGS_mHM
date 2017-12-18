 #ifdef CHEMAPP
/* ----------------------------------------------------------------------
 *  System      : ChemApp
 * ----------------------------------------------------------------------
 *  Module      : cacint.c (ChemApp C/C++ interface)
 *
 * ----------------------------------------------------------------------
 *  File        : $RCSfile: cacint.c,v $ 	
 *  Revision    : $Revision: 1.32 $	
 *  Last Change : $Date: 2006/06/02 09:21:53 $ by $Author: sp $  
 *
 *  Language    : C
 * ----------------------------------------------------------------------
 *  Subject     : This file contains the C/C++ interface to ChemApp
 * ----------------------------------------------------------------------
 */

/*
This file is Copyright (C) GTT-Technologies, Herzogenrath, Germany.
It may only be used together with GTT-Technologies´ ChemApp software.
Unauthorised duplication and distribution both in printed and online
form, also in parts, is prohibited.
*/
// #include "stdafx.h" /* MFC */

#include <string.h>
#include "cacint.h"

#ifdef __cplusplus
extern "C" {
#endif

static char *cacint_RevisionID="$Id: cacint.c,v 1.32 2006/06/02 09:21:53 sp Exp $";


#ifdef UNIX

CMT tqini_ (LIP NOERR);
CMT tqopen_(CHP FILE,LIP LUN, LIP NOERR,ftnlen);
CMT tqclos_(LIP LUN, LIP NOERR);
CMT tqgio_ (CHP OPTION, LIP IVAL, LIP NOERR,ftnlen);
CMT tqcio_ (CHP OPTION, LIP IVAL, LIP NOERR,ftnlen);
CMT tqrfil_(LIP NOERR);
CMT tqgsu_ (CHP OPTION, CHP UNIT, LIP NOERR,ftnlen,ftnlen);
CMT tqcsu_ (CHP OPTION, CHP UNIT, LIP NOERR,ftnlen,ftnlen);
CMT tqinsc_(CHP NAME, LIP INDEXS, LIP NOERR,ftnlen);
CMT tqgnsc_(LIP INDEXS, CHP NAME, LIP NOERR,ftnlen);
CMT tqnosc_(LIP NSCOM, LIP NOERR);
CMT tqstsc_(LIP INDEXS, DBP STOI, DBP WMASS, LIP NOERR);
CMT tqcsc_ (CHP NAME, LIP NOERR,ftnlen);
CMT tqinp_ (CHP NAME, LIP INDEXP, LIP NOERR,ftnlen);
CMT tqgnp_ (LIP INDEXP, CHP NAME, LIP NOERR,ftnlen);
CMT tqnop_ (LIP NPHASE, LIP NOERR);
CMT tqinpc_(CHP NAME, LIP INDEXP, LIP INDEXC, LIP NOERR,ftnlen);
CMT tqgnpc_(LIP INDEXP, LIP INDEXC, CHP NAME, LIP NOERR,ftnlen);
CMT tqnopc_(LIP INDEXP, LIP NPCON, LIP NOERR);
CMT tqstpc_(LIP INDEXP, LIP INDEXC, DBP STOI, DBP WMASS, LIP NOERR);
CMT tqgsp_ (LIP INDEXP, CHP OPTION, LIP NOERR,ftnlen);
CMT tqcsp_ (LIP INDEXP, CHP OPTION, LIP NOERR,ftnlen);
CMT tqgspc_(LIP INDEXP, LIP INDEXC, CHP OPTION, LIP NOERR,ftnlen);
CMT tqcspc_(LIP INDEXP, LIP INDEXC, CHP OPTION, LIP NOERR,ftnlen);
CMT tqsetc_(CHP OPTION, LIP INDEXP, LIP INDEX, DBP VAL, LIP NUMCON, LIP NOERR,ftnlen);
CMT tqremc_(LIP NUMCON, LIP NOERR);
CMT tqsttp_(CHP IDENTS, DBP VALS, LIP NOERR,ftnlen);
CMT tqstca_(CHP IDENTS, LIP INDEXP, LIP INDEXC, DBP VAL, LIP NOERR,ftnlen);
CMT tqstec_(CHP OPTION, LIP INDEXP, DBP VAL, LIP NOERR,ftnlen);
CMT tqstrm_(CHP IDENTS, LIP NOERR,ftnlen);
CMT tqce_  (CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR,ftnlen);
CMT tqcel_ (CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR,ftnlen);
CMT tqclim_(CHP OPTION, DBP VAL, LIP NOERR,ftnlen);
CMT tqgetr_(CHP OPTION, LIP INDEXP, LIP INDEX, DBP VAL, LIP NOERR,ftnlen);
CMT tqgdpc_(CHP OPTION, LIP INDEXP, LIP INDEXC,DBP VAL, LIP NOERR,ftnlen);
CMT tqshow_(LIP NOERR);
CMT tqerr_ (CHP MESS, LIP NOERR,ftnlen);

/* Added for ChemApp V113 */

CMT tqcprt_(LIP NOERR);
CMT tqvers_(LIP NVERS, LIP NOERR);
CMT tqsize_(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,LIP NJ,LIP NK,LIP NOERR);
CMT tqmodl_(LIP INDEXP, CHP NAME, LIP NOERR,ftnlen);
CMT tqstxp_(CHP IDENTS,CHP OPTION, DBP VAL, LIP NOERR,ftnlen,ftnlen);

/* Added for ChemApp V201 */

CMT tqlite_(LIP LITE, LIP NOERR);

/* Added for ChemApp V214 */

CMT tqrbin_(LIP NOERR);

/* Added for ChemApp V300 */

CMT tqmap_  (CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP ICONT, LIP NOERR,ftnlen);
CMT tqmapl_ (CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP ICONT, LIP NOERR,ftnlen);

/* Added for ChemApp V310 */

CMT tqpcis_(LIP INDEXP, LIP INDEXC, LIP ISPERM, LIP NOERR);

/* Added for ChemApp V330 */

CMT tqopna_(CHP FILE,LIP LUN, LIP NOERR,ftnlen);
CMT tqopnb_(CHP FILE,LIP LUN, LIP NOERR,ftnlen);

/* Added for ChemApp V331 */

CMT tqnosl_(LIP INDEXP, LIP NSUBL, LIP NOERR);
CMT tqnolc_(LIP INDEXP, LIP INDEXL, LIP NSLCON, LIP NOERR);
CMT tqinlc_(CHP NAME, LIP INDEXP, LIP INDEXL, LIP INDEXC, LIP NOERR,ftnlen);
CMT tqgnlc_(LIP INDEXP, LIP INDEXL, LIP INDEXC, CHP NAME, LIP NOERR,ftnlen);
CMT tqgtlc_(LIP INDEXP, LIP INDEXL, LIP INDEXC, DBP VAL, LIP NOERR);

/* Added for ChemApp V339 */
/*
CMT tqgopn_(CHP FILE,LIP LUN,CHP FFORM,CHP FSTAT,CHP FACC,LIP RECL,
	    LIP IOSTAT,LIP NOERR,ftnlen,ftnlen,ftnlen,ftnlen);
*/

CMT tqbond_(LIP INDEXP, LIP INDEXA, LIP INDEXB, LIP INDEXC, LIP INDEXD, 
	    DBP VAL, LIP NOERR);


/* Added for ChemApp V340 */
CMT tqused_(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,
	    LIP NJ,LIP NK,LIP NOERR);
CMT tqgtrh_(LIP TFHVER,
	    CHP TFHNWP,
	    LIP TFHVNW,
	    CHP TFHNRP,
	    LIP TFHVNR,
	    LIP TFHDTC,
	    LIP TFHDTE,
	    CHP TFHID,
	    CHP TFHUSR,
	    CHP TFHREM,
	    LIP NOERR,
	    ftnlen,
	    ftnlen,
	    ftnlen,
	    ftnlen,
	    ftnlen);
CMT tqopnt_(CHP FILE, LIP LUN, LIP NOERR, ftnlen);
CMT tqrcst_(LIP NOERR);
CMT tqgtid_(CHP ID, LIP NOERR, ftnlen);
CMT tqgtnm_(CHP NAME, LIP NOERR, ftnlen);
CMT tqgtpi_(CHP PID, LIP NOERR, ftnlen);
CMT tqwstr_(CHP OPTION, CHP CTXT, LIP NOERR,ftnlen,ftnlen);

/* Added for ChemApp V414 */

CMT tqgted_(LIP EDMON, LIP EDYEAR, LIP NOERR);
CMT tqgthi_(CHP HASPT, LIP HASPID, LIP NOERR,ftnlen);
CMT tqcen_ (CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR,ftnlen);
CMT tqcenl_(CHP OPTION, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR,ftnlen);

/* Added for ChemApp V540 */

CMT tqgdat_(LIP INDEXP, LIP INDEXC, CHP OPTION, LIP INDEXR, LIP NVALV, DBP VALV, LIP NOERR,ftnlen);


/*
********************************
Still missing:
TQLPAR
TQGPAR
********************************
*/
/*CMT tqlpar_(LIP INDEXP, CHP OPTION, LIP NOPAR, CHP CHRPAR, LIP NOERR,ftnlen,ftnlen);*/
CMT tqcdat_(LIP I1, LIP I2, LIP I3, LIP I4, LIP I5, DBP VAL, LIP NOERR);
CMT tqwasc_(CHP FILE, LIP NOERR,ftnlen);
CMT tqchar_(LIP INDEXP, LIP INDEXC, DBP VAL, LIP NOERR);

/* Added for ChemApp V544 */

CMT tqcnsc_ (LIP INDEXS, CHP NAME, LIP NOERR,ftnlen);


#else

CMT TQINI (LIP NOERR);
CMT TQOPEN(CHP FILE, LNT FILELEN, LIP LUN, LIP NOERR);
CMT TQCLOS(LIP LUN, LIP NOERR);
CMT TQGIO (CHP OPTION, LNT OPTIONLEN, LIP IVAL, LIP NOERR);
CMT TQCIO (CHP OPTION, LNT OPTIONLEN, LIP IVAL, LIP NOERR);
CMT TQRFIL(LIP NOERR);
CMT TQGSU (CHP OPTION, LNT OPTIONLEN, CHP UNIT, LNT UNITLEN, LIP NOERR);
CMT TQCSU (CHP OPTION, LNT OPTIONLEN, CHP UNIT, LNT UNITLEN, LIP NOERR);
CMT TQINSC(CHP NAME, LNT NAMELEN, LIP INDEXS, LIP NOERR);
CMT TQGNSC(LIP INDEXS, CHP NAME, LNT NAMELEN, LIP NOERR);
CMT TQNOSC(LIP NSCOM, LIP NOERR);
CMT TQSTSC(LIP INDEXS, DBP STOI, DBP WMASS, LIP NOERR);
CMT TQCSC (CHP NAME, LNT NAMELEN, LIP NOERR);
CMT TQINP (CHP NAME, LNT NAMELEN, LIP INDEXP, LIP NOERR);
CMT TQGNP (LIP INDEXP, CHP NAME, LNT NAMELEN, LIP NOERR);
CMT TQNOP (LIP NPHASE, LIP NOERR);
CMT TQINPC(CHP NAME, LNT NAMELEN, LIP INDEXP, LIP INDEXC, LIP NOERR);
CMT TQGNPC(LIP INDEXP, LIP INDEXC, CHP NAME, LNT NAMELEN, LIP NOERR);
CMT TQNOPC(LIP INDEXP, LIP NPCON, LIP NOERR);
CMT TQSTPC(LIP INDEXP, LIP INDEXC, DBP STOI, DBP WMASS, LIP NOERR);
CMT TQGSP (LIP INDEXP, CHP OPTION, LNT OPTIONLEN, LIP NOERR);
CMT TQCSP (LIP INDEXP, CHP OPTION, LNT OPTIONLEN, LIP NOERR);
CMT TQGSPC(LIP INDEXP, LIP INDEXC, CHP OPTION, LNT OPTIONLEN, LIP NOERR);
CMT TQCSPC(LIP INDEXP, LIP INDEXC, CHP OPTION, LNT OPTIONLEN, LIP NOERR);
CMT TQSETC(CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEX, DBP VAL, LIP NUMCON, LIP NOERR);
CMT TQREMC(LIP NUMCON, LIP NOERR);
CMT TQSTTP(CHP IDENTS, LNT IDENTSLEN, DBP VALS, LIP NOERR);
CMT TQSTCA(CHP IDENTS, LNT IDENTSLEN, LIP INDEXP, LIP INDEXC, DBP VAL, LIP NOERR);
CMT TQSTEC(CHP OPTION, LNT OPTIONLEN, LIP INDEXP, DBP VAL, LIP NOERR);
CMT TQSTRM(CHP IDENTS, LNT IDENTSLEN, LIP NOERR);
CMT TQCE  (CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR);
CMT TQCEL (CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR);
CMT TQCLIM(CHP OPTION, LNT OPTIONLEN, DBP VAL, LIP NOERR);
CMT TQGETR(CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEX, DBP VAL, LIP NOERR);
CMT TQGDPC(CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC,DBP VAL, LIP NOERR);
CMT TQSHOW(LIP NOERR);
CMT TQERR (CHP MESS, LNT MESSLEN, LIP NOERR);

/* Added for ChemApp V113 */

CMT TQCPRT(LIP NOERR);
CMT TQVERS(LIP NVERS, LIP NOERR);
CMT TQSIZE(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,LIP NJ,LIP NK,LIP NOERR);
CMT TQMODL(LIP INDEXP, CHP NAME, LNT NAMELEN, LIP NOERR);
CMT TQSTXP(CHP IDENTS,LNT IDENTSLEN,CHP OPTION,LNT OPTIONLEN, DBP VAL, LIP NOERR);

/* Added for ChemApp V201 */

CMT TQLITE(LIP LITE, LIP NOERR);

/* Added for ChemApp V214 */

CMT TQRBIN(LIP NOERR);

/* Added for ChemApp V300 */

CMT TQMAP  (CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP ICONT, LIP NOERR);
CMT TQMAPL (CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP ICONT, LIP NOERR);

/* Added for ChemApp V310 */

CMT TQPCIS (LIP INDEXP, LIP INDEXC, LIP ISPERM, LIP NOERR);

/* Added for ChemApp V330 */

CMT TQOPNA(CHP FILE, LNT FILELEN, LIP LUN, LIP NOERR);
CMT TQOPNB(CHP FILE, LNT FILELEN, LIP LUN, LIP NOERR);

/* Added for ChemApp V331 */

CMT TQNOSL(LIP INDEXP, LIP NSUBL, LIP NOERR);
CMT TQNOLC(LIP INDEXP, LIP INDEXL, LIP NSLCON, LIP NOERR);
CMT TQINLC(CHP NAME, LNT OPTIONLEN, LIP INDEXP, LIP INDEXL, LIP INDEXC, LIP NOERR);
CMT TQGNLC(LIP INDEXP, LIP INDEXL, LIP INDEXC, CHP NAME, LNT OPTIONLEN, LIP NOERR);
CMT TQGTLC(LIP INDEXP, LIP INDEXL, LIP INDEXC,DBP VAL, LIP NOERR);

/* Added for ChemApp V339 */
/*
CMT TQGOPN(CHP FILE,LNT FILELEN,LIP LUN,CHP FFORM,LNT FFORMLEN,
	   CHP FSTAT,LNT FSTATLEN,CHP FACC,LNT FACCLEN,LIP RECL,
	   LIP IOSTAT,LIP NOERR);
*/

CMT TQBOND(LIP INDEXP, LIP INDEXA, LIP INDEXB, LIP INDEXC, LIP INDEXD, 
	   DBP VAL, LIP NOERR);


/* Added for ChemApp V340 */
CMT TQUSED(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,
	   LIP NJ,LIP NK,LIP NOERR);
CMT TQGTRH(LIP TFHVER,
	   CHP TFHNWP,
	   LNT TFHNWPLEN,
	   LIP TFHVNW,
	   CHP TFHNRP,
	   LNT TFHNRPLEN,
	   LIP TFHVNR,
	   LIP TFHDTC,
	   LIP TFHDTE,
	   CHP TFHID,
	   LNT TFHIDLEN,
	   CHP TFHUSR,
	   LNT TFHUSRLEN,
	   CHP TFHREM,
	   LNT TFHREMLEN,
	   LIP NOERR);
CMT TQOPNT(CHP FILE, LNT FTNLEN, LIP LUN, LIP NOERR);
CMT TQRCST(LIP NOERR);
CMT TQGTID(CHP ID, LNT IDLEN, LIP NOERR);
CMT TQGTNM(CHP NAME, LNT FTNLEN, LIP NOERR);
CMT TQGTPI(CHP PID, LNT PIDLEN, LIP NOERR);
CMT TQWSTR(CHP OPTION, LNT OPTIONLEN, CHP CTXT, LNT CTXTLEN, LIP NOERR);

/* Added for ChemApp V414 */

CMT TQGTED(LIP EDMON, LIP EDYEAR, LIP NOERR);
CMT TQGTHI(CHP HASPT, LNT OPTIONLEN, LIP HASPID, LIP NOERR);
CMT TQCEN (CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR);
CMT TQCENL(CHP OPTION, LNT OPTIONLEN, LIP INDEXP, LIP INDEXC, DBP VALS, LIP NOERR);

/* Added for ChemApp V540 */

CMT TQGDAT(LIP INDEXP, LIP INDEXC, CHP OPTION, LNT OPTIONLEN, LIP INDEXR, LIP NVALV, DBP VALV, LIP NOERR);
/*
********************************
Still missing:
TQLPAR
TQGPAR
********************************
*/
/*CMT TQLPAR(LIP INDEXP, CHP OPTION, LNT OPTIONLEN, LIP NOPAR, CHP CHRPAR, LNT CHRPARLEN, LIP NOERR);*/
CMT TQCDAT(LIP I1, LIP I2, LIP I3, LIP I4, LIP I5, DBP VAL, LIP NOERR);
CMT TQWASC(CHP FILE, LNT FILELEN, LIP NOERR);
CMT TQCHAR(LIP INDEXP, LIP INDEXC, DBP VAL, LIP NOERR);

/* Added for ChemApp V544 */

CMT TQCNSC (LIP INDEXS, CHP NAME, LNT NAMELEN, LIP NOERR);

#endif

#ifdef __cplusplus
};
#endif      /* __cplusplus        */

/* Default global buffer for the Error Message */

char TQERRMSG[3][80];

/* internal procedure 
   removes trailing spaces from string t (length l) */

void TQremspaces(char *t,int l)
{
      int i;
      
      i=l-1;
      /* Warning, this deletes last char in string!!! */
      t[i--]=0;      
      
      while( (i>=0) && (t[i]==' ') ) t[i--]=0;
}

/* ChemApp Interface Routines */

int tqini (LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqini_(NOERR);
#else
  TQINI(NOERR);
#endif
  return(*NOERR);
}

int tqopen(CHP FILE, LI LUN, LIP NOERR)
{
  LI i1=LUN;
#if defined(UNIX) && !defined(CRAY)
  tqopen_(FILE,&i1,NOERR,(ftnlen)strlen(FILE));
#else
  TQOPEN(FILE,strlen(FILE),&i1,NOERR);
#endif
  return(*NOERR);
}

int tqclos(LI LUN, LIP NOERR)
{
  LI i1=LUN;
#if defined(UNIX) && !defined(CRAY)
  tqclos_(&i1,NOERR);
#else
  TQCLOS(&i1,NOERR);
#endif
  return(*NOERR);
}

int tqgio (CHP OPTION, LIP IVAL, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgio_(OPTION,IVAL,NOERR,(ftnlen)strlen(OPTION));
#else
  TQGIO(OPTION,strlen(OPTION),IVAL,NOERR);
#endif
  return(*NOERR);
}

int tqcio (CHP OPTION, LI IVAL, LIP NOERR)
{
  LI i1=IVAL;
#if defined(UNIX) && !defined(CRAY)
  tqcio_(OPTION,&i1,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCIO(OPTION,strlen(OPTION),&i1,NOERR);
#endif
  return(*NOERR);
}

int tqrfil(LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqrfil_(NOERR);
#else
  TQRFIL(NOERR);
#endif
  return(*NOERR);
}

int tqgsu (CHP OPTION, CHP UNIT, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgsu_ (OPTION,UNIT,NOERR,(ftnlen)strlen(OPTION),TQSTRLEN);
#else
  TQGSU (OPTION,strlen(OPTION),UNIT,TQSTRLEN,NOERR);
#endif
  TQremspaces(UNIT,TQSTRLEN); 
  return(*NOERR);
}

int tqcsu (CHP OPTION, CHP UNIT, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqcsu_ (OPTION,UNIT,NOERR,(ftnlen)strlen(OPTION),(ftnlen)strlen(UNIT));
#else
  TQCSU (OPTION,strlen(OPTION),UNIT,strlen(UNIT),NOERR);
#endif
  return(*NOERR);
}

int tqinsc(CHP NAME, LIP INDEXS, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqinsc_(NAME,INDEXS,NOERR,(ftnlen)strlen(NAME));
#else
  TQINSC(NAME,strlen(NAME),INDEXS,NOERR);
#endif
  return(*NOERR);
}

int tqgnsc(LI INDEXS, CHP NAME, LIP NOERR)
{            
  LI i1=INDEXS;
#if defined(UNIX) && !defined(CRAY)
  tqgnsc_(&i1,NAME,NOERR,TQSTRLEN);
#else
  TQGNSC(&i1,NAME,TQSTRLEN,NOERR);
#endif
  TQremspaces(NAME,TQSTRLEN);
  return(*NOERR);
}

int tqnosc(LIP NSCOM, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqnosc_(NSCOM,NOERR);
#else
  TQNOSC(NSCOM,NOERR);
#endif
  return(*NOERR);
}

int tqstsc(LI INDEXS, DBP STOI, DBP WMASS, LIP NOERR)
{            
  LI i1=INDEXS;
#if defined(UNIX) && !defined(CRAY)
  tqstsc_(&i1,STOI,WMASS,NOERR);
#else
  TQSTSC(&i1,STOI,WMASS,NOERR);
#endif
  return(*NOERR);
}

int tqcsc (CHP NAME, LIP NOERR)
{
  char buf[30*TQSTRLEN];
  int i,j; 
  i=0;
  while(*(NAME+i*TQSTRLEN)!=0) {                    
     strcpy(buf+i*(TQSTRLEN-1),NAME+i*TQSTRLEN);
     for(j=strlen(NAME+i*TQSTRLEN);j<TQSTRLEN-1;j++)
     {
    *(buf+i*(TQSTRLEN-1)+j)=' ';
     }  
     i++;
  }

#if defined(UNIX) && !defined(CRAY)
  tqcsc_ (buf,NOERR,(ftnlen)TQSTRLEN-1);
#else
  TQCSC (buf,TQSTRLEN-1,NOERR);
#endif
  return(*NOERR);
}

int tqinp (CHP NAME, LIP INDEXP, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqinp_ (NAME,INDEXP,NOERR,(ftnlen)strlen(NAME));
#else
  TQINP (NAME,strlen(NAME),INDEXP,NOERR);
#endif
  return(*NOERR);
}

int tqgnp (LI INDEXP, CHP NAME, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqgnp_ (&i1,NAME,NOERR,TQSTRLEN);
#else
  TQGNP (&i1,NAME,TQSTRLEN,NOERR);
#endif
  TQremspaces(NAME,TQSTRLEN);
  return(*NOERR);
}

int tqnop (LIP NPHASE, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqnop_ (NPHASE,NOERR);
#else
  TQNOP (NPHASE,NOERR);
#endif
  return(*NOERR);
}

int tqinpc(CHP NAME, LI INDEXP, LIP INDEXC, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqinpc_(NAME,&i1,INDEXC,NOERR,(ftnlen)strlen(NAME));
#else
  TQINPC(NAME,strlen(NAME),&i1,INDEXC,NOERR);
#endif
  return(*NOERR);
}

int tqgnpc(LI INDEXP, LI INDEXC, CHP NAME, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqgnpc_(&i1,&i2,NAME,NOERR,TQSTRLEN);
#else
  TQGNPC(&i1,&i2,NAME,TQSTRLEN,NOERR);
#endif
  TQremspaces(NAME,TQSTRLEN);
  return(*NOERR);
}

int tqnopc(LI INDEXP, LIP NPCON, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqnopc_(&i1,NPCON,NOERR);
#else
  TQNOPC(&i1,NPCON,NOERR);
#endif
  return(*NOERR);
}

int tqstpc(LI INDEXP, LI INDEXC, DBP STOI, DBP WMASS, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqstpc_(&i1,&i2,STOI,WMASS,NOERR);
#else
  TQSTPC(&i1,&i2,STOI,WMASS,NOERR);
#endif
  return(*NOERR);
}

int tqgsp (LI INDEXP, CHP OPTION, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqgsp_ (&i1,OPTION,NOERR,TQSTRLEN);
#else
  TQGSP (&i1,OPTION,TQSTRLEN,NOERR);
#endif
  TQremspaces(OPTION,TQSTRLEN);
  return(*NOERR);
}

int tqcsp (LI INDEXP, CHP OPTION, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqcsp_ (&i1,OPTION,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCSP (&i1,OPTION,strlen(OPTION),NOERR);
#endif
  return(*NOERR);
}

int tqgspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqgspc_(&i1,&i2,OPTION,NOERR,TQSTRLEN);
#else
  TQGSPC(&i1,&i2,OPTION,TQSTRLEN,NOERR);
#endif
  TQremspaces(OPTION,TQSTRLEN);
  return(*NOERR);
}

int tqcspc(LI INDEXP, LI INDEXC, CHP OPTION, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqcspc_(&i1,&i2,OPTION,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCSPC(&i1,&i2,OPTION,strlen(OPTION),NOERR);
#endif
  return(*NOERR);
}

int tqsetc(CHP OPTION, LI INDEXP, LI INDEX, DB VAL, LIP NUMCON, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEX;
  DB d1=VAL;     
#if defined(UNIX) && !defined(CRAY)
  tqsetc_(OPTION,&i1,&i2,&d1,NUMCON,NOERR,(ftnlen)strlen(OPTION));
#else
  TQSETC(OPTION,strlen(OPTION),&i1,&i2,&d1,NUMCON,NOERR);
#endif
  return(*NOERR);
}

int tqremc(LI NUMCON, LIP NOERR)
{
  LI i1=NUMCON;
#if defined(UNIX) && !defined(CRAY)
  tqremc_(&i1,NOERR);
#else
  TQREMC(&i1,NOERR);
#endif
  return(*NOERR);
}

int tqsttp(CHP IDENTS, DBP VALS, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqsttp_(IDENTS,VALS,NOERR,(ftnlen)strlen(IDENTS));
#else
  TQSTTP(IDENTS,strlen(IDENTS),VALS,NOERR);
#endif
  return(*NOERR);
}

int tqstca(CHP IDENTS, LI INDEXP, LI INDEXC, DB VAL, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
  DB d1=VAL;     
#if defined(UNIX) && !defined(CRAY)
  tqstca_(IDENTS,&i1,&i2,&d1,NOERR,(ftnlen)strlen(IDENTS));
#else
  TQSTCA(IDENTS,strlen(IDENTS),&i1,&i2,&d1,NOERR);
#endif
  return(*NOERR);
}

int tqstec(CHP OPTION, LI INDEXP, DB VAL, LIP NOERR)
{
  LI i1=INDEXP;
  DB d1=VAL;
#if defined(UNIX) && !defined(CRAY)
  tqstec_(OPTION,&i1,&d1,NOERR,(ftnlen)strlen(OPTION));
#else
  TQSTEC(OPTION,strlen(OPTION),&i1,&d1,NOERR);
#endif
  return(*NOERR);
}

int tqstrm(CHP IDENTS, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqstrm_(IDENTS,NOERR,(ftnlen)strlen(IDENTS));
#else
  TQSTRM(IDENTS,strlen(IDENTS),NOERR);
#endif
  return(*NOERR);
}

int tqce (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqce_  (OPTION,&i1,&i2,VALS,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCE  (OPTION,strlen(OPTION),&i1,&i2,VALS,NOERR);
#endif
  return(*NOERR);
}

int tqcel (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqcel_ (OPTION,&i1,&i2,VALS,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCEL (OPTION,strlen(OPTION),&i1,&i2,VALS,NOERR);
#endif
  return(*NOERR);
}

int tqclim(CHP OPTION, DB VAL, LIP NOERR)
{
  DB d1=VAL;
#if defined(UNIX) && !defined(CRAY)
  tqclim_(OPTION,&d1,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCLIM(OPTION,strlen(OPTION),&d1,NOERR);
#endif
  return(*NOERR);
}

int tqgetr(CHP OPTION, LI INDEXP, LI INDEX, DBP VAL, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEX;
#if defined(UNIX) && !defined(CRAY)
  tqgetr_(OPTION,&i1,&i2,VAL,NOERR,(ftnlen)strlen(OPTION));
#else
  TQGETR(OPTION,strlen(OPTION),&i1,&i2,VAL,NOERR);
#endif
  return(*NOERR);
}

int tqgdpc(CHP OPTION, LI INDEXP, LI INDEXC,DBP VAL, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqgdpc_(OPTION,&i1,&i2,VAL,NOERR,(ftnlen)strlen(OPTION));
#else
  TQGDPC(OPTION,strlen(OPTION),&i1,&i2,VAL,NOERR);
#endif
  return(*NOERR);
}

int tqshow(LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqshow_(NOERR);
#else
  TQSHOW(NOERR);
#endif
  return(*NOERR);
}

int tqerr(CHP MESS, LIP NOERR)
{ 
  char buf[250];
#if defined(UNIX) && !defined(CRAY)
  tqerr_ (buf,NOERR,80);
#else
  TQERR (buf,80,NOERR);
#endif
  buf[79]=0;buf[159]=0;buf[239]=0;
  TQremspaces(buf,78);
  TQremspaces(buf+80,78);
  TQremspaces(buf+160,78);
  strcpy(MESS,buf);
  strcpy(MESS+80,buf+80);
  strcpy(MESS+160,buf+160);
  return(*NOERR);
}

/* --------------------------------------------------------
     Added for ChemApp V113 , 24/07/97 cat
   -------------------------------------------------------- */

int tqcprt(LIP NOERR)
{                        
#if defined(UNIX) && !defined(CRAY)
  tqcprt_(NOERR);
#else
  TQCPRT(NOERR);
#endif
  return(*NOERR);
}

int tqvers(LIP NVERS, LIP NOERR)
{                                 
#if defined(UNIX) && !defined(CRAY)
  tqvers_(NVERS,NOERR);
#else
  TQVERS(NVERS,NOERR);
#endif
  return(*NOERR);
}

int tqsize(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,LIP NJ,LIP NK,LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqsize_(NA,NB,NC,ND,NE,NF,NG,NH,NI,NJ,NK,NOERR);
#else
  TQSIZE(NA,NB,NC,ND,NE,NF,NG,NH,NI,NJ,NK,NOERR);
#endif
  return(*NOERR);
}

int tqmodl(LI INDEXP, CHP NAME, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqmodl_(&i1,NAME,NOERR,TQSTRLEN);
#else
  TQMODL(&i1,NAME,TQSTRLEN,NOERR);
#endif
  TQremspaces(NAME,TQSTRLEN);
  return(*NOERR);
}

int tqstxp(CHP IDENTS,CHP OPTION, DBP VAL, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqstxp_(IDENTS,OPTION,VAL,NOERR,(ftnlen)strlen(IDENTS),(ftnlen)strlen(OPTION));
#else
  TQSTXP(IDENTS,strlen(IDENTS),OPTION,strlen(OPTION),VAL,NOERR);
#endif
  return(*NOERR);
}

int tqlite (LIP LITE, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqlite_ (LITE,NOERR);
#else
  TQLITE (LITE,NOERR);
#endif
  return(*NOERR);
}

int tqrbin(LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqrbin_(NOERR);
#else
  TQRBIN(NOERR);
#endif
  return(*NOERR);
}

int tqmap (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqmap_  (OPTION,&i1,&i2,VALS,ICONT,NOERR,(ftnlen)strlen(OPTION));
#else
  TQMAP  (OPTION,strlen(OPTION),&i1,&i2,VALS,ICONT,NOERR);
#endif
  return(*NOERR);
}

int tqmapl (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP ICONT, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqmapl_  (OPTION,&i1,&i2,VALS,ICONT,NOERR,(ftnlen)strlen(OPTION));
#else
  TQMAPL  (OPTION,strlen(OPTION),&i1,&i2,VALS,ICONT,NOERR);
#endif
  return(*NOERR);
}

int tqpcis(LI INDEXP, LI INDEXC, LIP ISPERM, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqpcis_(&i1,&i2,ISPERM,NOERR);
#else
  TQPCIS(&i1,&i2,ISPERM,NOERR);
#endif
  return(*NOERR);
}

int tqopna(CHP FILE, LI LUN, LIP NOERR)
{
  LI i1=LUN;
#if defined(UNIX) && !defined(CRAY)
  tqopna_(FILE,&i1,NOERR,(ftnlen)strlen(FILE));
#else
  TQOPNA(FILE,strlen(FILE),&i1,NOERR);
#endif
  return(*NOERR);
}

int tqopnb(CHP FILE, LI LUN, LIP NOERR)
{
  LI i1=LUN;
#if defined(UNIX) && !defined(CRAY)
  tqopnb_(FILE,&i1,NOERR,(ftnlen)strlen(FILE));
#else
  TQOPNB(FILE,strlen(FILE),&i1,NOERR);
#endif
  return(*NOERR);
}

int tqnosl(LI INDEXP, LIP NSUBL, LIP NOERR)
{
  LI i1=INDEXP;
#if defined(UNIX) && !defined(CRAY)
  tqnosl_(&i1,NSUBL,NOERR);
#else
  TQNOSL(&i1,NSUBL,NOERR);
#endif
  return(*NOERR);
}

int tqnolc(LI INDEXP, LI INDEXL, LIP NSLCON, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXL;
#if defined(UNIX) && !defined(CRAY)
  tqnolc_(&i1,&i2,NSLCON,NOERR);
#else
  TQNOLC(&i1,&i2,NSLCON,NOERR);
#endif
  return(*NOERR);
}

int tqinlc(CHP NAME, LI INDEXP, LI INDEXL, LIP INDEXC, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXL;
#if defined(UNIX) && !defined(CRAY)
  tqinlc_(NAME,&i1,&i2,INDEXC,NOERR,(ftnlen)strlen(NAME));
#else
  TQINLC(NAME,strlen(NAME),&i1,&i2,INDEXC,NOERR);
#endif
  return(*NOERR);
}

int tqgnlc(LI INDEXP, LI INDEXL, LI INDEXC, CHP NAME, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXL;    
  LI i3=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqgnlc_(&i1,&i2,&i3,NAME,NOERR,TQSTRLEN);
#else
  TQGNLC(&i1,&i2,&i3,NAME,TQSTRLEN,NOERR);
#endif
  TQremspaces(NAME,TQSTRLEN);
  return(*NOERR);
}

int tqgtlc(LI INDEXP, LI INDEXL, LI INDEXC, DBP VAL, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXL;    
  LI i3=INDEXC;    
#if defined(UNIX) && !defined(CRAY)
  tqgtlc_(&i1,&i2,&i3,VAL,NOERR);
#else
  TQGTLC(&i1,&i2,&i3,VAL,NOERR);
#endif
  return(*NOERR);
}

/*
int tqgopn (CHP FILE,LI LUN,CHP FFORM,CHP FSTAT,CHP FACC,LI RECL,
	    LIP IOSTAT,LIP NOERR)
{
  LI i1=LUN;
  LI i2=RECL;
#if defined(UNIX) && !defined(CRAY)
  tqgopn_ (FILE,&i1,FFORM,FSTAT,FACC,&i2,IOSTAT,NOERR,
	   (ftnlen)strlen(FILE),(ftnlen)strlen(FFORM),
	   (ftnlen)strlen(FSTAT),(ftnlen)strlen(FACC));
#else
  TQCSU (FILE,strlen(FILE),&i1,FFORM,strlen(FFORM),
	 FSTAT,strlen(FSTAT),FACC,strlen(FACC),&i2,IOSTAT,NOERR);
#endif
  return(*NOERR);
}

*/

int tqbond(LI INDEXP, LI INDEXA, LI INDEXB, LI INDEXC, LI INDEXD, 
	   DBP VAL, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXA;    
  LI i3=INDEXB;    
  LI i4=INDEXC;    
  LI i5=INDEXD;    
#if defined(UNIX) && !defined(CRAY)
  tqbond_(&i1,&i2,&i3,&i4,&i5,VAL,NOERR);
#else
  TQBOND(&i1,&i2,&i3,&i4,&i5,VAL,NOERR);
#endif
  return(*NOERR);
}

int tqused(LIP NA,LIP NB,LIP NC,LIP ND,LIP NE,LIP NF,LIP NG,LIP NH,LIP NI,
	   LIP NJ,LIP NK,LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqused_(NA,NB,NC,ND,NE,NF,NG,NH,NI,NJ,NK,NOERR);
#else
  TQUSED(NA,NB,NC,ND,NE,NF,NG,NH,NI,NJ,NK,NOERR);
#endif
  return(*NOERR);
}

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
	   LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgtrh_(TFHVER,
	  TFHNWP,
	  TFHVNW,
	  TFHNRP,
	  TFHVNR,
	  TFHDTC,
	  TFHDTE,
	  TFHID,
	  TFHUSR,
	  TFHREM,
	  NOERR,
	  40,
	  40,
	  255,
	  80,
	  80);
#else
  TQGTRH(TFHVER,
	 TFHNWP,
	 40,
	 TFHVNW,
	 TFHNRP,
	 40,
	 TFHVNR,
	 TFHDTC,
	 TFHDTE,
	 TFHID,
	 255,
	 TFHUSR,
	 80,
	 TFHREM,
	 80,
	 NOERR);
#endif
  TQremspaces(TFHNWP,40);
  TQremspaces(TFHNRP,40);
  TQremspaces(TFHID,255);
  TQremspaces(TFHUSR,80);
  TQremspaces(TFHREM,80);
  return(*NOERR);
}


int tqopnt(CHP FILE, LI LUN, LIP NOERR)
{
  LI i1=LUN;
#if defined(UNIX) && !defined(CRAY)
  tqopnt_(FILE,&i1,NOERR,(ftnlen)strlen(FILE));
#else
  TQOPNT(FILE,strlen(FILE),&i1,NOERR);
#endif
  return(*NOERR);
}

int tqrcst(LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqrcst_(NOERR);
#else
  TQRCST(NOERR);
#endif
  return(*NOERR);
}

int tqgtid(CHP ID, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgtid_(ID,NOERR,255);
#else
  TQGTID(ID,255,NOERR);
#endif
  TQremspaces(ID,255);
  return(*NOERR);
}

int tqgtnm(CHP NAME, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgtnm_(NAME,NOERR,80);
#else
  TQGTNM(NAME,80,NOERR);
#endif
  TQremspaces(NAME,80);
  return(*NOERR);
}

int tqgtpi(CHP PID, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgtpi_(PID,NOERR,TQSTRLEN);
#else
  TQGTPI(PID,TQSTRLEN,NOERR);
#endif
  TQremspaces(PID,TQSTRLEN);
  return(*NOERR);
}

int tqwstr(CHP OPTION, CHP CTXT, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqwstr_(OPTION,CTXT,NOERR,(ftnlen)strlen(OPTION),(ftnlen)strlen(CTXT));
#else
  TQWSTR(OPTION,strlen(OPTION),CTXT,strlen(CTXT),NOERR);
#endif
  return(*NOERR);
}


/* Added for ChemApp V414 */

int tqgted(LIP EDMON, LIP EDYEAR, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgted_(EDMON, EDYEAR, NOERR);
#else
  TQGTED(EDMON, EDYEAR, NOERR);
#endif
  return(*NOERR);
}

int tqgthi(CHP HASPT, LIP HASPID, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqgthi_ (HASPT,HASPID,NOERR,TQSTRLEN);
#else
  TQGTHI (HASPT,TQSTRLEN,HASPID,NOERR);
#endif
  TQremspaces(HASPT,TQSTRLEN);
  return(*NOERR);
}

int tqcen (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqcen_  (OPTION,&i1,&i2,VALS,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCEN  (OPTION,strlen(OPTION),&i1,&i2,VALS,NOERR);
#endif
  return(*NOERR);
}

int tqcenl (CHP OPTION, LI INDEXP, LI INDEXC, DBP VALS, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqcenl_  (OPTION,&i1,&i2,VALS,NOERR,(ftnlen)strlen(OPTION));
#else
  TQCENL  (OPTION,strlen(OPTION),&i1,&i2,VALS,NOERR);
#endif
  return(*NOERR);
}


int tqwasc(CHP FILE, LIP NOERR)
{
#if defined(UNIX) && !defined(CRAY)
  tqwasc_(FILE,NOERR,(ftnlen)strlen(FILE));
#else
  TQWASC(FILE,strlen(FILE),NOERR);
#endif
  return(*NOERR);
}

int tqgdat(LI INDEXP, LI INDEXC, CHP OPTION, LI INDEXR, LIP NVALV, DBP VALV, LIP NOERR)
{
  LI i1=INDEXP,i2=INDEXC,i3=INDEXR;
#if defined(UNIX) && !defined(CRAY)
  tqgdat_(&i1,&i2,OPTION,&i3,NVALV,VALV,NOERR,(ftnlen)strlen(OPTION));
#else
  TQGDAT(&i1,&i2,OPTION,strlen(OPTION),&i3,NVALV,VALV,NOERR);
#endif
  return(*NOERR);
}

int tqcdat(LI I1, LI I2, LI I3, LI I4, LI I5, DB VAL, LIP NOERR)
{
  LI i1=I1,i2=I2,i3=I3,i4=I4,i5=I5;
  DB d1=VAL;     
#if defined(UNIX) && !defined(CRAY)
  tqcdat_(&i1,&i2,&i3,&i4,&i5,&d1,NOERR);
#else
  TQCDAT(&i1,&i2,&i3,&i4,&i5,&d1,NOERR);
#endif
  return(*NOERR);
}

int tqchar(LI INDEXP, LI INDEXC, DBP VAL, LIP NOERR)
{
  LI i1=INDEXP;
  LI i2=INDEXC;
#if defined(UNIX) && !defined(CRAY)
  tqchar_(&i1,&i2,VAL,NOERR);
#else
  TQCHAR(&i1,&i2,VAL,NOERR);
#endif
  return(*NOERR);
}

int tqcnsc (LI INDEXS, CHP NAME, LIP NOERR)
{
  LI i1=INDEXS;
#if defined(UNIX) && !defined(CRAY)
  tqcnsc_ (&i1,NAME,NOERR,(ftnlen)strlen(NAME));
#else
  TQCNSC (&i1,NAME,strlen(NAME),NOERR);
#endif
  return(*NOERR);
}
 #endif

