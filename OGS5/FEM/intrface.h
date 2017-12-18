/**************************************************************************/
/* ROCKFLOW - Modul: intrface.h
 */
/* Aufgabe:
   Stellt allgemeine Funktionszeiger bereit.

   letzte Aenderung:  RK   03/2003
                                                                          */
/**************************************************************************/

#ifndef intrface_INC

#define intrface_INC
/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */

/* Deklarationen */

/* Elementdaten-Schnittstelle */

extern VoidFuncLDXD SetElementJacobiMatrix;
/* Setzt die inverse Jakobi-Matrix und deren Determinante
   im Zentrum des Elements number */

extern DoubleXFuncLDX GetElementJacobiMatrix;
/* Liefert die inverse Jakobi-Matrix und deren Determinatne
im Zentrum des Elements number */

extern VoidFuncLong DeleteElementJacobiMatrix;
/* Erdet den Jacobi-Matrix-Zeiger des Elements number */

extern DoubleFuncLong GetElementPecletNumber;

/* Weitere externe Objekte */
#endif
