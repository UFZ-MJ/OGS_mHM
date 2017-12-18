/*
 * CRSIO.h
 *
 *  Created on: Feb 7, 2011
 *      Author: TF
 */

#ifndef CRSIO_H_
#define CRSIO_H_

namespace FileIO {

template<class T> void
writeCompressedStorageFmt (std::ostream &os, unsigned n, unsigned* iA, unsigned* jA, T* A)
{
  os.write((char*) &n, sizeof(unsigned));
  os.write((char*) iA, (n+1)*sizeof(unsigned));
  os.write((char*) jA, iA[n]*sizeof(unsigned));
  os.write((char*) A, iA[n]*sizeof(T));
}

template<class T> void
readCompressedStorageFmt (std::istream &is, unsigned &n, unsigned* &iA, unsigned* &jA, T* &A)
{
  is.read((char*) &n, sizeof(unsigned));
  if (iA != NULL) {
    delete [] iA;
    delete [] jA;
    delete [] A;
  }
  iA = new unsigned[n+1];
  assert(iA!=NULL);
  is.read((char*) iA, (n+1)*sizeof(unsigned));

  jA = new unsigned[iA[n]];
  assert(jA!=NULL);
  is.read((char*) jA, iA[n]*sizeof(unsigned));

  A = new T[iA[n]];
  assert(A!=NULL);
  is.read((char*) A, iA[n]*sizeof(T));

#ifndef NDEBUG
  // do simple checks
  if (iA[0]!=0)
    std::cerr << std::endl
              << "CRS matrix: array iA doesn't start with 0" << std::endl;

  unsigned i = 0;
  while (i<iA[n] && jA[i]<n) ++i;
  if (i<iA[n])
    std::cerr << std::endl
              << "CRS matrix: the " << i << "th entry of jA has the value "
              << jA[i] << ", which is out of bounds." << std::endl;
#endif
}

} // end namespace FileIO

#endif /* CRSIO_H_ */
