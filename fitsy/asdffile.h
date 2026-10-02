// Copyright (C) 1999-2021
// Smithsonian Astrophysical Observatory, Cambridge, MA, USA
// For conditions of distribution and use, see copyright notice in "copyright"

#ifndef __fitsasdffile_h__
#define __fitsasdffile_h__

#include <string>
#include <vector>

#include "file.h"

// ASDF support, backed by libasdf. Named asdffile.h rather than asdf.h so
// it does not shadow libasdf's own <asdf.h>.
//
//   FitsAsdfFile  an open ASDF file: its tree text, its ndarrays, and raw
//                 access to its blocks (for the Tcl side, via tclfitsy)
//
//   FitsAsdf      one ndarray presented to DS9 as an image (a FitsFile)
//                 with a synthetic header

// libasdf's file handle, forward declared so that this header does not
// pull libasdf's headers into everything that includes it
struct asdf_file;

// deepest ndarray rank reported; only rank 2 and 3 are loadable
#define ASDF_MAXDIM 8

// why an ndarray cannot be loaded as pixels
#define ASDF_OK         0x00
#define ASDF_BYTEORDER  0x01
#define ASDF_STRIDES    0x02
#define ASDF_OFFSET     0x04
#define ASDF_MASK       0x08

// one block backed ndarray in the tree
class FitsAsdfNode {
 public:
  std::string path_;         // root relative, e.g. roman/data
  int source_;               // binary block index
  std::string datatype_;     // asdf name, e.g. float32
  int big_;                  // 1 big endian, 0 little, -1 absent
  int ndim_;
  long long dim_[ASDF_MAXDIM];  // row major: slowest axis first
  unsigned unsupported_;
  int hasMaskScalar_;
  double maskScalar_;

 public:
  FitsAsdfNode();
};

class FitsAsdfFile {
 protected:
  struct asdf_file* file_;
  std::string fn_;
  std::string tree_;
  int treeRead_;
  std::vector<FitsAsdfNode> nodes_;

  int valid_;
  std::string error_;

 protected:
  void enumNdarrays();
  void setError(const char*, ...);

 public:
  FitsAsdfFile(const char*);
  virtual ~FitsAsdfFile();

  int isValid() {return valid_;}
  const char* error() {return error_.empty() ? NULL : error_.c_str();}

  // libasdf's handle, for FitsAsdf
  struct asdf_file* file() {return file_;}

  // the YAML tree, verbatim
  const char* tree();
  size_t treeSize();

  int nnode() {return (int)nodes_.size();}
  const FitsAsdfNode* node(int ii);
  const FitsAsdfNode* find(const char*);

  int nblock();

  // Decompressed bytes of one block. The caller owns the result and
  // releases it with delete[]. NULL on failure, with error() set.
  char* block(int, size_t*);
};

class FitsAsdf : public FitsFile {
 protected:
  FitsAsdfFile* asdf_;
  std::string asdfError_;

 protected:
  int build(const char*);

 public:
  FitsAsdf(const char*, const char*);
  virtual ~FitsAsdf();

  virtual const char* errorMessage()
  {return asdfError_.empty() ? NULL : asdfError_.c_str();}

  FitsAsdfFile* asdfFile() {return asdf_;}
};

// One further plane of a rank 3 array. DS9 presents a cube as a chain of
// single plane images; this shares the buffer and header FitsAsdf owns
// (hence manageHead_ = 0) and points at the next plane in it.
class FitsAsdfNext : public FitsFile {
 public:
  FitsAsdfNext(FitsFile*);
};

// asdf datatype -> the BITPIX fitsy holds it in, or 0. fitsy accepts
// {8,16,-16,32,64,-32,-64}, where -16 is its own (non FITS) code for
// unsigned 16 bit.
extern int asdfBitpix(const char*);

// datatypes with no BITPIX of their own, paired with the narrowest one
// that holds every value exactly, or NULL
extern const char* asdfWiden(const char*);

#endif
