// Copyright (C) 1999-2021
// Smithsonian Astrophysical Observatory, Cambridge, MA, USA
// For conditions of distribution and use, see copyright notice in "copyright"

// libasdf first: card.h, pulled in by the fitsy headers, carries a
// `using namespace std'
#include <asdf.h>

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <new>
#include <string>
#include <vector>
using namespace std;

#include "asdffile.h"
#include "head.h"
#include "util.h"

// ASDF binary block magic (asdf-standard file_layout)
static const char ASDF_BLK[4] = {(char)0xd3, 'B', 'L', 'K'};

// ---------------------------------------------------------------- datatypes

int asdfBitpix(const char* datatype)
{
  if (!datatype)
    return 0;
  if (!strcmp(datatype,"float32")) return -32;
  if (!strcmp(datatype,"float64")) return -64;
  if (!strcmp(datatype,"int16")) return 16;
  if (!strcmp(datatype,"uint16")) return -16;
  if (!strcmp(datatype,"int32")) return 32;
  if (!strcmp(datatype,"int64")) return 64;
  if (!strcmp(datatype,"uint8")) return 8;
  if (!strcmp(datatype,"bool8")) return 8;
  return 0;
}

const char* asdfWiden(const char* datatype)
{
  // lossless only; uint64 has no lossless target and is refused
  if (!datatype)
    return NULL;
  if (!strcmp(datatype,"float16")) return "float32";
  if (!strcmp(datatype,"int8")) return "int16";
  if (!strcmp(datatype,"uint32")) return "int64";
  return NULL;
}

// ---------------------------------------------------------- FitsAsdfNode

FitsAsdfNode::FitsAsdfNode()
{
  source_ = -1;
  big_ = -1;
  ndim_ = 0;
  for (int ii=0; ii<ASDF_MAXDIM; ii++)
    dim_[ii] = 0;
  unsupported_ = ASDF_OK;
  hasMaskScalar_ = 0;
  maskScalar_ = 0;
}

// ---------------------------------------------------------- FitsAsdfFile

FitsAsdfFile::FitsAsdfFile(const char* fn)
{
  file_ = NULL;
  treeRead_ = 0;
  valid_ = 0;

  if (!fn || !*fn) {
    setError("no file name");
    return;
  }
  fn_ = fn;

  file_ = asdf_open_file(fn, "r");
  if (!file_) {
    setError("unable to open %s as ASDF", fn);
    return;
  }

  enumNdarrays();
  valid_ = 1;
}

FitsAsdfFile::~FitsAsdfFile()
{
  if (file_)
    asdf_close(file_);
}

void FitsAsdfFile::setError(const char* fmt, ...)
{
  char buf[1024];
  va_list args;
  va_start(args, fmt);
  vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  error_ = buf;
}

// The tree text verbatim: everything before the first binary block.
// libasdf has no public accessor for this yet. Only the Tcl side uses it
// (the header viewer, and the GWCS subtree handed to AST), so it is read
// on demand.
const char* FitsAsdfFile::tree()
{
  if (treeRead_)
    return tree_.c_str();
  treeRead_ = 1;

  FILE* ff = fopen(fn_.c_str(), "rb");
  if (!ff)
    return tree_.c_str();

  char buf[65536];
  size_t got;
  while ((got = fread(buf, 1, sizeof(buf), ff)) > 0) {
    // the magic may straddle two reads
    size_t from = tree_.size() >= 3 ? tree_.size()-3 : 0;
    tree_.append(buf, got);
    size_t at = tree_.find(string(ASDF_BLK, 4), from);
    if (at != string::npos) {
      tree_.resize(at);
      break;
    }
  }
  fclose(ff);

  return tree_.c_str();
}

size_t FitsAsdfFile::treeSize()
{
  tree();
  return tree_.size();
}

const FitsAsdfNode* FitsAsdfFile::node(int ii)
{
  if (ii < 0 || ii >= (int)nodes_.size())
    return NULL;
  return &nodes_[ii];
}

const FitsAsdfNode* FitsAsdfFile::find(const char* path)
{
  if (!path)
    return NULL;
  for (size_t ii=0; ii<nodes_.size(); ii++)
    if (nodes_[ii].path_ == path)
      return &nodes_[ii];
  return NULL;
}

// Every block backed ndarray anywhere in the tree. Inline ndarrays are
// skipped, as before: they are small, and only the GWCS path uses them.
void FitsAsdfFile::enumNdarrays()
{
  asdf_value_t* root = asdf_get_value(file_, "");
  if (!root)
    return;

  asdf_find_iter_t* iter = asdf_find_iter_init(root, asdf_value_is_ndarray);
  while (iter && asdf_value_find_iter_next(&iter)) {
    asdf_ndarray_t* nd = NULL;
    if (asdf_value_as_ndarray(iter->value, &nd) != ASDF_VALUE_OK || !nd)
      continue;

    if (asdf_ndarray_storage(nd) != ASDF_ARRAY_STORAGE_INTERNAL) {
      asdf_ndarray_destroy(nd);
      continue;
    }

    // a YAML alias to an ndarray is reached again, under the anchor's path
    const char* pp = asdf_value_path(iter->value);
    string path = pp ? (*pp == '/' ? pp+1 : pp) : "";
    if (find(path.c_str())) {
      asdf_ndarray_destroy(nd);
      continue;
    }

    FitsAsdfNode nn;
    nn.path_ = path;
    nn.source_ = (int)nd->source;
    const char* dt = asdf_scalar_datatype_to_string(nd->datatype.type);
    nn.datatype_ = dt ? dt : "";

    switch (nd->byteorder) {
    case ASDF_BYTEORDER_BIG:
      nn.big_ = 1;
      break;
    case ASDF_BYTEORDER_LITTLE:
      nn.big_ = 0;
      break;
    default:
      // required whenever source is present; refuse rather than guess
      nn.unsupported_ |= ASDF_BYTEORDER;
      break;
    }

    nn.ndim_ = (int)nd->ndim;
    for (int ii=0; ii<nn.ndim_ && ii<ASDF_MAXDIM; ii++)
      nn.dim_[ii] = (long long)nd->shape[ii];

    if (nd->strides)
      nn.unsupported_ |= ASDF_STRIDES;
    if (nd->offset)
      nn.unsupported_ |= ASDF_OFFSET;

    // A scalar mask names the missing value outright. An ndarray mask is
    // enumerated as its own <path>/mask node; anything else is unusable.
    string mpath = nn.path_ + "/mask";
    asdf_value_t* mv = asdf_get_value(file_, mpath.c_str());
    if (mv) {
      double dd;
      if (asdf_value_is_ndarray(mv))
	;
      else if (asdf_get_double(file_, mpath.c_str(), &dd) == ASDF_VALUE_OK) {
	nn.hasMaskScalar_ = 1;
	nn.maskScalar_ = dd;
      }
      else
	nn.unsupported_ |= ASDF_MASK;
      asdf_value_destroy(mv);
    }

    nodes_.push_back(nn);
    asdf_ndarray_destroy(nd);
  }
  if (iter)
    asdf_find_iter_destroy(iter);

  asdf_value_destroy(root);
}

int FitsAsdfFile::nblock()
{
  return (int)asdf_block_count(file_);
}

char* FitsAsdfFile::block(int source, size_t* size)
{
  *size = 0;

  if (source < 0 || source >= nblock()) {
    setError("block index %d out of range", source);
    return NULL;
  }

  asdf_block_t* bb = asdf_block_open(file_, (size_t)source);
  if (!bb) {
    setError("unable to open block %d", source);
    return NULL;
  }

  size_t nn = 0;
  const void* data = asdf_block_data(bb, &nn);
  if (!data) {
    asdf_block_close(bb);
    setError("unable to read block %d", source);
    return NULL;
  }

  char* dst = new (std::nothrow) char[nn ? nn : 1];
  if (!dst) {
    asdf_block_close(bb);
    setError("unable to allocate %zu bytes for block %d", nn, source);
    return NULL;
  }
  memcpy(dst, data, nn);
  asdf_block_close(bb);

  *size = nn;
  return dst;
}

// ------------------------------------------------------------------ mask
//
// asdf represents a numpy masked array as an explicit boolean `mask'
// ndarray beside the data; the core/ndarray schema also allows a scalar
// `mask', which "is used to represent missing values" - the FITS BLANK
// convention exactly. Both are mapped onto what FITS already does, so
// DS9's existing machinery does the work with nothing promoted to float:
//
//   integer data -> a BLANK sentinel value. FitsData keeps the native
//                   integer storage and substitutes NAN only at the
//                   getValueFloat() boundary, while minmax skips blank
//                   pixels (tksao/frame/fitsdata.C).
//   float data   -> NaN written into the masked pixels, which is the FITS
//                   convention for float nulls and costs nothing since
//                   the data is already floating point.
//
// Scalar matching is exact, deliberately unlike numpy's
// ma.masked_values() (rtol=1e-5): a fill value is stored exactly in
// practice, and exact equality is both predictable and what FITS BLANK
// means. A scalar mask of NaN means "NaN values are missing", so the test
// is isnan() rather than equality, the same way asdf's own reader does it.

static int64_t asdfLoadInt(const unsigned char* pp, int sz, int isSigned,
			   int big)
{
  uint64_t uu = 0;
  for (int ii=0; ii<sz; ii++) {
    unsigned char bb = big ? pp[ii] : pp[sz-1-ii];
    uu = (uu << 8) | bb;
  }

  if (isSigned && sz < 8) {
    uint64_t sign = (uint64_t)1 << (sz*8-1);
    if (uu & sign)
      uu |= ~((sign << 1) - 1);
  }

  return (int64_t)uu;
}

static void asdfStoreInt(unsigned char* dd, uint64_t uu, int sz, int big)
{
  for (int bb=0; bb<sz; bb++) {
    unsigned char cc = (unsigned char)(uu >> (8*bb));
    if (big)
      dd[sz-1-bb] = cc;
    else
      dd[bb] = cc;
  }
}

// Resolves a mask, rewriting `data' in place and in its own byteorder.
// Returns 1 when something was masked, 0 when nothing was (or the mask
// cannot be represented - the caller then loads the array unchanged).
// *hasBlank/*blank carry the BLANK keyword for an integer type.
static int asdfApplyMask(char* data, size_t dataLen, const char* datatype,
			 int big, const char* mask, size_t maskLen,
			 int haveScalar, double scalar,
			 int* hasBlank, long long* blank, string& err)
{
  *hasBlank = 0;
  *blank = 0;

  int sz, isFloat, isSigned = 0;
  int64_t tmin = 0, tmax = 0;

  // Only the datatypes with a real FITS BITPIX are eligible. That
  // excludes uint16 (it would need FITS's BZERO unsigned offset
  // convention) and the widened types - those load without nulls rather
  // than as if they had none.
  if (!strcmp(datatype,"uint8")) {
    sz = 1; isFloat = 0; isSigned = 0; tmin = 0; tmax = 255;
  }
  else if (!strcmp(datatype,"int16")) {
    sz = 2; isFloat = 0; isSigned = 1; tmin = -32768; tmax = 32767;
  }
  else if (!strcmp(datatype,"int32")) {
    sz = 4; isFloat = 0; isSigned = 1;
    tmin = -2147483648LL; tmax = 2147483647LL;
  }
  else if (!strcmp(datatype,"int64")) {
    // Not INT64_MIN/MAX. tmin/tmax are only ever used as BLANK candidates
    // (and to ask whether an unmasked pixel already holds one), and FITS
    // carries BLANK as an int the whole way down - FitsImageHDU::blank_
    // and FitsData::blank_ are both int - so a wider sentinel would be
    // truncated on the way in and then stop comparing equal to the very
    // pixels it was chosen to mark. For every other type here the full
    // range already fits, so this is the only one that has to be clamped.
    sz = 8; isFloat = 0; isSigned = 1; tmin = INT_MIN; tmax = INT_MAX;
  }
  else if (!strcmp(datatype,"float32")) {
    sz = 4; isFloat = 1;
  }
  else if (!strcmp(datatype,"float64")) {
    sz = 8; isFloat = 1;
  }
  else
    return 0;

  if (dataLen % sz)
    return 0;
  size_t nn = dataLen/sz;

  // The mask must be broadcastable to the data's shape. Only the
  // degenerate cases are accepted: one element per pixel, or a single
  // element covering everything. Anything else is refused rather than
  // guessed at, since a wrong partial broadcast would mismark real
  // pixels.
  if (mask && maskLen != nn && maskLen != 1) {
    err = "mask is not broadcastable to the array's shape";
    return 0;
  }

  int scalarIsNan = haveScalar && scalar != scalar;
  int64_t iscalar = haveScalar ? (int64_t)scalar : 0;
  int rewrite = 0;

  if (!isFloat) {
    int64_t common = 0;
    size_t nmasked = 0;
    int uniq = 1, minUsed = 0, maxUsed = 0;

    for (size_t ii=0; ii<nn; ii++) {
      int64_t vv = asdfLoadInt((const unsigned char*)data+ii*sz, sz,
			       isSigned, big);
      int mm = mask ? (maskLen == 1 ? mask[0] : mask[ii]) != 0
	: (vv == iscalar);
      if (mm) {
	if (!nmasked)
	  common = vv;
	else if (vv != common)
	  uniq = 0;
	nmasked++;
      }
      else {
	if (vv == tmin) minUsed = 1;
	if (vv == tmax) maxUsed = 1;
      }
    }

    if (!nmasked)
      return 0;

    if (haveScalar) {
      // a scalar mask names the sentinel outright, so there is nothing to
      // choose and nothing to rewrite - but it still has to be a value
      // FITS can carry, and there is no second choice available if it is
      // not, since rewriting the pixels would discard what the file said
      if (iscalar < INT_MIN || iscalar > INT_MAX) {
	err = "scalar mask value is outside the FITS BLANK range";
	return 0;
      }
      *blank = iscalar;
      *hasBlank = 1;
    }
    else {
      // Can the masked pixels' own value serve as BLANK? Only if they all
      // share one value and no unmasked pixel also carries it.
      int commonUsed = 0;
      if (uniq && mask) {
	for (size_t ii=0; ii<nn; ii++) {
	  if ((maskLen == 1 ? mask[0] : mask[ii]) != 0)
	    continue;
	  if (asdfLoadInt((const unsigned char*)data+ii*sz, sz,
			  isSigned, big) == common) {
	    commonUsed = 1;
	    break;
	  }
	}
      }

      // The masked pixels' own value is the first choice - it needs no
      // rewrite - but only if FITS can hold it; otherwise fall through to
      // the candidates, which are in range by construction.
      if (uniq && !commonUsed && common >= INT_MIN && common <= INT_MAX)
	*blank = common;
      else {
	int64_t c0 = isSigned ? tmin : tmax;
	int64_t c1 = isSigned ? tmax : tmin;
	int u0 = isSigned ? minUsed : maxUsed;
	int u1 = isSigned ? maxUsed : minUsed;

	if (!u0)
	  *blank = c0;
	else if (!u1)
	  *blank = c1;
	else {
	  err = "no unused value available for BLANK";
	  return 0;
	}
	rewrite = 1;
      }
      *hasBlank = 1;
    }
  }

  // second pass: mark the masked pixels
  size_t nout = 0;
  for (size_t ii=0; ii<nn; ii++) {
    unsigned char* dd = (unsigned char*)data + ii*sz;

    int mm;
    if (mask)
      mm = (maskLen == 1 ? mask[0] : mask[ii]) != 0;
    else if (isFloat) {
      uint64_t raw = 0;
      for (int bb=0; bb<sz; bb++) {
	unsigned char cc = big ? dd[bb] : dd[sz-1-bb];
	raw = (raw << 8) | cc;
      }
      if (sz == 4) {
	float ff;
	uint32_t r32 = (uint32_t)raw;
	memcpy(&ff, &r32, 4);
	mm = scalarIsNan ? (ff != ff) : (ff == (float)scalar);
      }
      else {
	double ff;
	memcpy(&ff, &raw, 8);
	mm = scalarIsNan ? (ff != ff) : (ff == scalar);
      }
    }
    else
      mm = asdfLoadInt(dd, sz, isSigned, big) == iscalar;

    if (!mm)
      continue;
    nout++;

    if (isFloat)
      asdfStoreInt(dd, sz == 4 ? 0x7fc00000ULL : 0x7ff8000000000000ULL,
		   sz, big);
    else if (rewrite)
      asdfStoreInt(dd, (uint64_t)*blank, sz, big);
  }

  if (!nout) {
    *hasBlank = 0;
    *blank = 0;
    return 0;
  }

  return 1;
}

// -------------------------------------------------------------- FitsAsdf

FitsAsdf::FitsAsdf(const char* fn, const char* path)
{
  asdf_ = NULL;

  // the ordinary fitsy name/extension/section parse
  parse(fn);
  if (!valid_)
    return;
  valid_ = 0;

  if (!pName_) {
    asdfError_ = "no file name";
    return;
  }

  asdf_ = new FitsAsdfFile(pName_);
  if (!asdf_->isValid()) {
    asdfError_ = asdf_->error() ? asdf_->error() : "unable to read";
    return;
  }

  if (!build(path && *path ? path : "roman/data"))
    return;

  valid_ = 1;
}

FitsAsdf::~FitsAsdf()
{
  // allocated by asdf_ndarray_read_all()
  if (data_)
    asdf_free(data_);

  if (asdf_)
    delete asdf_;
}

// ------------------------------------------------------------ FitsAsdfNext

FitsAsdfNext::FitsAsdfNext(FitsFile* prev)
{
  primary_ = prev->primary();
  managePrimary_ = 0;

  head_ = prev->head();
  manageHead_ = 0;

  FitsImageHDU* hdu = head_ ? (FitsImageHDU*)head_->hdu() : NULL;
  if (!hdu)
    return;

  data_ = (char*)prev->data() + hdu->imgbytes();
  dataSize_ = 0;
  dataSkip_ = 0;

  ext_ = prev->ext();
  inherit_ = prev->inherit();
  byteswap_ = prev->byteswap();
  endian_ = prev->endian();

  pcoord_ = prev->pcoord();
  pxvalid_ = prev->pxvalid();
  pxmin_ = prev->pxmin();
  pxmax_ = prev->pxmax();
  pyvalid_ = prev->pyvalid();
  pymin_ = prev->pymin();
  pymax_ = prev->pymax();
  pzvalid_ = prev->pzvalid();
  pzmin_ = prev->pzmin();
  pzmax_ = prev->pzmax();
  pbvalid_ = prev->pbvalid();
  pblock_ = prev->pblock();

  valid_ = 1;
}

// Reads the mask of a 2-d array, if it has one: either an ndarray, whose
// elements are returned as bytes, or a scalar missing-value sentinel.
// Returns 0 if there is no usable mask.
static int asdfReadMask(asdf_file_t* file, const string& path,
			vector<char>& mask,
			int* haveScalar, double* scalar, string& err)
{
  mask.clear();
  *haveScalar = 0;

  string mpath = path + "/mask";

  asdf_ndarray_t* mnd = NULL;
  if (asdf_get_ndarray(file, mpath.c_str(), &mnd) == ASDF_VALUE_OK && mnd) {
    size_t nn = (size_t)asdf_ndarray_size(mnd);
    asdf_scalar_datatype_t mt = mnd->datatype.type;

    if (mt == ASDF_DATATYPE_BOOL8 || mt == ASDF_DATATYPE_UINT8) {
      // one byte per element, so no conversion or byte order to handle
      // (and libasdf's read_all does not convert bool8 at all yet)
      size_t len = 0;
      const char* data = (const char*)asdf_ndarray_data(mnd, &len);
      if (data && len >= nn)
	mask.assign(data, data+nn);
    }
    else {
      void* buf = NULL;
      if (asdf_ndarray_read_all(mnd, ASDF_DATATYPE_UINT8, &buf) ==
	  ASDF_NDARRAY_OK && buf) {
	mask.assign((char*)buf, (char*)buf+nn);
	asdf_free(buf);
      }
    }
    asdf_ndarray_destroy(mnd);

    if (mask.empty()) {
      err = "unreadable mask";
      return 0;
    }
    return 1;
  }

  if (asdf_get_double(file, mpath.c_str(), scalar) == ASDF_VALUE_OK) {
    *haveScalar = 1;
    return 1;
  }

  return 0;
}

int FitsAsdf::build(const char* path)
{
  asdf_file_t* file = asdf_->file();

  asdf_ndarray_t* nd = NULL;
  if (asdf_get_ndarray(file, path, &nd) != ASDF_VALUE_OK || !nd) {
    asdfError_ = string("could not find ndarray ") + path;
    return 0;
  }

  if (nd->byteorder != ASDF_BYTEORDER_BIG &&
      nd->byteorder != ASDF_BYTEORDER_LITTLE) {
    asdf_ndarray_destroy(nd);
    asdfError_ = string("ndarray has no byteorder ") + path;
    return 0;
  }
  if (nd->strides || nd->offset) {
    asdf_ndarray_destroy(nd);
    asdfError_ = string("unsupported ndarray view ") + path;
    return 0;
  }

  // ASDF shape is row major, slowest axis first; FITS wants NAXIS1 first
  if (nd->ndim != 2 && nd->ndim != 3) {
    asdf_ndarray_destroy(nd);
    asdfError_ = string("unsupported ndarray rank ") + path;
    return 0;
  }
  for (uint32_t ii=0; ii<nd->ndim; ii++) {
    if (nd->shape[ii] == 0 || nd->shape[ii] > (uint64_t)INT_MAX) {
      asdf_ndarray_destroy(nd);
      asdfError_ = string("unsupported ndarray shape ") + path;
      return 0;
    }
  }

  int xdim, ydim, zdim = 1;
  if (nd->ndim == 2) {
    ydim = (int)nd->shape[0];
    xdim = (int)nd->shape[1];
  }
  else {
    zdim = (int)nd->shape[0];
    ydim = (int)nd->shape[1];
    xdim = (int)nd->shape[2];
  }

  const char* dt = asdf_scalar_datatype_to_string(nd->datatype.type);
  string datatype = dt ? dt : "";
  const char* widen = NULL;
  int bitpix = asdfBitpix(datatype.c_str());
  if (!bitpix) {
    widen = asdfWiden(datatype.c_str());
    if (!widen) {
      asdf_ndarray_destroy(nd);
      asdfError_ = string("unsupported ndarray datatype ") +
	(datatype.empty() ? "(non-scalar)" : datatype);
      return 0;
    }
    bitpix = asdfBitpix(widen);
  }

  // libasdf decompresses, widens, and converts to native byte order
  asdf_scalar_datatype_t target = widen ?
    asdf_scalar_datatype_from_string(widen) : ASDF_DATATYPE_SOURCE;
  void* buf = NULL;
  asdf_ndarray_err_t rr = asdf_ndarray_read_all(nd, target, &buf);
  asdf_ndarray_destroy(nd);
  if (rr != ASDF_NDARRAY_OK || !buf) {
    const char* ee = asdf_error(file);
    asdfError_ = string("unable to read ndarray ") + path +
      (ee ? string(": ") + ee : string());
    return 0;
  }

  size_t want = (size_t)(abs(bitpix)/8) * xdim * ydim * zdim;

  // Masks are resolved for a plain 2-d array of a directly representable
  // datatype only: a widened array's BITPIX is no longer the datatype the
  // mask was written against, and FITS BLANK is per file, not per plane.
  int hasBlank = 0;
  long long blank = 0;

  if (zdim == 1 && !widen) {
    vector<char> mask;
    int haveScalar = 0;
    double scalar = 0;
    string err;

    if (asdfReadMask(file, path, mask, &haveScalar, &scalar, err))
      asdfApplyMask((char*)buf, want, datatype.c_str(), !lsb(),
		    mask.empty() ? NULL : &mask[0], mask.size(),
		    haveScalar, scalar, &hasBlank, &blank, err);

    if (!err.empty()) {
      string mm = string("Warning: ASDF: ") + err +
	", loading without null values";
      internalError(mm.c_str());
    }
  }

  pWidth_ = xdim;
  pHeight_ = ydim;
  pDepth_ = zdim;
  pBitpix_ = bitpix;
  pSkip_ = 0;
  pArch_ = NATIVE;
  pEncoding_ = RAW;

  if (!validParams()) {
    asdf_free(buf);
    asdfError_ = string("unsupported ndarray shape or datatype ") + path;
    return 0;
  }

  head_ = new FitsHead(pWidth_, pHeight_, pDepth_, pBitpix_);
  if (!head_->isValid()) {
    asdf_free(buf);
    delete head_;
    head_ = NULL;
    asdfError_ = "unable to build a header";
    return 0;
  }
  manageHead_ = 1;

  if (hasBlank && pBitpix_ > 0) {
    head_->appendInteger("BLANK", (int)blank, "Null pixel value");
    head_->updateHDU();
  }

  data_ = buf;
  dataSize_ = want;
  dataSkip_ = 0;

  setByteSwap();

  return 1;
}
