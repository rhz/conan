/////////////////////////////////
// Generated header: TypeManip.h
// Forwards to the appropriate code
// that works on the detected compiler
// Generated on Mon Sep 30 23:14:48 2002
/////////////////////////////////

#ifdef LOKI_USE_REFERENCE
#	include "Reference/TypeManip.h"
#else
#	if (__BORLANDC__ >= 0x560)
#		include "Borland/TypeManip.h"
#	else
#		include "Reference/TypeManip.h"
#	endif
#endif
